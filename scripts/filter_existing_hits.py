import argparse
import csv
import hashlib
import io
import json
import re
import sys
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path
VERSION = '0.1.0'
INPUT_FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
MAPPING_FIELDS = ['reference_id', 'arg_subtype', 'arg_type']
ANNOTATION_FIELDS = ['orf_id', 'arg_subtype', 'arg_type']
EVIDENCE_FIELDS = ['pident', 'alignment_length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'query_coverage_percent', 'source_line']

class InputError(ValueError):
    pass

def error(path, number, message):
    raise InputError(f'{path}:{number}: {message}')

def number(value, label):
    if isinstance(value, bool) or value is None:
        raise InputError(f'{label} must be a finite number')
    try:
        result = Decimal(str(value))
    except (InvalidOperation, ValueError):
        raise InputError(f'{label} must be a finite number') from None
    if not result.is_finite():
        raise InputError(f'{label} must be a finite number')
    return result

def percentage(value, label):
    result = number(value, label)
    if not Decimal(0) <= result <= Decimal(100):
        raise InputError(f'{label} must be between 0 and 100')
    return result

def integer(value, label, minimum=1):
    if not re.fullmatch('[0-9]+', value):
        raise InputError(f'{label} must be an integer >= {minimum}')
    result = int(value)
    if result < minimum:
        raise InputError(f'{label} must be an integer >= {minimum}')
    return result

def fingerprint(path, content):
    return {'path': str(path.resolve()), 'sha256': hashlib.sha256(content).hexdigest()}

def reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise InputError(f'duplicate JSON key: {key}')
        result[key] = value
    return result

def read_settings(path, max_evalue):
    content = path.read_bytes()
    try:
        config = json.loads(content.decode('utf-8-sig'), object_pairs_hook=reject_duplicate_keys)
        settings = config['study_records_not_executed_by_this_wrapper']['diamond']
        identity = percentage(settings['reported_final_identity_percent'], 'final identity threshold')
        coverage = percentage(settings['query_coverage_percent'], 'query coverage threshold')
    except (UnicodeDecodeError, json.JSONDecodeError, KeyError, TypeError) as exc:
        raise InputError(f'{path}: invalid or missing central DIAMOND settings: {exc}') from None
    evalue = number(settings.get('final_max_evalue', '1e-5') if max_evalue is None else max_evalue, 'maximum E-value')
    if evalue < 0:
        raise InputError('maximum E-value must be nonnegative')
    return ((identity, coverage, evalue), fingerprint(path, content))

def read_mapping(path):
    if path is None:
        return (None, None)
    content = path.read_bytes()
    try:
        stream = io.StringIO(content.decode('utf-8-sig'), newline='')
    except UnicodeDecodeError:
        raise InputError(f'{path}: mapping must be UTF-8') from None
    reader = csv.reader(stream, delimiter='\t')
    if next(reader, None) != MAPPING_FIELDS:
        raise InputError(f"{path}: expected exact header: {' '.join(MAPPING_FIELDS)}")
    mapping = {}
    for row in reader:
        if len(row) != 3 or any((not item or item != item.strip() for item in row)):
            error(path, reader.line_num, 'expected three nonempty fields without surrounding whitespace')
        reference_id, subtype, arg_type = row
        if any((character.isspace() for character in reference_id)):
            error(path, reader.line_num, 'reference_id must be an exact nonempty token')
        if reference_id in mapping:
            error(path, reader.line_num, f'duplicate reference_id {reference_id!r}')
        mapping[reference_id] = {'arg_subtype': subtype, 'arg_type': arg_type}
    return (mapping, fingerprint(path, content))

def parse_hit(row, path, line):
    if len(row) != len(INPUT_FIELDS):
        error(path, line, 'expected exactly 14 tab-separated columns, without a header')
    if any((value != value.strip() for value in row)):
        error(path, line, 'leading or trailing whitespace in a field')
    raw = dict(zip(INPUT_FIELDS, row))
    for field in ('qseqid', 'sseqid'):
        if not raw[field] or any((character.isspace() for character in raw[field])):
            error(path, line, f'{field} must be a nonempty token')
    try:
        values = {field: integer(raw[field], field) for field in ('length', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen')}
        for field in ('mismatch', 'gapopen'):
            values[field] = integer(raw[field], field, minimum=0)
        values['pident'] = percentage(raw['pident'], 'pident')
        values['evalue'] = number(raw['evalue'], 'evalue')
        values['bitscore'] = number(raw['bitscore'], 'bitscore')
        if values['evalue'] < 0:
            raise InputError('evalue must be nonnegative')
        if max(values['qstart'], values['qend']) > values['qlen']:
            raise InputError('query coordinates exceed qlen')
        if max(values['sstart'], values['send']) > values['slen']:
            raise InputError('subject coordinates exceed slen')
    except InputError as exc:
        error(path, line, str(exc))
    span = abs(values['qend'] - values['qstart']) + 1
    with localcontext() as context:
        context.prec = 50
        query_coverage = Decimal(span) * 100 / Decimal(values['qlen'])
    return {'raw': raw, 'values': values, 'query_span': span, 'query_coverage': query_coverage, 'line': line}

def rank(hit):
    values = hit['values']
    return (values['evalue'], -values['bitscore'], -values['length'], hit['line'])

def read_hits(path, thresholds):
    identity, coverage, evalue = thresholds
    digest = hashlib.sha256()
    queries = set()
    best = {}
    counts = {key: 0 for key in ('input_rows', 'qualifying_rows', 'rejected_rows', 'failed_identity_rows', 'failed_query_coverage_rows', 'failed_evalue_rows')}
    with path.open('rb') as handle:
        for line, raw_bytes in enumerate(handle, 1):
            digest.update(raw_bytes)
            try:
                text = raw_bytes.decode('utf-8-sig' if line == 1 else 'utf-8').rstrip('\r\n')
            except UnicodeDecodeError:
                error(path, line, 'input must be UTF-8')
            hit = parse_hit(text.split('\t'), path, line)
            values, raw = (hit['values'], hit['raw'])
            counts['input_rows'] += 1
            queries.add(raw['qseqid'])
            failures = (('failed_identity_rows', values['pident'] < identity), ('failed_query_coverage_rows', Decimal(hit['query_span']) * 100 < coverage * values['qlen']), ('failed_evalue_rows', values['evalue'] > evalue))
            for name, failed in failures:
                counts[name] += int(failed)
            if any((failed for _, failed in failures)):
                counts['rejected_rows'] += 1
                continue
            counts['qualifying_rows'] += 1
            previous = best.get(raw['qseqid'])
            if previous is None or rank(hit) < rank(previous):
                best[raw['qseqid']] = hit
    counts.update(input_queries=len(queries), qualifying_queries=len(best), selected_rows=len(best), queries_without_qualifying_rows=len(queries) - len(best), nonselected_qualifying_rows=counts['qualifying_rows'] - len(best))
    source = {'path': str(path.resolve()), 'sha256': digest.hexdigest()}
    return (sorted(best.values(), key=lambda item: item['line']), counts, source)

def coverage_text(value):
    with localcontext() as context:
        context.prec = 50
        text = format(value.quantize(Decimal('0.000000000001')), 'f')
    return text.rstrip('0').rstrip('.') or '0'

def write_tsv(path, fields, rows):
    with path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def evidence_row(hit, mapping):
    raw = hit['raw']
    row = {'orf_id': raw['qseqid'], 'reference_id': raw['sseqid']}
    if mapping is not None:
        row.update(mapping[raw['sseqid']])
    row.update({field: raw[field] for field in INPUT_FIELDS if field not in {'qseqid', 'sseqid', 'length'}})
    row['alignment_length'] = raw['length']
    row['query_coverage_percent'] = coverage_text(hit['query_coverage'])
    row['source_line'] = hit['line']
    return row

def run(args):
    requested_output = Path(args.output_dir)
    if requested_output.exists() or requested_output.is_symlink():
        raise InputError(f'output directory already exists; refusing overwrite: {requested_output}')
    output = requested_output.resolve()
    if output.exists():
        raise InputError(f'output directory already exists; refusing overwrite: {output}')
    config_path, input_path = (Path(args.config), Path(args.input))
    mapping_path = Path(args.mapping) if args.mapping is not None else None
    thresholds, config_source = read_settings(config_path, args.max_evalue)
    mapping, mapping_source = read_mapping(mapping_path)
    selected, counts, hit_source = read_hits(input_path, thresholds)
    if mapping is not None:
        missing = sorted({hit['raw']['sseqid'] for hit in selected} - mapping.keys())
        if missing:
            raise InputError(f'mapping lacks {len(missing)} selected reference IDs; first: {missing[0]!r}')
    fields = ['orf_id', 'reference_id'] + ([] if mapping is None else ['arg_subtype', 'arg_type']) + EVIDENCE_FIELDS
    audit = {'completion_state': 'complete', 'tool': 'filter_existing_hits.py', 'version': VERSION, 'status': 'New submission-preparation bookkeeping; historical execution and study results not verified', 'input_columns': INPUT_FIELDS, 'inputs': {'alignments': hit_source, 'config': config_source, 'reference_mapping': mapping_source}, 'thresholds': {'identity_min_percent': str(thresholds[0]), 'query_coverage_min_percent': str(thresholds[1]), 'evalue_max': str(thresholds[2])}, 'threshold_sources': {'identity_and_coverage': 'central config study_records_not_executed_by_this_wrapper.diamond', 'evalue': '--max-evalue override' if args.max_evalue is not None else 'central config final_max_evalue; fallback 1e-5 from Methods 4.7'}, 'query_coverage_formula': '100 * (abs(qend - qstart) + 1) / qlen; qlen is required', 'coverage_display_decimal_places': 12, 'selection_order': ['evalue ascending', 'bitscore descending', 'alignment length descending', 'original input line ascending'], 'counts': counts, 'mapping_rows': None if mapping is None else len(mapping), 'count_note': 'Failure counts overlap; rejected_rows counts each rejected row once', 'limitation': 'Ranks available reported rows only. Cannot recover or rank unreported targets excluded by -k 1.', 'outputs': {'selected_hits': 'selected_hits.tsv', 'arg_annotations': None if mapping is None else 'arg_annotations.tsv'}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(exist_ok=False)
    write_tsv(output / 'selected_hits.tsv', fields, (evidence_row(hit, mapping) for hit in selected))
    if mapping is not None:
        write_tsv(output / 'arg_annotations.tsv', ANNOTATION_FIELDS, ({'orf_id': hit['raw']['qseqid'], **mapping[hit['raw']['sseqid']]} for hit in selected))
    with (output / 'audit.json').open('x', encoding='utf-8') as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    return audit

def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('--config', required=True, type=Path, help='Central package config.json')
    result.add_argument('--input', required=True, type=Path, help='Existing headerless 14-column alignment TSV for one sample/run')
    result.add_argument('--mapping', type=Path, help='Optional exact reference_id,arg_subtype,arg_type TSV')
    result.add_argument('--max-evalue', default=None, help='Optional override of config final_max_evalue (fallback: 1e-5)')
    result.add_argument('--output-dir', required=True, type=Path, help='New directory; existing paths are refused')
    return result

def main():
    try:
        audit = run(parser().parse_args())
    except (InputError, OSError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2
    print(json.dumps(audit['counts'], sort_keys=True))
    return 0
if __name__ == '__main__':
    sys.exit(main())
