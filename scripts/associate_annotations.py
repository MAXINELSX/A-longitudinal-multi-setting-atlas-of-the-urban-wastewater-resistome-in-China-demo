import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote
VERSION = '0.1.0'
COMPARTMENTS = {'chromosome', 'plasmid', 'virus', 'unclassified'}
OUTPUT_FIELDS = ['sample_id', 'orf_id', 'arg_subtype', 'arg_type', 'contig_id', 'compartment', 'compartment_status', 'contig_taxid', 'contig_taxon_rank', 'contig_taxon_name', 'host_species', 'host_status']

class InputError(ValueError):
    pass

def fail(path, line, message):
    raise InputError(f'{path}:{line}: {message}')

def parse_attributes(text, path, line):
    result = {}
    for item in text.strip().split(';'):
        if not item:
            continue
        if '=' not in item:
            fail(path, line, 'attribute must be key=value')
        key, value = item.split('=', 1)
        key, value = (unquote(key), unquote(value))
        if key in result:
            fail(path, line, f'duplicate attribute {key!r}')
        result[key] = value
    return result

def positive_integer(value, path, line, label):
    if not re.fullmatch('[0-9]+', value) or int(value) < 1:
        fail(path, line, f'{label} must be a positive integer')
    return int(value)

def read_gff(path):
    genes = {}
    contigs = set()
    with path.open(encoding='utf-8-sig') as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip('\r\n')
            if line == '##FASTA':
                break
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) != 9:
                fail(path, number, 'expected 9 GFF columns')
            if cols[2] != 'CDS':
                continue
            contig = cols[0]
            if not contig or any((character.isspace() for character in contig)):
                fail(path, number, 'CDS seqid must be a nonempty token')
            start = positive_integer(cols[3], path, number, 'start')
            end = positive_integer(cols[4], path, number, 'end')
            if end < start or cols[6] not in {'+', '-'}:
                fail(path, number, 'invalid CDS interval or strand')
            attrs = parse_attributes(cols[8], path, number)
            identifier = attrs.get('ID', '')
            if not identifier:
                fail(path, number, 'CDS requires an ID attribute')
            if identifier in genes:
                fail(path, number, f'duplicate CDS ID {identifier!r}')
            genes[identifier] = (contig, start, end, cols[6])
            contigs.add(contig)
    if not genes:
        raise InputError(f'{path}: no CDS records')
    return (genes, contigs)

def read_protein_headers(path, genes):
    proteins = {}
    seen_gff = set()
    have_header = False
    with path.open(encoding='utf-8-sig') as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip('\r\n')
            if not line:
                continue
            if not line.startswith('>'):
                if not have_header:
                    fail(path, number, 'content precedes first FASTA header')
                continue
            have_header = True
            parts = [part.strip() for part in line[1:].split('#')]
            if not parts[0]:
                fail(path, number, 'empty FASTA identifier')
            orf_id = parts[0].split()[0]
            if orf_id in proteins:
                fail(path, number, f'duplicate FASTA identifier {orf_id!r}')
            if len(parts) == 5:
                attrs = parse_attributes(parts[4], path, number)
                header_coordinates = (positive_integer(parts[1], path, number, 'header start'), positive_integer(parts[2], path, number, 'header end'), {'1': '+', '-1': '-'}.get(parts[3]))
                if header_coordinates[2] is None:
                    fail(path, number, 'header strand must be 1 or -1')
            elif len(parts) == 1:
                tokens = parts[0].split(maxsplit=1)
                if len(tokens) != 2:
                    fail(path, number, 'FASTA header requires explicit ID= attribute')
                attrs = parse_attributes(tokens[1], path, number)
                header_coordinates = None
            else:
                fail(path, number, 'expected Prodigal header or token ID=... header')
            gff_id = attrs.get('ID', '')
            if not gff_id or gff_id not in genes:
                fail(path, number, f'header ID {gff_id!r} is absent from GFF CDS records')
            if gff_id in seen_gff:
                fail(path, number, f'multiple FASTA headers reference GFF ID {gff_id!r}')
            if header_coordinates is not None and header_coordinates != genes[gff_id][1:]:
                fail(path, number, f'coordinates/strand disagree with GFF for {gff_id!r}')
            proteins[orf_id] = genes[gff_id][0]
            seen_gff.add(gff_id)
    missing = set(genes) - seen_gff
    if missing:
        raise InputError(f'{path}: {len(missing)} GFF CDS IDs lack a protein header')
    return proteins

def read_tsv(path, required, optional, sample_id):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        if fields is None or len(fields) != len(set(fields)):
            raise InputError(f'{path}: missing or duplicate TSV column names')
        if not required.issubset(fields):
            raise InputError(f"{path}: required columns: {', '.join(sorted(required))}")
        unknown = set(fields) - required - optional - {'sample_id'}
        if unknown:
            raise InputError(f"{path}: unsupported columns: {', '.join(sorted(unknown))}")
        rows = []
        for row in reader:
            number = reader.line_num
            if None in row or any((value is None for value in row.values())):
                fail(path, number, 'ragged TSV row')
            if any((row[key] == '' for key in required)):
                fail(path, number, 'required field is empty')
            if any((value != value.strip() for value in row.values())):
                fail(path, number, 'leading/trailing whitespace in a field')
            if 'sample_id' in row and row['sample_id'] != sample_id:
                fail(path, number, f'sample_id differs from selected sample {sample_id!r}')
            rows.append((number, row))
    return rows

def read_args(path, sample_id, proteins):
    annotations = []
    seen = set()
    for number, row in read_tsv(path, {'orf_id', 'arg_subtype'}, {'arg_type'}, sample_id):
        identifier = row['orf_id']
        if identifier in seen:
            fail(path, number, f'duplicate final annotation for ORF {identifier!r}')
        if identifier not in proteins:
            fail(path, number, f'ORF {identifier!r} has no protein/GFF mapping')
        seen.add(identifier)
        annotations.append(row)
    return annotations

def read_context(path, sample_id):
    labels = {}
    if path is None:
        return labels
    for number, row in read_tsv(path, {'contig_id', 'compartment'}, set(), sample_id):
        contig = row['contig_id']
        if contig in labels:
            fail(path, number, f'duplicate context for contig {contig!r}')
        if row['compartment'] not in COMPARTMENTS:
            fail(path, number, 'compartment must be chromosome/plasmid/virus/unclassified')
        labels[contig] = row['compartment']
    return labels

def read_kraken_report(path):
    records = {}
    with path.open(encoding='utf-8-sig') as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            cols = raw.rstrip('\r\n').split('\t')
            if len(cols) != 6:
                fail(path, number, 'expected standard 6-column Kraken report')
            taxid = cols[4].strip()
            if not re.fullmatch('[0-9]+', taxid):
                fail(path, number, 'report taxid must be numeric')
            taxid = str(int(taxid))
            rank, name = (cols[3].strip(), cols[5].strip())
            if not rank or not name:
                fail(path, number, 'report rank and name must be nonempty')
            if taxid in records:
                fail(path, number, f'duplicate report taxid {taxid!r}')
            records[taxid] = (rank, name)
    return records

def read_kraken_output(path):
    records = {}
    with path.open(encoding='utf-8-sig') as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            cols = raw.rstrip('\r\n').split('\t')
            if len(cols) != 5:
                fail(path, number, 'expected default 5-column Kraken classification output')
            status, contig, taxid = cols[:3]
            if status not in {'C', 'U'} or not contig:
                fail(path, number, 'invalid Kraken status or empty contig ID')
            if not re.fullmatch('[0-9]+', taxid):
                fail(path, number, 'Kraken taxid must be numeric; use default output format')
            taxid = str(int(taxid))
            if status == 'U' and taxid != '0' or (status == 'C' and taxid == '0'):
                fail(path, number, 'Kraken status and taxid are inconsistent')
            if contig in records:
                fail(path, number, f'duplicate Kraken contig ID {contig!r}')
            records[contig] = (status, taxid)
    return records

def existing_taxonomy(contig, kraken, report):
    empty = {'contig_taxid': '', 'contig_taxon_rank': '', 'contig_taxon_name': '', 'host_species': '', 'host_status': 'no_kraken_annotation'}
    if contig not in kraken:
        return empty
    status, taxid = kraken[contig]
    empty['contig_taxid'] = taxid
    if status == 'U':
        empty['host_status'] = 'unclassified'
        return empty
    if taxid not in report:
        empty['host_status'] = 'taxid_missing_from_report'
        return empty
    rank, name = report[taxid]
    empty.update(contig_taxon_rank=rank, contig_taxon_name=name)
    if rank == 'S':
        empty.update(host_species=name, host_status='species_exact_rank')
    else:
        empty['host_status'] = 'classified_non_species_rank'
    return empty

def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def write_tsv(path, fields, rows):
    with path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def run(args):
    if not args.sample_id or args.sample_id != args.sample_id.strip() or any((character in args.sample_id for character in '\t\r\n')):
        raise InputError('sample-id must be a nonempty single-line identifier')
    if bool(args.kraken_output) != bool(args.kraken_report):
        raise InputError('--kraken-output and --kraken-report must be supplied together')
    inputs = {name: getattr(args, name) for name in ('gff', 'proteins', 'arg_table', 'context_table', 'kraken_output', 'kraken_report') if getattr(args, name) is not None}
    resolved_inputs = [path.resolve() for path in inputs.values()]
    if len(resolved_inputs) != len(set(resolved_inputs)):
        raise InputError('each input role must have a distinct file')
    audit_path = args.audit or args.output.with_suffix('.audit.json')
    outputs = [args.output, audit_path] + ([args.contig_counts] if args.contig_counts else [])
    resolved_outputs = [path.resolve() for path in outputs]
    if len(resolved_outputs) != len(set(resolved_outputs)):
        raise InputError('output paths must be distinct')
    if set(resolved_outputs) & set(resolved_inputs):
        raise InputError('an output path would overwrite an input')
    for path in outputs:
        if path.exists():
            raise InputError(f'output already exists: {path}')
        if not path.parent.is_dir():
            raise InputError(f'output parent directory does not exist: {path.parent}')
    genes, gff_contigs = read_gff(args.gff)
    proteins = read_protein_headers(args.proteins, genes)
    annotations = read_args(args.arg_table, args.sample_id, proteins)
    labels = read_context(args.context_table, args.sample_id)
    kraken = read_kraken_output(args.kraken_output) if args.kraken_output else {}
    report = read_kraken_report(args.kraken_report) if args.kraken_report else {}
    rows = []
    for annotation in annotations:
        orf = annotation['orf_id']
        contig = proteins[orf]
        rows.append({'sample_id': args.sample_id, 'orf_id': orf, 'arg_subtype': annotation['arg_subtype'], 'arg_type': annotation.get('arg_type', ''), 'contig_id': contig, 'compartment': labels.get(contig, 'unclassified'), 'compartment_status': 'provided_label' if contig in labels else 'no_annotation', **existing_taxonomy(contig, kraken, report)})
    groups = defaultdict(list)
    for row in rows:
        groups[row['sample_id'], row['contig_id']].append(row)
    contig_rows = [{'sample_id': sample, 'contig_id': contig, 'arg_orf_count': len(group), 'distinct_arg_subtype_count': len({row['arg_subtype'] for row in group})} for (sample, contig), group in sorted(groups.items())]
    annotated_contigs = {row['contig_id'] for row in rows}
    audit = {'implementation_version': VERSION, 'implementation_status': 'new implementation; unvalidated against study data', 'sample_id': args.sample_id, 'association_method': 'explicit identifiers from existing annotation records only', 'inputs': {name: {'path': str(path), 'sha256': sha256(path)} for name, path in inputs.items()}, 'counts': {'gff_cds': len(genes), 'gff_contigs_with_cds': len(gff_contigs), 'protein_headers': len(proteins), 'final_arg_annotations': len(annotations), 'output_arg_orf_rows': len(rows), 'unique_sample_contigs': len(groups), 'context_rows': len(labels), 'kraken_rows': len(kraken), 'context_rows_not_joined_to_arg_contigs': len(set(labels) - annotated_contigs), 'kraken_rows_not_joined_to_arg_contigs': len(set(kraken) - annotated_contigs), 'context_contigs_not_in_gff_cds': len(set(labels) - gff_contigs), 'kraken_contigs_not_in_gff_cds': len(set(kraken) - gff_contigs), 'arg_orfs_by_compartment_status': dict(Counter((row['compartment_status'] for row in rows))), 'arg_orfs_by_host_status': dict(Counter((row['host_status'] for row in rows)))}, 'interpretation': ['One output row per final ARG ORF annotation; no best-hit selection occurs.', 'Compartment is a supplied existing label; missing labels remain unclassified.', "host_species copies the contig's final Kraken assignment only when report rank is exactly S.", 'Same-contig bookkeeping does not establish biological host, transfer, infectivity, or mobility.', 'GFF/FAA have no sample field; the caller must supply files from one sample and one run.']}
    write_tsv(args.output, OUTPUT_FIELDS, rows)
    if args.contig_counts:
        write_tsv(args.contig_counts, ['sample_id', 'contig_id', 'arg_orf_count', 'distinct_arg_subtype_count'], contig_rows)
    with audit_path.open('x', encoding='utf-8') as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False)
        handle.write('\n')
    return audit

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-id', required=True)
    for name in ('gff', 'proteins', 'arg-table'):
        parser.add_argument('--' + name, required=True, type=Path)
    for name in ('context-table', 'kraken-output', 'kraken-report', 'contig-counts', 'audit'):
        parser.add_argument('--' + name, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    try:
        audit = run(args)
    except (InputError, OSError, csv.Error) as error:
        parser.exit(2, f'ERROR: {error}\n')
    print(json.dumps({'output': str(args.output), 'counts': audit['counts']}, ensure_ascii=False))
if __name__ == '__main__':
    main()
