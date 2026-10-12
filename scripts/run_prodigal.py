import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()

def fasta_records(path):
    opener = gzip.open if path.suffix.lower() == '.gz' else open
    with opener(path, 'rt', encoding='ascii') as handle:
        header, sequence = (None, [])
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield (header, ''.join(sequence))
                header, sequence = (line[1:], [])
                if not header.strip():
                    raise ValueError('Empty FASTA header')
            elif header is None:
                raise ValueError('FASTA sequence precedes first header')
            else:
                sequence.append(line)
        if header is not None:
            yield (header, ''.join(sequence))

def prepare_contigs(source, target, minimum):
    seen, retained = (set(), {})
    dropped = 0
    with target.open('w', encoding='ascii', newline='\n') as handle:
        for header, sequence in fasta_records(source):
            identifier = header.split()[0]
            if identifier in seen:
                raise ValueError(f'Duplicate contig ID: {identifier}')
            seen.add(identifier)
            if not sequence:
                raise ValueError(f'Empty contig: {identifier}')
            if any((c.isspace() for c in sequence)):
                raise ValueError(f'Whitespace inside sequence: {identifier}')
            if len(sequence) < minimum:
                dropped += 1
                continue
            retained[identifier] = len(sequence)
            handle.write('>' + header + '\n')
            for pos in range(0, len(sequence), 80):
                handle.write(sequence[pos:pos + 80] + '\n')
    if not retained:
        raise ValueError('No contigs meet the configured minimum length')
    return (retained, dropped)

def command_for(executable, folder, settings):
    command = [executable, '-i', str(folder / 'input.filtered.fna'), '-p', settings['mode'], '-g', str(settings['translation_table_argument'])]
    if settings['closed_ends']:
        command.append('-c')
    return command + ['-a', str(folder / 'proteins.faa'), '-d', str(folder / 'genes.fna'), '-f', 'gff', '-o', str(folder / 'genes.gff')]

def fasta_attribute_ids(path):
    identifiers = set()
    first_tokens = set()
    for header, sequence in fasta_records(path):
        match = re.search('(?:^|[;\\s])ID=([^;\\s]+)', header)
        if not match:
            raise ValueError(f'FASTA header lacks Prodigal ID attribute: {header}')
        if match.group(1) in identifiers or header.split()[0] in first_tokens:
            raise ValueError(f'Duplicate predicted-gene identifier in {path.name}')
        if not sequence:
            raise ValueError(f'Empty predicted-gene sequence in {path.name}')
        identifiers.add(match.group(1))
        first_tokens.add(header.split()[0])
    return identifiers

def validate_outputs(folder, contigs):
    ids, translation_tables = (set(), set())
    with (folder / 'genes.gff').open(encoding='utf-8') as handle:
        for raw in handle:
            if raw.startswith('# Model Data:'):
                match = re.search('(?:^|;)transl_table=(\\d+)', raw)
                if match:
                    translation_tables.add(int(match.group(1)))
            if raw.startswith('#') or not raw.strip():
                continue
            fields = raw.rstrip('\r\n').split('\t')
            if len(fields) != 9 or fields[2] != 'CDS':
                raise ValueError('Unexpected non-CDS or malformed Prodigal GFF record')
            if fields[0] not in contigs:
                raise ValueError(f'GFF references unknown input contig: {fields[0]}')
            start, end = (int(fields[3]), int(fields[4]))
            if not 1 <= start <= end <= contigs[fields[0]] or fields[6] not in ['+', '-']:
                raise ValueError('Invalid GFF coordinates or strand')
            attrs = dict((part.split('=', 1) for part in fields[8].split(';') if '=' in part))
            identifier = attrs.get('ID')
            if not identifier or identifier in ids:
                raise ValueError('Missing or repeated GFF ID')
            ids.add(identifier)
    if fasta_attribute_ids(folder / 'proteins.faa') != ids:
        raise ValueError('Protein FASTA IDs disagree with GFF IDs')
    if fasta_attribute_ids(folder / 'genes.fna') != ids:
        raise ValueError('Nucleotide FASTA IDs disagree with GFF IDs')
    return (len(ids), sorted(translation_tables))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output-root', type=Path)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', args.sample_id):
        raise ValueError('Unsafe sample ID')
    config = json.loads(args.config.read_text(encoding='utf-8'))
    settings = config['prodigal']
    if settings['mode'] != 'meta' or type(settings['closed_ends']) is not bool:
        raise ValueError('Expected explicit meta mode and boolean closed_ends')
    if type(settings['minimum_contig_length']) is not int or settings['minimum_contig_length'] < 1:
        raise ValueError('minimum_contig_length must be positive')
    source = args.input.resolve()
    if not source.is_file() or not source.stat().st_size:
        raise ValueError(f'Missing/empty input: {source}')
    outroot = args.output_root.resolve() if args.output_root else (args.config.resolve().parent / settings['output_root']).resolve()
    final = outroot / args.sample_id
    if args.dry_run:
        print(json.dumps({'status': 'DRY_RUN', 'source': str(source), 'preparation': {'explicit_gzip_decompression': source.suffix.lower() == '.gz', 'minimum_contig_length': settings['minimum_contig_length']}, 'command_argv': command_for(settings['executable'], final, settings), 'note': 'New reconstruction; no Prodigal process executed'}, indent=2))
        return
    executable = shutil.which(settings['executable'])
    if not executable:
        raise ValueError('Prodigal executable not found')
    version = subprocess.run([executable, '-v'], capture_output=True, text=True)
    version_text = (version.stdout + version.stderr).strip()
    match = re.search('(\\d+\\.\\d+\\.\\d+)', version_text)
    if not match or match.group(1) != settings['expected_version']:
        raise ValueError(f'Prodigal version mismatch: {version_text!r}')
    outroot.mkdir(parents=True, exist_ok=True)
    lock = outroot / ('.' + args.sample_id + '.lock')
    try:
        lock.mkdir()
    except FileExistsError:
        raise ValueError(f'Existing lock: {lock}; inspect active jobs before removal')
    temporary = None
    try:
        if final.exists():
            raise ValueError(f'Refusing to overwrite existing output: {final}; choose a new output root')
        temporary = Path(tempfile.mkdtemp(prefix='.' + args.sample_id + '.incomplete.', dir=outroot))
        lengths, dropped = prepare_contigs(source, temporary / 'input.filtered.fna', settings['minimum_contig_length'])
        command = command_for(executable, temporary, settings)
        started = time.perf_counter()
        with (temporary / 'prodigal.log').open('w', encoding='utf-8') as log:
            log.write('COMMAND_ARGV=' + json.dumps(command) + '\n')
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        elapsed = time.perf_counter() - started
        gene_count, tables = validate_outputs(temporary, lengths)
        record = {'status': 'NEW_RECONSTRUCTED_WORKER', 'sample_id': args.sample_id, 'input_sha256': sha256(source), 'input_path': str(source), 'settings': settings, 'version_output': version_text, 'retained_contigs': len(lengths), 'dropped_short_contigs': dropped, 'gene_count': gene_count, 'translation_tables_reported_in_gff': tables, 'command_argv': command, 'elapsed_seconds': elapsed, 'output_sha256': {name: sha256(temporary / name) for name in ['input.filtered.fna', 'proteins.faa', 'genes.fna', 'genes.gff', 'prodigal.log']}}
        (temporary / 'completed.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
        temporary.rename(final)
        temporary = None
        print(f'COMPLETED {final}')
    finally:
        if temporary is not None:
            print(f'Incomplete output retained: {temporary}', file=sys.stderr)
        lock.rmdir()
if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
