import argparse
import csv
import gzip
import hashlib
import json
import os
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
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def resolved(root, value):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()

def read_manifest(path, root):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['sample_id', 'read1', 'read2']:
            raise ValueError('Manifest header must be: sample_id<TAB>read1<TAB>read2')
        rows = list(reader)
    if not rows:
        raise ValueError('Manifest contains no samples')
    seen_ids, seen_paths = (set(), set())
    for row in rows:
        if None in row or any((not row.get(k) for k in ['sample_id', 'read1', 'read2'])):
            raise ValueError('Manifest row has extra columns or an empty field')
        sample = row['sample_id']
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
            raise ValueError(f'Unsafe sample_id: {sample!r}')
        if sample in seen_ids:
            raise ValueError(f'Duplicate sample_id: {sample}')
        seen_ids.add(sample)
        for key in ['read1', 'read2']:
            path = resolved(root, row[key])
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f'Missing or empty input: {path}')
            path_key = os.path.normcase(str(path))
            if path_key in seen_paths:
                raise ValueError(f'Input path reused within manifest: {path}')
            seen_paths.add(path_key)
            row[key] = path
    return rows

def validate_settings(settings, threads):
    for key in ['cut_window_size', 'cut_mean_quality', 'length_required']:
        value = settings[key]
        if type(value) is not int or value < 1:
            raise ValueError(f'{key} must be a positive integer')
    if settings['cut_mean_quality'] > 93:
        raise ValueError('cut_mean_quality is outside the FASTQ Phred+33 range')
    if threads < 1:
        raise ValueError('threads must be positive')
    for key in ['detect_adapter_for_pe', 'cut_front', 'cut_tail']:
        if type(settings[key]) is not bool:
            raise ValueError(f'{key} must be true or false')
    if not isinstance(settings.get('expected_version'), str) or not settings['expected_version']:
        raise ValueError('fastp.expected_version must be explicitly set')

def command_for(executable, row, folder, settings, threads):
    command = [executable, '-i', str(row['read1']), '-I', str(row['read2']), '-o', str(folder / 'R1.clean2.fq.gz'), '-O', str(folder / 'R2.clean2.fq.gz'), '-w', str(threads)]
    for key in ['detect_adapter_for_pe', 'cut_front', 'cut_tail']:
        if settings[key]:
            command.append('--' + key)
    for key in ['cut_window_size', 'cut_mean_quality', 'length_required']:
        command.extend(['--' + key, str(settings[key])])
    command.extend(['-h', str(folder / 'fastp.html'), '-j', str(folder / 'fastp.json')])
    return command

def validate_fastq(path, mate):
    with gzip.open(path, 'rt', encoding='ascii', newline='') as handle:
        count = 0
        while True:
            header = handle.readline().rstrip('\r\n')
            if not header:
                if handle.read(1):
                    raise ValueError(f'Blank header in {path}')
                return count
            sequence = handle.readline().rstrip('\r\n')
            plus = handle.readline().rstrip('\r\n')
            quality = handle.readline().rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+') or (not sequence):
                raise ValueError(f'Invalid FASTQ record in {path}')
            if len(sequence) != len(quality):
                raise ValueError(f'Sequence/quality lengths differ in {path}')
            count += 1

def fingerprint_inputs(row, settings, threads, version):
    run_settings = {key: settings[key] for key in ['detect_adapter_for_pe', 'cut_front', 'cut_tail', 'cut_window_size', 'cut_mean_quality', 'length_required']}
    return {'sample_id': row['sample_id'], 'input_sha256': {key: sha256(row[key]) for key in ['read1', 'read2']}, 'settings': run_settings, 'threads': threads, 'fastp_version': version}

def check_existing(folder, fingerprint):
    marker = folder / 'completed.json'
    if not marker.is_file():
        raise ValueError(f'Existing result has no completion record: {folder}; use a new output root')
    saved = json.loads(marker.read_text(encoding='utf-8'))
    if saved.get('fingerprint') != fingerprint:
        raise ValueError('Existing result uses different inputs, settings, or fastp version; use a new output root')
    expected_names = {'R1.clean2.fq.gz', 'R2.clean2.fq.gz', 'fastp.html', 'fastp.json', 'fastp.log'}
    recorded = saved.get('output_sha256', {})
    if set(recorded) != expected_names:
        raise ValueError('Completion record is missing required output checksums')
    for name, checksum in recorded.items():
        path = folder / name
        if not path.is_file() or sha256(path) != checksum:
            raise ValueError(f'Existing output checksum failed: {name}')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--sample-id')
    selection.add_argument('--task-index', type=int)
    parser.add_argument('--threads', type=int)
    parser.add_argument('--output-root', type=Path)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    config_path = args.config.resolve()
    root = config_path.parent
    config = json.loads(config_path.read_text(encoding='utf-8'))
    if config.get('schema_version') != 1:
        raise ValueError('Unsupported config schema_version')
    settings = config['fastp']
    threads = args.threads if args.threads is not None else int(os.environ.get('SLURM_CPUS_PER_TASK', settings['threads']))
    validate_settings(settings, threads)
    rows = read_manifest(resolved(root, config['manifest']), root)
    if args.sample_id is not None:
        matches = [row for row in rows if row['sample_id'] == args.sample_id]
        if not matches:
            raise ValueError(f'sample_id not found: {args.sample_id}')
        row = matches[0]
    else:
        index = args.task_index if args.task_index is not None else int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
        if not 0 <= index < len(rows):
            raise ValueError(f'task-index must be between 0 and {len(rows) - 1}')
        row = rows[index]
    outroot = args.output_root.resolve() if args.output_root else resolved(root, config['output_root'])
    final = outroot / row['sample_id']
    if args.dry_run:
        print(json.dumps({'status': 'DRY_RUN_NO_FASTP_EXECUTION', 'sample_id': row['sample_id'], 'command_argv': command_for(settings['executable'], row, final, settings, threads), 'output_directory': str(final), 'expected_version_not_checked': settings['expected_version']}, indent=2))
        return
    executable = shutil.which(settings['executable'])
    if not executable:
        raise ValueError(f"Executable not found: {settings['executable']}")
    version_process = subprocess.run([executable, '--version'], capture_output=True, text=True, check=True)
    version_text = (version_process.stdout + version_process.stderr).strip()
    match = re.search('(?<![0-9.])(\\d+\\.\\d+\\.\\d+)(?![0-9.])', version_text)
    if not match or match.group(1) != settings['expected_version']:
        raise ValueError(f"fastp version mismatch: expected {settings['expected_version']}, received {version_text!r}")
    fingerprint = fingerprint_inputs(row, settings, threads, match.group(1))
    outroot.mkdir(parents=True, exist_ok=True)
    lock = outroot / ('.' + row['sample_id'] + '.lock')
    try:
        lock.mkdir()
    except FileExistsError:
        raise ValueError(f'Lock exists: {lock}; check for an active job before manually removing it')
    temporary = None
    try:
        if final.exists():
            check_existing(final, fingerprint)
            print(f'REUSE_VERIFIED {final}')
            return
        temporary = Path(tempfile.mkdtemp(prefix='.' + row['sample_id'] + '.incomplete.', dir=outroot))
        command = command_for(executable, row, temporary, settings, threads)
        started = time.perf_counter()
        with (temporary / 'fastp.log').open('w', encoding='utf-8') as log:
            log.write('COMMAND_ARGV=' + json.dumps(command) + '\n')
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        elapsed = time.perf_counter() - started
        counts = [validate_fastq(temporary / name, mate) for name, mate in [('R1.clean2.fq.gz', 1), ('R2.clean2.fq.gz', 2)]]
        if counts[0] != counts[1]:
            raise ValueError('Output mate counts differ')
        report = json.loads((temporary / 'fastp.json').read_text(encoding='utf-8'))
        if not isinstance(report.get('summary'), dict):
            raise ValueError('fastp JSON lacks a summary object')
        if report['summary']['after_filtering']['total_reads'] != counts[0] + counts[1]:
            raise ValueError('FASTQ read count disagrees with fastp JSON')
        if not (temporary / 'fastp.html').stat().st_size:
            raise ValueError('fastp HTML report is empty')
        output_files = ['R1.clean2.fq.gz', 'R2.clean2.fq.gz', 'fastp.html', 'fastp.json', 'fastp.log']
        completion = {'wrapper_status': 'NEW_REVIEW_IMPLEMENTATION', 'fingerprint': fingerprint, 'version_output': version_text, 'elapsed_seconds': elapsed, 'output_pairs': counts[0], 'command_argv': command, 'output_sha256': {name: sha256(temporary / name) for name in output_files}}
        (temporary / 'completed.json').write_text(json.dumps(completion, indent=2) + '\n', encoding='utf-8')
        temporary.rename(final)
        temporary = None
        print(f'COMPLETED {final} elapsed_seconds={elapsed:.3f}')
    finally:
        if temporary is not None:
            print(f'Incomplete outputs retained for inspection: {temporary}', file=sys.stderr)
        lock.rmdir()
if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError, subprocess.SubprocessError, json.JSONDecodeError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
