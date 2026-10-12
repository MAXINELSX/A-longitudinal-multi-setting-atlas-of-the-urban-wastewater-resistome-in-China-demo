import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
TOOL_KEYS = ('bowtie2_build', 'bowtie2', 'samtools', 'depth', 'metabat2', 'checkm2', 'drep', 'mash', 'fastani', 'gtdbtk')

def write_table(path, columns, rows, delimiter='\t'):
    with Path(path).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter, lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def table(path, required, delimiter='\t'):
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise ValueError(f'{path}: required columns {required}')
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f'{path}: duplicate column names')
        rows = list(reader)
        if any((None in row or any((v is None for v in row.values())) for row in rows)):
            raise ValueError(f'{path}: malformed row')
        return rows

def load_samples(path, check_files=True):
    path = Path(path).resolve()
    rows = table(path, ('sample_id', 'assembly', 'r1', 'r2'))
    if not rows:
        raise ValueError('Sample manifest is empty')
    seen = set()
    for row in rows:
        name = row['sample_id']
        if not name or any((c in name for c in '\t\r\n')) or name in seen:
            raise ValueError(f'Invalid or duplicate sample_id: {name!r}')
        seen.add(name)
        slug = re.sub('[^A-Za-z0-9]+', '_', name).strip('_')[:32] or 'sample'
        row['sample_key'] = 's_' + slug + '_' + hashlib.sha256(name.encode()).hexdigest()[:12]
        for field in ('assembly', 'r1', 'r2'):
            target = Path(row[field]).expanduser()
            target = (path.parent / target).resolve() if not target.is_absolute() else target.resolve()
            if check_files and (not target.is_file() or target.stat().st_size == 0):
                raise ValueError(f'Missing/empty {field}: {target}')
            row[field] = str(target)
        if row['r1'] == row['r2']:
            raise ValueError(f'Identical paired read files: {name}')
        if row['assembly'].lower().endswith('.gz'):
            raise ValueError('Assembly FASTA must be uncompressed; provide a separately decompressed file')
    rows.sort(key=lambda r: r['sample_id'])
    if len({r['sample_key'] for r in rows}) != len(rows):
        raise ValueError('Sample key collision')
    return rows

def load_config(path, executing=False):
    with Path(path).open(encoding='utf-8') as handle:
        cfg = json.load(handle)['mag']
    threads = cfg.get('threads')
    if type(threads) is not int or threads < 1:
        raise ValueError('mag.threads must be a positive integer')
    if cfg.get('metabat_min_contig') != 1500:
        raise ValueError('This Methods 4.10 workflow requires metabat_min_contig=1500')
    for key in TOOL_KEYS:
        spec = cfg['tools'][key]
        if not spec.get('executable') or not isinstance(spec.get('version_args'), list):
            raise ValueError(f'Tool {key} needs executable and version_args')
        if executing and (not spec.get('expected_version_text')):
            raise ValueError(f'Set the actual tested expected_version_text for {key}')
    seed = cfg.get('metabat_seed')
    if seed is not None and (type(seed) is not int or seed < 0):
        raise ValueError('Optional metabat_seed must be a nonnegative integer; null uses the tool default')
    scoring = cfg.get('drep_scoring')
    if scoring is not None:
        if not isinstance(scoring, dict):
            raise ValueError('Optional drep_scoring must be an object; null uses tool defaults')
        for key in ('completeness', 'contamination', 'n50', 'size', 'centrality'):
            value = scoring.get(key)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or (not math.isfinite(value)):
                raise ValueError(f'Set finite drep_scoring.{key}')
    if executing:
        if not Path(cfg['checkm2_database']).is_file():
            raise ValueError('checkm2_database must point to the compatible .dmnd file')
        if not Path(cfg['gtdbtk_database']).is_dir() or cfg.get('gtdb_release') != '226':
            raise ValueError('Set gtdbtk_database to extracted release 226; gtdb_release must be 226')
        metadata = Path(cfg['gtdbtk_database']) / 'metadata' / 'metadata.txt'
        lines = metadata.read_text(encoding='utf-8').splitlines()
        release = next((line.split('=', 1)[1].strip() for line in lines if line.startswith('VERSION_DATA=')), None)
        if release != 'r226':
            raise ValueError(f'GTDB metadata VERSION_DATA must be r226, observed {release!r}')
        cfg['checkm2_database'] = str(Path(cfg['checkm2_database']).resolve())
        cfg['gtdbtk_database'] = str(Path(cfg['gtdbtk_database']).resolve())
    return cfg

def executable(cfg, key):
    return cfg['tools'][key]['executable']

def sample_commands(cfg, sample, out):
    folder = Path(out) / 'samples' / sample['sample_key']
    index = folder / 'index' / 'assembly'
    sam, bam = (folder / 'alignment.sam', folder / 'alignment.sorted.bam')
    depth = folder / 'depth.tsv'
    threads = str(cfg['threads'])
    seed_args = [] if cfg.get('metabat_seed') is None else ['--seed', str(cfg['metabat_seed'])]
    return [('index', [executable(cfg, 'bowtie2_build'), sample['assembly'], str(index)]), ('map', [executable(cfg, 'bowtie2'), '-x', str(index), '-1', sample['r1'], '-2', sample['r2'], '-p', threads, '-S', str(sam)]), ('sort', [executable(cfg, 'samtools'), 'sort', '-@', threads, '-o', str(bam), str(sam)]), ('bam_check', [executable(cfg, 'samtools'), 'quickcheck', '-v', str(bam)]), ('depth', [executable(cfg, 'depth'), '--outputDepth', str(depth), str(bam)]), ('bin', [executable(cfg, 'metabat2'), '-i', sample['assembly'], '-a', str(depth), '-o', str(folder / 'bins' / 'bin'), '-m', '1500', '-t', threads] + seed_args)]

def downstream_commands(cfg, out):
    out = Path(out)
    threads = str(cfg['threads'])
    scoring = cfg.get('drep_scoring') or {}
    drepcmd = [executable(cfg, 'drep'), 'dereplicate', str(out / 'drep'), '-g', str(out / 'retained_genomes.txt'), '--genomeInfo', str(out / 'drep_genome_info.csv'), '-p', threads, '-l', '0', '-comp', '50', '-con', '10', '-pa', '0.90', '-sa', '0.95', '-nc', '0.30', '--S_algorithm', 'fastANI']
    if scoring:
        for key, flag in (('completeness', '-comW'), ('contamination', '-conW'), ('n50', '-N50W'), ('size', '-sizeW'), ('centrality', '-centW')):
            drepcmd += [flag, str(scoring[key])]
    return [('checkm2', [executable(cfg, 'checkm2'), 'predict', '--input', str(out / 'all_bins'), '--output-directory', str(out / 'checkm2'), '--threads', threads, '-x', 'fna', '--database_path', str(cfg.get('checkm2_database') or '<CHECKM2_DMND>')]), ('drep', drepcmd), ('gtdbtk', [executable(cfg, 'gtdbtk'), 'classify_wf', '--genome_dir', str(out / 'retained_bins'), '--out_dir', str(out / 'gtdbtk'), '--extension', 'fna', '--cpus', threads])]

def fasta_ids(path):
    ids = []
    with Path(path).open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f'Empty FASTA identifier: {path}')
                ids.append(fields[0])
    if not ids or len(ids) != len(set(ids)):
        raise ValueError(f'Missing or duplicate FASTA identifiers: {path}')
    return ids

def input_hashes(samples):
    cache, rows = ({}, [])
    for sample in samples:
        for role in ('assembly', 'r1', 'r2'):
            path = Path(sample[role])
            if str(path) not in cache:
                before = path.stat()
                digest = hashlib.sha256()
                with path.open('rb') as handle:
                    for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
                        digest.update(chunk)
                after = path.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise ValueError(f'Input changed while hashing: {path}')
                cache[str(path)] = (after.st_size, digest.hexdigest())
            size, digest = cache[str(path)]
            rows.append({'sample_id': sample['sample_id'], 'role': role, 'path': str(path), 'bytes': size, 'sha256': digest})
    return rows

def read_quality(path, bins):
    expected = {b['bin_id']: b for b in bins}
    qualities = {}
    for row in table(path, ('Name', 'Completeness', 'Contamination')):
        name = row['Name']
        if name not in expected and name.endswith('.fna'):
            name = name[:-4]
        if name not in expected or name in qualities:
            raise ValueError(f"Unknown/duplicate CheckM2 identifier: {row['Name']}")
        comp, con = (float(row['Completeness']), float(row['Contamination']))
        if not math.isfinite(comp) or not math.isfinite(con) or (not 0 <= comp <= 100) or (con < 0):
            raise ValueError(f'Invalid CheckM2 quality value: {name}')
        qualities[name] = {**expected[name], 'completeness': comp, 'contamination': con, 'retained': comp >= 50 and con <= 10, 'high_quality': comp >= 90 and con <= 5}
    if set(qualities) != set(expected):
        raise ValueError('CheckM2 report does not cover every bin; refusing dRep fallback to CheckM1')
    return [qualities[k] for k in sorted(qualities)]

def catalogue_rows(cdb_path, wdb_path, retained):
    expected = {r['bin_id'] + '.fna': r for r in retained}
    clusters = {}
    for row in table(cdb_path, ('genome', 'secondary_cluster'), ','):
        name = row['genome']
        if name not in expected or name in clusters:
            raise ValueError(f'Unknown/duplicate Cdb genome: {name}')
        clusters[name] = row['secondary_cluster']
    if set(clusters) != set(expected):
        raise ValueError('dRep clustered genome set differs from supplied quality-filtered set')
    winners = {}
    for row in table(wdb_path, ('genome', 'cluster'), ','):
        name, cluster = (row['genome'], row['cluster'])
        if cluster in winners or clusters.get(name) != cluster:
            raise ValueError('Inconsistent dRep winning genome table')
        winners[cluster] = name
    if set(winners) != set(clusters.values()):
        raise ValueError('dRep did not provide one representative per cluster')
    return [{'sample_id': expected[name]['sample_id'], 'bin_id': expected[name]['bin_id'], 'secondary_cluster': clusters[name], 'representative_id': winners[clusters[name]][:-4], 'is_representative': name == winners[clusters[name]]} for name in sorted(expected)]

class Runner:

    def __init__(self, out, cfg):
        self.out = Path(out)
        self.env = dict(os.environ)
        self.env['GTDBTK_DATA_PATH'] = str(Path(cfg['gtdbtk_database']).resolve())
        helper_dirs = []
        for key in ('mash', 'fastani'):
            path = shutil.which(executable(cfg, key))
            if path:
                helper_dirs.append(str(Path(path).resolve().parent))
        self.env['PATH'] = os.pathsep.join(helper_dirs + [self.env.get('PATH', '')])

    def run(self, stage, args):
        with (self.out / 'commands.jsonl').open('a', encoding='utf-8') as log:
            log.write(json.dumps({'stage': stage, 'argv': args}) + '\n')
        with (self.out / 'logs' / (stage + '.log')).open('w', encoding='utf-8') as log:
            result = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, env=self.env, cwd=self.out, shell=False)
        if result.returncode:
            raise RuntimeError(f'{stage} failed ({result.returncode}); see logs/{stage}.log')

def verify_versions(cfg, out, runner):
    records = {}
    for key in TOOL_KEYS:
        spec = cfg['tools'][key]
        path = shutil.which(spec['executable'])
        if not path:
            raise ValueError(f"Executable not found: {spec['executable']}")
        result = subprocess.run([path] + spec['version_args'], capture_output=True, text=True, env=runner.env, shell=False)
        text = result.stdout + result.stderr
        records[key] = {'path': str(Path(path).resolve()), 'argv': [path] + spec['version_args'], 'exit_code': result.returncode, 'output': text, 'expected_version_text': spec['expected_version_text']}
        (out / 'versions.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
        if spec['expected_version_text'] not in text:
            raise ValueError(f'Version probe for {key} lacks expected text; see versions.json')
        spec['executable'] = str(Path(path).resolve())
    for key, canonical in (('mash', 'mash'), ('fastani', 'fastANI')):
        resolved = shutil.which(canonical, path=runner.env['PATH'])
        if not resolved or Path(resolved).resolve() != Path(records[key]['path']).resolve():
            raise ValueError(f'dRep helper {canonical} on PATH differs from configured executable')

def execute(cfg, samples, out):
    out = Path(out).resolve()
    if ' ' in str(out):
        raise ValueError('GTDB-Tk 2.6.1 does not support genome paths containing spaces; use an output path without spaces')
    out.mkdir(parents=True, exist_ok=False)
    (out / 'logs').mkdir()
    (out / 'all_bins').mkdir()
    (out / 'retained_bins').mkdir()
    (out / 'configuration.json').write_text(json.dumps({'mag': cfg}, indent=2), encoding='utf-8')
    write_table(out / 'samples.tsv', ('sample_id', 'sample_key', 'assembly', 'r1', 'r2'), samples)
    (out / 'RUNNING').write_text('New reconstructed run, not historical source\n', encoding='utf-8')
    runner = Runner(out, cfg)
    try:
        verify_versions(cfg, out, runner)
        write_table(out / 'input_files.tsv', ('sample_id', 'role', 'path', 'bytes', 'sha256'), input_hashes(samples))
        bins, links = ([], [])
        for sample in samples:
            assembly_ids = set(fasta_ids(sample['assembly']))
            folder = out / 'samples' / sample['sample_key']
            (folder / 'index').mkdir(parents=True)
            (folder / 'bins').mkdir()
            for stage, args in sample_commands(cfg, sample, out):
                runner.run(sample['sample_key'] + '_' + stage, args)
            files = sorted((folder / 'bins').glob('bin.*.fa'))
            assigned = set()
            for index, path in enumerate(files, 1):
                bin_id = sample['sample_key'] + f'_b{index:06d}'
                ids = fasta_ids(path)
                if not set(ids).issubset(assembly_ids) or assigned.intersection(ids):
                    raise ValueError('Bin contains unknown contigs or a contig assigned to multiple bins')
                assigned.update(ids)
                dest = out / 'all_bins' / (bin_id + '.fna')
                shutil.copyfile(path, dest)
                bins.append({'sample_id': sample['sample_id'], 'bin_id': bin_id, 'original_bin': str(path), 'fasta': str(dest)})
                links.extend(({'sample_id': sample['sample_id'], 'bin_id': bin_id, 'contig_id': name} for name in sorted(ids)))
        if not bins:
            raise ValueError('No MAG bins produced; no quality/dereplication commands were run')
        write_table(out / 'all_bins.tsv', ('sample_id', 'bin_id', 'original_bin', 'fasta'), bins)
        write_table(out / 'contig_to_bin.tsv', ('sample_id', 'bin_id', 'contig_id'), links)
        stages = dict(downstream_commands(cfg, out))
        runner.run('checkm2', stages['checkm2'])
        quality = read_quality(out / 'checkm2' / 'quality_report.tsv', bins)
        write_table(out / 'bin_quality.tsv', ('sample_id', 'bin_id', 'completeness', 'contamination', 'retained', 'high_quality', 'fasta'), quality)
        retained = [r for r in quality if r['retained']]
        if not retained:
            raise ValueError('No bins passed completeness >=50 and contamination <=10')
        for row in retained:
            dest = out / 'retained_bins' / (row['bin_id'] + '.fna')
            shutil.copyfile(row['fasta'], dest)
            row['fasta'] = str(dest)
        write_table(out / 'retained_bins.tsv', ('sample_id', 'bin_id', 'completeness', 'contamination', 'high_quality', 'fasta'), retained)
        (out / 'retained_genomes.txt').write_text(''.join((r['fasta'] + '\n' for r in retained)), encoding='utf-8')
        external_quality = [{'genome': Path(r['fasta']).name, 'completeness': r['completeness'], 'contamination': r['contamination']} for r in retained]
        write_table(out / 'drep_genome_info.csv', ('genome', 'completeness', 'contamination'), external_quality, ',')
        if {r['genome'] for r in external_quality} != {Path(r['fasta']).name for r in retained}:
            raise ValueError('Incomplete dRep external quality table')
        runner.run('drep', stages['drep'])
        mapping = catalogue_rows(out / 'drep' / 'data_tables' / 'Cdb.csv', out / 'drep' / 'data_tables' / 'Wdb.csv', retained)
        for row in mapping:
            if not (out / 'drep' / 'dereplicated_genomes' / (row['representative_id'] + '.fna')).is_file():
                raise ValueError('dRep representative FASTA is missing')
        write_table(out / 'catalogue_mapping.tsv', ('sample_id', 'bin_id', 'secondary_cluster', 'representative_id', 'is_representative'), mapping)
        runner.run('gtdbtk', stages['gtdbtk'])
        tax = {}
        for domain in ('bac120', 'ar53'):
            path = out / 'gtdbtk' / f'gtdbtk.{domain}.summary.tsv'
            if path.is_file():
                for row in table(path, ('user_genome', 'classification')):
                    if row['user_genome'] in tax:
                        raise ValueError('Duplicate GTDB genome identifier')
                    tax[row['user_genome']] = row['classification']
        if set(tax) != {r['bin_id'] for r in retained}:
            raise ValueError('GTDB summary coverage differs from retained bins; inspect GTDB warnings/failed genomes')
        write_table(out / 'retained_taxonomy.tsv', ('sample_id', 'bin_id', 'classification'), [{'sample_id': r['sample_id'], 'bin_id': r['bin_id'], 'classification': tax[r['bin_id']]} for r in retained])
        (out / 'COMPLETE').write_text('All stages and output interface checks passed\n', encoding='utf-8')
        (out / 'RUNNING').unlink()
    except Exception as error:
        (out / 'FAILED.txt').write_text(str(error) + '\n', encoding='utf-8')
        raise

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--out', required=True, help='New output directory; must not already exist')
    parser.add_argument('--dry-run', action='store_true', help='Print command plan; no tools executed or outputs written')
    args = parser.parse_args(argv)
    try:
        cfg = load_config(args.config, executing=not args.dry_run)
        samples = load_samples(args.samples, check_files=not args.dry_run)
        out = Path(args.out).resolve()
        if args.dry_run:
            commands = [{'stage': sample['sample_key'] + '_' + stage, 'argv': cmd} for sample in samples for stage, cmd in sample_commands(cfg, sample, out)]
            commands += [{'stage': stage, 'argv': cmd} for stage, cmd in downstream_commands(cfg, out)]
            print(json.dumps({'status': 'plan_only_not_execution', 'GTDBTK_DATA_PATH': cfg.get('gtdbtk_database'), 'note': 'CheckM2 input bins, retained manifests and dRep mappings are generated only during execution; placeholders require author configuration.', 'commands': commands}, indent=2))
        else:
            execute(cfg, samples, out)
    except (KeyError, OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0
if __name__ == '__main__':
    sys.exit(main())
