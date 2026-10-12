import argparse
import csv
import gzip
from pathlib import Path
LABELS = ('plasmid', 'virus', 'chromosome', 'unclassified')

class InputError(ValueError):
    pass

def iter_header_ids(path):
    opener = gzip.open if path.suffix.lower() == '.gz' else Path.open
    seen_header = False
    with opener(path, mode='rt', encoding='utf-8-sig') as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip('\r\n')
            if not line.strip():
                continue
            if line.startswith('>'):
                tokens = line[1:].split()
                if not tokens:
                    raise InputError(f'{path}:{number}: empty FASTA identifier')
                seen_header = True
                yield tokens[0]
            elif not seen_header:
                raise InputError(f'{path}:{number}: nonblank content precedes first FASTA header')

def export_labels(inputs, output):
    if not inputs:
        raise InputError('at least one classified FASTA input is required')
    if set(inputs) - set(LABELS):
        raise InputError('unsupported existing class label')
    input_paths = [path.resolve() for path in inputs.values()]
    if len(input_paths) != len(set(input_paths)):
        raise InputError('each input label must have a distinct file')
    if output.resolve() in input_paths:
        raise InputError('output path would overwrite an input')
    if output.exists():
        raise InputError(f'output already exists: {output}')
    if not output.parent.is_dir():
        raise InputError(f'output parent directory does not exist: {output.parent}')
    records = {}
    for label, path in inputs.items():
        for identifier in iter_header_ids(path):
            if identifier in records:
                raise InputError(f'{path}: duplicate FASTA identifier {identifier!r} in supplied {records[identifier]!r} and {label!r} labels')
            records[identifier] = label
    with output.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['contig_id', 'compartment'])
        writer.writerows(records.items())
    return len(records)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for label in LABELS:
        parser.add_argument(f'--{label}-fasta', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    inputs = {label: getattr(args, label + '_fasta') for label in LABELS if getattr(args, label + '_fasta') is not None}
    try:
        count = export_labels(inputs, args.output)
    except (InputError, OSError, UnicodeError, EOFError) as error:
        parser.exit(2, f'ERROR: {error}\n')
    print(f'Exported {count} existing header labels to {args.output}')
if __name__ == '__main__':
    main()
