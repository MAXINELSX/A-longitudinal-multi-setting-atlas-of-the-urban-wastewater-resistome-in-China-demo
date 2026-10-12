import argparse
import hashlib
import json
import os
from pathlib import Path
import urllib.request


def checksum(path):
    digest = hashlib.md5()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def download(item, destination):
    expected_size = int(item['bytes'])
    expected_md5 = item['md5']
    if destination.exists():
        if destination.stat().st_size == expected_size and checksum(destination) == expected_md5:
            print('OK ' + str(destination))
            return
        raise FileExistsError('Existing file does not match the source checksum: ' + str(destination))
    temporary = destination.with_name(destination.name + '.part')
    if temporary.exists():
        raise FileExistsError('Remove or relocate the incomplete download before retrying: ' + str(temporary))
    digest = hashlib.md5()
    size = 0
    print('Downloading ' + item['url'], flush=True)
    request = urllib.request.Request(item['url'], headers={'User-Agent': 'nature-code-demo/1.0'})
    with urllib.request.urlopen(request, timeout=120) as response, temporary.open('xb') as handle:
        if response.status != 200:
            raise ValueError('Expected a complete HTTP response')
        for chunk in iter(lambda: response.read(8 * 1024 * 1024), b''):
            handle.write(chunk)
            digest.update(chunk)
            size += len(chunk)
            if size > expected_size:
                raise ValueError('Download exceeds the expected file size')
    if size != expected_size or digest.hexdigest() != expected_md5:
        raise ValueError('Source file size or MD5 mismatch: ' + str(temporary))
    if destination.exists():
        raise FileExistsError(str(destination))
    os.rename(temporary, destination)
    print('OK ' + str(destination), flush=True)


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata', type=Path, default=root / 'demo' / 'dataset.json')
    parser.add_argument('--output-dir', type=Path, default=root / 'data')
    parser.add_argument('--print-urls', action='store_true')
    args = parser.parse_args()
    metadata = json.loads(args.metadata.read_text(encoding='utf-8'))
    files = metadata['fastq_files']
    if sorted(item['mate'] for item in files) != [1, 2]:
        parser.error('The metadata must contain exactly one file per mate')
    run = metadata['run_accession']
    if not run.isalnum():
        parser.error('Invalid run accession')
    if args.print_urls:
        print('\n'.join(item['url'] for item in files))
        return
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for item in files:
        destination = args.output_dir / (run + '_' + str(item['mate']) + '.fastq.gz')
        download(item, destination)


if __name__ == '__main__':
    main()
