import argparse
import json
from pathlib import Path
import shlex
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('assembly', 'prodigal'))
    parser.add_argument('--config', type=Path, default=Path(__file__).resolve().parent / 'slurm.json')
    parser.add_argument('--array', required=True)
    parser.add_argument('--export', default='ALL')
    parser.add_argument('--print-command', action='store_true')
    args = parser.parse_args()
    config_path = args.config.resolve()
    job = json.loads(config_path.read_text(encoding='utf-8'))[args.stage]
    script = (config_path.parent / job['script']).resolve()
    if not script.is_file():
        parser.error(str(script))
    options = list(job['sbatch_arguments'])
    command = ['sbatch', *options, '--array=' + args.array, '--export=' + args.export, '--wrap', shlex.join(['bash', str(script)])]
    if args.print_command:
        print(json.dumps(command, indent=2))
        return
    for index, option in enumerate(options):
        if option in ('-o', '-e', '--output', '--error'):
            Path(options[index + 1]).parent.mkdir(parents=True, exist_ok=True)
        elif option.startswith(('--output=', '--error=')):
            Path(option.split('=', 1)[1]).parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, check=True)


if __name__ == '__main__':
    main()
