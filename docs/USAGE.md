# Command-line reference

Run examples from the package root. Python commands use the standard library. Every script supports `--help`.

## Download the complete example

```bash
python download_demo.py
python download_demo.py --print-urls
```

| Option | Meaning |
| --- | --- |
| `--metadata PATH` | JSON metadata; default `demo/dataset.json` beside the script |
| `--output-dir PATH` | FASTQ destination; default `data/` beside the script |
| `--print-urls` | Print the two source URLs without downloading |

Metadata requires `run_accession` and `fastq_files`, with exactly two file entries having `mate` 1 and 2. Each entry supplies `url`, `bytes` and `md5`. Output filenames are `<run_accession>_1.fastq.gz` and `<run_accession>_2.fastq.gz`.

The downloader verifies the complete file size and MD5 before renaming its `.part` file. Matching completed files are reused. A mismatching completed file or an existing `.part` file causes an error. Inspect and relocate an incomplete file before retrying; downloads are not resumed automatically. If you change the destination, update `samples.tsv` accordingly.

## Read QC

```bash
python scripts/fastp_one_sample.py --config config.json --sample-id SRR14234107 --threads 4
```

| Option | Meaning |
| --- | --- |
| `--config PATH` | Required JSON configuration |
| `--sample-id ID` | Select a manifest sample |
| `--task-index INTEGER` | Select a zero-based manifest row; mutually exclusive with `--sample-id` |
| `--threads INTEGER` | Override the thread count |
| `--output-root PATH` | Override the output root |
| `--dry-run` | Validate settings and input paths, then print the command |

Without explicit sample selection, the script uses `SLURM_ARRAY_TASK_ID`, or row zero if that variable is absent. Thread precedence is `--threads`, then `SLURM_CPUS_PER_TASK`, then `fastp.threads`.

The manifest must have exactly these tab-separated columns in this order:

```text
sample_id	read1	read2
SRR14234107	data/SRR14234107_1.fastq.gz	data/SRR14234107_2.fastq.gz
```

Sample identifiers must match `[A-Za-z0-9][A-Za-z0-9_.-]*` and be unique. Input paths must be nonempty files and cannot be reused across rows. The manifest path, relative read paths and configured output root resolve relative to `config.json`, not the manifest directory. An explicit `--output-root` resolves relative to the current directory.

Configuration keys are `schema_version: 1`, `manifest`, `output_root`, and `fastp`. The `fastp` object contains `executable`, `expected_version`, `threads`, `detect_adapter_for_pe`, `cut_front`, `cut_tail`, `cut_window_size`, `cut_mean_quality`, and `length_required`. The supplied settings enable adapter detection and front/tail trimming with a four-base window, mean quality 20 and minimum retained length 50 nt.

Outputs are placed under `<output_root>/<sample_id>/`. Matching completed results are reused only after settings, version and checksums have been verified. Incompatible or partial results require a new output root. The Bash launcher accepts the same arguments:

```bash
bash scripts/fastp_one_sample.sh --config config.json --sample-id SRR14234107
```

## Assembly

Run MEGAHIT directly on the paired QC files as shown in the [README](../README.md). The fresh-run arguments are paired inputs, output directory, thread count and `--min-contig-len 500`; the output directory must be new. The cluster alternative is described below.

## Prodigal

```bash
python scripts/run_prodigal.py --config config.json --sample-id SRR14234107 --input results/megahit/SRR14234107/final.contigs.fa
```

Required options are `--config`, `--sample-id` and `--input`. Optional options are `--output-root` and `--dry-run`. The sample identifier follows the QC identifier rule. Input may be uncompressed FASTA or `.gz` FASTA and resolves relative to the current directory. The configured output root resolves relative to the configuration file; an explicit output-root override resolves relative to the current directory.

The `prodigal` configuration contains `executable`, `expected_version`, `mode`, `translation_table_argument`, `closed_ends`, `minimum_contig_length` and `output_root`. Supplied values use version 2.6.3, `meta` mode, a translation-table argument of 11, no closed-ends override, and a 500 bp minimum. Actual model translation tables reported in the GFF are recorded in `completed.json`.

The wrapper prepares `input.filtered.fna`, rejects duplicate contig identifiers, and checks agreement among predicted GFF, protein and nucleotide identifiers. It creates `<output_root>/<sample_id>/` and refuses an existing result directory. Logs and `completed.json` accompany the FAA, FNA and GFF outputs.

## MAG workflow

```bash
python scripts/mag_workflow.py --config config.json --samples mag_samples.tsv --out results/mag_SRR14234107
```

`--config`, `--samples` and `--out` are required. `--dry-run` prints a command plan without running tools or writing results. The MAG manifest requires these columns:

```text
sample_id	assembly	r1	r2
SRR14234107	results/megahit/SRR14234107/final.contigs.fa	results/fastp/SRR14234107/R1.clean2.fq.gz	results/fastp/SRR14234107/R2.clean2.fq.gz
```

Sample IDs must be unique. Relative assembly/read paths resolve against the manifest directory. Assembly FASTA must be uncompressed. Read mates must have different paths. The output directory must be new and its absolute path must not contain spaces.

Under `mag.tools`, each tool needs `executable`, `version_args` and `expected_version_text`. Before execution, replace null version text with text reported by the installed tool's version command. The wrapper checks that text and saves the full probe output. Tool keys are `bowtie2_build`, `bowtie2`, `samtools`, `depth`, `metabat2`, `checkm2`, `drep`, `mash`, `fastani` and `gtdbtk`. dRep must resolve the configured `mash` and `fastANI` helpers through their canonical executable names. External environment dependencies must also be available.

Set `mag.checkm2_database` to the compatible CheckM2 `.dmnd` file and `mag.gtdbtk_database` to the extracted r226 directory. Absolute paths avoid dependence on the current directory. The wrapper requires `metadata/metadata.txt` to report `VERSION_DATA=r226` and sets `GTDBTK_DATA_PATH` for subprocesses.

Other keys are `threads`, `metabat_min_contig` (1500), `metabat_seed`, `drep_scoring` and `gtdb_release` (`226`). Null seed/scoring values omit overrides. A supplied scoring object must contain numeric `completeness`, `contamination`, `n50`, `size` and `centrality` values.

The workflow performs read mapping, coverage calculation, binning, quality filtering, dereplication and taxonomy assignment. Retention uses completeness ≥50% and contamination ≤10%; the quality table separately flags completeness ≥90% and contamination ≤5%. Dereplication uses primary ANI 0.90, secondary FastANI 0.95 and minimum coverage 0.30. Pre-dereplication quality-passing bins are retained and classified, with a separate catalogue mapping.

Output tables include `all_bins.tsv`, `contig_to_bin.tsv`, `bin_quality.tsv`, `retained_bins.tsv`, `catalogue_mapping.tsv` and `retained_taxonomy.tsv`. Configuration, input hashes, versions and executed commands are recorded alongside stage logs. `COMPLETE` indicates successful completion. Zero bins or zero retained genomes produce an error without a completion marker.

## Existing alignment tables

```bash
python scripts/filter_existing_hits.py --config config.json --input sample.hits.tsv --mapping reference_labels.tsv --output-dir results/sample_hits
```

Required options are `--config`, `--input` and `--output-dir`. Optional options are `--mapping` and `--max-evalue`. The output directory must be new. Input is an uncompressed, headerless TSV with exactly 14 columns:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

The column names above describe order; actual records must be separated by tabs. The optional mapping TSV must have exactly `reference_id`, `arg_subtype`, `arg_type`, with one unique reference ID per row. Labels are copied from this table.

Thresholds come from `study_records_not_executed_by_this_wrapper.diamond`: `reported_final_identity_percent`, `query_coverage_percent` and `final_max_evalue`. Supplied values are 80%, 80% and 1e-5. `--max-evalue` overrides the E-value ceiling. Query coverage is `100 × (abs(qend − qstart) + 1) / qlen`. Selection uses the reported rows only, sorting by E-value ascending, bit score descending, alignment length descending, then original row order. Targets absent from an earlier search output cannot be recovered by this module.

Outputs are `selected_hits.tsv` and `audit.json`; providing the mapping also creates `arg_annotations.tsv`. An empty input yields header-only tables. This module does not execute an aligner.

## Existing context labels

```bash
python scripts/context_labels_from_headers.py --plasmid-fasta existing_plasmids.fasta --virus-fasta existing_viruses.fasta --output contexts.tsv
```

Available inputs are `--chromosome-fasta`, `--plasmid-fasta`, `--virus-fasta` and `--unclassified-fasta`; provide at least one. `--output` is required. Inputs may be FASTA or `.gz` FASTA. The supplied option determines the copied existing label. The first whitespace-delimited FASTA header token becomes `contig_id`, producing a TSV with `contig_id` and `compartment`. Duplicate IDs are rejected. The output file must not already exist and its parent directory must exist.

## Existing annotation associations

```bash
mkdir -p results/sample_association
python scripts/associate_annotations.py --sample-id SAMPLE --gff genes.gff --proteins proteins.faa --arg-table arg_annotations.tsv --context-table contexts.tsv --kraken-output contigs.kraken.out --kraken-report contigs.kraken.report --output results/sample_association/associations.tsv --contig-counts results/sample_association/contig_counts.tsv
```

Required options are `--sample-id`, `--gff`, `--proteins`, `--arg-table` and `--output`. Optional options are `--context-table`, `--kraken-output`, `--kraken-report`, `--contig-counts` and `--audit`. Kraken output and report must be supplied together. All files must belong to the same sample and assembly. Relative paths resolve against the current directory. Output parents must exist; existing output files are refused.

| Input | Required format |
| --- | --- |
| GFF | Prodigal CDS records with original contig IDs and explicit `ID` attributes |
| Proteins | Matching Prodigal FAA headers; the first token is the ORF ID and the header's `ID=` attribute matches the GFF |
| ARG table | TSV with `orf_id`, `arg_subtype`; optional `arg_type` and `sample_id`; exactly one annotation per ORF |
| Context table | TSV with `contig_id`, `compartment`; optional `sample_id`; labels `chromosome`, `plasmid`, `virus`, `unclassified` |
| Kraken output | Default five-column contig-classification output |
| Kraken report | Standard six-column report from the same run |

Identifiers are matched exactly. Each ARG ORF remains a separate row; missing context labels remain unclassified. `host_species` is populated only for the contig's exact final Kraken assignment at report rank `S`. The output preserves supplied `virus` labels. The optional count table reports distinct sample/contig records. A JSON execution record defaults to the output filename with `.audit.json` replacing its suffix; `--audit` overrides this path. These associations do not establish biological transfer.

## Slurm submission

`submit_slurm.py` submits the included cluster workers. `slurm.json` stores `script` and `sbatch_arguments` separately for `assembly` and `prodigal`. Script paths resolve against the Slurm configuration directory. The cluster workers retain the original environment setup and data layout, and do not read `config.json`; adapt their paths, environment activation and `slurm.json` to your cluster before submission.

```bash
python submit_slurm.py assembly --array '0-9%2' --export 'ALL,TAG=batch01,DATA_DIR=/path/to/raw,USB_ROOT=/path/to/run' --print-command
python submit_slurm.py prodigal --array '0-9%2' --export 'ALL,LISTFILE=/path/to/contigs.list' --print-command
```

The positional stage is `assembly` or `prodigal`. `--array` is required; `--export` defaults to `ALL`; `--config` defaults to the package's `slurm.json`. `--print-command` only displays the submission command. Remove it to submit.

The assembly worker uses one subdirectory per sample and expects one `*_1.clean.fq.gz` and one `*_2.clean.fq.gz` raw-input file in each. Its launcher builds `sample_list.<TAG>.txt`; array indices start at zero. Its supplied resource setting is 30 CPUs. The Prodigal worker uses one FASTA path per line in `LISTFILE`, with zero-based task indices and eight CPUs in the supplied configuration. Use uncompressed FASTA for that original worker. The portable download/QC examples use a different input layout and are run directly with the Python commands above.
