# Wastewater metagenome processing code

This package provides read quality control, sample assembly, gene prediction, MAG reconstruction and processing of existing annotation tables. Configuration is stored in `config.json`; detailed command-line options and input formats are in [USAGE](docs/USAGE.md).

## System requirements

The Python scripts require Python 3.9 or later and use only its standard library. Scientific tools are intended to run on Linux; Slurm is needed only for the cluster submission scripts. No GPU is required. The scientific-tool Linux distribution and full-sample resource usage have not been measured. MAG classification requires substantial RAM and database storage; use an appropriately provisioned compute node.

| Component | Configured study version | Installation documentation |
| --- | --- | --- |
| fastp | 1.0.1 | [fastp](https://github.com/OpenGene/fastp) |
| MEGAHIT | 1.2.9 | [MEGAHIT](https://github.com/voutcn/megahit) |
| Prodigal | 2.6.3 | [Prodigal](https://github.com/hyattpd/Prodigal) |
| Bowtie2 and samtools | Set the installed versions in `config.json` | [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), [samtools](https://www.htslib.org/download/) |
| MetaBAT2 and `jgi_summarize_bam_contig_depths` | MetaBAT2 2.18; record the installed depth-tool version text | [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) |
| CheckM2 | 1.1.0, with its compatible DIAMOND database | [CheckM2](https://github.com/chklovski/CheckM2) |
| dRep | 3.6.2 | [dRep](https://drep.readthedocs.io/en/latest/installation.html) |
| Mash and FastANI | Set the installed versions in `config.json` | [Mash](https://github.com/marbl/Mash), [FastANI](https://github.com/ParBLiSS/FastANI) |
| GTDB-Tk | 2.6.1, reference release r226 | [GTDB-Tk 2.6.1](https://github.com/Ecogenomics/GTDBTk/tree/2.6.1) |

These are configured versions, not a claim that the external tools were tested in this package. The Python modules were tested with Python 3.12.14 on Windows 11 build 26200 using synthetic inputs and mocked scientific-tool calls. The complete public run and real scientific tools have not been processed here.

## Installation

1. Extract the package into a working directory. Run the commands below from that directory.
2. Install Python. No `pip install` step is required for these scripts.
3. Install the indicated scientific-tool releases using their official instructions. Use a QC/assembly environment for fastp, MEGAHIT and Prodigal, and separate compatible environments for CheckM2, dRep and GTDB-Tk. Install each tool's required dependencies in its own environment.
4. Place the active QC tools on `PATH`. For the MAG wrapper, set executable paths and observed version text under `mag.tools` in `config.json`; executable paths may point to separate environments. Make the dependencies required by each executable available through that environment or a launcher. Set `mag.checkm2_database` to the compatible `.dmnd` file and `mag.gtdbtk_database` to the extracted r226 directory. Use absolute database paths.

The supplied QC and MAG configurations use four threads. Adjust thread counts to the resources allocated to your job. Optional MAG seed and representative-scoring values of `null` retain the installed software defaults. Installation time, including reference database downloads, has not been measured on a normal desktop computer.

## Complete influent example

The example uses **all paired reads from [SRR14234107](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14234107)**, BioProject [PRJNA721748](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA721748), from the Tong Zhang group at The University of Hong Kong. Sample `ST_IN_0406_merged` is wastewater influent collected in Hong Kong on 6 April 2018. Its 33,648,454 PE150 read pairs contain 10,094,536,200 bases. The two compressed FASTQ files total 4,989,936,269 bytes, approximately 4.99 GB. Accession metadata, provenance links, download locations and MD5 values are in [dataset.json](demo/dataset.json).

### Download and quality control

```bash
python download_demo.py
python scripts/fastp_one_sample.py --config config.json --sample-id SRR14234107 --threads 4
```

The downloader retrieves the complete paired FASTQ files and checks their recorded byte counts and MD5 values. It writes `data/SRR14234107_1.fastq.gz` and `data/SRR14234107_2.fastq.gz`, matching `samples.tsv`. It does not extract a read subset. Existing files are reused only when their checksums and sizes match.

QC produces `results/fastp/SRR14234107/R1.clean2.fq.gz`, `R2.clean2.fq.gz`, `fastp.html`, `fastp.json`, `fastp.log` and `completed.json`. The completion record includes input/output checksums, settings and measured fastp elapsed time. An optional `--dry-run` prints the command after input validation without running fastp.

### Assembly and gene prediction

With MEGAHIT and Prodigal available, run a fresh assembly:

```bash
mkdir -p results/megahit
megahit -1 results/fastp/SRR14234107/R1.clean2.fq.gz -2 results/fastp/SRR14234107/R2.clean2.fq.gz -o results/megahit/SRR14234107 -t 4 --min-contig-len 500
python scripts/run_prodigal.py --config config.json --sample-id SRR14234107 --input results/megahit/SRR14234107/final.contigs.fa
```

Use a new MEGAHIT output directory. Expected outputs include `final.contigs.fa` and the Prodigal files `results/prodigal/SRR14234107/proteins.faa`, `genes.fna`, `genes.gff`, `input.filtered.fna`, `prodigal.log` and `completed.json`. Contig and gene counts depend on the actual execution; fixed expected counts for this complete run are not supplied.

### MAG reconstruction

After installing and configuring all MAG tools and databases, the supplied `mag_samples.tsv` links the assembly to the corresponding QC reads:

```bash
python scripts/mag_workflow.py --config config.json --samples mag_samples.tsv --out results/mag_SRR14234107 --dry-run
python scripts/mag_workflow.py --config config.json --samples mag_samples.tsv --out results/mag_SRR14234107
```

Use a new output directory without spaces. Expected outputs include `contig_to_bin.tsv`, `bin_quality.tsv`, `retained_bins.tsv`, `catalogue_mapping.tsv`, `retained_taxonomy.tsv`, retained genome FASTA files and stage logs. `COMPLETE` is written only after all stages and output checks succeed. The workflow retains quality-filtered genomes before dereplication as well as their mapping to catalogue representatives. MAG yield is data-dependent; a run with no bins or no quality-passing bins stops with a descriptive error.

### Runtime

Installation and complete-sample runtime have not been measured on a normal desktop. For an approximately 10 Gbp sample on an approximately 60-core node with sufficient RAM and installed databases, **8–24 hours** is a preliminary scheduling estimate for the broader processing workflow. This is not a benchmark of SRR14234107; it excludes downloads, installation and queue time. The four-thread examples above may take longer.

## Use with your data

Create a tab-separated QC manifest with `sample_id`, `read1`, `read2`, and set `manifest` in `config.json`. Use the same sample identifier for assembly and Prodigal outputs. For MAG processing, supply `sample_id`, `assembly`, `r1`, `r2` in a separate manifest. See [USAGE](docs/USAGE.md) for path rules, full options, existing annotation-table formats and Slurm submission.

The table-processing modules accept existing alignment results and existing contig labels. They do not run DIAMOND, Kraken2 or geNomad. ARGs-OAP v3.2.4 read profiling and normalization, those external annotation steps, and the paper's downstream statistical analyses are not implemented by this package. Consequently, these commands are not a complete reproduction of every quantitative result in the manuscript.

## License and release

The package license has not yet been assigned, and this package version has not been published as a repository release. External software and public data retain their respective licenses and terms.
