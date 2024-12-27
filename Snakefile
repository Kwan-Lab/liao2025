"""
Process Single nuclei reads with cellranger
"""
# snakemake-templates -  flag to run as modular pipeline
import os
import glob
# Pandas is only used to read a column of a csv. This can probably be dropped
# for the python csv module instead
import pandas

# Input and Output directories
FASTQ_DIR = "/nfs/turbo/umms-kykwan/raw_data"
CELLRANGER_DIR = "cellranger"
KALLISTO_DIR = "kb"
FASTQC_DIR = "fastqc"
MULTIQC_DIR = "multiqc"
COMBINED_10X_METRICS = "10x_metrics.csv"

# Programs
SEQTK = "/nfs/turbo/umms-kykwan/projects/programs/seqtk/seqtk"

# References
KALLISTO_INDEX = "/nfs/turbo/umms-kykwan/projects/reference/kallisto-introns/lamanno-sf"
CELLRANGER_REFERENCE = "/nfs/turbo/umms-kykwan/projects/reference/cellranger-mm10/refdata-gex-mm10-2020-A"

# Intermediate directories
SCRATCH_FASTQ_DIR = "/scratch/kykwan_root/kykwan/shared_data/alex_kwan/fastq"
SCRATCH_KALLISTO_DIR = "/scratch/kykwan_root/kykwan/shared_data/alex_kwan/kb"
SCRATCH_CELLRANGER_DIR = "/scratch/kykwan_root/kykwan/shared_data/alex_kwan/cellranger"

# Snake variables
# Read pattern once it is transferred out of raw_data folder to scratch
# This specific pattern is required by cellranger
BCL2_READS_PATTERN = "{{sample}}_S1_L001_R{read}_001.fastq.gz"
BCL2_FASTQC_PATTERN = "{{sample}}_S1_L001_R{read}_001_fastqc.{ext}"
BCL2_MULTIQC_IN_PATTERN = "{sample}_S1_L001_R{{read}}_001_fastqc.zip"

# Variables for subsetting of 250k reads/cell datasets
# 325M reads is the maximum number of reads for samples that are 25k reads/cell
MAX_READS = 325_000_000

# Raw data directory
SAMPLE_FILE = "PsiloKet snRNA Seq (2022-11-19).csv"
sample_df = pandas.read_csv(SAMPLE_FILE)
def get_unique_samples():
    # List of ids for use in read patterns
    samples =  sample_df["Seq ID"].dropna().tolist()

    # Maps sample id to list of fastq paths
    samples_map = {s: glob.glob(os.path.join(FASTQ_DIR, f"{s[:-2]}/fastqs_{s[:-2]}/{s}*.fastq.gz"))
                  for s in samples}

    # Maps sample id to its 10x cellranger analysis (outputed by rule cellranger)
    cellranger_map = {s: os.path.join(SCRATCH_CELLRANGER_DIR, f"{s}/outs/metrics_summary.csv")
                      for s in samples}

    return samples, samples_map, cellranger_map

SAMPLES, SAMPLE_TO_PATH, SAMPLE_TO_10X = get_unique_samples()

# seed is 11 because that is the default of SEQTK and I am making that explicit
def subsample(in_fagz, out_fagz, reads_or_percent, seed = "11"):
    import subprocess
    cmd = f"{SEQTK} sample -2 -s {seed} {in_fagz} {reads_or_percent} | gzip -c > {out_fagz}"
    subprocess.run(cmd, shell=True)

def out(*dirs):
    return expand(os.path.join(*dirs), sample=SAMPLES)


rule all:
    input:
        out(KALLISTO_DIR, "{sample}"),
        out(CELLRANGER_DIR, "{sample}"),
        COMBINED_10X_METRICS,
        expand(os.path.join(MULTIQC_DIR, "R{read}/R{read}_multiqc_report.html"), read=["1", "2"]),


# Great Lakes computing has a limitation with analyzing data stored on turbo, so copy it first to scratch
rule prepare_scratch:
    input: lambda wc: sorted(SAMPLE_TO_PATH.get(wc.sample))
    output: expand(os.path.join(SCRATCH_FASTQ_DIR, "{{sample}}", BCL2_READS_PATTERN), read=["1", "2"])
    resources:
        cpus_per_task = 2,
        mem_mb = '14GB'
    run:
        is_25k = (sample_df[sample_df["Seq ID"] == wildcards.sample]["reads/cell requested"] == "25K").any()
        if is_25k:
            print("Copying without subsampling")
            import shutil
            shutil.copyfile(input[0], output[0])
            shutil.copyfile(input[1], output[1])
        # Subsample if too large
        else:
            print("Subsampling")
            subsample(input[0], output[0], MAX_READS)
            subsample(input[1], output[1], MAX_READS)

rule kbtools:
    input: expand(os.path.join(SCRATCH_FASTQ_DIR, "{{sample}}", BCL2_READS_PATTERN), read=["1", "2"])
    output: directory(os.path.join(SCRATCH_KALLISTO_DIR, "{sample}"))
    # || true is there because kb info returns error by default. (See https://github.com/pachterlab/kb_python/issues/179) for issue tracking
    resources:
        cpus_per_task = 10,
        mem_mb = f"70GB",
        runtime = '60h'
    shell:
        """
        kb info || true
        kb count \
        -t 8 -m 50G \
        -i {KALLISTO_INDEX}/index.idx \
        -g {KALLISTO_INDEX}/t2g.txt \
        -c1 {KALLISTO_INDEX}/cdna_t2c.txt \
        -c2 {KALLISTO_INDEX}/intron_t2c.txt \
        -x 10xv3 --workflow nucleus --filter bustools \
        --h5ad -o {output} {input}
        """

# TODO rewrite with a general form
rule retrieve_scratch_kallisto:
    input: os.path.join("SCRATCH_KALLISTO_DIR", "{sample}")
    output: directory(os.path.join(KALLISTO_DIR, "{sample}"))
    shell: "cp -r {input} {KALLISTO_DIR}"

CELL_RANGER_MEM_GB = "70"
# First just be stupid and use the data ported to scratch
rule cellranger:
    input: 
        expand(os.path.join(SCRATCH_FASTQ_DIR, "{{sample}}", BCL2_READS_PATTERN), read=["1", "2"])
    output: 
        cellranger = directory(os.path.join(SCRATCH_CELLRANGER_DIR, "{sample}")),
        # metrics_summary links explicitly to metrics_summary rule since its inputs aren't created yet
        # Snakemake also creates the os.path.join(SCRATCH_CELLRANGER_DIR, '{sample}', 'outs') directory
        # Cell Ranger doesn't like that (google "NOT A PIPESTANCE" for details)
        # So we need to delete {output.cellranger} in the shell script before running cellranger
        metrics_summary = os.path.join(SCRATCH_CELLRANGER_DIR, "{sample}/outs/metrics_summary.csv")
    resources:
        cpus_per_task = 10,
        mem_mb = f"{CELL_RANGER_MEM_GB}GB",
        runtime = '60h'
    shell: # see above comments for why we rm -rf
        """
        set +u
        ml Bioinformatics cellranger/7.1.0
        set -u
        echo "CellRanger version: $(cellranger --version)"

        return_wd="$PWD"
        mkdir -p {SCRATCH_CELLRANGER_DIR}
        cd {SCRATCH_CELLRANGER_DIR}
        rm -rf {output.cellranger}
        cellranger count \
           --id {wildcards.sample} \
           --transcriptome {CELLRANGER_REFERENCE} \
           --fastqs {SCRATCH_FASTQ_DIR}/{wildcards.sample} \
           --include-introns true \
           --localcores {resources.cpus_per_task} \
           --localmem {CELL_RANGER_MEM_GB}
        cd $return_wd
        """

rule retrieve_scratch_cellranger:
    input: os.path.join(SCRATCH_CELLRANGER_DIR, "{sample}")
    output: directory(os.path.join(CELLRANGER_DIR, "{sample}"))
    shell: "cp -r {input}/ {CELLRANGER_DIR}"

# Get sequencing saturation from cellranger report
rule metrics_summary:
    input: SAMPLE_TO_10X.values()
    output: COMBINED_10X_METRICS
    run:
        import pandas as pd
        import re
        # e.g. 6887-CL-3
        regex = re.compile('[0-9]{4}-[A-Z]{2}-[0-9]')
        results = []
        for metrics_path in input:
            key = regex.search(metrics_path).group(0)
            val = pd.read_csv(metrics_path).set_axis([key]).rename_axis("Sample")
            results.append(val)
        pd.concat(results).to_csv(output[0])


rule fastqc:
    input: expand(os.path.join(SCRATCH_FASTQ_DIR, "{{sample}}", BCL2_READS_PATTERN), read=["1", "2"])
    output:
        expand(os.path.join(FASTQC_DIR, BCL2_FASTQC_PATTERN),
               read=["1", "2"], ext=["html", "zip"])
    log: "logs/{sample}_fastqc.log"
    resources:
        cpus_per_task = 1,
        mem_mb = f"30GB",
        runtime = '20h'
    shell:
        """
        ml Bioinformatics fastqc
        fastqc --version
        fastqc --outdir {FASTQC_DIR} {input} &> {log}
        """

rule multiqc:
    input: out(FASTQC_DIR, BCL2_MULTIQC_IN_PATTERN)
    # {read} is the wildcard
    params:
        outdir = os.path.join(MULTIQC_DIR, "R{read}")
    output:
        rep = os.path.join(MULTIQC_DIR, "R{read}", "R{read}_multiqc_report.html"),
        filelist = temp(os.path.join(MULTIQC_DIR, "R{read}", "filelist.txt"))
    shell:
        """
        rm -rf {params.outdir}
        mkdir {params.outdir}
        echo {input} | tr " " "\n" > {output.filelist}
        multiqc --version
        multiqc --file-list {output.filelist} \
            --filename R{wildcards.read}_multiqc_report \
            --outdir {params.outdir} --force
        """
