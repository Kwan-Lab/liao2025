"""
Process Single nuclei reads with cellranger
"""
import os
import glob
import pandas

# Input and Output directories
CELLRANGER_DIR = "cellranger"
FASTQC_DIR = "fastqc"
MULTIQC_DIR = "multiqc"
COMBINED_10X_METRICS = "10x_metrics.csv"

# References
CELLRANGER_REFERENCE = "cellranger-mm10/refdata-gex-mm10-2020-A"

# Snake variables
# Read pattern once it is transferred out of raw_data folder to scratch
# This specific pattern is required by cellranger
BCL2_READS_PATTERN = "{{sample}}_S1_L001_R{read}_001.fastq.gz"
BCL2_FASTQC_PATTERN = "{{sample}}_S1_L001_R{read}_001_fastqc.{ext}"
BCL2_MULTIQC_IN_PATTERN = "{sample}_S1_L001_R{{read}}_001_fastqc.zip"

# Variables for subsetting of 250k reads/cell datasets
# 325M reads is the maximum number of reads for samples that are 25k reads/cell
MAX_READS = 325_000_000

def out(*dirs):
    return expand(os.path.join(*dirs), sample=SAMPLES)

rule all:
    input:
        out(CELLRANGER_DIR, "{sample}"),
        COMBINED_10X_METRICS,
        expand(os.path.join(MULTIQC_DIR, "R{read}/R{read}_multiqc_report.html"), read=["1", "2"]),

CELL_RANGER_MEM_GB = "70"
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
        module load Bioinformatics cellranger/7.1.0
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
        module load Bioinformatics fastqc
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
