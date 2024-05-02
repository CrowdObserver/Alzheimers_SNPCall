# Snakemake SNP Calling pipeline that extracts and examines 3 genes (APP, SOD1, DYRK1A).

## Rules:
- download&emsp;&emsp;&emsp;&emsp; - Downloads fastq files from specified directory.
- fastqc&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;- Performing FastQC on fastq files.
- bwa_align&emsp;&emsp;&emsp;&emsp; - Aligning fastq files to provided genome.
- call_variants&emsp;&emsp;&emsp;&nbsp;- Calling variants with bcftools.
- snp_cleaning&emsp;&emsp;&nbsp;&nbsp;&nbsp;- Filtering and cleaning variants.
- snpeff&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp;- Annotate variants with snpEff.
- SNP_Extraction&emsp;&emsp;- Extracting SNPs belonging to APP, SOD1 or DYRK1A.
- bwa_check&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp;- Checking for primary reads %.
- heatmap&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;- Generating heatmap for every fastq file.

## Outputs
- genes.vcf&emsp;&emsp;&emsp;&emsp;- Extracted SNPs.
- heat.png&emsp;&emsp;&emsp;&emsp;- Heatmap.
- bwa_check.png&emsp;&nbsp;- Primary Reads % plot.

## TODO:
- [ ] Update Slurm job file
- [ ] Update Report
- [ ] Upgrade logging
- [ ] Make better tests

## DAG

![DAG](https://github.com/CrowdObserver/Assignment_1_Alzheimers/blob/master/dag.png?raw=true)
