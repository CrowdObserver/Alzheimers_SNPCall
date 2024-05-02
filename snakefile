# Variables and wildcards.
fastq_source = "</path/to/fastq>"
genome = "</path/to/genome>"
SNPEFF_JAR="</path/to/snpEff.jar>"
down, load, = glob_wildcards("</path/to/{down}.{load}.fq.gz")
sample = [d + "." + l for d, l in zip(down,load)]

# All rule
rule all:
    input:
        heat="heat.png",
        plot="bwa_check.png",

# Download every fq.gz file and unzip them.
rule download:
    input:
        fqs="/path/to/fastq/{sample}.fq.gz",
    output:
        fq="000.fastq/{sample}.fq",
    log:
        "logs/fq/{sample}_fq.log",
    shell:
        """
        ####    DOWNLOADING AND UNPACKING FASTQ FILES TO 000.fastq/    ####

        mkdir -p 000.fastq/
        cp -u {input.fqs} 000.fastq/ 
        gzip -dfk {output.fq}.gz

        ####    TESTING ALL FASTQ ARE GOOD    ####

        echo "ERROR: {output.fq} is broken, line count not divisible by 4." > {log}
        wc -l {output.fq} | awk '{{print ($1%4==0) ? true : false}}'        # This crashes the pipeline if one of the fastq is broken.
        echo "{output.fq} is in good health." > {log}

        """


# Downloading fq files, FastQC and Testing
rule fastqc:
    input:
        fq="000.fastq/{sample}.fq",
    output:
        fastqc_zip="010.fastqc/{sample}_fastqc.zip",
        html=report("010.fastqc/{sample}_fastqc.html", 
                    category="FastQC"),
        summarytxt="010.fastqc/{sample}_fastqc/summary.txt",
        summarydata="010.fastqc/{sample}_fastqc/fastqc_data.txt",
    log:
        "logs/FastQC/{sample}_fastqc.log",
    shell:
        """
        ####    RUNNING FASTQC    ####

        echo "Input Fastq: {input.fq}  \n" >> {log}
        fastqc -o 010.fastqc {input.fq} --extract

        ####    TESTING    ####

        # Check for fails in the summary

        if grep -q fail {output.summarydata}; then
            failed_tests=$(grep fail {output.summarydata})

            echo "FastQC FAILED the following tests in the fastqc quality check:" >> {log}
            echo "$failed_tests" >> {log}
            echo "" >> {log}

            if grep -q warn {output.summarydata}; then
                warnings=$(grep warn {output.summarydata})
                echo "FastQC WARNING:" >> {log}
                echo "$warnings" >> {log}
                echo "" >> {log}
            fi
            
            echo "Finished FastQC." >> {log}
        else    # Check for warnings in the pipeline
            if grep -q warn {output.summarydata}; then
                warnings=$(grep warn {output.summarydata})

                echo "FastQC PASSED with warnings:" >> {log}
                echo "$warnings" >> {log}
                echo "\n" >> {log}
                echo "Finished FastQC." >> {log}
            else
                echo "FastQC successful with no warnings." >> {log}
            fi
        fi
        """

# BWA alignment and Testing
rule bwa_align:
    input:
        fq="000.fastq/{sample}.fq",
        rep="010.fastqc/{sample}_fastqc.html",
    output:
        bam="020.alignment/{sample}.bam",
    log:
        report("logs/bwa/{sample}_bwa.log", category="Logs", subcategory="BWA"),
    shell:
        """
        ####    RUNNING BWA    ####

        echo "Running bwa on {input.fq}... \n" >> {log}
        bwa mem {genome} {input.fq} | samtools sort -o {output.bam} -
        samtools index {output.bam}

        samtools idxstats {output.bam} >> {log}
        
        ####    TESTING    ####
        # Checking % of primary mapped reads.

        if samtools flagstat {output.bam} | awk '/primary mapped/ {{gsub(/[(%)]/, "", $6); if ($6 > 95) exit true; else exit false}}'; then
            true
        else
            echo "WARNING: Primary mapped reads are less than 95% \n" >> {log}
        fi

        # Creating bwa_check.log for later plotting.
        echo -n '{output.bam};' >> 020.alignment/bwa_check.log && samtools flagstat {output.bam} | grep -u 'primary mapped' >> 020.alignment/bwa_check.log

        echo "BWA alignment completed." >> {log}
        """


# Calling variants and Testing
rule call_variants:
    input:
        bams=expand("020.alignment/{sample}.bam", sample=sample)
    output:
        raw_snps=report("030.vcf/Alz_raw_snps.vcf", category="VCFs", subcategory="Raw SNPs"),
        stats="030.vcf/vcf_stats/Alz_stats.vcfstats",
        plots=directory("030.vcf/vcf_stats/plots_Alz"),
    log:
        "logs/raw_variants.log",
    shell:
        """
        ####    RUNNING BCFTOOLS AND CALLING VARIANTS    ####

        echo "Calling variants for files in 020.alignment/ ... \n" >> {log}
        bcftools mpileup -Ou -f {genome} {input.bams} \
        | bcftools call -mv -Ov -o {output.raw_snps} 
        
        bcftools stats {output.raw_snps} > {output.stats}

        # Making plots with vcfstats
        plot-vcfstats -P -p {output.plots} {output.stats}

        ####    TESTING    ####

        # Testing for TS/TV ratio
        if tail -1 {output.plots}/tstv_by_af.0.dat | awk '{{if ($3 > 2.5 || $3 < 1.5) exit false; else exit true}}'; then
            echo -n "TS/TV ratio is: " >> {log} && tail -1 {output.plots}/tstv_by_af.0.dat | awk '{{print $3}}' >> {log} && echo "" >> {log}
        else
            echo "WARNING: TS/TV ratio is above 2.5 or below 1.5, this might mean that something is wrong. \n" >> {log}
        fi
        echo "VCF file generated." >> {log}
        """


# Cleaning SNPs and Primitive Testing
rule snp_cleaning:
    input:
        vcf="030.vcf/Alz_raw_snps.vcf",
        raw_vcf="030.vcf/vcf_stats/Alz_stats.vcfstats"
    output:
        clean_vcf=report("040.vcf_clean/Alz_clean_snps.vcf",
                        category="VCFs", subcategory="Cleaned, Filtered SNPs"),
        debug="040.vcf_clean/Alz_debug.txt",
    log:
        "logs/clean_variants.log",
    shell:
        """
        ####    RUNNING VT TO CLEAN AND FILTER    ####

        echo "Cleaning and filtering {input.vcf}... \n" >> {log}
        cat {input.vcf} \
            | vt view -f "QUAL>20" -h - \
            | vt decompose - 2> {output.debug}\
            | vt normalize -n -r {genome} - \
            | vt uniq - \
            > {output.clean_vcf}

        
        ####    TESTING    ####

        # Checking that the number of variants has changed after filtering.
        previous=$(awk '/number of records:/ {{print $6; exit}}' {input.raw_vcf})

        if awk -v prev="$previous" '/no. variants/ {{if (int($5) < int(prev)) exit 0; else exit 1}}' {output.debug}; then
            true
        else
            echo "-------------------------------------------" >> {log}
            echo "WARNING: Number of Variants has not changed." >> {log}
            echo "------------------------------------------- \n" >> {log}
        fi
        echo "Cleaning and filtering finished." >> {log}
        """


# Annotating SNPs 
rule snpeff:
    input:
        clean_vcf="040.vcf_clean/Alz_clean_snps.vcf",
        debug="040.vcf_clean/Alz_debug.txt",
    output:
        annotated_vcf=report("050.annotated_vcf/Alz_annotated_snps.vcf",
                            category="VCFs", subcategory="Annotated SNPs"),
        genes="050.annotated_vcf/Alz_snpEff_genes.txt",
        summary=report("050.annotated_vcf/Alz_snpEff_summary.html",
                        category="VCFs", subcategory="Annotated SNPs"),
    log:
        "logs/annotated_variants.log"
    shell:
        """
        ####    RUNNING SNPEFF    ####

        echo "Running snpEff on {input.clean_vcf}... \n" >> {log}

        java -Xmx3400m -jar {SNPEFF_JAR} eff hg38 \
            -dataDir /path/to/snpeff_db \
            {input.clean_vcf} \
            > {output.annotated_vcf}
        
        mv snpEff_genes.txt {output.genes}
        mv snpEff_summary.html {output.summary}


        ####    TESTING    ####

        # Checking that number of variants has not changed after annotating.
        prev=$(awk '/total no. of biallelics/ {{print $6; exit}}' {input.debug})
        if awk '/^#/{{next}}; {{count++}} END {{print count}}' {output.annotated_vcf} \
        | awk -v prev="$prev" '{{if (int($1) == int(prev)) exit 0; else exit 1}}'; then
            true
        else
            echo "-------------------------------------------" >> {log}
            echo "ERROR: Number of Variants has changed!" >> {log}
            echo "------------------------------------------- \n" >> {log}
            false
        fi
        echo "Variants annotated successfully." >> {log}
        """

# Extracting SNPs for APP, SOD1 and DYRK1A from the snpEff vcf and Testing.
rule SNP_Extraction:
    input:
        vcf="050.annotated_vcf/Alz_annotated_snps.vcf",
    output:
        genes=report("genes.vcf", category="VCFs", subcategory="Genes VCF")
    log:
        "logs/extracted_genes.log"
    shell:
        """
        ####    EXTRACTING APP, SOD1 AND DYRK1A FROM ANNOTATED VCF    ####

        echo "Beginning extraction of APP, SOD1 and DYRK1A...\n" >> {log}

        
        ### File quality check.

        if [ $(grep -c '#' {input.vcf}) -gt 0 ]; then
            grep -u '#' {input.vcf} >> {output.genes}
        else
            echo "ERROR: VCF might be corrupted, no headers found. \n" >> {log}
        fi

        
        ### Adding SNPs to the output files.

        if [ $(grep 'ANN=[^|]*|[^|]*|[^|]*|APP|' {input.vcf} | grep -v 'INDEL' | wc -l) -gt 0 ]; then
            grep -u 'ANN=[^|]*|[^|]*|[^|]*|APP|' {input.vcf} | grep -v 'INDEL' >> {output.genes}
        else
            echo "ERROR: APP was not found in the vcf. \n" >> {log}
        fi
        if [ $(grep 'ANN=[^|]*|[^|]*|[^|]*|SOD1|' {input.vcf} | grep -v 'INDEL' | wc -l) -gt 0 ]; then
            grep -u 'ANN=[^|]*|[^|]*|[^|]*|SOD1|' {input.vcf} | grep -v 'INDEL' >> {output.genes}
        else
            echo "ERROR: SOD1 was not found in the vcf. \n" >> {log}
        fi
        if [ $(grep 'ANN=[^|]*|[^|]*|[^|]*|DYRK1A|' {input.vcf} | grep -v 'INDEL' | wc -l) -gt 0 ]; then
            grep -u 'ANN=[^|]*|[^|]*|[^|]*|DYRK1A|' {input.vcf} | grep -v 'INDEL' >> {output.genes} 
        else
            echo "ERROR: DYRK1A was not found in the vcf. \n" >> {log}
        fi

        bcftools sort {output.genes} -o {output.genes}

        echo "Extraction completed in genes.vcf" >> {log}
        """

# Plotting Primary Mapped Reads in bwa_align.
rule bwa_check:
    input:
        bams=expand("020.alignment/{sample}.bam", sample=sample),
    output:
        plot=report("bwa_check.png", category="BWA Check"),
    run:
        # Importing libraries.
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        # Making a dictionary of samples and relative PrimaryMapped%.
        check = dict()
        for elem in sample:
            check[elem] = {"Primary_Mapped" : 0}
        
        # Read values from bwa_check.log .
        with open("020.alignment/bwa_check.log","r") as f:
            for line in f.readlines():
                name, value = line.split(";")
                name = name.split("/")[1][:-4]
                value = value.split(" ")[5][1:-1]
                check[name]["Primary_Mapped"] = float(value)                

        df = pd.DataFrame.from_dict(check).T
        df.reset_index(inplace=True)
        
        # Make an histogram.
        fig, ax = plt.subplots(figsize=(9, 6))
        sns.barplot(data=df, y="index", x="Primary_Mapped", color="lightgreen", ax=ax)
        ax.bar_label(ax.containers[0], fmt=lambda x: f'{x:.2f}%', position=(-45,0))
        plt.xlim(0,100)
        plt.xlabel("Primary Mapped %")
        plt.ylabel("")
        plt.title("Primary Mapped Values Check")
        plt.tight_layout()
        plt.savefig("bwa_check.png", dpi=300)



# Creating Heatmap of the three genes
rule heatmap:
    input:
        genes="genes.vcf",
    output:
        heat=report("heat.png", category="Heatmap"),
    log:
        "logs/heatmap.log",
    run:
        # Importing Libraries
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import vcfpy
        import sys

        # Redirecting output and error on log file
        sys.stdout = open("logs/heatmap.log","w")
        sys.stderr = open("logs/heatmap.log","w")

        # Reading VCF file from rule "SNP_Extraction"
        print("Creating Heatmap... \n")
        vcf = vcfpy.Reader.from_path(input.genes)

        # Creating a dictionary to start counting SNPs for every bam file
        table = dict()
        for elem in vcf.header.samples.names:
            table[elem.split("/")[1][:-4]] = {'APP' : 0, 'SOD1' : 0, 'DYRK1A' : 0}
        
        # Counting SNPs for every record and adding them to the dict
        print("- Counting SNPs...")
        for i, record in enumerate(vcf):
            for call in record.calls:
                sample = call.sample.split("/")[1][:-4]
                if "|APP|" in record.INFO.get('ANN')[0] and call.data.get('GT') not in ('./.', '0/0'):
                    table[sample]['APP'] +=1
                if "|SOD1|" in record.INFO.get('ANN')[0] and call.data.get('GT') not in ('./.', '0/0'):
                    table[sample]['SOD1'] +=1
                if "|DYRK1A|" in record.INFO.get('ANN')[0] and call.data.get('GT') not in ('./.', '0/0'):
                    table[sample]['DYRK1A'] +=1
        
        # Dictionary to dataframe
        df = pd.DataFrame.from_dict(table, orient="index")
        df.columns = ['APP', 'SOD1', 'DYRK1A']

        # Plotting the heatmap with seaborn
        print("- Plotting...\n")
        fig, ax = plt.subplots(figsize=(7, 6))
        sns.heatmap(df, annot=True, cmap="Spectral", linewidths=0.05, cbar=True, fmt='.0f', center=0, ax=ax)
        ax.set_title("Heatmap of SNPs for each gene", fontsize=18, pad=21)
        plt.tight_layout()
        plt.savefig(output.heat, dpi=300)
        print("Done.")
