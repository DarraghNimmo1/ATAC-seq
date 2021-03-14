####################################################
##Darragh Nimmo
## March 2021
#snakemake workflow for ATAC-seq preprocessing.
####################################################


###############################################################################################
###IMPORTS and variables
###############################################################################################

import os

import functools

configfile: 'config.yaml'

READ_1 = config['read_1']

READ_2 = config['read_2']

THREADS = config['threads']

GENOME = config['genome']

GENOME2 = config['genome2']

SAMPLE = config['Sample']

ATACQC = config['ATACseqQC_script']

CONDA_ENV_DIR = config['conda_env_dir']

directory_function = functools.partial(os.path.join, config['results'])
FASTQC_DIR = directory_function('Fastqc')
BAM_DIR = directory_function('Bam')
READS_DIR = directory_function('READS')
ATACseqQC_DIR = directory_function('ATACseqQC')

###############################################################################################
###Rules
###############################################################################################

rule all:
        input:
                SAMPLE+'_paired.sorted.bam',
                SAMPLE+'_paired.sorted.bam.complexity.txt'

rule fastqc:
        input:
                read1 = READ_1
                read2 = READ_2
        output:
                os.path.join(FASTQC_DIR, os.path.basename(READ_1).rstrip(".fq.gz")+"_fastqc.html"),
                os.path.join(FASTQC_DIR, os.path.basename(READ_2).rstrip(".fq.gz")+"_fastqc.html"),
                os.path.join(FASTQC_DIR, os.path.basename(READ_1).rstrip(".fq.gz")+"_fastqc.zip"),
                os.path.join(FASTQC_DIR, os.path.basename(READ_2).rstrip(".fq.gz")+"_fastqc.zip"),
                os.path.join(FASTQC_DIR,'multiqc_data')
                os.path.join(FASTQC_DIR, 'multiqc_report.html')
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Examining quality of sequencing reads using fastqc and multiqc."   
        shell:
                """
                fastqc -t 48 {input.read1} {input.read2}| multiqc -o FASTQC_DIR  .
                """
                
rule trim_adapters:
        input:
                read1 = READ_1,
                read2 = READ_2
        output:
                read1 = os.path.join(READS_DIR, os.path.basename(READ_1).rstrip(".fq.gz")+"_val_1.fq.gz"),
                read2 = os.path.join(READS_DIR, os.path.basename(READ_2).rstrip(".fq.gz")+"_val_2.fq.gz")
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Trimming reads from the sequencing reads using trim galore."
        shell:
                """
                trim_galore  --paired {input}
                """
                
rule alignment:
        input:
                read1 = rules.trim_adapters.output.read1,
                read2 = rules.trim_adapters.output.read2
        params:
                index = GENOME
        output:
                bam = protected(os.path.join(BAM_DIR, SAMPLE+'_mapped.bam'))
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Aligning the sequencing reads to hg38 using bowtie2. This may take a while..."
        shell:
                """
                bowtie2 --very-sensitive -X 1000 -x {params.index} --threads 48
                -1 {input.read1} -2 {input.read2}| samtools view -@ 48 -bS - > {output}
                """
                
rule coordinate_sort_index_1:
        input:
                bam = rules.alignment.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, SAMPLE+'_mapped.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, SAMPLE+'_mapped.sorted.bam.bai'))
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')   
        message:
            "Coordinate sorting and indexing alignment bam."
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index1
                -@ 48 {output.bam}
                """
                
rule remove_chrM:
        input:
                bam = rules.coordinate_sort_index_1.output.bam,
                index = rules.coordinate_sort_index_1.output.index
        output:
                bam = temp(os.path.join(BAM_DIR,SAMPLE+'_noMT.bam'))
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Removing mitochondrial alignments from the bam file."
        shell:
                """
                samtools view -@ 48 -h -f 3 {input.bam}| python
                 /home/darragh/ATAC-Seq/Python_scripts/removeChrom.py  - - chrM
                  | samtools view -@ 48 -bh - > {output}
                """

rule coordinate_sort_index2:
        input:
                bam = rules.remove_chrM.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, SAMPLE+'_noMT.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, SAMPLE+'_noMT.sorted.bam.bai'))
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Coordinate sorting and indexing -chrM bamfile"
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam};
                 samtools index -@ 48 {output.bam}
                """

rule paired_reads_only:
        input:
                bam = rules.coordinate_sort_index2.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, SAMPLE+'_paired.bam'))
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Keeping just paired reads in the bam file."
        shell:
                """
                samtools view -@ 48 -bh -f 3 {input.bam} > {output.bam}
                """

rule coordinate_sort_index3:
        input:
                bam = rules.paired_reads_only.output.bam
        output:
                bam = protected(os.path.join(BAM_DIR, SAMPLE+'_paired.sorted.bam')),
                index = protected(os.path.join(BAM_DIR, SAMPLE+'_paired.sorted.bam.bai'))
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Coordinate sorting and indexing the bam file of just paired reads"
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index {output.bam}
                """

rule ATACseqQC_quality_checks:
    input:
        bam = rules.coordinate_sort_index3.output.bam
    output:
        Lib_complex_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_library_complexity.pdf'),
        Lib_complexity_table = os.path.join(ATACseqQC_DIR, SAMPLE+'_library_complexity.txt'),
        Frag_size_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_fragment_size_distribution.pdf'),
        TSS_enrichment_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_transcription_start_site.pdf'),
        Heat_map_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_heat_map.pdf'),
        Feat_align_dist_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_heat_map.pdf')
    threads: THREADS
    conda:
        os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
    message:
        "Creating Plots to assess the quality of the ATAC-seq data with ATACseqQC"
    shell:
        """
        Rscript ATACQC {input.bam}
        """

rule subsample_by_lib_complexity:
        input:
                bam = rules.coordinate_sort_index3.output.bam,
                index = rules.coordinate_sort_index3.output.index,
        output:
                bam = temp(SAMPLE+'_subsampled.bam')
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Subsampling the bam file according to the level of library complexity"
        shell:
                """
                samtools view -@ 48 -h -b -s 1.33 {input.bam} > {output.bam}
                """

rule coordinate_sort_index4:
        input:
                bam = rules.subsample_by_lib_complexity.output.bam
        output:
                bam = temp(SAMPLE+'_subsampled.sorted.bam'),
                index = temp(SAMPLE+'_subsampled.sorted.bam.bai')
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Coordinate sorting and indexing the bam file following subsampling."
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index {output.bam}
                """
                
rule mark_duplicates:
    input:
        bam=rules.coordinate_sort_index4.output.bam
    output:
        bam= temp(SAMPLE+'_markedDuplicates.bam'),
        metrics= temp(SAMPLE+"_markedDuplicates.txt")
    threads: THREADS
    conda:
        os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
    message:
        "Removing duplicates from the bam file."
    shell:
        """
        "java -jar /home/darragh/picard.jar MarkDuplicates -I {input.bam} 
        -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true
        """
        
rule coordinate_sort_index5:
        input:
                bam = rules.mark_duplicates.output.bam
        output:
                bam = protected(SAMPLE+'_md.sorted.bam'),
                index = protected(SAMPLE+'_md.sorted.bam.bai')
        threads: THREADS
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Coordinate sorting and indexing the final bam file."
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index {output.bam}
                """
