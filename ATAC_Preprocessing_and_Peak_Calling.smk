####################################################
## Darragh Nimmo
## Trinity College Dublin
## March 2021
# Snakemake workflow for ATAC-seq preprocessing and peak calling.
####################################################


###############################################################################################
###IMPORTS and variables
###############################################################################################

import os

import functools

configfile: 'config.yaml'

READ_1 = config['read_1']

READ_2 = config['read_2']

GENOME = config['genome']

SAMPLE = config['Sample']

SUBSAMPLE = config['subsample']

BLACK_LIST = config['black_list'] 

PICARD_JAR = config['picard_jar']

ATAC_script = config['ATACseqQC_script']

chrM_SCRIPT = config['chrM_script']

bedpe_SCRIPT = config['bedpe_script']

CONDA_ENV_DIR = config['conda_env_dir']

directory_function = functools.partial(os.path.join, config['results'])
FASTQC_DIR = directory_function('Fastqc')
BAM_DIR = directory_function('Bam')
READS_DIR = directory_function('READS')
ATACseqQC_DIR = directory_function('ATACseqQC')
PEAK_DIR = directory_function('Peak')

###############################################################################################
###Rules
###############################################################################################

rule all:
        input:
                os.path.join(PEAK_DIR, SAMPLE+'_processed_peaks.bed')
rule fastqc:
    input:
        read1 = READ_1,
        read2 = READ_2
    output:
        os.path.join(FASTQC_DIR, os.path.basename(READ_1).rstrip(".fq.gz")+"_fastqc.html"),
        os.path.join(FASTQC_DIR, os.path.basename(READ_2).rstrip(".fq.gz")+"_fastqc.html"),
        os.path.join(FASTQC_DIR, os.path.basename(READ_1).rstrip(".fq.gz")+"_fastqc.zip"),
        os.path.join(FASTQC_DIR, os.path.basename(READ_2).rstrip(".fq.gz")+"_fastqc.zip"),
    params:
        dir = FASTQC_DIR
    conda:
        os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
    message:
        "Examining quality of sequencing reads using fastqc."
    shell:
        """
        fastqc -t 48 {input.read1} {input.read2} -o {params.dir}
        """

rule trim_adapters:
    input:
        read1 = READ_1,
        read2 = READ_2
    output:
        read1 = os.path.join(READS_DIR, os.path.basename(READ_1).rstrip(".fq.gz")+"_val_1.fq.gz"),
        read2 = os.path.join(READS_DIR, os.path.basename(READ_2).rstrip(".fq.gz")+"_val_2.fq.gz")
    params:
        dir = READS_DIR
conda:
        os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
    message:
        "Trimming adapters from the sequencing reads using trim galore."
    shell:
        """
        trim_galore --paired {input} -o {params.dir}
        """


rule alignment:
        input:
                read1 = rules.trim_adapters.output.read1,
                read2 = rules.trim_adapters.output.read2
        params:
                index = GENOME
        output:
                bam = protected(os.path.join(BAM_DIR, SAMPLE+'_mapped.bam'))
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
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
  message:
            "Coordinate sorting and indexing alignment bam."
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index -@ 48 {output.bam}
                """

rule remove_chrM:
        input:
                bam = rules.coordinate_sort_index_1.output.bam,
                index = rules.coordinate_sort_index_1.output.index
        output:
                bam = temp(os.path.join(BAM_DIR,SAMPLE+'_noMT.bam'))
	params:
		script = chrM_SCRIPT
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Removing mitochondrial alignments from the bam file."
        shell:
                """
                bash {params.script} {input.bam} {output.bam}
                """

rule coordinate_sort_index2:
        input:
                bam = rules.remove_chrM.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, SAMPLE+'_noMT.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, SAMPLE+'_noMT.sorted.bam.bai'))
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
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Coordinate sorting and indexing the bam file of just paired reads"
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index {output.bam}
                """

#rule ATACseqQC_quality_checks:
#    input:
#        bam = rules.coordinate_sort_index3.output.bam
#    output:
#        Lib_complex_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_library_complexity.pdf'),
#        Lib_complexity_table = os.path.join(ATACseqQC_DIR, SAMPLE+'_library_complexity.txt'),
#        Frag_size_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_fragment_size_distribution.pdf'),
#        TSS_enrichment_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_transcription_start_site.pdf'),
#        Heat_map_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_heat_map.pdf'),
#        Feat_align_dist_plot = os.path.join(ATACseqQC_DIR, SAMPLE+'_feature_aligned_distribution.pdf')
#    params:
#        script = ATAC_script
#    conda:
#        os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
#    message:
#        "Creating Plots to assess the quality of the ATAC-seq data with ATACseqQC"
#
#    shell:
#        """
#        Rscript {params.script} {input.bam}
#        """

rule subsample_by_lib_complexity:
        input:
                bam = rules.coordinate_sort_index3.output.bam,
                index = rules.coordinate_sort_index3.output.index,
        output:
                bam = temp(os.path.join(BAM_DIR, SAMPLE+'_subsampled.bam'))
	params:
		sub = SUBSAMPLE
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Subsampling the bam file according to the level of library complexity"
        shell:
                """
                samtools view -@ 48 -h -b -s {params.sub} {input.bam} > {output.bam}
                """
rule coordinate_sort_index4:
        input:
                bam = rules.subsample_by_lib_complexity.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, SAMPLE+'_subsampled.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, SAMPLE+'_subsampled.sorted.bam.bai'))
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
        bam= temp(os.path.join(BAM_DIR, SAMPLE+'_markedDuplicates.bam')),
        metrics= temp(os.path.join(BAM_DIR, SAMPLE+"_markedDuplicates.txt"))
    params:
	picard = PICARD_JAR
    conda:
        os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
    message:
        "Removing duplicates from the bam file."
    shell:
        """
        java -jar {params.picard} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true
        """

rule coordinate_sort_5:
        input:
                bam = rules.Mark_Duplicates.output.bam
        output:
                bam = protected(os.path.join(BAM_DIR, SAMPLE+'_md.sorted.bam'))
        conda:
            os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
        message:
            "Coordinate sorting bam file by read name."
        shell:
                """
                samtools sort -@ 48 -n -o {output.bam} {input.bam}
                """
		
rule bam_to_bed_and_Tn5_shift:
    input:
        bam = rules.coordinate_sort_5.output.bam
    output:
        bedpe = temp(os.path.join(BAM_DIR, SAMPLE+'_md.sorted.bedpe'))
    conda:
        os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
    message:
        "Converting bam file to a bed file and performing Tn5 shift."
    shell:
        """
        samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 {input.bam} |bamToBed -i stdin -bedpe | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($9 == "+") {$2 = $2 + 4; $6 = $6 - 5} \
	else if ($9 == "-") {$3 = $3 - 5; $5 = $5 + 4} print $0}' > {output.bedpe}
        """
rule Convert_bedpe_macs:
	input:
		bedpe = rules.bam_to_bed_and_Tn5_shift.output.bedpe
	output:
		bedpe = protected(os.path.join(BAM_DIR, SAMPLE+'_macs_input.bedpe'))
	params:
		script = bedpe_SCRIPT
	conda:
		os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
	message:
        	"Converting bedpe to MACS2 compatible format."
	shell:
		"""
		bash {params.script} {input.bedpe} > {output.bedpe}
		"""
	
	
rule peak_calling:
    input:
        bed = rules.Convert_bedpe_macs.output.bed
    output:
        BroadPeak = protected(os.path.join(PEAK_DIR, SAMPLE+'_peaks.broadPeak'))
    params:
        dir = PEAK_DIR,
        sample = SAMPLE
    conda:
        os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
    message:
        "Peak calling."
    shell:
        """
        macs2 callpeak -t {input.bed} --nomodel --shift -100 --extsize 200 --broad -f BED --name {params.sample} -g hs --outdir {params.dir}
        """

rule preprocess_peaks:
	input:
		BroadPeak = rules.peak_calling.output.BroadPeak,
	output:
		bed = protected(os.path.join(PEAK_DIR, SAMPLE+'_processed_peaks.bed'))
	params:
		BlackList = BLACK_LIST
	conda:
		os.path.join(CONDA_ENV_DIR, 'ATACseq_Preprocessing_env.yaml')
	message:
		"Preprocessing peaks."
	shell:
		"""
		bedtools subtract -a {input.BroadPeak} -b {params.BlackList} > {output.bed}
		"""
