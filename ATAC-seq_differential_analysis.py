####################################################
## Darragh Nimmo
## Trinity College Dublin
## April 2021
## Snakemake workflow for ATAC-seq differential accessibility analysis.
## An adaptation of the workfow designed by Reske et al 2020 with modifications and written as a snakemake workflow.
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

SUBSAMPLE = config['subsample']

PICARD_JAR = config['picard_jar']

SAMPLES = config['samples']

SUBSAMPLE = config['subsample']

SAMPLE1 = config['sample1']

SAMPLE2 = config['sample2']

BEDPE_SCRIPT = config['bedpe_script'] 

BLACK_LIST = config['black_list']

NAIVE_OVERLAP = CONFIG['naive_overlap']

directory_function = functools.partial(os.path.join, config['results'])
READS_DIR = directory_function('READS')
BAM_DIR = directory_fucntion('Bam')
PEAK_DIR = directory_fucntion('Peak')
BED_DIR = directory_function('Bed')

###############################################################################################
###Rules
###############################################################################################

rule all:
  input:
     
     
rule trim_adapters:
    input:
        read1 = {cell}_rep_{replicate}_READ_1,
        read2 = {cell}_rep_{replicate}_READ_2
    output:
        read1 = os.path.join(READS_DIR, os.path.basename({cell}_rep_{replicate}_READ_1).rstrip(".fq.gz")+"_val_1.fq.gz"),
        read2 = os.path.join(READS_DIR, os.path.basename({cell}_rep_{replicate}_READ_2).rstrip(".fq.gz")+"_val_2.fq.gz")
    params:
        dir = READS_DIR
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
                bam = protected(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_mapped.bam'))
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
                bam = temp(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_mapped.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_mapped.sorted.bam.bai'))
        message:
            "Coordinate sorting and indexing alignment bam files."
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index -@ 48 {output.bam}
                """
            
rule remove_chrM:
        input:
                bam = rules.coordinate_sort_index_1.output.bam,
                index = rules.coordinate_sort_index_1.output.index
        output:
                bam = os.path.join(BAM_DIR, '{cell}_rep_{replicate}_noMT.bam')
	      params:
		            script = chrM_SCRIPT
        message:
            "Removing mitochondrial alignments from the bam files."
        shell:
                """
                bash {params.script} {input.bam} {output.bam}
                """
            
rule coordinate_sort_index2:
        input:
                bam = rules.remove_chrM.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_noMT.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_noMT.sorted.bam.bai'))
        message:
            "Coordinate sorting and indexing chrM bamfiles"
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam};
                 samtools index -@ 48 {output.bam}
                """     
            
rule paired_reads_only:
        input:
                bam = rules.coordinate_sort_index2.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_paired.bam'))
        message:
            "Keeping just paired reads in the bam files."
        shell:
                """
                samtools view -@ 48 -bh -f 3 {input.bam} > {output.bam}
                """
            
rule coordinate_sort_index3:
        input:
                bam = rules.paired_reads_only.output.bam
        output:
                bam = protected(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_paired.sorted.bam')),
                index = protected(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_paired.sorted.bam.bai'))
        message:
            "Coordinate sorting and indexing the bam files of just paired reads"
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index {output.bam}
                """
     
rule ATACseqQC_quality_checks:
    input:
        bam = rules.coordinate_sort_index3.output.bam
    output:
        Lib_complex_plot = os.path.join(ATACseqQC_DIR, '{cell}_rep_{replicate}_library_complexity.pdf'),
        Lib_complexity_table = os.path.join(ATACseqQC_DIR, '{cell}_rep_{replicate}_library_complexity.txt'),
        Frag_size_plot = os.path.join(ATACseqQC_DIR, '{cell}_rep_{replicate}_fragment_size_distribution.pdf'),
    params:
        script = ATAC_script
    message:
        "Creating Plots to assess the quality of the ATAC-seq data with ATACseqQC"
    shell:
        """
        Rscript {params.script} {input.bam}
        """
#stop here and check the level you need to subsample libraries - check the Lib_complexity_table.
#Need to add the input file and subsample parameters

rule subsample_by_lib_complexity:
        input:
                bam = rules.coordinate_sort_index3.output.bam,
                index = rules.coordinate_sort_index3.output.index,
        output:
                bam = temp(os.path.join(BAM_DIR, '_subsampled.bam'))
	      params:
		            sub = SUBSAMPLE
        message:
            "Subsampling the bam files according to the level of library complexity"
        shell:
                """
                samtools view -@ 48 -h -b -s {params.sub} {input.bam} > {output.bam}
                """
 
rule coordinate_sort_index4:
        input:
                bam = rules.subsample_by_lib_complexity.output.bam
        output:
                bam1 = protected(os.path.join(BAM_DIR, SAMPLE1+'_md_input.bam')),
                index = protected(os.path.join(BAM_DIR, SAMPLE1+'_md_input.bam.bai'))
        message:
            "Coordinate sorting and indexing the bam files following subsampling."
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index {output.bam}
                """
rule change_name:
    input:
        not_subsampled_bam = os.path.join(BAM_DIR, SAMPLE2+'_paired.sorted.bam'))
        not_subsampled_index = os.path.join(BAM_DIR, SAMPLE2+'_paired.sorted.bam.bai'))
    output:
        not_subsampled_bam = os.path.join(BAM_DIR, SAMPLE2+'_md_input.bam'))
        not_subsampled_index = os.path.join(BAM_DIR, SAMPLE2+'_md_input.bam.bai'))
    shell:
        """
        mv {input.not_subsampled_bam} {output.not_subsampled_bam}; mv {input.not_subsampled_index} {output.not_subsampled_index}
        """

            
rule mark_duplicates:
    input:
        bam = os.path.join(BAM_DIR, '{cell}_rep_{replicate}_md_input.bam')
        index = os.path.join(BAM_DIR, '{cell}_rep_{replicate}_md_input.bam.bai')
    output:
        bam= temp(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_markedDuplicates.bam')),
        metrics= temp(os.path.join(BAM_DIR, '{cell}_rep_{replicate}_markedDuplicates.txt'))
    params:
	      picard = PICARD_JAR
    message:
        "Removing duplicates from the bam files."
    shell:
        """
        java -jar {params.picard} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true
        """  

rule coordinate_sort_5:
        input:
                bam = rules.mark_duplicates.output.bam
        output:
                bam = os.path.join(BAM_DIR, '{cell}_rep_{replicate}_md.sorted.bam')
        message:
            "Coordinate sorting bam files by read name."
        shell:
                """
                samtools sort -@ 48 -n -o {output.bam} {input.bam}
                """  
rule fixmate:
    input:
        bam = rules.coordinate_sort_5.output.bam
    output:
        bam = os.path.join(BAM_DIR, '{cell}_rep_{replicate}_fixmate.bam'
    shell:
        """
        samtools fixmate {input.bam} {output.bam}
        """

rule bam_to_bed_and_Tn5_shift:
    input:
        bam = rules.fixmate.output.bam
    output:
        bedpe = os.path.join(BED_DIR, '{cell}_rep_{replicate}_md.sorted.bedpe')
    message:
        "Converting bam files to a bedpe files and performing the Tn5 shift."
    shell:
        """
        samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 {input.bam} |bamToBed -i stdin -bedpe | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($9 == "+") {$2 = $2 + 4; $6 = $6 - 5} \
	else if ($9 == "-") {$3 = $3 - 5; $5 = $5 + 4} print $0}' > {output.bedpe}
        """
                           
rule convert_bedpe_macs:
    input:
         bedpe = rules.bam_to_bed_and_Tn5_shift.output.output.bedpe
    output:
         bedpe = os.path.join(BED_DIR, '{cell}_rep_{replicate}_macs_input.bedpe')
    params:
         script = BEDPE_SCRIPT
    message:
        	"Converting bedpe files to a MACS2 compatible format."
    shell:
        """
        bash {params.script} {input.bedpe} > {output.bedpe}
        """
                           
rule peak_calling_all:
    input:
        bedpe = os.path.join(BED_DIR, '{cell}_rep_{replicate}_macs_input.bedpe'),
    output:
        BroadPeak = os.path.join(PEAK_DIR, '{cell}_rep_{replicate}_peaks.broadPeak')
    params:
        dir = PEAK_DIR,
        sample = {cell}_rep_{replicate}
    message:
        "Peak calling all."
    shell:
        """
        macs2 callpeak -t {input.bedpe} --nomodel --shift -100 --extsize 200 --broad -f BED --name {params.sample} -g hs --outdir {params.dir}
        """
                          
preprocess_peaks_all:
       input:
            peaks = rules.peak_calling_all.output.BroadPeak
       output:
            peaks = os.path.join(PEAK_DIR, '{cell}_rep_{replicate}_processed_peaks.broadPeak')
       params:
            black_list = BLACK_LIST
       shell:
           """
           bedtools intersect -v -a {input.peaks} -b {params.black_list} | grep -P 'chr[\tdXY]+[\t]' > {output.peaks}
           """
                           
rule peak_calling_pooled:
       input:
            rep1 = os.path.join(BED_DIR, '{cell}_rep_1_macs_input.bedpe'),
            rep2 = os.path.join(BED_DIR, '{cell}_rep_2_macs_input.bedpe') 
        output:
            BroadPeak = os.path.join(PEAK_DIR, '{cell}_pooled_reps.broadPeak')   
        params:
            dir = PEAK_DIR,
            sample = {cell}
        message:
            "Peak calling pooled"
        shell:
            """
            macs2 callpeak -t {input.rep1} {input.rep2} --nomodel --shift -100 --extsize 200 --broad -f BED --name {params.sample} -g hs --outdir {params.dir}
            """
                           
rule preprocess_peaks_pooled:
       input:
            peaks = rules.peak_calling_pooled.output.BroadPeak
       output:
            peaks = os.path.join(PEAK_DIR, '{cell}_pooled_reps_processed_peaks.broadPeak')
       params:
            black_list = BLACK_LIST
       shell:
           """
           bedtools intersect -v -a {input.peaks} -b {params.black_list} | grep -P 'chr[\tdXY]+[\t]' > {output.peaks}
           """
                           
rule naive_overlap:
               input:
                    rep1 = os.path.join(PEAK_DIR, '{cell}_rep_1_processed_peaks.broadPeak'),
                    rep2 = os.path.join(PEAK_DIR, '{cell}_rep_2_processed_peaks.broadPeak')
                    pooled = os.path.join(PEAK_DIR, '{cell}_pooled_reps_processed_peaks.broadPeak')
               output:
                    overlap = os.path.join(PEAK_DIR, '{cell}_pooled_reps_processed_peaks.broadPeak')
               params:
                    script = NAIVE_OVERLAP
               shell:
                    """
                    bash {params.script} {input.rep1} {input.rep2} {input.pooled} > {output.overlap}
                    """
                           
                           
                           
                           
                           
                           
                           
    
    
