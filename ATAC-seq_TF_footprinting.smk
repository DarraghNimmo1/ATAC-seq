####################################################
## Darragh Nimmo
## Trinity College Dublin
## March 2021
# Snakemake workflow for ATAC-seq footprinting.
####################################################

###############################################################################################
###IMPORTS and variables
###############################################################################################

import os

import functools

configfile: 'config.yaml'

CELL_TYPE = config['cell_type']

REPLICATE = config['replicate']

IN_BAM = config['in_bam']

BLACK_LIST = config['black_list']

PICARD_JAR = config['picard_jar']

GENOME = config['genome']

MOTIFS = config['motifs']

directory_function = functools.partial(os.path.join, config['results'])
PEAK_DIR  = directory_function('Peaks')
OUT_BAM = directory_function('Out_Bam')
ATAC_CORRECT = directory_function('ATAC_Correct')
FOOT_PRINT_SCORE = directory_function('Foot_print_score')
BIN_DETECT = directory_function('Bin_Detect')

###############################################################################################
###Rules
###############################################################################################

rule all:
    input:
        expand(os.path.join(PEAK_DIR, '{cell_type}_rep_{replicate}_peaks.broadPeak'), cell_type = CELL_TYPE, replicate = REPLICATE),
        expand(os.path.join(PEAK_DIR, '{cell_type}_rep_{replicate}_preprocessed.broadPeak'), cell_type = CELL_TYPE, replicate = REPLICATE),
        expand(os.path.join(OUT_BAM, '{cell_type}_merged.bam'), cell_type = CELL_TYPE),
        expand(os.path.join(PEAK_DIR, '{cell_type}_merged.bed'), cell_type = CELL_TYPE),
        expand(os.path.join(ATAC_CORRECT, '{cell_type}_corrected.bw'), cell_type = CELL_TYPE),
        expand(os.path.join(FOOT_PRINT_SCORE, '{cell_type}_FootPrintScore.bw'), cell_type = CELL_TYPE),
        expand(os.path.join(BIN_DETECT, '{cell_type}', 'bindetect_figures.pdf'), cell_type = CELL_TYPE)

rule peak_calling:
    input:
        bam = os.path.join(IN_BAM, '{cell_type}_rep_{replicate}.bam')
    output:
        BroadPeak = os.path.join(PEAK_DIR, '{cell_type}_rep_{replicate}_peaks.broadPeak')
    params:
        name = '{cell_type}_rep_{replicate}',
        dir = PEAK_DIR
    message:
        "Peak calling."
    shell:
        """
        macs2 callpeak -t {input.bam} --nomodel --shift -100 --extsize 200 --broad -f BAM --name {params.name} -g hs --outdir {params.dir}
        """
rule preprocess_peaks:
        input:
                BroadPeak = os.path.join(PEAK_DIR, '{cell_type}_rep_{replicate}_peaks.broadPeak')
        output:
                BroadPeak = os.path.join(PEAK_DIR, '{cell_type}_rep_{replicate}_preprocessed.broadPeak')
        params:
                BlackList = BLACK_LIST
        message:
                "Preprocessing peaks."
        shell:
                """
                bedtools subtract -a {input.BroadPeak} -b {params.BlackList} > {output.BroadPeak}
                """

rule merge_bams:
    input:
        rep1 = os.path.join(IN_BAM, '{cell_type}_rep_one.bam'),
        rep2 = os.path.join(IN_BAM, '{cell_type}_rep_two.bam')
    output:
        bams = os.path.join(OUT_BAM, '{cell_type}_merged.bam'),
        index = os.path.join(OUT_BAM, '{cell_type}_merged.bam.bai')
    params:
        jar = PICARD_JAR
    message:
        "Merging replicate bam files and indexing."
    shell:
        """
        java -jar {params.jar} MergeSamFiles I= {input.rep1} I= {input.rep2} O= {output.bam};
        samtools index -@ 48 {output.bam}
        """



rule merge_peaks:
    input:
        rep1 = os.path.join(PEAK_DIR, '{cell_type}_rep_one_preprocessed.broadPeak'),
        rep2 = os.path.join(PEAK_DIR, '{cell_type}_rep_two_preprocessed.broadPeak')
    output:
        peaks = os.path.join(PEAK_DIR, '{cell_type}_merged.bed')
    message:
        "Merging replicate peak files."
    shell:
        """
        cat {input.rep1} {input.rep2}| sort -k1,1 -k2,2n | bedtools merge -d 5 -c 4 -o distinct > {output.peaks}
        """

rule ATAC_correct:
    input:
        peak = rules.merge_peaks.output.peaks,
        bam = rules.merge_bams.output.bams
    output:
        uncorrected = os.path.join(ATAC_CORRECT, '{cell_type}_uncorrected.bw'),
        bias = os.path.join(ATAC_CORRECT, '{cell_type}_bias.bw'),
        expected = os.path.join(ATAC_CORRECT, '{cell_type}_expected.bw'),
        corrected = os.path.join(ATAC_CORRECT, '{cell_type}_corrected.bw'),
        pdf = os.path.join(ATAC_CORRECT, '{cell_type}_atacorrect.pdf'),
        pickle = os.path.join(ATAC_CORRECT, '{cell_type}_AtacBias.pickle')
    params:
        genome = GENOME,
        black_list = BLACK_LIST,
        dir = ATAC_CORRECT,
        name = '{cell_type}'
    message:
        "Correcting for Tn5 insertion bias and shifting reads to correct for Tn5 insertions."
    shell:
        """
        TOBIAS ATACorrect --bam {input.bam} --genome {params.genome} --peaks {input.peak} --blacklist {params.black_list} --outdir {params.dir} --prefix {params.name}  --cores 48
        """

rule FootPrintScore:
    input:
        signal = rules.ATAC_correct.output.corrected,
        region = rules.merge_peaks.output.peaks
    output:
        bigwig = os.path.join(FOOT_PRINT_SCORE, '{cell_type}_FootPrintScore.bw')
    message:
        "Calculating footprint score."
    shell:
        """
        TOBIAS ScoreBigwig --signal {input.signal} --regions {input.region} --output {output.bigwig} --cores 48
        """

rule BinDetect:
    input:
        signal = rules.FootPrintScore.output.bigwig,
        region = rules.merge_peaks.output.peaks
    output:
        TF = os.path.join(BIN_DETECT, '{cell_type}', 'bindetect_figures.pdf')
    params:
        genome = GENOME,
        motifs = MOTIFS,
        dir = os.path.join(BIN_DETECT, '{cell_type}')
    message:
        "Running bin detect on to make predictions on specific transcription factor binding at peaks"
    shell:
        """
        TOBIAS BINDetect --motifs {params.motifs} --signals {input.signal} --genome {params.genome} --peaks {input.region} --outdir {params.dir} --cores 48
        """
