
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

configfile: 'config.yaml'

BAMS = config['bams']

CELL_TYPE = config['cell_type']

REPLICATE = config['replicate']

GENOME = config['genome']

BLACK_LIST = config['black_list']


directory_function = functools.partial(os.path.join, config['results'])
IN_BAM = directory_function('In_Bam')
OUT_BAM = directory_function('Out_Bam')
PEAK_DIR  = directory_function('Peaks')
ATAC_CORRECT = directory_function('ATAC_Correct')
FOOT_PRINT_SCORE = directory_function('Foot_Print_Score')

###############################################################################################
###Rules
###############################################################################################

rule all:
        input:


rule Peak_calling:
    input:
        bam = expand(os.path.join(IN_BAM, '{cell_type}_rep_{replicate}.bam'), cell_type = CELL_TYPE, replicate = REPLICATE)
    output:
        BroadPeak = expand(os.path.join(PEAK_DIR, '{cell_type}_rep_{replicate}.broadPeak'), cell_type = CELL_TYPE, replicate = REPLICATE)
    params:
	name = expand('{cell_type}_rep_{replicate}', cell_type = CELL_TYPE, replicate = REPLICATE)
        dir = PEAK_DIR
    message:
        "Peak calling."
    shell:
        """
        macs2 callpeak -t {input.bam} --nomodel --shift -100 --extsize 200 --broad -f BAM --name {params.name} -g hs --outdir {params.dir}
        """
		     
rule preprocess_peaks:
	input:
		BroadPeak = rules.peak_calling.output.BroadPeak,
	output:
		BroadPeak = expand(os.path.join(PEAK_DIR, '{cell_type}_rep_{replicate}_preprocessed.broadPeak'), cell_type = CELL_TYPE, replicate = REPLICATE)
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
                rep1 = expand(os.path.join(IN_BAM, '{cell_type}_rep_one.bam'), cell_type = CELL_TYPE),
                rep2 = expand(os.path.join(IN_BAM, '{cell_type}_rep_two.bam'), cell_type = CELL_TYPE)
        output:
                bam = expand(os.path.join(IN_BAM, '{cell_type}_merged.bam'), cell_type = CELL_TYPE),
                index = expand(os.path.join(OUT_BAM, '{cell_type}_merged.bam.bai'), cell_type = CELL_TYPE)
        message:
            "Merging replicate bam files and indexing."
        shell:
            """
            java -jar PICARD_JAR MergeSamFiles I= {input.rep1} I= {input.rep2} O= {output.bam};
            samtools index -@ 48 {output.bam}
            """

rule merge_peaks:
	input:
		rep1 = expand(os.path.join(PEAK_DIR, '{cell_type}_rep_one_preprocessed.broadPeak'), cell_type = CELL_TYPE),
		rep2 = expand(os.path.join(PEAK_DIR, '{cell_type}_rep_two_preprocessed.broadPeak'), cell_type = CELL_TYPE)
	output:
		peaks = expand(os.path.join(PEAK_DIR, '{cell_type}_merged.bed'), cell_type = CELL_TYPE)
	message:
		"Merging replicate peak files."
	shell:
		"""
		cat {input.rep1} {input.rep2}| sort -k1,1 -k2,2n | bedtools merge -d 5 -c 4 -o distinct > {output.peaks}
		"""
rule ATAC_correct:
	input:
		peak = rules.merge_peaks.ouput.peaks,
		bam = rules.merge_bams.output.peaks
	output:
		uncorrected = expand(os.path.join(ATAC_CORRECT, '{cell_type}_uncorrected.bw'), cell_type = CELL_TYPE),
		bias = expand(os.path.join(ATAC_CORRECT, '{cell_type}_bias.bw'), cell_type = CELL_TYPE),	     
		expected = expand(os.path.join(ATAC_CORRECT, '{cell_type}_expected.bw'), cell_type = CELL_TYPE),
		corrected = expand(os.path.join(ATAC_CORRECT, '{cell_type}_corrected.bw'), cell_type = CELL_TYPE),	      
		pdf = expand(os.path.join(ATAC_CORRECT, '{cell_type}_atacorrect.pdf'), cell_type = CELL_TYPE),
		pickle = expand(os.path.join(ATAC_CORRECT, '{cell_type}_AtacBias.pickle'), cell_type = CELL_TYPE)
	params:
		genome = GENOME,
		black_list = BLACK_LIST,
		dir = ATAC_CORRECT,
		name = expand('{cell_type}', cell_type = CELL_TYPE)
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
		bigwig = expand(os.path.join(FOOT_PRINT_SCORE, '{cell_type}_FootPrintScore.bw'), cell_type = CELL_TYPE)
	message:
		"Calculating footprint score."
	shell:
		"""
		TOBIAS FootprintScores --signal {input.signal} --regions {input.region} --output {output.bigwig} --cores 48 
		"""
		
		
		
		
		
