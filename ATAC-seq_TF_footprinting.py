
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

###############################################################################################
###Rules
###############################################################################################

rule all:
        input:


rule Peak_calling:
    input:
        bam = expand(os.path.join(IN_BAM, '{cell_type}_rep_{replicate}.bam', cell_type = CELL_TYPE, replicate = REPLICATE)
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
                rep1 = expand(os.path.join(IN_BAM, '{cell_type}_rep_one.bam', cell_type = CELL_TYPE),
                rep2 = expand(os.path.join(IN_BAM, '{cell_type}_rep_two.bam', cell_type = CELL_TYPE)
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
		rep1 = expand(os.path.join(PEAK_DIR, '{cell_type}_rep_one_preprocessed.broadPeak', cell_type = CELL_TYPE),
		rep2 = expand(os.path.join(PEAK_DIR, '{cell_type}_rep_two_preprocessed.broadPeak', cell_type = CELL_TYPE)
	output:
		peaks = expand(os.path.join(PEAK_DIR, '{cell_type}_merged.broadPeak', cell_type = CELL_TYPE)
	message:
		"Merging replicate peak files."
	shell:
		"""
		"""
			     
			      
			      
