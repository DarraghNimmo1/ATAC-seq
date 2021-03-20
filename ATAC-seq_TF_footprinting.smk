
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

SAMPLES = config['samples']

GENOME = config['genome']


directory_function = functools.partial(os.path.join, config['results'])
IN_BAM = directory_function('In_Bam')
OUT_BAM = directory_function('Out_Bam')

###############################################################################################
###Rules
###############################################################################################

rule all:
        input:

		
rule Tn5_shift:
    input:
        bam = expand(os.path.join(IN_BAM, '{sample}_md.sorted.bam', sample=SAMPLES)
    output:
        bam = expand(os.path.join(OUT_BAM, '{sample}_Tn5_shifted.bam', sample=SAMPLES)
    message:
        "Converting bam file to bed file, performing Tn5 shift, converting back to bam file."
    shell:
        """
        samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 {input.bam} |bamToBed -i stdin | awk -F $'\t' 'BEGIN {{OFS = FS}}{{ if ($9 == "+") {{$2 = $2 + 4; $6 = $6 - 5}} \
	else if ($9 == "-") {{$3 = $3 - 5; $5 = $5 + 4}} print $0}}'| bedToBam -i stdin -g GENOME > {output.bam}
        """
