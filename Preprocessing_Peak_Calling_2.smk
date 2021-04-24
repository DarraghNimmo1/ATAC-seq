####################################################
## Darragh Nimmo
## Trinity College Dublin
## April 2021
# Snakemake workflow for ATAC-seq preprocessing and peak calling for single sample with replicates.
####################################################

###############################################################################################
###IMPORTS and variables
###############################################################################################

import os

import functools

configfile: 'config.yaml'

IN_DIR = config['in_dir']

GENOME = config['genome']

  
directory_function = functools.partial(os.path.join, config['results'])
BAM_DIR = directory_function('Bam')

###############################################################################################
###Rules
###############################################################################################

rule alignment:
        input:
                read1 = os.path.join(IN_DIR, '{sample}_1.fq.gz'),
                read2 = os.path.join(IN_DIR, '{sample}_2.fq.gz')
        output:
                bam = os.path.join(BAM_DIR, '{sample}_mapped.bam'))
        params:
                index = GENOME
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
        message:
            "Coordinate sorting and indexing alignment bam."
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index -@ 48 {output.bam}
                """            
