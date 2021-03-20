
####################################################
## Darragh Nimmo
## Trinity College Dublin
## March 2021
# Snakemake workflow for ATAC-seq footprinting.
####################################################

###############################################################################################
###IMPORTS and variables
###############################################################################################

configfile: 'config.yaml'

BAMS = config['bams']


directory_function = functools.partial(os.path.join, config['results'])


rule bam_to_bed_and_Tn5_shift:
    input:
        bam = BAMS
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
