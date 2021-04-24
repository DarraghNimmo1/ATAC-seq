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

chrM_SCRIPT = config['chrM_script']

bedpe_SCRIPT = config['bedpe_script']

BLACK_LIST = config['black_list']

  
directory_function = functools.partial(os.path.join, config['results'])
BAM_DIR = directory_function('Bam')
BED_DIR = directory_function('Bed')
PEAK_DIR = directory_function('Peak')

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
                bam = temp(os.path.join(BAM_DIR, '{sample}_mapped.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, '{sample}_mapped.sorted.bam.bai'))
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
                bam = os.path.join(BAM_DIR, '{sample}_noMT.bam'))
	      params:
		            script = chrM_SCRIPT
        message:
                "Removing mitochondrial alignments from the bam file."
        shell:
                """
                bash {params.script} {input.bam} {output.bam}
                """
            
            
rule coordinate_sort_index_2:
        input:
                bam = rules.remove_chrM.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, '{sample}_noMT.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, '{sample}_noMT.sorted.bam.bai'))
        message:
                "Coordinate sorting and indexing -chrM bamfile"
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam};
                 samtools index -@ 48 {output.bam}
                """  

rule paired_reads_only:
        input:
                bam = rules.coordinate_sort_index_2.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, '{sample}_paired.bam'))
        message:
            "Keeping just paired reads in the bam file."
        shell:
                """
                samtools view -@ 48 -bh -f 3 {input.bam} > {output.bam}
                """
            
rule name_sort:
        input:
                bam = rules.paired_reads_only.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, '{sample}_namesorted.bam'))
        message:
                "Name sorting bam file"
        shell:
                """
                samtools sort -n -@ 48 {input.bam} -o {output.bam}
                """
                             
rule fixmate:
       input:
               bam = rules.name_sort.output.bam
       output:
               bam = temp(os.path.join(BAM_DIR, '{sample}_fixmate.bam'))
       message:
              "Performing samtools fixmate"
       shell:
              """
              samtools fixmate -@48 -m {input.bam} {output.bam}
              """
 
rule coordinate_sort_index_3:
        input:
                bam = rules.fixmate.output.bam
        output:
                bam = temp(os.path.join(BAM_DIR, '{sample}_fixmate.sorted.bam')),
                index = temp(os.path.join(BAM_DIR, '{sample}_fixmate.sorted.bam.bai'))
        message:
                "Coordinate sorting and indexing fixmate bamfile"
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam};
                 samtools index -@ 48 {output.bam}
                """  
                      
rule mark_duplicates:
        input:
              bam = rules.coordinate_sort_index_3.output.bam
        output:
              bam = temp(os.path.join(BAM_DIR, '{sample}_md.bam'))
        message:
              "Removing duplicate fragments from the bam file"
        shell:
              """
              samtools markdup -r {input.bam} {output.bam} -@ 48
              """
                                                                
rule coordinate_sort_index_4:
        input:
                bam = rules.mark_duplicates.output.bam
        output:
                bam = os.path.join(BAM_DIR, '{sample}_md.sorted.bam')),
                index = os.path.join(BAM_DIR, '{sample}_md.sorted.bam.bai'))
        message:
                "Coordinate sorting and indexing md bamfile"
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam};
                 samtools index -@ 48 {output.bam}
                """                                        
                                        
rule bam_to_bed_and_Tn5_shift:
    input:
        bam = rules.coordinate_sort_4.output.bam
    output:
        bedpe = os.path.join(BED_DIR, '{sample}_Tn5_shift.bedpe')
    message:
        "Converting bam file to a bed file and performing Tn5 shift."
    shell:
        """
        bamToBed -i {input.bam} -bedpe | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($9 == "+") {$2 = $2 + 4; $6 = $6 - 5} 
        else if ($9 == "-") {$3 = $3 - 5; $5 = $5 + 4} print $0}' > {output.bedpe}
        """             
                                        
rule convert_bedpe_bam:
     input:
        bedpe = rules.bam_to_bed_and_Tn5_shift.output.bedpe
     output:
        bam = os.path.join(BAM_DIR, '{sample}_Tn5_shift.bam')
     message:
        "Converting Tn5 bedpe back to bam"
     shell:
        """
        bedpeToBam -i {input.bedpe} -g {params.genome} > {output.bam}
        """
        
rule convert_bedpe_macs:
	input:
		bedpe = rules.bam_to_bed_and_Tn5_shift.output.bedpe
	output:
		bedpe = os.path.join(BED_DIR, '{sample}_macs_input.bedpe'))
	params:
		script = bedpe_SCRIPT
	message:
        	"Converting bedpe to MACS2 compatible format."
	shell:
		"""
		bash {params.script} {input.bedpe} > {output.bedpe}
		""" 
          
rule peak_calling:
    input:
        bam = rules.convert_bedpe_macs.output.bam
    output:
        BroadPeak = os.path.join(PEAK_DIR, '{sample}_peaks.broadPeak')
    params:
        dir = PEAK_DIR
    message:
        "Peak calling."
    shell:
        """
        macs2 callpeak -t {input.bam} --nomodel --shift -100 --extsize 200 --broad -f BED --name {params.sample} -g hs --outdir {params.dir}
        """

rule preprocess_peaks:
	input:
		BroadPeak = rules.peak_calling.output.BroadPeak,
	output:
		bed = os.path.join(PEAK_DIR, SAMPLE+'_processed_peaks.bed')
	params:
		BlackList = BLACK_LIST
	message:
		"Preprocessing peaks."
	shell:
		"""
		bedtools subtract -a {input.BroadPeak} -b {params.BlackList} > {output.bed}
		"""
                                        
