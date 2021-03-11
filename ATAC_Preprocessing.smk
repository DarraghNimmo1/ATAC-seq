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

SAMPLE = config['Sample']

ATACQC = config['ATACseqQC_script']

directory_function = functools.partial(os.path.join, config['results'])
FASTQC2_DIR = prefix_results('fastqc_2')


###############################################################################################
###Rules
###############################################################################################

rule all:
        input:
                SAMPLE+'_paired.sorted.bam',
                SAMPLE+'_paired.sorted.bam.complexity.txt'

#check read quality using fastqc - which will output a html which can be viewed for a quick check
rule fastqc:
        input:
                read1 = READ_1
                read2 = READ_2
        output:
                os.path.basename(READ_1).rstrip(".fq.gz")+"_fastqc.html",
                os.path.basename(READ_2).rstrip(".fq.gz")+"_fastqc.html",
                os.path.basename(READ_1).rstrip(".fq.gz")+"_fastqc.zip",
                os.path.basename(READ_2).rstrip(".fq.gz")+"_fastqc.zip"
               
        shell:
                """
                fastqc -t 48 {input.read1} {input.read2}| multiqc .
                """
#If fastqc indicates adapter contamination of reads then you need to trim these
rule Trim_Adapters:
        input:
                read1 = READ_1,
                read2 = READ_2
        output:
                read1 = os.path.basename(READ_1).rstrip(".fq.gz")+"_val_1.fq.gz",
                read2 = os.path.basename(READ_2).rstrip(".fq.gz")+"_val_2.fq.gz"
        threads: THREADS
        shell:
                """
                trim_galore  --paired {input}

                """
#align reads to the reference genome
rule Alignment:
        input:
                read1 = rules.Trim_Adapters.output.read1,
                read2 = rules.Trim_Adapters.output.read2
        params:
                index = GENOME
        output:
                bam = SAMPLE+'_mapped.bam'
        shell:
                """
                bowtie2 --very-sensitive -X 1000 -x {params.index} --threads 48 -1 {input.read1} -2 {input.read2}| samtools view -@ 48 -bS - > {output}
                """
#corrdinate sort bam file
rule coordinate_sort1:
        input:
                bam = rules.Alignment.output.bam
        output:
                bam = temp(SAMPLE+'_mapped.sorted.bam'),
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}
                """
#index bam file
rule index1:
        input:
                bam = rules.coordinate_sort1.output.bam
        output:
                index = temp(SAMPLE+'_mapped.sorted.bam.bai')
        shell:
                """
                samtools index -@ 48 {input.bam}
                """
#remove the mitochondrial chromosome - this is common practice for ATAC-seq analysis as the M chromosome is free of nucleosomes so there will be a concentration of reads here
rule remove_chrM:
        input:
                bam = rules.coordinate_sort1.output.bam,
                index = rules.index1.output.index
        output:
                bam = temp(SAMPLE+'_noMT.bam')

        shell:
                """
                samtools view -@ 48 -h -f 3 {input.bam}| python /home/darragh/ATAC-Seq/Python_scripts/removeChrom.py  - - chrM | samtools view -@ 48 -bh - > {output}
                """
#coordinate sort and index the bam file
rule coordinate_sort_index2:
        input:
                bam = rules.remove_chrM.output.bam
        output:
                bam = SAMPLE+'_noMT.sorted.bam',
                index = SAMPLE+'_noMT.sorted.bam.bai'
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index -@ 48 {output.bam}
                """

#only keep paired reads in the bam file
rule paired_reads_only:
        input:
                bam = rules.coordinate_sort_index2.output.bam
        output:
                bam = temp(SAMPLE+'_paired.bam')
        shell:
                """
                samtools view -@ 48 -bh -f 3 {input.bam} > {output.bam}

                """
#coordinate sort bam file
rule Coordinate_sort3:
        input:
                bam = rules.paired_reads_only.output.bam
        output:
                bam = SAMPLE+'_paired.sorted.bam'
        shell:
                """
                samtools sort -@ 48 -n -o {output.bam} {input.bam}; samtools index {output.bam}
                """
#index the bam file
rule index3:
        input:
                bam = rules.Coordinate_sort3.output.bam
        output:
                index = SAMPLE+'_paired.sorted.bam.bai'
        shell:
                """
                samtools index  {input.bam}
                """

#file in size coordinates
rule samtools_fixmate:
        input:
                bam = rules.Coordinate_sort3.output.bam,
                index = rules.index3.output.index
        output:
                bam = temp(SAMPLE+'_fixmate.bam')
        shell:
                """
                samtools fixmate -m {input.bam} {output.bam}
                """
#coordinate sort and index bam file
rule Coordinate_sort_index4:
        input:
                bam = rules.samtools_fixmate.output.bam
        output:
                bam = temp(SAMPLE+'_fixmate.sorted.bam'),
                index = temp(SAMPLE+'_fixmate.sorted.bam.bai)
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index {output.bam}
                """
#remove duplicate reads from the bam file
rule remove_duplicates:
        input:
                bam = rules.Coordinate_sort_index4.bam
        output:
                bam = temp(SAMPLE+'_nodups.bam')
        shell:
                """
                samtools markdup -@ 48 -r -s {input.bam} {output.bam}
                """
#coordinate sort and index the bam file               
rule Coordinate_sort_index5:
        input:
                bam = rules.remove_duplicates.output.bam
        output:
                bam = SAMPLE+'_nodups.sorted.bam',
                index = SAMPLE+'_nodups.sorted.bam.bai
        shell:
                """
                samtools sort -@ 48 -n {output.bam} {input.bam}; samtools index {output.bam}
                """
#convert the bam file back into fastqs
rule convert_bam_fastq:
        input:
                bam = rules.Coordinate_sort_index5.bam
        output:
                read1 = SAMPLE+'_processed_read_1.fq
                read2 = SAMPLE+'_processed_read_2.fq
        shell:
                """
                bedtools bamtofastq -i {input.bam} -fq {output.read1} -fq2 {output.read2}
                """
#run fastqc analysis again after performing all of the quality checks     
rule fastqc_2:
        input:
                read1 = rules.convert_bam_fastq.output.read1
                read2 = rules.convert_bam_fastq.output.read2
        output:
                os.path.join(FASTQC2_DIR, os.path.basename(READ_1).rstrip(".fq.gz")+"_fastqc.html"),
                os.path.join(FASTQC2_DIR, os.path.basename(READ_2).rstrip(".fq.gz")+"_fastqc.html"),
                os.path.join(FASTQC2_DIR, os.path.basename(READ_1).rstrip(".fq.gz")+"_fastqc.zip"),
                os.path.join(FASTQC2_DIR, os.path.basename(READ_2).rstrip(".fq.gz")+"_fastqc.zip")
            
        params:
                outdir=FASTQC2_DIR
               
        shell:
                """
                fastqc -t 48 -o {params.outdir} {input.read1} {input.read2}
                """
               

        
#rule chrom_length_info:
#        input:
#                bam = rules.paired_reads_only.output.bam
#        output:
#                chrom = SAMPLE+'_auto.x.chrom.human.txt'
#        shell:
#                """
#                samtools view -@ 48 -H {input.bam} | grep -P 'SN:' | cut -f2,3  | perl -p -e 's/[SL]N://g' | awk 'BEGIN{FS=OFS="\t"} {print $1, 1, $2, "hg38"}' > {output.chroms}
#                """
#rule ATACseqQC:
#        input:
#                script = ATACQC,
#                bam = rules.Coordinate_sort_index3.output.bam,
#                chroms = rules.chrom_length_info.output.chrom
#        output:
#                libcomplex_plot = SAMPLE+'_paired.sorted.bam.library.complexity.pdf',
#                libcomplex_table = SAMPLE+'_paired.sorted.bam.complexity.txt',
#                fragsize_plot = SAMPLE+'_paired.sorted.fragment.size.distribution.pdf',
#                heatmap_plot = SAMPLE+'_paired.sorted.heatmap and averaged coverage.pdf',
#                cumulatpercent_plot = SAMPLE+'_paired.sorted.Cumulative_Percentage.pdf',
#                feature_align_dist_plot = SAMPLE+'__paired.sorted.Feature_Aligned_Distribution.pdf'#

#        shell:
#                """
#                Rscript {input.script} {input.bam} {input.chroms}#

#                """


              
