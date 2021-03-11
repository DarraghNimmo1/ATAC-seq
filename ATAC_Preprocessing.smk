###############################################################################################
###IMPORTS and variables
###############################################################################################
import os

configfile: 'config.yaml'

READ_1 = config['read_1']

READ_2 = config['read_2']

THREADS = config['threads']

GENOME = config['genome']

SAMPLE = config['Sample']

ATACQC = config['ATACseqQC_script']

###############################################################################################
###Functions
###############################################################################################


###############################################################################################
###Rules
###############################################################################################

rule all:
        input:
                SAMPLE+'_paired.sorted.bam',
                SAMPLE+'_paired.sorted.bam.complexity.txt'


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

rule coordinate_sort1:
        input:
                bam = rules.Alignment.output.bam
        output:
                bam = temp(SAMPLE+'_mapped.sorted.bam'),
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}
                """

rule index1:
        input:
                bam = rules.coordinate_sort1.output.bam
        output:
                index = temp(SAMPLE+'_mapped.sorted.bam.bai')
        shell:
                """
                samtools index -@ 48 {input.bam}
                """
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

rule paired_reads_only:
        input:
                bam = rules.coordinate_sort_index2.output.bam
        output:
                bam = temp(SAMPLE+'_paired.bam')
        shell:
                """
                samtools view -@ 48 -bh -f 3 {input.bam} > {output.bam}

                """
rule Coordinate_sort_index3:
        input:
                bam = rules.paired_reads_only.output.bam
        output:
                bam = SAMPLE+'_paired.sorted.bam',
                index = SAMPLE+'_paired.sorted.bam.bai'
        shell:
                """
                samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index {output.bam}
                """

rule samtools_fixmate:
        input:
                bam = rules.Coordinate_sort_index3.ouput.bam
        output:
                bam = SAMPLE+'_fixmate.bam'
        shell:
                """
                samtools fixmate -m {input.bam} {output.bam}
                """

rule Coordinate_sort_index4:
        input:
                bam = rules.Coordinate_sort_index4.bam
        output:
                bam = SAMPLE+'_fixmate.sorted.bam'

rule chrom_length_info:
        input:
                bam = rules.paired_reads_only.output.bam
        output:
                chrom = SAMPLE+'_auto.x.chrom.human.txt'
        shell:
                """
                samtools view -@ 48 -H {input.bam} | grep -P 'SN:' | cut -f2,3  | perl -p -e 's/[SL]N://g' | awk 'BEGIN{FS=OFS="\t"} {print $1, 1, $2, "hg38"}' > {output.chroms}
                """
rule ATACseqQC:
        input:
                script = ATACQC,
                bam = rules.Coordinate_sort_index3.output.bam,
                chroms = rules.chrom_length_info.output.chrom
        output:
                libcomplex_plot = SAMPLE+'_paired.sorted.bam.library.complexity.pdf',
                libcomplex_table = SAMPLE+'_paired.sorted.bam.complexity.txt',
                fragsize_plot = SAMPLE+'_paired.sorted.fragment.size.distribution.pdf',
                heatmap_plot = SAMPLE+'_paired.sorted.heatmap and averaged coverage.pdf',
                cumulatpercent_plot = SAMPLE+'_paired.sorted.Cumulative_Percentage.pdf',
                feature_align_dist_plot = SAMPLE+'__paired.sorted.Feature_Aligned_Distribution.pdf'

        shell:
                """
                Rscript {input.script} {input.bam} {input.chroms}

                """


              
