import pandas as pd
import pybedtools as pybed
import numpy as np
import pysam
import os
import glob
import  csv


#input file paths Breast - MCF7
gtf = "/home/darragh/h38_protein_coding_genes.minimal.gtf"
gtf2 = "/home/darragh/hg38_minimal.gtf"
EED_peak = pybed.BedTool("/home/darragh/ChIP-seq_2/data/concat/Breast/MDAMB468/EED/Peak/EED_q_0.05_peaks.narrowPeak")
SUZ12_peak = pybed.BedTool("/home/darragh/ChIP-seq_2/data/concat/Breast/MDAMB468/SUZ12/Peak/SUZ12_q_0.05_peaks.narrowPeak")
REST_peak = pybed.BedTool("/home/darragh/ChIP-seq_2/data/concat/Breast/MDAMB468/REST/Peak/REST_q_0.05_peaks.narrowPeak")
H3K27ac_peak = pybed.BedTool("/home/darragh/ChIP-seq_2/data/concat/Breast/MDAMB468/H3K27ac/Peak/H3K27ac_q_0.01_peaks.narrowPeak")
H3K27me3_peak = pybed.BedTool("/home/darragh/ChIP-seq_2/data/concat/Breast/MDAMB468/H3K27me3/Peak/H3K27me3_q_0.01_peaks.broadPeak")
H3K4me1_peak = pybed.BedTool("/home/darragh/ChIP-seq_2/data/concat/Breast/MDAMB468/H3K4me1/Peak/H3K4me1_q_0.01_peaks.broadPeak")
TAC_peak = pybed.BedTool("/home/darragh/ATAC-Seq/Round2/naive_2/HMEC/HMEC_pool_peaks_preprocessed.bed")
Annotation_dir = "/home/darragh/Annotation_reg_elements/MDAMB468"


def Identify_coding_promoters(gtf):
    hg38_gtf = pd.read_table(gtf)
    hg38_gtf.columns = ["seqname","start","end","strand", "gene"]
    seqnames_we_like = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    hg38_gtf = hg38_gtf[hg38_gtf['seqname'].isin(seqnames_we_like)]
    hg38_gtf.loc[hg38_gtf['strand'] == "+", ['start']] = hg38_gtf['start']-2000
    hg38_gtf.loc[hg38_gtf['strand'] == "+", ['end']] = hg38_gtf['start']+5000
    hg38_gtf.loc[hg38_gtf['strand'] == "-", ['start']] = hg38_gtf['end']-3000
    hg38_gtf.loc[hg38_gtf['strand'] == "-", ['end']] = hg38_gtf['start']+5000
    
    return(hg38_gtf)    
  
def Identify_promoters(gtf2):
    hg38_gtf = pd.read_table(gtf2, delimiter="\t",header=None)
    hg38_gtf.columns = ["seqname","start","end","strand", "feature"]
    hg38_gtf = hg38_gtf[hg38_gtf["feature"] == "gene"]
    seqnames_we_like = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    hg38_gtf = hg38_gtf[hg38_gtf['seqname'].isin(seqnames_we_like)]
    hg38_gtf.loc[hg38_gtf['strand'] == "+", ['start']] = hg38_gtf['start']-2000
    hg38_gtf.loc[hg38_gtf['strand'] == "+", ['end']] = hg38_gtf['start']+5000
    hg38_gtf.loc[hg38_gtf['strand'] == "-", ['start']] = hg38_gtf['end']-3000
    hg38_gtf.loc[hg38_gtf['strand'] == "-", ['end']] = hg38_gtf['start']+5000

    return hg38_gtf

  
  def annotate_ATAC_seq_peaks(ATAC_peak, H3K27ac_peak, H3K27me3_peak, H3K4me1_peak, EED_peak, SUZ12_peak, REST_peak, gtf2):
  ATAC_peak_df = pd.read_table(ATAC_peak.fn,header=None)
  ATAC_peak_df = ATAC_peak_df[ATAC_peak_df.columns[:4]]
  ATAC_peak_df.columns = ["chr", "start", "end", "peak"]
  
  #H3K27ac
  H3K27ac = ATAC_peak.intersect(H3K27ac_peak, wa=True, u=True, f=0.2)
  H3K27ac_df = pd.read_table(H3K27ac.fn,header=None)
  H3K27ac_df = H3K27ac_df[H3K27ac_df.columns[:4]]
  H3K27ac_df.columns = ["chr", "start", "end", "peak"]
  H3K27ac_df['H3K27ac'] = "Yes"
  ATAC_peak_df = pd.merge(ATAC_peak_df, H3K27ac_df, how="outer")
  
  #H3K27me3
  H3K27me3 = ATAC_peak.intersect(H3K27me3_peak, wa=True, u=True, f=0.2)
  H3K27me3_df = pd.read_table(H3K27me3.fn,header=None)
  H3K27me3_df = H3K27me3_df[H3K27me3_df.columns[:4]]
  H3K27me3_df.columns = ["chr", "start", "end", "peak"]
  H3K27me3_df['H3K27me3'] = "Yes"
  ATAC_peak_df = pd.merge(ATAC_peak_df, H3K27me3_df, how="outer")
  
  #H3K4me1
  H3K4me1 = ATAC_peak.intersect(H3K4me1_peak, wa=True, u=True, f=0.2)
  H3K4me1_df = pd.read_table(H3K4me1.fn,header=None)
  H3K4me1_df = H3K4me1_df[H3K4me1_df.columns[:4]]
  H3K4me1_df.columns = ["chr", "start", "end", "peak"]
  H3K4me1_df['H3K4me1'] = "Yes"
  ATAC_peak_df = pd.merge(ATAC_peak_df, H3K4me1_df, how="outer")
  
  #EED
  EED = ATAC_peak.intersect(EED_peak, wa=True, u=True)
  EED_df = pd.read_table(EED.fn,header=None)
  EED_df = EED_df[EED_df.columns[:4]]
  EED_df.columns = ["chr", "start", "end", "peak"]
  EED_df['EED'] = "Yes"
  ATAC_peak_df = pd.merge(ATAC_peak_df, EED_df, how="outer")
  
  #SUZ12
  SUZ12 = ATAC_peak.intersect(SUZ12_peak, wa=True, u=True)
  SUZ12_df = pd.read_table(SUZ12.fn,header=None)
  SUZ12_df = SUZ12_df[SUZ12_df.columns[:4]]
  SUZ12_df.columns = ["chr", "start", "end", "peak"]
  SUZ12_df['SUZ12'] = "Yes"
  ATAC_peak_df = pd.merge(ATAC_peak_df, SUZ12_df, how="outer")
  
  #REST
  REST = ATAC_peak.intersect(REST_peak, wa=True, u=True)
  REST_df = pd.read_table(REST.fn,header=None)
  REST_df = REST_df[REST_df.columns[:4]]
  REST_df.columns = ["chr", "start", "end", "peak"]
  REST_df['REST'] = "Yes"
  ATAC_peak_df = pd.merge(ATAC_peak_df, REST_df, how="outer")
  
  #Proximal
  Promoter = pybed.BedTool.from_dataframe(Identify_promoters(gtf2))
  Proximal = ATAC_peak.intersect(Promoter, wa=True, u=True, f=0.2)
  Proximal_df = pd.read_table(Proximal.fn,header=None)
  Proximal_df = Proximal_df[Proximal_df.columns[:4]]
  Proximal_df.columns = ["chr", "start", "end", "peak"]
  Proximal_df['Proximal'] = "Yes"
  ATAC_peak_df = pd.merge(ATAC_peak_df, Proximal_df, how="outer")
  
  ATAC_peak_df = ATAC_peak_df.fillna('No')
  
  #Active promoters
  active_promoter = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'Yes') & (ATAC_peak_df['H3K27ac'] == 'Yes') & (ATAC_peak_df['SUZ12'] == 'No') & (ATAC_peak_df['EED'] == 'No') & (ATAC_peak_df['H3K27me3'] == 'No')]
  active_promoter.to_csv(os.path.join(Annotation_dir, "Active_promoter") ,sep = "\t")
  
  #Poised promoters
  poised_promoter_suz12 = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'Yes') & (ATAC_peak_df['H3K27me3'] == 'Yes') &  (ATAC_peak_df['SUZ12'] == 'Yes')]
  poised_promoter_eed = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'Yes') & (ATAC_peak_df['H3K27me3'] == 'Yes') &  (ATAC_peak_df['EED'] == 'Yes')]
  poised_promoter = pd.concat([poised_promoter_suz12, poised_promoter_eed])
  poised_promoter = poised_promoter.drop_duplicates()
  poised_promoter.to_csv(os.path.join(Annotation_dir, "Poised_promoter") ,sep = "\t")

  #Silencer proximal
  rest_proximal = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'Yes') & (ATAC_peak_df['H3K27ac'] == 'No') &  (ATAC_peak_df['EED'] == 'No') & (ATAC_peak_df['SUZ12'] == 'No') & (ATAC_peak_df['REST'] == 'Yes')]
  rest_proximal.to_csv(os.path.join(Annotation_dir, "Rest_silencer_proximal") ,sep = "\t")
  
  #Uncharacterised proximal
  Proximal = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'Yes')]
  not_active_promoter = pd.merge(Proximal,active_promoter, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
  Uncharacterised_proximal = pd.merge(not_active_promoter,poised_promoter, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
  Uncharacterised_proximal.to_csv(os.path.join(Annotation_dir, "Uncharacterised_proximal") ,sep = "\t")
  
  #Active enhancer
  active_enhancer = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'No') & (ATAC_peak_df['H3K27ac'] == 'Yes') & (ATAC_peak_df['H3K4me1'] == 'Yes') & (ATAC_peak_df['SUZ12'] == 'No') & (ATAC_peak_df['EED'] == 'No') & (ATAC_peak_df['H3K27me3'] == 'No')]
  active_enhancer.to_csv(os.path.join(Annotation_dir, "Active_enhancer") ,sep = "\t")
  
  #Poised enhancer
  poised_enhancer_suz12 = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'No') & (ATAC_peak_df['H3K27me3'] == 'Yes') &  (ATAC_peak_df['SUZ12'] == 'Yes') & (ATAC_peak_df['H3K4me1'] == 'Yes')]
  poised_enhancer_eed = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'No') & (ATAC_peak_df['H3K27me3'] == 'Yes') &  (ATAC_peak_df['EED'] == 'Yes') & (ATAC_peak_df['H3K4me1'] == 'Yes')]
  poised_enhancer = pd.concat([poised_enhancer_suz12, poised_enhancer_eed])
  poised_enhancer = poised_enhancer.drop_duplicates()
  poised_enhancer.to_csv(os.path.join(Annotation_dir, "Poised_enhancer") ,sep = "\t")

  #Uncharacterised distal
  Distal = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'No')]
  not_active_enhancer = pd.merge(Distal,active_enhancer, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
  Uncharacterised_distal = pd.merge(not_active_enhancer,poised_enhancer, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
  Uncharacterised_distal.to_csv(os.path.join(Annotation_dir, "Uncharacterised_distal") ,sep = "\t")

  #Silencer distal
  rest_distal = ATAC_peak_df.loc[(ATAC_peak_df['Proximal'] == 'No') & (ATAC_peak_df['H3K27ac'] == 'No') &  (ATAC_peak_df['EED'] == 'No') & (ATAC_peak_df['SUZ12'] == 'No') & (ATAC_peak_df['REST'] == 'Yes')]
  rest_distal.to_csv(os.path.join(Annotation_dir, "Rest_silencer_distal") ,sep = "\t")

  return(active_promoter,poised_promoter,Uncharacterised_proximal,active_enhancer,poised_enhancer,Uncharacterised_distal)
