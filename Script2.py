#--- standard library imports
#
import os, psutil
import pandas as pd
import time
import gc
import resource
#--- project specific imports
import pysam


#---------------------##Getting base and read support inforamtion --step3
##read the barcode merged file 
pd.set_option("display.max_rows", None, "display.max_columns", None)
bar=pd.read_csv('Barcode_merged_file.txt',sep="\t",iterator=True, chunksize=20000)
barcoded_file = pd.concat(bar,ignore_index=True)
del bar
gc.collect()
##sorting the entries on position
#barcoded_file.POS = pd.to_numeric(barcoded_file.POS, errors='coerce')
barcoded_file1=barcoded_file.sort_values('POS')
#del barcoded_file
gc.collect()
alt_bam=barcoded_file1.groupby(["POS"])['READ-BASE'].value_counts()  #count of base support in alt bam
alt_vcf=barcoded_file1.groupby(["POS"])['ALT'].value_counts()  ##count of base support in alt vcf
a3=alt_vcf.to_frame().to_csv('Base_support_at_snp_position_vcf.txt')
a2=alt_bam.to_frame().to_csv('Base_support_at_snp_position_bam.txt')
