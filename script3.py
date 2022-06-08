#--- standard library imports
#
import os, psutil
import pandas as pd
import time
import gc
import resource
#--- project specific imports
import pysam
pd.set_option("display.max_rows", None, "display.max_columns", None)

b_s=pd.read_csv('Base_support_at_snp_position_bam.txt',sep=",")
b_s.columns=['POS','ALT_bam','Count']
v_s=pd.read_csv('Base_support_at_snp_position_vcf.txt',sep=",")
v_s.columns=['POS','ALT_vcf','Count']
result=v_s.merge(b_s, how='left', on="POS")

#now removing the bases if the highest base support variant isn't equal to alt vcf
##keep only first vaLUE in bams and if it is equal to results alt_vcf keep else remove
# Using DataFrame.drop_duplicates() to keep first duplicate row
df2 = b_s.drop_duplicates(subset=["POS"],keep='first')
df2.columns=['POS','ALT_bam1','Count']
#if df2 in result ok, if not remove
result2=df2.merge(result, how='left', on="POS")
print(result2.head(10))

df2.columns=['POS','ALT_bam1','Count']
rfinal=result2[result2['ALT_bam1']==result2['ALT_vcf']]
##now getting the % of base support
rfinal2=rfinal[['POS','Count_x','ALT_vcf','ALT_bam','Count_y']].copy()
rfinal2.columns=['POS','Count_vcf','ALT_vcf','ALT_bam','Count_bam']

##filtering reads:base support should be atleast 0.15x of total read support
rfinal2['0.15x_alt_base'] = (15/100)*(rfinal2['Count_vcf'])
rfinal2=rfinal2[rfinal2['Count_bam']>=rfinal2['0.15x_alt_base']]
print(rfinal2.head(20))
#now reading barcoded file
df_bar=pd.read_csv('Barcode_merged_file.txt',sep="\t",iterator=True, chunksize=20000)
df_bar2= pd.concat(df_bar,ignore_index=True)
#del df_bar, b_s, v_s ,result_base
gc.collect()
##merging barcode and base support file, for that we need to change column names
df_bar2.rename(columns={'READ-BASE': 'ALT_bam'}, inplace=True)
df_f=df_bar2.merge(rfinal2, on=['POS','ALT_bam'], how='left',copy=False) #merging the two files based on a common position and Chromosome
df_f0 = df_f.dropna()
del df_bar2, df_f
gc.collect()

##Now sort the data frame on position
df_f00=df_f0.sort_values('POS')
del df_f0

##copying specific columns
t=df_f00[['CHROM','POS','Barcode','#Read-Name','REF','ALT_vcf','REF-BASE','ALT_bam']].copy()
t.to_csv('Final_file.txt',index=False,sep="\t")
del df_f00
gc.collect()
