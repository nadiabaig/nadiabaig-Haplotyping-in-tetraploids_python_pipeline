#--- standard library imports
#
import os, psutil
import pandas as pd
import time
import gc
import resource
#--- project specific imports
import pysam


with open ('Modified_chr1.vcf') as f1:
    filt_list=[]
    
    for i in f1:
        j=i.strip().split('\t')
        ref=j[2]
        alt=j[3]
        if len(ref)==1 and len(alt)==1:
            filt_list.append(j)
df=pd.DataFrame(filt_list)
df.to_csv('02_Filtered_chr1.vcf',sep="\t",index=False)  

time.sleep(0.1)


bam=pd.read_csv('Chr1_tsv_bam.txt',sep="\t",chunksize=200000,low_memory=False,dtype={"#Read-Name": "category","READ-BASE": "category","CHROM": "category","REF-BASE": "category"})#bam in tsv format--1000000 lines
t = time.process_time()
#do some stuff
dfList = []
for df in bam:
        #df1 = df[df['REF-POS1']!='.']   #some positions showed . (missing info)
        dfList.append(df)
del bam
gc.collect()

df_bam = pd.concat(dfList,sort=False)
elapsed_time = time.process_time() - t
print('dataframe filtered in',elapsed_time)
df_bam=df_bam.rename(columns={"REF-POS1": 'POS'})
df_bam.info(verbose=False, memory_usage="deep")

#--------------------------------------------Reading the VCF file- 1 Chromosome at a time
vcf=pd.read_csv('Modified_chr1.vcf',sep="\t")
print('vcf read')

vcf.columns=['CHROM','REF-POS1','REF','ALT','GT']
vcf2=vcf[['CHROM','REF-POS1','REF','ALT']].copy()
vcf2=vcf2.rename(columns={"REF-POS1": 'POS'})
print('vcf copied')
del vcf
gc.collect()

##type conversions
vcf2.POS=vcf2.POS.astype(int)
df_bam.POS=df_bam.POS.astype(int)

e_t = time.process_time() - t
print('VCF read in secs:',e_t)
vcf2.info(verbose=False, memory_usage="deep")
#-------------------------------------- Convert to categorical data types (if every value is unique, don't bother!)
for df_temp in [vcf2, df_bam]:
    for col in ['CHROM']:
        df_temp.loc[:, col] = df_temp[col].astype('category')
# Merge using less memory
result = pd.merge(vcf2, df_bam, on=["CHROM","POS"], how='left')
del vcf2,df_bam
gc.collect()
result.info(verbose=False, memory_usage="deep")

#----------------------------------------- Combining barcode information into the result dataframe containg info from vcf and bam
samfile = pysam.AlignmentFile("Barcode_filtered.bam", check_sq=False)
barcodes=[]
for read in samfile.fetch(until_eof=True):
    
        t=read.query_name,read.get_tag('BX') ## getting read and barcode name
        barcodes.append(t)
df5 = pd.DataFrame(barcodes)
df5.columns=['#Read-Name','Barcode'] 

del samfile
del barcodes
gc.collect()
## now combine this inforamtion with result
for df_temp1 in [df5, result]:
    for col in ["CHROM","#Read-Name"]:
        df_temp.loc[:, col] = df_temp[col].astype('category')
result1 = pd.merge(df5, result, on=["#Read-Name"], how='left')
##Now removing NA values
df_final1=result1.dropna()
df_final1.POS=df_final1.POS.astype(int)
df_final1.info(verbose=False, memory_usage="deep")
df_final1.to_csv('Barcode_merged_file.txt',sep="\t",index=False)




