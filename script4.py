#--- standard library imports
#
import os, psutil
import pandas as pd
import time
import gc
import resource
import csv
from ast import literal_eval
#--- project specific imports


f1=pd.read_csv('Final_file.txt',sep="\t",usecols=['CHROM','POS','Barcode','ALT_bam'])
f1["Pos_alt"] = f1["POS"].astype(str)+"_"+f1["ALT_bam"]
f2=f1[['CHROM','Pos_alt','Barcode','ALT_bam']].copy()
del f1
gc.collect()
##grouping
grouped_df = f2.groupby("Pos_alt")
grouped_lists = grouped_df["Barcode"].apply(list)  #groping on same pos and alt then on same barcode
grouped_lists = grouped_lists.reset_index()
grouped_lists[['POS', 'ALT_bam']] = grouped_lists['Pos_alt'].str.split('_', expand=True)
del grouped_df
gc.collect()
grouped_lists.to_csv('Final.txt',sep="\t",index=False)
#-------------------reading file
df1=pd.read_csv('Final.txt',sep='\t',low_memory=False,quoting=csv.QUOTE_NONE)
df1['Barcode'] = df1['Barcode'].apply(literal_eval) #convert to list type
t=df1.explode('Barcode')
g_allel = t.groupby("Barcode")
g_lists = g_allel["Pos_alt"].apply(list) #I see same barcode repeating twice in some cases for one of the allele of snp position-- can't figure it out why but I removed those positions and kept 1
g_lists = g_lists.reset_index()
g_lists['Haps'] = [','.join(map(str, l)) for l in g_lists['Pos_alt']]
del df1,t,g_allel
gc.collect()
def Haplotype(Haps):
    return ','.join(set(Haps.split(',')))

g_lists['Haplotype'] = g_lists.Haps.apply(Haplotype)
g_lists["Haplotype"] = g_lists["Haplotype"] .str.split(",").apply(lambda x: sorted(x))

flist=g_lists[['Barcode','Haplotype']].copy()
del g_lists
gc.collect()
flist.to_csv('Final_haps.txt',sep="\t",index=False)
