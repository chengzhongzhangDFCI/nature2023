import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pkl
import os
chrom_list = ['chr'+str(i) for i in range(1,23)]
chrom_list.append('chrX')

par_samples = ['par_'+str(i) for i in range(1,11)]
norm_counts = pd.read_csv('peakmatrixNORM220815.txt', sep=' ')
#top = norm_counts[par_samples].median(axis=1).sort_values()[0:75000]

cn_fragments = pd.read_csv('ATAC_frag_counts_CN.txt', sep=' ', index_col=0).dropna()
cn_fragments = cn_fragments[cn_fragments.min(axis=1)>0]
print(cn_fragments.shape)
#istop = cn_fragments.copy().index.isin(top.index)
#cn_fragments = cn_fragments[istop]

print(cn_fragments.shape)

peaks = cn_fragments.index
cn_fragments = cn_fragments.reset_index(drop=True)
parse_ranges = lambda x: [i.split('-')[x] for i in peaks]
chroms = pd.Series(parse_ranges(0))
starts = pd.Series(parse_ranges(1)).astype('int')


plot_samples = ['par_'+str(i) for i in range(1,11)]
plot_samples.extend(['F9_'+str(i) for i in range(1,11)])
plot_samples.extend(['F11_'+str(i) for i in range(1,11)])
plot_samples.extend(['parent', 'B4', 'E5', 'E8','F2', 'F3', 'F9', 'K11', 'N6', 'N11', 'O16', 'F11', 'K2'])
for sample_id in plot_samples:
    pklfile = './PermTest/'+sample_id+'_fc_background_withGC.pkl'
    fc_background = pkl.load(open(pklfile, 'rb'))
    #fc_background = fc_background[istop]
    fc_background = fc_background.reset_index(drop=True)
    fc_obs = cn_fragments.copy()[sample_id]/cn_fragments.copy()['parent']
    scaled = fc_obs/fc_background.mean(axis=1)
    scaled.index = [str(i)+'-'+str(j) for i,j in zip(chroms, starts)]
    #scaled.to_csv(sample_id+'_perm-scaled_top75000.txt', sep='\t', header=None)
