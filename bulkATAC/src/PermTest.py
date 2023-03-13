import numpy as np
import pandas as pd
import pickle as pkl
import sys
cn_fragments = pd.read_csv('/singlecellcenter/greg/ATAC_22-2-24/ATAC_frag_counts_CN.txt', sep=' ', index_col=0).dropna()
cn_fragments = cn_fragments[cn_fragments.min(axis=1)>0]
peaks = cn_fragments.index
background_file = '/singlecellcenter/greg/ATAC_22-2-24/bg_gc.dat'
bg = pd.read_csv(background_file, sep='\t').astype('int')
bg = bg.copy()-1
bg = bg.reset_index(drop=True)
cn_fragments = cn_fragments.reset_index(drop=True)
sample_id = sys.argv[1]
sample_index = list(cn_fragments.columns).index(sample_id)
nperms = 10000
bg_mat = bg.to_numpy()
sample = np.random.choice(50, [bg.shape[0], nperms])
rows = np.tile([range(bg.shape[0])], (nperms, 1)).T
sample_peaks = bg_mat[[rows], [sample]]
perm_peaks = sample_peaks[0]
frag_mat = cn_fragments.to_numpy()
sample_cols = np.full((frag_mat.shape[0], nperms), sample_index).astype('int')
par_cols = np.full((frag_mat.shape[0], nperms), 11).astype('int')
sample_perms = frag_mat[[perm_peaks],[sample_cols]]
par_perms = frag_mat[[perm_peaks],[par_cols]]
fc = pd.DataFrame(sample_perms[0].copy()/par_perms[0].copy())
fc_obs = cn_fragments.copy()[sample_id]/cn_fragments.copy()['parent']
pkl.dump(fc, open('/singlecellcenter/greg/ATAC_22-2-24/PermTest/'+sample_id+'_fc_background_withGC.pkl', 'wb'), protocol=4)