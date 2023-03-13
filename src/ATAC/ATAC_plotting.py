import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm
import matplotlib.patches as mpatches
import seaborn as sns
chrom_list = ['chr'+str(i) for i in range(1,23)]
chrom_list.append('chrX')
bridge_clones = ['B4', 'K2', 'F2', 'K11', 'F9', 'F3', 'E8', 'F11', 'N11', 'E5', 'O16','N6' ]
par_list = ['par_'+str(i) for i in range(1,11)]

def getFC_dict(sample_list, window_size, cutoff):
    fc_dict = {}
    m_dict = {}
    for sample in sample_list:
        sample_means = []
        sample_intervals = []
        sample_fc = pd.read_csv('../../data/ATAC/FC_tables/'+sample+'_FC_table_1Mb.txt', sep='\t', index_col=0)
        sample_fc['Perm-scaled_FC'] = sample_fc['Perm-scaled_FC']/sample_fc['Perm-scaled_FC'].median()
        peaks = sample_fc.index
        sample_fc = sample_fc.reset_index(drop=True)
        chroms = pd.Series([i.split('-')[0] for i in peaks])
        starts = pd.Series([i.split('-')[1] for i in peaks]).astype('int')
        
        for chrom in chrom_list:
            counter = 0
            chrom_fc = sample_fc.copy()[chroms==chrom]
            chrom_starts = starts[chroms==chrom]
            
            while counter < window_size:
                chrom_fc = chrom_fc[chrom_starts>counter]
                chrom_starts = chrom_starts[chrom_starts>counter]
                chrom_fc['z_signal'] = chrom_fc['Perm-scaled_FC']*chrom_fc['Peak_count']
                
                chrom_bins = pd.Series(range(counter, max(chrom_starts), window_size))
                bin_ids = np.digitize(chrom_starts, chrom_bins)
                in_bins = pd.Series([chrom_bins[i-1] for i in np.unique(bin_ids)])

                bin_counts = chrom_fc.groupby(bin_ids).sum()['Peak_count']
                chrom_sample_means = chrom_fc.groupby(bin_ids).mean()['z_signal']/chrom_fc.groupby(bin_ids).mean()['Peak_count']
                
                chrom_sample_means = chrom_sample_means[bin_counts>=cutoff]
                chrom_intervals = pd.Series([chrom+'-'+str(i) for i in in_bins])[(bin_counts>=cutoff).values]
                
                sample_means.extend(chrom_sample_means)
                sample_intervals.extend(chrom_intervals)
                
                if chrom=='chr4' and counter==int(8e6):
                    m_dict[sample] = chrom_sample_means.values[np.where(chrom_intervals=='chr4-28000000')[0]]

                counter+=int(1e6)

            
            fc_dict[sample] = sample_means
        
    return fc_dict, m_dict, sample_intervals

def plotSignal(window):
    plt.rcParams["figure.figsize"] = (4,3)
    intervals_all = {}
    fc_all = {}
    fc_all_par = {}
    
    fc_dict_bridge, m_bridge, intervals = getFC_dict(bridge_clones, window, window*1e-5)
    fc_dict_par, m_par, intervals = getFC_dict(par_list, window, window*1e-5)
    plt.rcParams["figure.figsize"] = (8,4)
    bridge_df = pd.DataFrame(fc_dict_bridge)
    par_df = pd.DataFrame(fc_dict_par)

    chroms = pd.Series([i.split('-')[0] for i in intervals])
    starts = pd.Series([i.split('-')[1] for i in intervals]).astype('int')

    df = pd.DataFrame()
    for sample in bridge_df.columns:
        sample_df = pd.DataFrame()
        sample_df['FC'] = bridge_df[sample]
        sample_df['sample'] = sample
        df = df.append(sample_df)
    for sample in par_df.columns:
        sample_df = pd.DataFrame()
        sample_df['FC'] = par_df[sample]
        sample_df['sample'] = sample
        df = df.append(sample_df)
    df = df.reset_index(drop=True)
    fig, ax = plt.subplots()
    sns.boxplot(data=df, x=df['sample'], y=df['FC'],
                   ax=ax, color='skyblue')
    ax.set_ylabel('Normalized ATAC signal', fontsize=14)
    ax.set_ylim([0, 5])

    interval_points_bridge = bridge_df.copy()[((chroms=='chr4') & (starts>=int(27e6)) & (starts<=int(37e6)))]
    interval_points_bridge = interval_points_bridge.reset_index(drop=True)

    interval_points_par = par_df.copy()[((chroms=='chr4') & (starts>=int(27e6)) & (starts<=int(37e6)))]
    interval_points_par = interval_points_par.reset_index(drop=True)

    for i in range(bridge_df.shape[1]):
        for j in range(interval_points_bridge.shape[0]):
            ax.scatter(i, interval_points_bridge[bridge_clones[i]][j], c='maroon', s=3)

    for i in range(par_df.shape[1]):
        for j in range(interval_points_par.shape[0]):
            ax.scatter(i+12, interval_points_par[par_list[i]][j], c='maroon', s=3)


    #plt.savefig('1Mb.pdf')
    plt.show()