from multiprocessing import Pool
from itertools import chain
import pandas as pd
import numpy as np
import ateam.data.shiny as shiny
import sys
from scipy.stats import pearsonr


ephys_path = '../data/human_mouse_ephys_all_0127.csv'
ephys_df = pd.read_csv(ephys_path, index_col=0)

md_path = "../human_IVSCC_excitatory_L23_consolidated_0131.csv"
md_df = pd.read_csv(md_path, index_col=0)
human_df = (md_df.join(ephys_df)
            .loc[lambda df: df["SeuratMapping"].str.contains("FREM")]
           )

genes = pd.read_csv('../data/gene_names.csv', index_col=0).iloc[:,0].values
efeatures = list(ephys_df.columns)
join_on = 'sample_id'

def calc(genes):
    genes_df = (pd.read_feather('../dataH.feather', columns=np.append(genes, [join_on]))
                .set_index(join_on)
                .apply(lambda x: np.log2(x+1))
               )
    data = human_df.join(genes_df, on='sample_id')

    results = []
    efeature = 'depth'
    for gene in genes:
        if data[gene].pipe(lambda x: sum(x>1)) < 5:
            continue
        df = data.dropna(subset=[gene, efeature])
        r, p = pearsonr(df[gene], df[efeature])
        out = {
        "gene":gene,
        "feature":efeature,
        'r_corr':r,
        'rsquared_corr':r**2,
        'p_corr':p,
        }
        results.append(out)
    return results
# n=1
# out = list(map(calc, np.array_split(genes[:50], 10)))

n = sys.argv[1]
step = 2500
nproc = 24
start = int(n)*step
end = min(start+step, len(genes))
pool = Pool()
out = pool.map(calc, np.array_split(genes[start:end], nproc))

metrics = pd.DataFrame.from_records(chain(*out))
metrics.to_csv(f'../data/frem_genes_depth/frem_gene_depthonly_ephys_{n}.csv')
