import pandas as pd
import numpy as np
from pathlib import Path

ephys_path = Path('../data/human_mouse_ephys_all_0127.csv')
mouse_cells_path = Path("../data/mouse_IVSCC_excitatory_L23_consolidated_0129.csv")
human_cells_path = Path("../data/human_IVSCC_excitatory_L23_consolidated_0131.csv")
mouse_morph_path = Path("../data/All_Mouse_Cells_Lockdown_All_raw_features.csv")
human_morph_path = Path("../data/All_L23_Lockdown_all_raw_features.csv")
feather_path = Path('../data/feather/dataH.feather')

def load_data():
    ephys_df = pd.read_csv(ephys_path, index_col=0).dropna(how='all')
    morph_df = pd.concat([
        pd.read_csv(mouse_morph_path, index_col=0),
        pd.read_csv(human_morph_path, index_col=0),
    ])

    human_df = (pd.read_csv(human_cells_path, index_col=0)
                .join(ephys_df).join(morph_df)
                .assign(species='human')
                .pipe(fix_df)
            )

    mouse_df = (pd.read_csv(mouse_cells_path, index_col=0)
                .join(ephys_df).join(morph_df)
                .assign(species='mouse')
                .pipe(fix_df)
            )
    drop_neg = [
        'input_resistance',
        'fi_fit_slope',
        'sag'
    ]
    for feat in drop_neg:
        mouse_df[feat] = mouse_df[feat].apply(lambda x: x if (x>0) else np.nan)
        human_df[feat] = human_df[feat].apply(lambda x: x if (x>0) else np.nan)
    return human_df, mouse_df, ephys_df, morph_df

depth = "L23_cell_depth"
cluster = "SeuratMapping"
types_ordered = [ 'LTK', 'GLP2R', 'FREM3', 'CARM1P1', 'COL22A1',
      'Adamts2', 'Rrad', 'Agmat', ]

from pandas.api.types import CategoricalDtype
ttype_categorical = CategoricalDtype(categories=types_ordered, ordered=True)

def fix_df(df):
    df = df.assign(cluster=lambda df: df[cluster]
                        .apply(lambda name: name.split(' ')[-1])
                        .astype(ttype_categorical),
                    depth=lambda df: df[depth]).sample(frac=1, random_state=42)
    df.cluster.cat.remove_unused_categories(inplace=True)
    return df

def join_gene_data(df, genes):
    join_on = 'sample_id'
    genes_df = (pd.read_feather(feather_path, columns=genes+[join_on])
                .set_index(join_on)
                .apply(lambda x: np.log2(x+1))
               )
    return df.join(genes_df, on='sample_id')