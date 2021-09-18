import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from itertools import combinations
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.anova import anova_lm

def ols_model(data, formula_rhs, feature, anova_type=2, cov_type='HC3'):
    metrics = ['aic', 'bic', 'fvalue', 'f_pvalue', 'llf', 'rsquared', 'rsquared_adj', 'nobs']
    formula = f"{feature} ~ {formula_rhs}"
    res = smf.ols(formula=formula, data=data).fit(cov_type=cov_type)
    fit_dict = {name: getattr(res, name) for name in metrics}

    if anova_type:
        anova = anova_lm(res, typ=anova_type)
        pvals = anova["PR(>F)"].dropna().rename(lambda x: f"pval_{x}")
        fvals = anova["F"].dropna().rename(lambda x: f"fval_{x}")
        eta = anova["sum_sq"].dropna().apply(lambda x: x/(x+res.ssr)).rename(lambda x: f"eta_p_{x}")
        fit_dict.update(pvals.to_dict())
        fit_dict.update(fvals.to_dict())
        fit_dict.update(eta.to_dict())
    return fit_dict, res
    
def fit_models(data, formulas, features, formula_names=None, feature_names=None, cov_type='HC3'):
    feature_names = feature_names or [str(x) for x in features]
    formula_names = formula_names or [str(x) for x in formulas]
    all_fits = []
    for feature, feature_name in zip(features, feature_names):
        for formula, formula_name in zip(formulas, formula_names):
            fit_dict, results = ols_model(data, formula, feature, cov_type=cov_type)
            all_fits.append(dict(fit_dict, model=formula_name, feature=feature_name))
        
    fits_df = pd.DataFrame(all_fits)
    return fits_df

def plot_fit(data, feature, formula, x=None, cluster='cluster', ax=None, legend=False,
            print_attr=None, print_pvals=True, **sns_args):
    if not ax:
        fig, ax = plt.subplots()
    # factors = [c for c in data.columns if c in formula]
    # data = data.dropna(subset=factors+[feature])
    out_dict, res = ols_model(data, formula, feature)

    x = x or formula.replace('*','+').split('+')[0].strip()
    sns_args['s'] = sns_args.get('s', 25)
    # if cluster is not None and cluster not in formula:
    #     data[cluster] = data[cluster].fillna('none')
    sns.scatterplot(data=data, y=feature, x=x, hue=cluster, ax=ax, legend=legend, **sns_args)
    
    if cluster is not None and cluster in formula:
        hue = cluster 
        c = None
    else:
        hue = None
        c = 'k'
    y_fit = res.fittedvalues.reindex(data.index)
    sns.lineplot(data=data, y=y_fit, x=x, hue=hue, color=c, legend=False, ax=ax)
    ax.set_xlabel(getattr(x, "label", None) or x)
    ax.set_ylabel(getattr(feature, "label", None) or feature)

    summary = ''
    if print_attr:
        if not isinstance(print_attr, list):
            print_attr = [print_attr]
        for attr in print_attr:
            value = out_dict.get(attr)
            attr_name = getattr(attr, "label", None) or attr
            summary += f"{attr_name} = {value:.2g}\n"
    if print_pvals:
        anova = anova_lm(res, typ=2)
        pvals = anova["PR(>F)"].dropna()
        summary += ", ".join(f"p_{key}={pvals[key]:.2g}" for key in pvals.index)
    ax.text(0.5, 0.99, summary, transform=ax.transAxes,
        verticalalignment='top', horizontalalignment='center')
    sns.despine()
    return out_dict

def plot_mw_bars(data, var, group, group_vals=None, pairs='all', cutoff=0.05, label='stars', ax=None, y0=None):
    group_vals = group_vals or data[group].unique().tolist()
    pairs_list, pairs_idx, pvals = pairwise_mw(data, var, group, group_vals, pairs)
    pvals = multipletests(pvals, method='fdr_bh')[1]
    y0 = data[ data[group].isin(set.union(*map(set,pairs_list))) ][var].max()
    plot_sig_bars(pvals, pairs_idx, cutoff, label=label, ax=ax, y0=y0)

def plot_sig_bars(pvals, pairs_idx, cutoff=0.05, label='stars', ax=None, y0=None):
    ax = ax or plt.gca()
    ylim = ax.get_ylim()
    pairs_sig = np.flatnonzero(np.array(pvals)<cutoff)
    
    y0 = y0 or ylim[1]
    n = len(pairs_sig)
    dy = 0.04*(ylim[1]-ylim[0]) # use 4% of y-axis range
    yvals = y0 + dy*np.arange(1, n+1)
    for i, pair in enumerate(pairs_sig):
        plot_sig_bar(pvals[pair], yvals[i], pairs_idx[pair], label=label, ax=ax)

def plot_sig_bar(pval, y, x_pair, label='stars', ax=None):
    ax = ax or plt.gca()
    ax.plot(x_pair, [y, y], 'grey')
    if label=='stars':
        text = np.choose(np.searchsorted([1e-3, 1e-2, 5e-2], pval), ['***','**','*',''])
    elif label=='pval':
        text = "p={p:.2}".format(p=pval)
    else:
        text = ''
    ax.text(np.mean(x_pair), y, text, horizontalalignment='center', verticalalignment='bottom')

def pairwise_mw(data, var, group, group_vals=None, pairs='all'):
    data = data[~data[var].isna()]
    group_vals = group_vals or data[group].unique().tolist()
    if pairs is 'all':
        pairs = list(combinations(group_vals, 2))
    pvals = []
    pairs_idx = []
    groups = data.groupby(group)[var]
    for pair in pairs:
        u, p = mannwhitneyu(groups.get_group(pair[0]), groups.get_group(pair[1]), alternative='two-sided')
        pvals.append(p)
        pairs_idx.append([group_vals.index(pair[0]), group_vals.index(pair[1])])
    return pairs, pairs_idx, pvals

def outline_boxplot(ax):
    for i,artist in enumerate(ax.artists):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        col = artist.get_facecolor()
        artist.set_edgecolor(col)
        artist.set_facecolor('None')

        # Each box has 5(?) associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        n=5
        for j in range(i*n,i*n+n):
            line = ax.lines[j]
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)
            
def plot_boxplot_multiple(data, features, x='cluster', labels=None, horizontal=False, figsize=(4,8),
                            palette=None, strip_width=0.2, pairs_sets=[], cutoff=0.05, 
                             invert_y=False, label_yaxis=False, pad_title=0, title_loc='right', label_counts=True,
                             highlight=None
                            ):
    n = len(features)
    labels = labels or features
    if horizontal:
        fig, axes = plt.subplots(1,n, figsize=figsize, sharex=True)
    else:    
        fig, axes = plt.subplots(n,1, figsize=figsize, sharex=True)
        
    for i, ax in enumerate(axes):
        plot_box_cluster_feature(data, y=features[i], x=x, label=labels[i], ax=ax,
                    palette=palette, strip_width=strip_width, pairs_sets=pairs_sets, cutoff=cutoff, 
                    invert_y=invert_y, label_yaxis=label_yaxis, pad_title=pad_title, title_loc=title_loc,
                    label_counts=label_counts, highlight=highlight
                    )

def plot_box_cluster_feature(data, y, x='cluster', label=None, ax=None,
                    palette=None, strip_width=0.2, pairs_sets=[], cutoff=0.05, 
                    invert_y=False, label_yaxis=False, pad_title=0, title_loc='right',
                    label_counts=True, label_color=False, size=3, highlight=None, **kwargs
                    ):
    data = data.dropna(subset=[x,y]).copy()
    if not hasattr(data[x], 'cat'):
        # make column ordered categorical
        data[x] = data[x].astype('category')
        
    label = label or y
    if ax is None:
        fig, ax = plt.subplots()
    if ax=='gca':
        ax = plt.gca()
    sns.stripplot(data=data, x=x, y=y, palette=palette, ax=ax, jitter=strip_width, size=size, alpha=0.7, **kwargs)
    sns.boxplot(data=data, x=x, y=y, palette=palette, ax=ax, showfliers=False, **kwargs)
    if highlight is not None:
        sns.stripplot(data=data.loc[data.index.intersection(highlight)], 
                      x=x, y=y, palette=palette, ax=ax, jitter=strip_width, size=size*4, marker='X', **kwargs)
    
    outline_boxplot(ax)
    ax.set_xlabel(None)
    sns.despine()
    if label_yaxis:
        ax.set_ylabel(label)
    else:
        ax.set_ylabel(None)
        ax.set_title(label, loc=title_loc, pad=pad_title)
    xlabels = [text.get_text() for text in ax.get_xmajorticklabels()]
    if label_counts:
        counts = data[x].value_counts()
        newlabels = [f"$\it{{{name}}}$\n(N={counts.loc[name]})" 
                     for name in xlabels]
        newlabels = [name.replace(' ','\\ ') for name in newlabels]
        ax.set_xticklabels(newlabels, rotation=90, ha='center',)
        if label_color:
            for xtick, name in zip(ax.get_xticklabels(), xlabels):
                xtick.set_color(palette[name])
    else:
        ax.set_xticklabels(xlabels, rotation=45, fontstyle='italic', ha='right')
    if data[y].min()>0:
        ax.set_ylim(0, None, auto=True)
        
    for pairs in pairs_sets:
        group_vals=data[x].cat.categories.to_list()
        plot_mw_bars(data, y, x, group_vals=group_vals, pairs=pairs, ax=ax, cutoff=cutoff, label=None)
        
    if data[y].min()>0:
        ax.set_ylim(0, None, auto=True)
    if invert_y:
        ax.invert_yaxis()

def run_cluster_anova(df, features, cluster_var='cluster', pval='pval_cluster', cov_type='HC3', fdr_method='fdr_bh',):
    df = df.copy()
    df['cluster'] = df[cluster_var]
    fdr_method='fdr_bh'
    
    results = (fit_models(df, ['cluster'], features, cov_type=cov_type).set_index('feature')
        .sort_values('rsquared', ascending=False)
    )
    results[pval] = (
        results[pval].pipe(pd.DataFrame)#in case this is Series
        .fillna(1)
        .apply(lambda col: multipletests(col, method=fdr_method)[1]).astype(float)
    )
    return results

def plot_cluster_anova_bar(results, pval='pval_cluster', val='rsquared', cov_type='HC3', ylabels=None, 
                               figsize=(1.5,8), nshow=20):

    data = results.loc[:,pval]
    stars = pd.cut(data.iloc[:nshow], [0, 0.001, 0.01, 0.05, 1], labels=['***','**','*',''])


    fig, ax = plt.subplots(figsize=figsize)
    bardata = results.iloc[:nshow].loc[:,val].reset_index().melt(id_vars=['feature'])
    sns.barplot(data=bardata, y='feature', x='value', hue='variable')
    sns.despine()
    if ylabels is not None:
        ax.set_yticklabels([ylabels[label.get_text()] for label in ax.get_yticklabels()])
    ax.set_ylabel(None)
    ax.set_xlabel(None)
    ax.set_title('Cluster ANOVA $\eta^2$')
    ax.get_legend().remove()

    nfeat = min(nshow, len(ylabels)) if ylabels is not None else nshow
    for i, p in enumerate(ax.patches):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        annot = stars[i]
        if not annot:
            col = p.get_facecolor()
            p.set_edgecolor(col)
            p.set_facecolor('None')
        else:
            space = 0.005
            _x = p.get_x() + p.get_width() + float(space)
            _y = p.get_y() + p.get_height()
            ax.text(_x, _y, annot, ha="left")


def run_anova_pairs(mouse_df, human_df, features, cluster='cluster', cov_type='HC3', fdr_method='fdr_bh',):
    pval = f'pval_{cluster}'
    label = 'rsquared'
    fdr_method='fdr_bh'
    
    results = (pd.concat(
        [
        fit_models(human_df, [cluster], features, cov_type=cov_type).set_index('feature'),
        fit_models(mouse_df, [cluster], features, cov_type=cov_type).set_index('feature'),
        ],
        keys=[
            'human',
            'mouse',
             ],
        axis=1)
        .sort_values(('human','rsquared'), ascending=False)
        .swaplevel(axis=1)
    )
    results = results.dropna(subset=[(pval,'human')])
    if fdr_method is not None:
        results[pval] = (
            results[pval].pipe(pd.DataFrame)#in case this is Series
            .fillna(1)
            .apply(lambda col: multipletests(col, method=fdr_method)[1]).astype(float)
        )
    return results

import warnings
def plot_cluster_anova_species_bar(mouse_df, human_df, features, pval='pval_cluster', cov_type='HC3', ylabels=None, 
                               figsize=(1.5,8), nshow=20, palette=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = run_anova_pairs(mouse_df, human_df, features, cov_type=cov_type)

    data = results.loc[:,pval]
    stars = data.iloc[:nshow].apply(lambda x: pd.cut(x, [0, 0.001, 0.01, 0.05, 1], labels=['***','**','*',''])).astype(str).values


    fig, ax = plt.subplots(figsize=figsize)
    bardata = results.iloc[:nshow].loc[:,'rsquared'].reset_index().melt(id_vars=['feature'])
    
    palette = palette or {
        # 'human': '#4d8000',
        # 'mouse': '#a84920'
        'human': '#6f2af7',
        'mouse': '#3db516'
    }
    sns.barplot(data=bardata, y='feature', x='value', hue='variable', palette=palette)
    sns.despine()
    if ylabels is not None:
        ax.set_yticklabels([ylabels[label.get_text()] for label in ax.get_yticklabels()])
    ax.set_ylabel(None)
    ax.set_xlabel(None)
    ax.set_title('Cluster ANOVA $\eta^2$')
    ax.get_legend().remove()
    
    nfeat = min(nshow, len(ylabels))
    for i, p in enumerate(ax.patches):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        ifeat = i % nfeat
        ispecies = i // nfeat
        annot = stars[ifeat, ispecies]
        if not annot:
            col = p.get_facecolor()
            p.set_edgecolor(col)
            p.set_facecolor('None')
        else:
            space = 0.005
            _x = p.get_x() + p.get_width() + float(space)
            _y = p.get_y() + p.get_height()
            ax.text(_x, _y, annot, ha="left")

