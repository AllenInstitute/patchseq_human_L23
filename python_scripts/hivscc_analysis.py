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
#     data = data.dropna(subset=variables+[feature])
    out_dict, res = ols_model(data, formula, feature)

    x = x or formula.replace('*','+').split('+')[0].strip()
    sns_args['s'] = sns_args.get('s', 25)
    sns.scatterplot(data=data, y=feature, x=x, hue=cluster, ax=ax, legend=legend, **sns_args)
    
    hue = cluster if cluster in formula else None
    c = None if cluster in formula else 'k'
    y_fit = res.fittedvalues.reindex(data.index)
    sns.lineplot(data=data, y=y_fit, x=x, hue=hue, color=c, legend=False, ax=ax)
    ax.set_xlabel(getattr(x, "label", None) or x)
    ax.set_ylabel(getattr(feature, "label", None) or feature)

    summary = ''
    if print_attr:
        value = out_dict.get(print_attr)
        attr_name = getattr(print_attr, "label", None) or print_attr
        summary = f"{attr_name} = {value:.2g}\n"
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
    
    y0 = y0 or ylim[0]
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