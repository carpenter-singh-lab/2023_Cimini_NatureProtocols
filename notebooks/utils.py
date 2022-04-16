#adapted with appreciation from https://github.com/jump-cellpainting/pilot-cpjump1-analysis

import os
import pandas as pd
import numpy as np
import scipy
import random
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_array, as_float_array
from sklearn.base import TransformerMixin, BaseEstimator
import kneed
import seaborn as sns

random.seed(9000)

def get_metacols(df):
    """return a list of metadata columns"""
    return [c for c in df.columns if c.startswith("Metadata_")]


def get_featurecols(df):
    """returna  list of featuredata columns"""
    return [c for c in df.columns if not c.startswith("Metadata")]


def get_metadata(df):
    """return dataframe of just metadata columns"""
    return df[get_metacols(df)]


def get_featuredata(df):
    """return dataframe of just featuredata columns"""
    return df[get_featurecols(df)]

def remove_negcon_empty_wells(df):
    """return dataframe of non-negative-control wells"""
    df = (
        df.query('Metadata_control_type!="negcon"')
        .dropna(subset=['Metadata_broad_sample'])
        .reset_index(drop=True)
    )
    return df


def remove_empty_wells(df):
    """return dataframe of non-empty wells"""
    df = (
        df.dropna(subset=['Metadata_broad_sample'])
        .reset_index(drop=True)
    )
    return df

def percent_score(null_dist, corr_dist, how='right'):
    """
    Calculates the Percent replicating
    :param null_dist: Null distribution
    :param corr_dist: Correlation distribution
    :param how: "left", "right" or "both" for using the 5th percentile, 95th percentile or both thresholds
    :return: proportion of correlation distribution beyond the threshold
    """
    if how == 'right':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        return np.mean(above_threshold.astype(float)), perc_95
    if how == 'left':
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return np.mean(below_threshold.astype(float)), perc_5
    if how == 'both':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return np.mean(above_threshold.astype(float)) + np.mean(below_threshold.astype(float)), perc_95, perc_5
    
def corr_between_replicates(df, group_by_feature):
    """
    Correlation between replicates
    :param df: pd.DataFrame
    :param group_by_feature: Feature name to group the data frame by
    :return: list-like of correlation values
    """
    replicate_corr = []
    replicate_grouped = df.groupby(group_by_feature)
    for name, group in replicate_grouped:
        group_features = get_featuredata(group)
        corr = np.corrcoef(group_features)
        if len(group_features) == 1:  # If there is only one replicate on a plate
            replicate_corr.append(np.nan)
        else:
            np.fill_diagonal(corr, np.nan)
            replicate_corr.append(np.nanmedian(corr))  # median replicate correlation
    return replicate_corr


def corr_between_non_replicates(df, n_samples, n_replicates, metadata_compound_name):
    """
    Null distribution between random "replicates".
    :param df: pandas.DataFrame
    :param n_samples: int
    :param n_replicates: int
    :param metadata_compound_name: Compound name feature
    :return: list-like of correlation values, with a  length of `n_samples`
    """
    df.reset_index(drop=True, inplace=True)
    null_corr = []
    while len(null_corr) < n_samples:
        compounds = random.choices([_ for _ in range(len(df))], k=n_replicates)
        sample = df.loc[compounds].copy()
        if len(sample[metadata_compound_name].unique()) == n_replicates:
            sample_features = get_featuredata(sample)
            corr = np.corrcoef(sample_features)
            np.fill_diagonal(corr, np.nan)
            null_corr.append(np.nanmedian(corr))  # median replicate correlation
    return null_corr

class ZCA_corr(BaseEstimator, TransformerMixin):
    def __init__(self, copy=False):
        self.copy = copy

    def estimate_regularization(self, eigenvalue):
        x = [_ for _ in range(len(eigenvalue))]
        kneedle = kneed.KneeLocator(x, eigenvalue, S=1.0, curve='convex', direction='decreasing')
        reg = eigenvalue[kneedle.elbow]/10.0
        return reg # The complex part of the eigenvalue is ignored

    def fit(self, X, y=None):
        """
        Compute the mean, sphereing and desphereing matrices.
        Parameters
        ----------
        X : array-like with shape [n_samples, n_features]
            The data used to compute the mean, sphereing and desphereing
            matrices.
        """
        X = check_array(X, accept_sparse=False, copy=self.copy, ensure_2d=True)
        X = as_float_array(X, copy=self.copy)
        self.mean_ = X.mean(axis=0)
        X_ = X - self.mean_
        cov = np.dot(X_.T, X_) / (X_.shape[0] - 1)
        V = np.diag(cov)
        df = pd.DataFrame(X_)
        corr = np.nan_to_num(df.corr()) # replacing nan with 0 and inf with large values
        G, T, _ = scipy.linalg.svd(corr)
        regularization = self.estimate_regularization(T.real)
        t = np.sqrt(T.clip(regularization))
        t_inv = np.diag(1.0 / t)
        v_inv = np.diag(1.0/np.sqrt(V.clip(1e-3)))
        self.sphere_ = np.dot(np.dot(np.dot(G, t_inv), G.T), v_inv)
        return self

    def transform(self, X, y=None, copy=None):
        """
        Parameters
        ----------
        X : array-like with shape [n_samples, n_features]
            The data to sphere along the features axis.
        """
        check_is_fitted(self, "mean_")
        X = as_float_array(X, copy=self.copy)
        return np.dot(X - self.mean_, self.sphere_.T)


def sphere_plate_zca_corr(plate):
    """
    sphere each plate to the DMSO negative control values
    Parameters:
    -----------
    plate: pandas.DataFrame
        dataframe of a single plate's featuredata and metadata
    Returns:
    -------
    pandas.DataFrame of the same shape as `plate`
    """
    # sphere featuredata to DMSO sphereing matrix
    sphereer = ZCA_corr()
    dmso_df = plate.loc[plate.Metadata_control_type=="negcon"]
    # dmso_df = plate.query("Metadata_pert_type == 'control'")
    dmso_vals = get_featuredata(dmso_df).to_numpy()
    all_vals = get_featuredata(plate).to_numpy()
    sphereer.fit(dmso_vals)
    sphereed_vals = sphereer.transform(all_vals)
    # concat with metadata columns
    feature_df = pd.DataFrame(
        sphereed_vals, columns=get_featurecols(plate), index=plate.index
    )
    metadata = get_metadata(plate)
    combined = pd.concat([feature_df, metadata], axis=1)
    assert combined.shape == plate.shape
    return combined

def calculate_percent_replicating(batch_path,plate):
    metadata_plate = 'Metadata_Plate'
    metadata_compound_name = 'Metadata_pert_iname'
    metadata_moa = 'Metadata_moa'
    corr_strong_df= pd.DataFrame()
    corr_recall_df= pd.DataFrame()
    compound_names = []
    n_samples_strong = 10000
    data_df = pd.read_csv(os.path.join(batch_path, plate,
                                           plate+'_normalized_feature_select_negcon.csv.gz'))

    data_df = sphere_plate_zca_corr(data_df)

    data_df = remove_negcon_empty_wells(data_df)

    replicate_corr = list(corr_between_replicates(data_df, metadata_compound_name))
    null_corr = list(corr_between_non_replicates(data_df, n_samples=n_samples_strong, n_replicates=4, metadata_compound_name = metadata_compound_name))

    prop_95, value_95 = percent_score(null_corr, replicate_corr)

    return(prop_95)

def plot_simple_comparison(df,x,hue):
    sns.set_style("ticks")
    sns.set_context("paper",font_scale=1.5)
    g = sns.catplot(data=df, x = x ,y='Percent Replicating',hue=hue,palette='Set1',s=8,linewidth=1)
    g.set(ylim=([0,1]))
    plotname = f"../figures/{x}-{hue}.png"
    g.savefig(plotname,dpi=300)
    print(f'Saved to {plotname}')