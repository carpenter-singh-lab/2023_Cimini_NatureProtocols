#adapted with appreciation from https://github.com/jump-cellpainting/pilot-cpjump1-analysis

import os
import random
import textwrap

import pandas as pd
import numpy as np
import scipy
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_array, as_float_array
from sklearn.base import TransformerMixin, BaseEstimator
import kneed
import matplotlib.pyplot as plt
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
        return 100 * np.mean(above_threshold.astype(float)), perc_95
    if how == 'left':
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return 100 * np.mean(below_threshold.astype(float)), perc_5
    if how == 'both':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return 100 * (np.mean(above_threshold.astype(float)) + np.mean(below_threshold.astype(float))), perc_95, perc_5
    
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

def corr_between_replicates_across_plates(df, reference_df, pertcol = 'Metadata_pert_iname'):
    items = list(df[pertcol].unique())
    common_columns = [x for x in df.columns if x in reference_df.columns]
    df = df[common_columns]
    reference_df = reference_df[common_columns]
    replicate_corr = []
    for item in items:
        compound_df = df.query(pertcol+' == @item')
        compound_reference_df = reference_df.query(pertcol+' == @item')

        compound_df_profiles = get_featuredata(compound_df).values
        compound_reference_df_profiles = get_featuredata(compound_reference_df).values
        try:
            corr = np.corrcoef(compound_df_profiles, compound_reference_df_profiles)
            corr = corr[0:len(compound_df_profiles), len(compound_df_profiles):]

            corr_median_value = np.nanmedian(corr, axis=1)
            corr_median_value = np.nanmedian(corr_median_value)

            replicate_corr.append(corr_median_value)
        except:
            print(item,corr)
    return replicate_corr

def corr_between_non_replicates_across_plates(df, reference_df, n_samples, pertcol = 'Metadata_pert_iname'):
    np.random.seed(9000)
    common_columns = [x for x in df.columns if x in reference_df.columns]
    df = df[common_columns]
    reference_df = reference_df[common_columns]
    null_corr = []
    compounds = list(df[pertcol].unique())  
    while len(null_corr) < n_samples:
        both_compounds = np.random.choice(compounds, size=2, replace=False)
        compound1 = both_compounds[0]
        compound2 = both_compounds[1]

        compound1_df = df.query(pertcol+' == @compound1')
        compound2_df = reference_df.query(pertcol+' == @compound2')

        compound1_df_profiles = get_featuredata(compound1_df).values
        compound2_df_profiles = get_featuredata(compound2_df).values
        try:
            corr = np.corrcoef(compound1_df_profiles, compound2_df_profiles)
            corr = corr[0:len(compound1_df_profiles), len(compound1_df_profiles):]

            corr_median_value = np.nanmedian(corr, axis=1)
            corr_median_value = np.nanmedian(corr_median_value)

            null_corr.append(corr_median_value)
        except:
            pass

    return null_corr

def corr_between_compound_moa(df, metadata_moa, metadata_compound_name):
    """
        Correlation between compounds with the same MOA
        Parameters:
        -----------
        df: pd.DataFrame
        metadata_moa: MOA feature
        metadata_compound_name: Compound name feature
        Returns:
        --------
        list-like of correlation values
     """
    replicate_corr = []

    profile_df = (
        get_metadata(df)
        .assign(profiles=list(get_featuredata(df).values))
    )

    replicate_grouped = (
        profile_df.groupby([metadata_moa, metadata_compound_name]).profiles
            .apply(list)
            .reset_index()
    )

    moa_grouped = (
        replicate_grouped.groupby([metadata_moa]).profiles
            .apply(list)
            .reset_index()
    )

    for i in range(len(moa_grouped)):
        if len(moa_grouped.iloc[i].profiles) > 1:
            compound1_profiles = moa_grouped.iloc[i].profiles[0]
            compound2_profiles = moa_grouped.iloc[i].profiles[1]

            corr = np.corrcoef(compound1_profiles, compound2_profiles)
            corr = corr[0:len(moa_grouped.iloc[i].profiles[0]), len(moa_grouped.iloc[i].profiles[0]):]
            # np.fill_diagonal(corr, np.nan)
            replicate_corr.append(np.nanmedian(corr))

    return replicate_corr

def null_corr_between_compound_moa(df, n_samples, metadata_moa, metadata_compound_name):
    """
        Null distribution between random pairs of compounds.
        Parameters:
        ------------
        df: pandas.DataFrame
        n_samples: int
        metadata_moa: MOA feature
        metadata_compound_name: Compound name feature
        Returns:
        --------
        list-like of correlation values, with a  length of `n_samples`
    """
    df.reset_index(drop=True, inplace=True)
    null_corr = []

    profile_df = (
        get_metadata(df)
        .assign(profiles=list(get_featuredata(df).values))
    )

    replicate_grouped = (
        profile_df.groupby([metadata_moa, metadata_compound_name]).profiles
            .apply(list)
            .reset_index()
    )

    while len(null_corr) < n_samples:
        compounds = random.choices([_ for _ in range(len(replicate_grouped))], k=2)
        compound1_moa = replicate_grouped.iloc[compounds[0]].Metadata_moa
        compound2_moa = replicate_grouped.iloc[compounds[1]].Metadata_moa
        if compound1_moa != compound2_moa:
            compound1_profiles = replicate_grouped.iloc[compounds[0]].profiles
            compound2_profiles = replicate_grouped.iloc[compounds[1]].profiles
            corr = np.corrcoef(compound1_profiles, compound2_profiles)
            corr = corr[0:len(replicate_grouped.iloc[0].profiles), len(replicate_grouped.iloc[0].profiles):]
            # np.fill_diagonal(corr, np.nan)
            null_corr.append(np.nanmedian(corr))  # median replicate correlation
    return null_corr

def correlation_between_modalities(modality_1_df, modality_2_df, modality_1, modality_2, metadata_common, metadata_perturbation):
    """
    Compute the correlation between two different modalities.
    :param modality_1_df: Profiles of the first modality
    :param modality_2_df: Profiles of the second modality
    :param modality_1: feature that identifies perturbation pairs
    :param modality_2: perturbation name feature
    :param metadata_common: perturbation name feature
    :param metadata_perturbation: perturbation name feature
    :return: list-like of correlation values
    """
    list_common_perturbation_groups = list(np.intersect1d(list(modality_1_df[metadata_common]), list(modality_2_df[metadata_common])))

    merged_df = pd.concat([modality_1_df, modality_2_df], ignore_index=False, join='inner')

    modality_1_df = merged_df.query('Metadata_modality==@modality_1')
    modality_2_df = merged_df.query('Metadata_modality==@modality_2')

    corr_modalities = []

    for group in list_common_perturbation_groups:
        modality_1_perturbation_df = modality_1_df.loc[modality_1_df[metadata_common] == group]
        modality_2_perturbation_df = modality_2_df.loc[modality_2_df[metadata_common] == group]

        for sample_1 in modality_1_perturbation_df[metadata_perturbation].unique():
            for sample_2 in modality_2_perturbation_df[metadata_perturbation].unique():
                modality_1_perturbation_sample_df = modality_1_perturbation_df.loc[modality_1_perturbation_df[metadata_perturbation] == sample_1]
                modality_2_perturbation_sample_df = modality_2_perturbation_df.loc[modality_2_perturbation_df[metadata_perturbation] == sample_2]

                modality_1_perturbation_profiles = get_featuredata(modality_1_perturbation_sample_df)
                modality_2_perturbation_profiles = get_featuredata(modality_2_perturbation_sample_df)

                corr = np.corrcoef(modality_1_perturbation_profiles, modality_2_perturbation_profiles)
                corr = corr[0:len(modality_1_perturbation_profiles), len(modality_1_perturbation_profiles):]
                corr_modalities.append(np.nanmedian(corr))  # median replicate correlation

    return corr_modalities

def null_correlation_between_modalities(modality_1_df, modality_2_df, modality_1, modality_2, metadata_common, metadata_perturbation, n_samples):
    """
    Compute the correlation between two different modalities.
    :param modality_1_df: Profiles of the first modality
    :param modality_2_df: Profiles of the second modality
    :param modality_1: "Compound", "ORF" or "CRISPR"
    :param modality_2: "Compound", "ORF" or "CRISPR"
    :param metadata_common: feature that identifies perturbation pairs
    :param metadata_perturbation: perturbation name feature
    :param n_samples: int
    :return:
    """
    list_common_perturbation_groups = list(np.intersect1d(list(modality_1_df[metadata_common]), list(modality_2_df[metadata_common])))

    merged_df = pd.concat([modality_1_df, modality_2_df], ignore_index=False, join='inner')

    modality_1_df = merged_df.query('Metadata_modality==@modality_1')
    modality_2_df = merged_df.query('Metadata_modality==@modality_2')

    null_modalities = []

    count = 0

    while count < n_samples:
        perturbations = random.choices(list_common_perturbation_groups, k=2)
        modality_1_perturbation_df = modality_1_df.loc[modality_1_df[metadata_common] == perturbations[0]]
        modality_2_perturbation_df = modality_2_df.loc[modality_2_df[metadata_common] == perturbations[1]]

        for sample_1 in modality_1_perturbation_df[metadata_perturbation].unique():
            for sample_2 in modality_2_perturbation_df[metadata_perturbation].unique():
                modality_1_perturbation_sample_df = modality_1_perturbation_df.loc[modality_1_perturbation_df[metadata_perturbation] == sample_1]
                modality_2_perturbation_sample_df = modality_2_perturbation_df.loc[modality_2_perturbation_df[metadata_perturbation] == sample_2]

                modality_1_perturbation_profiles = get_featuredata(modality_1_perturbation_sample_df)
                modality_2_perturbation_profiles = get_featuredata(modality_2_perturbation_sample_df)

                corr = np.corrcoef(modality_1_perturbation_profiles, modality_2_perturbation_profiles)
                corr = corr[0:len(modality_1_perturbation_profiles), len(modality_1_perturbation_profiles):]
                null_modalities.append(np.nanmedian(corr))  # median replicate correlation
        count += 1

    return null_modalities

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

def calculate_percent_replicating_MOA(batch_path,plate,data_df=None):
    """
    For plates treated with the JUMP-MOA source plates, at least 
    4 copies of each perturbation are present on each plate.
    Percent replicating is therefore calculated per plate.
    """
    metadata_compound_name = 'Metadata_pert_iname'
    n_samples_strong = 10000
    if type(data_df)!=pd.DataFrame:
        data_df = pd.read_csv(os.path.join(batch_path, plate,
                                            plate+'_normalized_feature_select_negcon.csv.gz'))

    data_df = sphere_plate_zca_corr(data_df)

    data_df = remove_negcon_empty_wells(data_df)

    replicate_corr = list(corr_between_replicates(data_df, metadata_compound_name))
    null_corr = list(corr_between_non_replicates(data_df, n_samples=n_samples_strong, n_replicates=4, metadata_compound_name = metadata_compound_name))

    prop_95, _ = percent_score(null_corr, replicate_corr)

    return(prop_95)

def calculate_percent_replicating_across_plates_MOA(batch_path1,plate1,batch_path2,plate2 ):
    """
    For plates treated with the JUMP-MOA source plates, at least 
    4 copies of each perturbation are present on each plate.
    Percent replicating is therefore calculated per plate.
    """
    metadata_compound_name = 'Metadata_pert_iname'
    n_samples_strong = 10000

    data_df1 = pd.read_csv(os.path.join(batch_path1, plate1,
                                           plate1+'_normalized_feature_select_negcon.csv.gz'))
    data_df1 = sphere_plate_zca_corr(data_df1)
    data_df1 = remove_negcon_empty_wells(data_df1)

    data_df2 = pd.read_csv(os.path.join(batch_path2, plate2,
                                           plate2+'_normalized_feature_select_negcon.csv.gz'))
    data_df2 = sphere_plate_zca_corr(data_df2)
    data_df2 = remove_negcon_empty_wells(data_df2)

    replicate_corr = corr_between_replicates_across_plates(data_df1, data_df2)
    null_corr = corr_between_non_replicates_across_plates(data_df1, data_df2, n_samples=n_samples_strong)

    prop_95, _ = percent_score(null_corr, replicate_corr)

    return(prop_95)

def calculate_percent_matching_MOA(batch_path,plate):
    """
    For plates treated with the JUMP-MOA source plates, at least 
    4 copies of each perturbation are present on each plate.
    Percent replicating is therefore calculated per plate.
    """
    metadata_moa_name = 'Metadata_moa'
    metadata_compound_name = 'Metadata_pert_iname'
    n_samples_strong = 10000
    data_df = pd.read_csv(os.path.join(batch_path, plate,
                                           plate+'_normalized_feature_select_negcon.csv.gz'))

    data_df = sphere_plate_zca_corr(data_df)

    data_df = remove_negcon_empty_wells(data_df)

    replicate_corr = list(corr_between_compound_moa(data_df, metadata_moa_name, metadata_compound_name))
    null_corr = list(null_corr_between_compound_moa(data_df, n_samples_strong, metadata_moa_name, metadata_compound_name))

    prop_95, _ = percent_score(null_corr, replicate_corr)

    return(prop_95)

def calculate_percent_replicating_Target(batch_path,platelist,sphere=None,
suffix = '_normalized_feature_select_negcon.csv.gz',n_replicates=4):
    """
    For plates treated with the JUMP-Target source plates, most 
    perturbations are only present in one or two 2 copies per plate. 
    Percent replicating is therefore calculated per group of replicate plates.

    Since feature selection happens on a per-plate level, an inner join
    is performed across all plates in the replicate, meaning only common
    features are used in calculation of percent replicating.

    It doesn't look like sphering was done consistently in previous 
    analysis of these plates, therefore it is configurable here; either 
    not done, done at the plate level by passing 'sphere=plate', or 
    done at the batch level by passing 'sphere=batch'.
    """
    metadata_compound_name = 'Metadata_broad_sample'
    n_samples_strong = 10000

    data_dict = {}

    for plate in platelist:
        plate_df = pd.read_csv(os.path.join(batch_path, plate,
                                            plate+suffix))
        
        if sphere == 'plate':
            plate_df = sphere_plate_zca_corr(plate_df)

        data_dict[plate] = plate_df
    
    data_df = pd.concat(data_dict, join='inner', ignore_index=True)

    if sphere == 'batch':
        data_df = sphere_plate_zca_corr(data_df)

    data_df = remove_negcon_empty_wells(data_df)

    replicate_corr = list(corr_between_replicates(data_df, metadata_compound_name))
    null_corr = list(corr_between_non_replicates(data_df, n_samples=n_samples_strong, n_replicates=n_replicates, metadata_compound_name = metadata_compound_name))

    prop_95, _ = percent_score(null_corr, replicate_corr)

    return(prop_95)

def calculate_percent_replicating_across_plates_Target(batch_path1,plate1,batch_path2,plate2 ):
    """
    For Target 1 vs Target 2
    """
    metadata_compound_name = 'Metadata_broad_sample'
    n_samples_strong = 10000

    data_df1 = pd.read_csv(os.path.join(batch_path1, plate1,
                                           plate1+'_normalized_feature_select_negcon.csv.gz'))
    data_df1 = remove_negcon_empty_wells(data_df1)

    data_df2 = pd.read_csv(os.path.join(batch_path2, plate2,
                                           plate2+'_normalized_feature_select_negcon.csv.gz'))
    data_df2 = remove_negcon_empty_wells(data_df2)

    replicate_corr = corr_between_replicates_across_plates(data_df1, data_df2,pertcol=metadata_compound_name)
    null_corr = corr_between_non_replicates_across_plates(data_df1, data_df2, n_samples=n_samples_strong,pertcol=metadata_compound_name)

    prop_95, _ = percent_score(null_corr, replicate_corr)

    return(prop_95)

def calculate_percent_matching_Target(batch_path_1,platelist_1,modality_1, batch_path_2,platelist_2, modality_2,
sphere=None,suffix = '_normalized_feature_select_negcon.csv.gz'):
    """

    It doesn't look like sphering was done consistently in previous 
    analysis of these plates, therefore it is configurable here; either 
    not done, done at the plate level by passing 'sphere=plate', or 
    done at the batch level by passing 'sphere=batch'.
    """
    n_samples = 10000

    data_dict_1 = {}
    for plate in platelist_1:
        plate_df = pd.read_csv(os.path.join(batch_path_1, plate,
                                            plate+suffix))
        if sphere == 'plate':
            plate_df = sphere_plate_zca_corr(plate_df)

        data_dict_1[plate] = plate_df   
    data_df_1 = pd.concat(data_dict_1, join='inner', ignore_index=True)
    if modality_1 =='Compounds':
        data_df_1.rename(columns={'Metadata_target':'Metadata_genes'},inplace=True)
    data_df_1['Metadata_modality'] = modality_1
    if sphere == 'batch':
        data_df_1 = sphere_plate_zca_corr(data_df_1)
    data_df_1 = remove_negcon_empty_wells(data_df_1)

    data_dict_2 = {}
    for plate in platelist_2:
        plate_df = pd.read_csv(os.path.join(batch_path_2, plate,
                                            plate+suffix))
        if sphere == 'plate':
            plate_df = sphere_plate_zca_corr(plate_df)

        data_dict_2[plate] = plate_df   
    data_df_2 = pd.concat(data_dict_2, join='inner', ignore_index=True)
    if modality_2 =='Compounds':
        data_df_2.rename(columns={'Metadata_target':'Metadata_genes'},inplace=True)
    data_df_2['Metadata_modality'] = modality_2
    if sphere == 'batch':
        data_df_2 = sphere_plate_zca_corr(data_df_2)
    data_df_2 = remove_negcon_empty_wells(data_df_2)

    replicate_corr = list(correlation_between_modalities(data_df_1, data_df_2, modality_1, modality_2, 'Metadata_genes', 'Metadata_broad_sample'))
    null_corr = list(null_correlation_between_modalities(data_df_1, data_df_2, modality_1, modality_2, 'Metadata_genes', 'Metadata_broad_sample', n_samples))

    prop_95, _, _ = percent_score(null_corr, replicate_corr, how='both')

    return(prop_95)

def calculate_percent_matching_with_feature_dropout_Target(batch_path_1,platelist_1,modality_1, batch_path_2,platelist_2, modality_2,
sphere=None,suffix = '_normalized_negcon.csv.gz',drop='AGP'):
    """

    It doesn't look like sphering was done consistently in previous 
    analysis of these plates, therefore it is configurable here; either 
    not done, done at the plate level by passing 'sphere=plate', or 
    done at the batch level by passing 'sphere=batch'.
    """
    import pycytominer
    n_samples = 10000

    data_dict_1 = {}
    for plate in platelist_1:
        plate_df = pd.read_csv(os.path.join(batch_path_1, plate,
                                            plate+suffix))
        cols_to_drop = [x for x in plate_df.columns if drop in x]
        plate_df.drop(columns=cols_to_drop,inplace=True)
        feature_select_features = pycytominer.cyto_utils.infer_cp_features(
        plate_df
        )
        plate_df = pycytominer.feature_select(
        profiles=plate_df,
        features=feature_select_features,
        operation=['variance_threshold','correlation_threshold',
        'drop_na_columns','blocklist']
        )
        if sphere == 'plate':
            plate_df = sphere_plate_zca_corr(plate_df)

        data_dict_1[plate] = plate_df   
    data_df_1 = pd.concat(data_dict_1, join='inner', ignore_index=True)
    if modality_1 =='Compounds':
        data_df_1.rename(columns={'Metadata_target':'Metadata_genes'},inplace=True)
    data_df_1['Metadata_modality'] = modality_1
    if sphere == 'batch':
        data_df_1 = sphere_plate_zca_corr(data_df_1)
    data_df_1 = remove_negcon_empty_wells(data_df_1)

    data_dict_2 = {}
    for plate in platelist_2:
        plate_df = pd.read_csv(os.path.join(batch_path_1, plate,
                                            plate+suffix))
        cols_to_drop = [x for x in plate_df.columns if drop in x]
        plate_df.drop(columns=cols_to_drop,inplace=True)
        feature_select_features = pycytominer.cyto_utils.infer_cp_features(
        plate_df
        )
        plate_df = pycytominer.feature_select(
        profiles=plate_df,
        features=feature_select_features,
        operation=['variance_threshold','correlation_threshold',
        'drop_na_columns','blocklist']
        )
        if sphere == 'plate':
            plate_df = sphere_plate_zca_corr(plate_df)

        data_dict_2[plate] = plate_df   
    data_df_2 = pd.concat(data_dict_2, join='inner', ignore_index=True)
    if modality_2 =='Compounds':
        data_df_2.rename(columns={'Metadata_target':'Metadata_genes'},inplace=True)
    data_df_2['Metadata_modality'] = modality_2
    if sphere == 'batch':
        data_df_2 = sphere_plate_zca_corr(data_df_2)
    data_df_2 = remove_negcon_empty_wells(data_df_2)

    replicate_corr = list(correlation_between_modalities(data_df_1, data_df_2, modality_1, modality_2, 'Metadata_genes', 'Metadata_broad_sample'))
    null_corr = list(null_correlation_between_modalities(data_df_1, data_df_2, modality_1, modality_2, 'Metadata_genes', 'Metadata_broad_sample', n_samples))

    prop_95, _, _ = percent_score(null_corr, replicate_corr, how='both')

    return(prop_95)

def plot_simple_comparison(df,x,hue,y='Percent Replicating',order=None,hue_order=None,
col=None, col_order=None, col_wrap=None,row=None,row_order=None,jitter=0.25,dodge=True,plotname=None,
ylim=None, title=None,aspect=1,sharex=True,facet_kws={}):
    plt.rcParams["legend.markerscale"] =1.5
    sns.set_style("ticks")
    sns.set_context("paper",font_scale=1.5)
    g = sns.catplot(data=df, x = x ,y = y, order=order,
    hue=hue, hue_order=hue_order, col=col, col_order = col_order, col_wrap=col_wrap,row=row,
    row_order = row_order, palette='Set1',s=80,linewidth=1,jitter=jitter,
    alpha=0.8,dodge=dodge,aspect=aspect,sharex=sharex,facet_kws=facet_kws)
    if sharex:
        labels = []
        if not order:
            orig_labels = list(dict.fromkeys(df[x].values).keys())
        else:
            orig_labels = order
        for label in orig_labels:
            if type(label)!= str:
                label = str(int(label))
            labels.append(textwrap.fill(label, width=45/len(orig_labels),break_long_words=False))
        g.set_xticklabels(labels=labels,rotation=0)
    if ylim:
        ymin,ymax=ylim
    else:
        ymin = 50
        ymax = 80 
    if df[y].min()<ymin:
        ymin = df[y].min()-2
    if df[y].max()>ymax:
        ymax = df[y].max()+2
    g.set(ylim=([ymin,ymax]))    
    if plotname:
        plotname = f"../figures/{plotname}"
    else:
        plotname = f"../figures/{x}-{y}-{hue}-{col}-{row}.svg"
    if not col:
        if not row:
            if title:
                g.set(title=title)
            else:
                g.set(title=f"{x}-{y}")
    g.savefig(plotname,dpi=300)
    print(f'Saved to {plotname}')

def plot_two_comparisons(df,x='Percent Replicating',y='Percent Matching',hue = None, hue_order=None,
col=None, col_order=None,col_wrap=None,row=None,row_order=None,style=None,xlim=None,ylim=None,title=None,
title_variable = None, facet_kws={'sharex':True}):
    plt.rcParams["legend.markerscale"] =1.5
    sns.set_style("ticks")
    sns.set_context("paper",font_scale=1.5)
    g = sns.relplot(data=df, x = x ,y= y, hue=hue, hue_order=hue_order, col=col, col_order = col_order, 
    col_wrap=col_wrap, row=row, row_order = row_order, style = style, palette='Set1',edgecolor='k',alpha=0.8,s=80,
    facet_kws=facet_kws)
    if xlim:
        xmin,xmax=xlim
    else:
        xmin = 50
        xmax = 80
    if df[x].min()<xmin:
        xmin = df[x].min()-2
    if df[x].max()>xmax:
        xmax = df[x].max()+2
    g.set(xlim=([xmin,xmax]))  
    if ylim:
        ymin,ymax=ylim
    else:
        ymin = 5
        ymax = 40 
    if df[y].min()<ymin:
        ymin = df[y].min()-2
    if df[y].max()>ymax:
        ymax = df[y].max()+2
    g.set(ylim=([ymin,ymax]))
    if title:
        g.set(title=title)
    elif title_variable:
        g.set(title=title_variable)
    plotname = f"../figures/{x}-{y}-{hue}-{col}-{row}-{style}.svg"
    g.savefig(plotname,dpi=300)
    print(f'Saved to {plotname}')

def enforce_modality_match_order(modality1,modality2):
    modality_dict_forward = {'Compounds':1,'ORF':2,'CRISPR':3}
    modality_dict_reverse = {v:k for k,v in modality_dict_forward.items()}
    modality_list = [modality_dict_forward[modality1],modality_dict_forward[modality2]]
    modality_list.sort()
    return f"{modality_dict_reverse[modality_list[0]]} - {modality_dict_reverse[modality_list[1]]}"

def enforce_timepoint_order(timepoint1,timepoint2):
    timepoint_list = [int(timepoint1),int(timepoint2)]
    timepoint_list.sort()
    return f"{timepoint_list[0]}-{timepoint_list[1]}"

def enforce_timepoint_order_in_plot(timepointlist):
    timepointlist=list(set(timepointlist))
    intlist = []
    timepoint_dict = {}
    for x in timepointlist:
        first,second=x.split('-')
        intlist.append([int(first),int(second)])
    intlist.sort()
    outlist = []
    for x in intlist:
        outlist.append(f"{x[0]}-{x[1]}")
    return outlist

def enforce_modality_match_order_in_plot(modalitylist):
    modalitylist=list(set(modalitylist))
    modality_dict_forward = {'Compounds':1,'ORF':2,'CRISPR':3}
    modality_dict_reverse = {v:k for k,v in modality_dict_forward.items()}
    master_modalitylist = []
    for eachpair in modalitylist:
        modality1,modality2 = eachpair.split(' - ')
        sublist = [modality_dict_forward[modality1],modality_dict_forward[modality2]]
        sublist.sort()
        master_modalitylist.append(sublist)
    master_modalitylist.sort()
    outlist = [f"{modality_dict_reverse[x[0]]} - {modality_dict_reverse[x[1]]}" for x in master_modalitylist]    
    return outlist
