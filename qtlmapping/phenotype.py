# -*- UTF-8 -*-

# TODO: Finish the document
import os
import sys
import copy
import math
import json

from os.path import join as pjoin

# Third-party
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sbn
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

from sklearn.decomposition import PCA

# Local libraries
from preproc import PreProc
from utls import log_me, not_runnable
from utls import UnknownCorrelationMethod


class Phenotype(PreProc):
    """An phenotype object
    Attributes:
    Examples:
    To-dos:
    Notes:
        1. If there are less than three phenotypes, remove outliers by PCA isn't plausible.
    """

    def __init__(self, ipt, wk_dir):
        """A class of phenotype.

        Args: 
            ipt (str): REQUIRED the path to the phenotype files. It requires in specific format, please refer the "Misc" -> "Phenotype file format" section in README for more information
            wk_dir (str): REQUIRED
        """
        if not (os.path.isfile(ipt) or os.path.isdir(ipt)):
            raise TypeError("ipt should be an readable file or direcotry")

        if not os.path.isdir(wk_dir):
            raise TypeError("wk_dir should be a direcotry")

        super().__init__(ipt)
        self.wk_dir = wk_dir
        self.norm_plots = None          # Histgram to check the normality of each phenotype measurement
        self.pca = None                 # The result of PCA
        self.pca_plot = None            # PCA plot to show the outliers. Phenotype vector in new space could be informative
        self.safe_zone = None           # [mean - pad * std, mean + pad * std]
        self.pw_corr = None             # The pairwise correlation of phenotype values
        self.pw_corr_plot = None           # Pair-wise correlation plot of phenotypes
        self.pheno_signals = None       # Phenotypes with too many repeated values
        self._pad = 0

    def check_norm(self):
        """Check the normality of each phenotype.

        Args: None

        Raises: None

        Returns:
            self: the instance itself

        To-dos: this should be achieve as interactive.
        """
        mnfl_dtfm = self.wk_dtfm.fillna(self.wk_dtfm.mean())

        norm_plots = sbn.PairGrid(mnfl_dtfm)
        norm_plots.map_upper(sbn.regplot, scatter_kws={"s": 1}, line_kws={"linewidth": 1})
        norm_plots.map_diag(sbn.distplot)
        norm_plots.map_lower(sbn.kdeplot, linewidth=0.5)

        self.norm_plots = norm_plots
        
        return self
        
    def transform_by(self, func=math.log10, target=None):
        """Transform the phenotype measurement using a function (e.g. log10).

        Args:
            func (callable): [math.log10] Function used to transform the phenotype data.
            target (None, list, tuple): [None] Columns will be transformed by given `func`.

        Raises: None

        Returns:
            Self: the instance itself
        """
        if target:
            log_me.info("Transform the following phenotypes by {}: {}".format(func.__name__, ", ".join(target)))
            self.wk_dtfm.loc[:, target] = self.wk_dtfm.loc[:, target].applymap(func)

            if not isinstance(target, (list, tuple)):
                target = [target]
        else:
            log_me.info("Transform all phenotypes by {}".format(func.__name__))
            self.wk_dtfm = self.wk_dtfm.applymap(func)

            target = list(self.wk_dtfm.columns)

        name_map = dict(zip(target, [x + "_" + func.__name__ for x in target]))

        self.wk_dtfm = self.wk_dtfm.rename(columns=name_map)

        return self

    def calc_pca(self, target=None, pad=3):
        """Calulate PCA to pin-point outliers.

        Args:
            target (list, tuple, None): [None] Columns will be used to do PCA analysis. If it's `None`, the program will take all columns to do the analysis.
            pad (int): [3] Times of standard deviation.

        Raises: None

        Returns:
            self: the instance itself
        """
        def calc_comp_safe_zone(vals, pad):
            _mean = np.mean(vals)
            _std = np.std(vals)
            return [_mean - pad * _std, _mean + pad * _std]

        pca = PCA(n_components="mle")  # Automatic choosing dimentionality
        if target is None:
            log_me.info("Will use all features in `wk_dtfm` to do PCA analysis")
            target = self.wk_dtfm.columns

            if len(target) < 3:
                log_me.warn("PCA requires > 2 features, but found {}".format(len(target)))
                return self

            self.pca = pca.fit_transform(self.wk_dtfm)
        elif isinstance(target, (list, tuple)):
            if len(target) < 3:
                log_me.warn("PCA requires > 2 features, but found {}".format(len(target)))
                return self

            self.pca = pca.fit_transform(self.wk_dtfm.loc[:, target])
        else:
            raise TypeError("Required list, tuple or None. If it's None, take all features.")

        comp1_ratio, comp2_ratio, *_ = pca.explained_variance_ratio_  # Explained variance
        comp1_data, comp2_data = self.pca[:, 0], self.pca[:, 1]

        comp1_safe = calc_comp_safe_zone(comp1_data, pad)
        comp2_safe = calc_comp_safe_zone(comp2_data, pad)

        sample_id = self.wk_dtfm.index

        outlier_dtfm = pd.DataFrame(
            [
                [sample_id[idx], x, y, False]  # Non-outliers
                if (comp1_safe[0] < x < comp1_safe[1]) \
                    or (comp2_safe[0] < y < comp2_safe[1])
                else 
                    [sample_id[idx], x, y, True]  # Outliers
                for idx, (x, y) in enumerate(zip(comp1_data, comp2_data))
            ],
            columns=["sample_id", "PC1", "PC2", "outlier"],
        ) # Outliers defined by "pad" times of standard devition

        grid = sbn.JointGrid(x="PC1", y="PC2", data=outlier_dtfm)
        grid.plot_joint(sbn.scatterplot, data=outlier_dtfm,
            hue="outlier", hue_order=[True, False],
            size="outlier", size_order=[True, False]
        )
        grid.plot_marginals(sbn.distplot, kde=False, color="0.5")

        grid.set_axis_labels(
            "PC1 ({:.2f}%)".format(comp1_ratio * 100),
            "PC2 ({:.2f}%)".format(comp2_ratio * 100)
        )

        for text, x, y, is_outlier in outlier_dtfm.values:
            if is_outlier == "Outlier":
                grid.ax_joint.annotate(text, (x, y))

        self.pca_plot = grid

        return self

    def mask_outliers(self, by=np.NaN, pad=3):
        """Mask outliers by given value.

        Args:
            by (numpy.NaN, int, float): [numpy.NaN] Will replace the measure of outliers by numpy.NaN
            pad (int): [3] Times of standard deviation
        
        Raises: None

        Returns:
            Self: the instance itself
        """
        def is_outlier(pheno_vals, safe_zone):
            pheno_name = pheno_vals.name
            lower = safe_zone.loc[pheno_name, "lower"]
            upper = safe_zone.loc[pheno_name, "upper"]
            return ~pheno_vals.between(lower, upper)

        self._pad = 3
        means = self.wk_dtfm.apply(np.mean, axis=0)
        stds = self.wk_dtfm.apply(np.std, axis=0)
        lower, upper = means - pad * stds, means + pad * stds
        self.safe_zone = pd.DataFrame(dict(lower=lower, upper=upper))

        outlier_matrix = self.wk_dtfm.apply(is_outlier, safe_zone=self.safe_zone)
        self.wk_dtfm = self.wk_dtfm.mask(outlier_matrix)  

        return self

    def calc_corr(self, target=None, method="spearman", dist_mthd="ward"):
        """Caluculate the correlation between phenotypes (Pairwise).

        Args:
            target (list, tuple, None): [None] All the combination or subset of all phenotypes, the later can be specified by and list
            method (str, callable): ['spearman'] The method use to calculate the correlation. Options: [spearman, pearson, kendall]. For more info, refer to `pandas.DataFrame.coor()`
            dist_mthd (str): ['ward'] Method used to calculate the distance between correlations.

        Raises:
            UknownCorrelationMethod: Unknown type of correlation method>
            TypeError: When the target augument given value type beyound None, tuple, and list

        Returns:
            self: The instance itself
        """
        if method in ["spearman", "pearson", "kendall"]:
            log_me.info("Will use {}".format(method))
        elif callable(method):
            log_me.info("User custom callable {}".format(method))
        else:
            raise UnknownCorrelationMethod

        if target is None:
            self.pw_corr = self.wk_dtfm.corr(method=method)  # Pairwise
        elif isinstance(target, (tuple, list)):
            self.pw_corr = self.wk_dtfm.loc[:, target].corr(method=method)
        else:
            raise TypeError("`target` should be a list or tuple, otherwise leave target as None using all")

        self.pw_corr_plot = sbn.clustermap(self.pw_corr)
        
        return self

    def pop_bad_phenos(self, threshold=0.5, target=None):
        """Remove phenotypes with low signal.
        Args:
            threshold (float): [0.5] The ratio of the mode value should not be higher than the threshold
            target (list, tuple, None): [None]. Only check given phenotypes by target

        Raises: None

        Return:
            self: the instance itself
        """
        if isinstance(target, (list, tuple)):
            log_me.info("Will only test: {}".format(", ".join(target)))
            pheno_mode_n = self.wk_dtfm.loc[:, target].apply(pd.value_counts).max()
        else:
            log_me.info("Will check all the phenotypes")
            pheno_mode_n = self.wk_dtfm.apply(pd.value_counts).max()
        
        n_samples, _ = self.wk_dtfm.shape
        cast_phenos = pheno_mode_n[pheno_mode_n >= threshold * n_samples]
        log_me.info("Phenotypes with low signal: {}".format(", ".join(cast_phenos.index)))

        keep_phenos = [
            idx for idx in self.wk_dtfm.columns 
            if idx in set(self.wk_dtfm.columns) - set(cast_phenos.index)
        ]

        self.wk_dtfm = self.wk_dtfm.loc[:, keep_phenos]
        self.pheno_signals = pheno_mode_n / float(n_samples)

        return self
        
    def into_report(self, dump_phenos=True, dump_ppr=True, dump_pw_corr=False):
        """Write the preprocessed result and report to the disk.

        Args:
            dump_phenos (bool): [False] Whether dump processed phenotype dataframe to the disk
            dump_ppr (bool): [True] Whether dump preprocessing report into disk
            dump_pw_corr (bool): [False] Whether dump pair-wise correlation between phenotypes to disk

        Raises: None

        Return: None
        """
        ppr_path = pjoin(self.wk_dir, "Preprocessing")

        if self.pca_plot is not None:
            path = pjoin(ppr_path, "phenotypes_pca.pdf")
            log_me.info("Dump PCA plot to file: {}".format(path))
            self.pca_plot.savefig(path)

        if self.norm_plots is not None:
            path = pjoin(ppr_path, "phenotypes_normality.pdf")
            log_me.info("Dump histgrams of normality of each phenotype to file: {}".format(path))
            self.norm_plots.savefig(path)
        
        if dump_phenos:
            path = pjoin(ppr_path, "phenotypes_preprocessed.tsv")
            log_me.info("Dump preprocessed phenotypes to file: {}".format(path))
            self.dump_wkdtfm(ppr_path, "phenotypes_preprocessed.tsv")

        if self.pw_corr_plot is not None:
            path = pjoin(ppr_path, "phenotypes_pairwise_correlation.pdf")
            log_me.info("Dump phenotypes pairwise correlation plot to file: {}".format(path))
            self.pw_corr_plot.savefig(path)

        if dump_pw_corr:
            path = pjoin(ppr_path, "phenotypes_pairwise_correlations.tsv")
            log_me.info("Dump phenotype pairwise correlations to file: {}".format(path))
            self.pw_corr.to_csv(path, sep="\t", header=True, index=True)

        if dump_ppr:
            report_dict = {
                "The percentage of modes in each phenotypes": self.pheno_signals.to_list(),
                "The processed dataframe has multiple index": self.mltp_index,
                "The {} times interval of each phenotypes".format(self._pad): self.safe_zone.to_dict(),
                "The matrix's transposing times": self.trps_counts,
                "The removed 'columns'": self.rmd_cols,
                "The removed 'rows'": self.rmd_rows,
            }

            path = pjoin(ppr_path, "phenotypes_preprocessing_report.json")
            log_me.info("Dump preprocessing report in JSON format to file: {}".format(path))
            with open(path, "w") as opth:
                json.dump(report_dict, opth, indent=4)

        return self


if __name__ == "__main__":
    not_runnable()
