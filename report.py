# -*- UTF-8 -*-
"""I'm meant to be empty."""

#
## report.py
#

# Standard libraries
import os
import sys
import copy

from os.path import join as pjoin

# Third party libraries
import scipy as sp
import numpy as np
import pandas as pd
import gffutils as gfu

import matplotlib.pyplot as plt

from numpy import log10, linspace, logspace
from scipy.stats import linregress
from scipy.stats import mstats

# Customer libraries
from preproc import PreProc
from genotype import Genotype
from phenotype import Phenotype
from utls import not_runnable, log_me


class Reporter(PreProc):
    """A class to generate QTL mapping reports.

    Attributes:
    Methods:
    InherenteFrom:
    Examples:
    Todos:
    Notes:
    """

    def __init__(self, wk_dir, annos_path=None):
        """Initilize a instance of Reporter.
        Args:
            qtls:
            annos:
        Return:
            self: the instance itself
        """
        self.wk_dir = wk_dir
        self.qtls_path = os.path.join(wk_dir, "QTLInformation", "qtls.tsv")
        self.annos_path = annos_path

        # This only read one file
        super(Reporter, self).__init__(self.qtls_path)

        self.qqplot = None

    def load_files(self):
        """Load files"""
        super().load_files(dtfm_sig="qtls", index_col=False, sep=",")

        if self.annos_path:
            super().load_files(self.annos_path, dtfm_sig="annotations")

        return self

    def plot_qq(self, all_in_one=True, pval_btm=0.5, n_quant=100):
        """Draw Q-Q plot for the QTLs.
        Args:
            all_in_one (bool; optional): Whether put all Q-Q plot on one figure. Default: True
            threshold (float; optional): The threshold of pvalue of SNPs that will annotated by (pvalue and beta)
        Raises:
        Returns:
            self
        Todos:
            1. Clean up non-informative comments.
            2. Get options functional.

        Reference:
            1. The algorithm was reffered to [assocplots](https://github.com/khramts/assocplots).
        """

        # Fetch the dataframe including QTLs information
        dtfm = copy.deepcopy(self.get_dtfm(0, "qtls"))
        # Group the QTLs information per "gene" (a.k.a. phenotype)
        dtfm_grouped = dtfm.groupby("gene")

        fig = plt.figure(figsize=(10, 10))
        ax = fig.subplots()

        th_line = [0, -log10(dtfm.loc[:, "pvalue"].min())]
        ax.plot(th_line, th_line, linestyle="--", linewidth=0.1)

        # By default the qtls are sorted by pvalue
        for pntp_name, pntp_idx in dtfm_grouped.groups.items():
            # Pick up subdataframe grouped by gene
            pntp_dtfm = dtfm.loc[pntp_idx, :]
            # Length of the `linspace`
            pval_n = len(pntp_idx)
            pvalue = copy.deepcopy(pntp_dtfm.loc[:, "pvalue"])

            # FIXME: The pvalue_tho impact the pvalue_obs a lot!! BUT don't know why!!!
            # FIXME: There is infinite value, which covers the most significant pvalue with high chance
            th_pval = np.concatenate([linspace(0, pval_btm, n_quant), logspace(-log10(pval_n), 0, n_quant)])
            ob_pval = mstats.mquantiles(pvalue, prob=th_pval, alphap=0, betap=1, limit=(0, 1))

            th_pval_log10, ob_pval_log10 = -log10(th_pval), -log10(ob_pval)
            ax.scatter(th_pval_log10, ob_pval_log10, s=5)

            slope, interc, r_val, p_val, stderr = linregress(ob_pval, th_pval)
            xy = sorted(list(zip(th_pval_log10, ob_pval_log10)), key=lambda x: x[0])[-2]
            ax.annotate("{}\n(R: {:0.3f}, p: {:0.3f})".format(pntp_name, r_val, p_val), xy)

        fig.savefig(os.path.join(self.wk_dir, "ManhattanPlots", "QQ_plot.pdf"))

        self.qqplot = ax

        return self

    def plot_manhattan(self, gntp: Genotype):
        """Draw a manhattan plot.
        Args:
        Raises:
        Returns:
            self
        """
        log_me.info("Plot Manhattan...")

        CHR, POS, PVAL, PNTP, SNPID = ["SequenceName", "Position", "pvalue", "gene", "snps"]
        TOP_ONLY = False
        GWSIGN = -log10(5e-8)  # Genome-wide significant line
        SGSIGN = -log10(5e-6)  # Genome-wide suggesstive line

        qtls_dtfm = copy.deepcopy(self.get_dtfm(0, "qtls").dropna(axis=1, how="all"))
        gntp_info_dtfm = copy.deepcopy(gntp.get_dtfm(1, "info").dropna(axis=1, how="all"))

        # In the qtl mapping results, the SNP id could be duplicated, if there're more than one phenotypes, while in genotype and genotype infomation files, SNP id is unique.
        gntp_info_dtfm[SNPID] = gntp_info_dtfm.index  
        qtls_dtfm = qtls_dtfm.merge(gntp_info_dtfm, how="inner", on=SNPID)
        qtls_dtfm_grouped = qtls_dtfm.groupby(PNTP)

        groups = qtls_dtfm_grouped.groups
        # snps gene statistic pvalue FDR beta SequenceName Position EffectAllele AlternativeAllele
        for name, sub_dtfm_index in groups.items():
            sub_cols = [CHR, POS, PVAL, SNPID] # The SequenceName could be a function argument
            sub_dtfm = copy.deepcopy(qtls_dtfm.loc[sub_dtfm_index, sub_cols])

            sub_dtfm.sort_values(by=[CHR, POS], inplace=True)
            sub_dtfm.loc[:, POS] = range(sub_dtfm.shape[0])
            sub_dtfm.loc[:, PVAL] = -log10(sub_dtfm.loc[:, PVAL])

            sub_dtfm_grouped = sub_dtfm.groupby(CHR)
            sub_dtfm_groups = sub_dtfm_grouped.groups

            fig = plt.figure(figsize=(16, 9))
            ax = fig.subplots()

            xticks, xlabels = [], []
            for _chrom, _chrom_dtfm_index in sub_dtfm_groups.items():
                chrom_dtfm = sub_dtfm.loc[_chrom_dtfm_index, ]

                # Non-genome-wide significant SNPs
                chr_ngws_dtfm = chrom_dtfm.loc[chrom_dtfm.loc[:, PVAL]<GWSIGN-2, ]
                chr_ngws_pos = chr_ngws_dtfm.loc[:, POS]
                chr_ngws_pval = chr_ngws_dtfm.loc[:, PVAL]
                ax.scatter(chr_ngws_pos, chr_ngws_pval, s=2)

                # Genome-wide significant SNPs
                chr_gws_dtfm = chrom_dtfm.loc[chrom_dtfm.loc[:, PVAL]>=GWSIGN-2, ]
                if not chr_gws_dtfm.empty:
                    chr_gws_pos = chr_gws_dtfm.loc[:, POS]
                    chr_gws_pval = chr_gws_dtfm.loc[:, PVAL]
                    chr_gws_snpid = chr_gws_dtfm.loc[:, SNPID]
                    ax.scatter(chr_gws_pos, chr_gws_pval, s=10)

                    if TOP_ONLY:
                        topsnp_idx = chr_gws_pval.idxmax()
                        _text = chr_gws_snpid.loc[topsnp_idx]
                        _x = chr_gws_pos.loc[topsnp_idx]
                        _y = chr_gws_pval.loc[topsnp_idx]
                        ax.annotate(_text, xy=(_x, _y), xytext=(_x*1.01, _y*1.01))
                    else:
                        for _, _x, _y, _text in chr_gws_dtfm.values:
                            ax.annotate(_text, xy=(_x, _y), xytext=[_x*1.01, _y*1.01])
                else:
                    log_me.warn("No genome-wide significant SNPs on chromosome {}".format(_chrom))

                xticks.append(chr_ngws_pos.mean())
                xlabels.append(_chrom)

            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabels)

            ax.set_ylabel("-log10(p value)")

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            pntp_pval_max = sub_dtfm.loc[:, PVAL].max()
            if pntp_pval_max >= GWSIGN:
                ax.axhline(GWSIGN, color="red", linewidth=1)  # Draw the genome wide significant line 5e-8

            if pntp_pval_max >= SGSIGN:
                ax.axhline(SGSIGN, color="blue", linewidth=1)  # Draw the genome wide suggesstive line 5e-6

            fig.savefig(pjoin(self.wk_dir, "ManhattanPlots", "ManhattanPlot_{}.pdf".format(name)))

        return self

    def _drawer(self, qtls_row, gntp_obj, pntp_obj, fig_name=None):
        """The exace painter of `plot_pntp_per_gntp()` method.

        Notes:
            1. The phenotype columns name could be with suffix like
            '_log10' if it's been transformed in preprocessing steps
        """
        pntp = qtls_row.name
        chr, pos, snpid, pval, beta = qtls_row

        gntp_row = gntp_obj \
            .encode(snpid) \
            .gntp_dtfm \
            .loc[snpid, :] # the `name` is the snpid

        pntp_row = pntp_obj \
            .wk_dtfm \
            .loc[pntp, :] # the `name` is pntp

        gntp_pntp_dtfm = pd \
            .concat([gntp_row, pntp_row], axis=1, sort=False) \
            .dropna(axis=0)

        fig = plt.figure(figsize=(10, 10))
        ax = fig.subplots()

        _labels, _values = [], []
        for _key, _idx in gntp_pntp_dtfm.groupby(snpid, axis=0).groups.items():
            _labels.append("{}({})".format(_key, len(_idx)))
            _values.append(gntp_pntp_dtfm.loc[_idx, pntp])

        ax.boxplot(_values, positions=[0, 1, 2], labels=_labels)

        for _x, _y in enumerate(_values):
            ax.scatter([_x] * len(_y), _y, s=5)

        ax.set_title("Phenotype level per genotype for {} in {}".format(snpid, pntp))
        ax.set_xlabel("Genotypes(n), p-value: {:.3E}, beta: {:.3f}".format(pval, beta))
        ax.set_ylabel("Phenotype levels: {}".format(pntp))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if fig_name is None:  # By default it will check the 
            fig_name = os.path.join(
                self.wk_dir, "PhenotypeLevelPlots",
                "phenotypeLevelPerGenotypePlot_{}_{}.pdf".format(pntp, snpid)
            )

        fig.savefig(fig_name)

        return qtls_row


    def plot_pntp_per_gntp(self, gntp: Genotype, pntp: Phenotype, gwth: float =5e-8, corr_for: list =None):
        """Plot dosage for each genotype, currently not possible to specify the SNP id.
        Args:
            gntp (Genptype): Genptype object. Required
            pntp (Phenotype): Phenotype object. Required
            gwth (float): The genome-wide significant threshold. Default: 5e-8
            corr_for (list): Correction for any given covariates. Default: None
        Raises:
        Returns:
            self
        """
        log_me.info("Plot phenotype level per genotype")

        CHR, POS, PVAL, PNTP, SNPID, BETA = "SequenceName", "Position", "pvalue", "gene", "snps", "beta"

        # Fetch and merge DataFrame (qtls, info), its duplicated in `plot_manhattan()`
        qtls_dtfm = copy.deepcopy(self.get_dtfm(0, "qtls").dropna(axis=1, how="all"))
        gntp_info_dtfm = copy.deepcopy(gntp.get_dtfm(1, "info").dropna(axis=1, how="all"))
        gntp_info_dtfm[SNPID] = gntp_info_dtfm.index  
        qtls_dtfm = qtls_dtfm.merge(gntp_info_dtfm, how="inner", on=SNPID)

        _ = qtls_dtfm.groupby(PNTP) \
            .min() \
            .query('{} <= @gwth'.format(PVAL)) \
            .loc[:, [CHR, POS, SNPID, PVAL, BETA]] \
            .apply(self._drawer, axis=1, gntp_obj=gntp, pntp_obj=pntp)

    def dump_affected_genes(self, qtl_pos=-1, distance=2.5e5):
        """Fetch genes around a loci by distance
        Args:
            qtl_pos: the position of candidate QTL(s)
            distance: Maximum distance to a gene that will be fetched
        Raises:
        Returns:
        """
        # gene annotations, top QTLs, genotype information (Genotype: C)

    def pathway_enrichment(self, qtls=None):
        """Pathway enrichment analysis for given QTLs.
        Args:
        Raises:
        Returns:
            self
        """
        not_implemented()

    def plot_locuszoom(self):
        """Plot the LocusZoom.
        Args:
        Raises:
        Returns:
            self
        TODO: 1. Not yet implemented. It could be implemented locally, however,
            the official scipt from LocusZoom looks doesn't work. Perhaps pull
            result from the official site helps.
        """
        # QTLs, genome annotations

    def dump_reports(self, path):
        """Dump out reports into given path.
        Args:
            path (str): required. Where the report will be dumpped into.
        Raises: None
        Returns: None
        """


if __name__ == "__main__":
    not_runnable()
