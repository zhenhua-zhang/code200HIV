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
from utls import not_runnable, log_me, not_implemented


class Reporter(PreProc):
    """A class to generate QTL mapping reports.

    Attributes:
    Methods:
    InherenteFrom:
    Examples:
    Todos:
        1. How to properly override parent's methods.
    Notes:
    """

    def __init__(self, wk_dir, qtls, gntp=None, gntp_ptrn=None, gtif=None,
                 gtif_ptrn=None, pntp=None, cvrt=None, anot=None):
        """Initilize a instance of Reporter.
        Args:
        Return:
            self: the instance itself
        """
        self.wk_dir = wk_dir
        self.qtls = qtls
        self.gntp = gntp
        self.gntp_ptrn = gntp_ptrn
        self.gtif = gtif
        self.gtif_ptrn = gtif_ptrn
        self.pntp = pntp
        self.cvrt = cvrt
        self.anot = anot

        super(Reporter, self).__init__(self.qtls)

    def create_meta_db(self, gntp: Genotype):
        """Construct meta database"""
        CHR, POS, PVAL, PNTP, SNPID = ["SequenceName", "Position", "pvalue", "gene", "snps"]

        qtls_dtfm = copy.deepcopy(self.get_dtfm("qtls").dropna(axis=1, how="all"))
        gtif_dtfm = copy.deepcopy(gntp.get_dtfm("info").dropna(axis=1, how="all"))
        gtif_dtfm[SNPID] = gtif_dtfm.index
        qtls_dtfm = qtls_dtfm.merge(gtif_dtfm, how="inner", on=SNPID)
        self.qtls_gntp_info_dtfm = qtls_dtfm.merge(gtif_dtfm, how="inner", on=SNPID)

        return self

    def plot_qq(self, pval_btm: float = 0.5, n_quant: int = 100):
        """Draw Q-Q plot for the QTLs.

        Args:
            pval_btm (float): The top boundary of p-value for the interpolation. Default: 0.5
            n_quant (int): Number of quantiles will be used to construct the interpolated space.  Default: 100
            fig_name (str): The path (dir/file-name) plot to be written into. Default: None
        Returns:
            self
        Todos:
            1. Clean up non-informative comments.
            2. Get options functional.
        """
        log_me.info("Plotting Q-Q plot")

        self.load_files(self.qtls, dtfm_sig="qtls", sep=",", index_col=None)
        dtfm = copy.deepcopy(self.get_dtfm("qtls"))
        dtfm_grouped = dtfm.groupby("gene")

        for pntp_name, pntp_idx in dtfm_grouped.groups.items():
            fig = plt.figure(figsize=(10, 10))
            axis = fig.subplots()

            th_line = [0, -log10(dtfm.loc[:, "pvalue"].min())]
            axis.plot(th_line, th_line, linestyle="--", linewidth=0.1)

            pntp_dtfm = dtfm.loc[pntp_idx, :]
            pval_n = len(pntp_idx)
            pvalue = copy.deepcopy(pntp_dtfm.loc[:, "pvalue"])

            # FIXME: The pvalue_tho impact the pvalue_obs a lot!! BUT don't know why!!!
            # FIXME: There are infinite values, which covers the most significant pvalue probablly
            th_pval = np.concatenate(
                [linspace(0, pval_btm, n_quant/2), logspace(-log10(pval_n), 0, n_quant/2)]
            )
            ob_pval = mstats.mquantiles(pvalue, prob=th_pval, alphap=0, betap=1, limit=(0, 1))

            th_pval_log10, ob_pval_log10 = -log10(th_pval), -log10(ob_pval)
            axis.scatter(th_pval_log10, ob_pval_log10, s=5)

            slope, interc, r_val, p_val, stderr = linregress(ob_pval, th_pval)
            xy = sorted(list(zip(th_pval_log10, ob_pval_log10)), key=lambda x: x[0])[-2]
            axis.annotate("{}\n(R: {:0.3f}, p: {:0.3f})".format(pntp_name, r_val, p_val), xy)

            fig_name = pjoin(self.wk_dir, "ManhattanPlots", "QQPlot_{}.pdf".format(pntp_name))

            log_me.info("Check {} for Q-Q plot".format(fig_name))
            fig.savefig(fig_name)

        return self

    def plot_manhattan(self, minpval: float = 0.05):
        """Draw a manhattan plot.
        Args:
            minpval (float; optional): The minmum pvalue of SNPs will be painted on the Manhattan plot. Default: 0.05
        Raises:
        Returns:
            self
        """
        log_me.info("Plot Manhattan plot")

        CHR, POS, PVAL, PNTP, SNPID = ["SequenceName", "Position", "pvalue", "gene", "snps"]
        TOP_ONLY = False
        GWSIGN = -log10(5e-8)  # Genome-wide significant line
        SGSIGN = -log10(5e-6)  # Genome-wide suggesstive line

        if "qtls" not in self.dtfm_bin:
            self.load_files(self.qtls, dtfm_sig="qtls", sep=",", index_col=None)

        if "gtif" not in self.dtfm_bin:
            self.load_files(self.gtif, self.gtif_ptrn, dtfm_sig="gtif", sep=" ")

        qtls_dtfm = self.get_dtfm("qtls").query("pvalue <= {}".format(minpval))
        gtif_dtfm = self.get_dtfm("gtif")

        gtif_dtfm[SNPID] = gtif_dtfm.index
        qtls_gtif_dtfm = pd.merge(qtls_dtfm, gtif_dtfm, how="inner", on=SNPID)
        qtls_gtif_dtfm_grouped = qtls_gtif_dtfm.groupby(PNTP)

        groups = qtls_gtif_dtfm_grouped.groups
        # snps gene statistic pvalue FDR beta SequenceName Position EffectAllele AlternativeAllele
        for name, sub_dtfm_index in groups.items():
            sub_cols = [CHR, POS, PVAL, SNPID] # The SequenceName could be a function argument
            sub_dtfm = copy.deepcopy(qtls_gtif_dtfm.loc[sub_dtfm_index, sub_cols])

            sub_dtfm.sort_values(by=[CHR, POS], inplace=True)
            sub_dtfm.loc[:, POS] = range(sub_dtfm.shape[0])
            sub_dtfm.loc[:, PVAL] = -log10(sub_dtfm.loc[:, PVAL])

            sub_dtfm_grouped = sub_dtfm.groupby(CHR)
            sub_dtfm_groups = sub_dtfm_grouped.groups

            fig = plt.figure(figsize=(16, 9))
            axis = fig.subplots()

            xticks, xlabels = [], []
            for _chrom, _chrom_dtfm_index in sub_dtfm_groups.items():
                chrom_dtfm = sub_dtfm.loc[_chrom_dtfm_index, ]

                # Non-genome-wide significant SNPs
                chr_ngws_dtfm = chrom_dtfm.loc[chrom_dtfm[PVAL] < GWSIGN, ]
                chr_ngws_pos = chr_ngws_dtfm.loc[:, POS]
                chr_ngws_pval = chr_ngws_dtfm.loc[:, PVAL]
                axis.scatter(chr_ngws_pos, chr_ngws_pval, s=2)

                # Genome-wide significant SNPs
                chr_gws_dtfm = chrom_dtfm.loc[chrom_dtfm[PVAL] >= GWSIGN, ]
                if not chr_gws_dtfm.empty:
                    chr_gws_pos = chr_gws_dtfm.loc[:, POS]
                    chr_gws_pval = chr_gws_dtfm.loc[:, PVAL]
                    chr_gws_snpid = chr_gws_dtfm.loc[:, SNPID]
                    axis.scatter(chr_gws_pos, chr_gws_pval, s=10)

                    if TOP_ONLY:
                        topsnp_idx = chr_gws_pval.idxmax()
                        _text = chr_gws_snpid.loc[topsnp_idx]
                        _x = chr_gws_pos.loc[topsnp_idx]
                        _y = chr_gws_pval.loc[topsnp_idx]
                        axis.annotate(_text, xy=(_x, _y), xytext=(_x*1.01, _y*1.01))
                    else:
                        for _, _x, _y, _text in chr_gws_dtfm.values:
                            axis.annotate(_text, xy=(_x, _y), xytext=[_x*1.01, _y*1.01])
                else:
                    log_me.warn("No genome-wide significant SNPs on chromosome {} in phenotype {}".format(_chrom, name))

                xticks.append(chr_ngws_pos.mean())
                xlabels.append(_chrom)

            axis.set_xticks(xticks)
            axis.set_xticklabels(xlabels)

            axis.set_ylabel("-log10(p value)")

            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)

            pntp_pval_max = sub_dtfm.loc[:, PVAL].max()
            if pntp_pval_max + 1 >= GWSIGN: # Draw the genome wide significant line 5e-8
                axis.axhline(GWSIGN, color="red", linewidth=1, linestyle="--")

            if pntp_pval_max + 1 >= SGSIGN: # Draw the genome wide suggesstive line 5e-6
                axis.axhline(SGSIGN, color="blue", linewidth=1, linestyle="--")

            fig_name = pjoin(self.wk_dir, "ManhattanPlots", "ManhattanPlot_{}.pdf".format(name))

            fig.savefig(fig_name)

        return self

    def pntp_per_gntp_plotter(self, qtls_row, gntp_obj, pntp_obj, fig_name=None):
        """The exace painter of `plot_pntp_per_gntp()` method.

        Notes:
            1. The phenotype columns name could be with suffix like
            '_log10' if it's been transformed in preprocessing steps
        """
        pntp = qtls_row.name
        snpid, pval, beta = qtls_row

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

        print(gntp_pntp_dtfm)

        fig = plt.figure(figsize=(10, 10))
        axis = fig.subplots()

        _labels, _values = [], []
        for _key, _idx in gntp_pntp_dtfm.groupby(snpid, axis=0).groups.items():
            _labels.append("{}({})".format(_key, len(_idx)))
            _values.append(gntp_pntp_dtfm.loc[_idx, pntp])

        axis.boxplot(_values, positions=[0, 1, 2], labels=_labels)

        for _x, _y in enumerate(_values):
            axis.scatter([_x] * len(_y), _y, s=5)

        axis.set_title("Phenotype level per genotype for {} in {}".format(snpid, pntp))
        axis.set_xlabel("Genotypes(n), p-value: {:.3E}, beta: {:.3f}".format(pval, beta))
        axis.set_ylabel("Phenotype levels: {}".format(pntp))

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

        if fig_name is None:  # By default it will check the 
            fig_name = os.path.join(
                self.wk_dir, "PhenotypeLevelPlots",
                "phenotypeLevelPerGenotypePlot_{}_{}.pdf".format(pntp, snpid)
            )

        fig.savefig(fig_name)

        return qtls_row

    def plot_pntp_per_gntp(self, snpid=None, pntp=None, gwth: float = 5e-8):
        """Plot dosage for each genotype, currently not possible to specify the SNP id.

        Args:
            gwth (float): The genome-wide significant threshold. Default: 5e-8

        Returns:
            self
        """
        log_me.info("Plotting phenotype level per genotype")

        PVAL, PNTP, SNPID, BETA, CHROM = "pvalue", "gene", "snps", "beta", "SequenceName"

        if "qtls" not in self.dtfm_bin:
            self.load_files(self.qtls, dtfm_sig="qtls", sep=",", index_col=None)

        qtls_dtfm = copy.deepcopy(self.get_dtfm("qtls").dropna(axis=1, how="all"))
        sign_qtls_dtfm = qtls_dtfm.query('{} <= {}'.format(PVAL, gwth))
        sign_qtls_snpids = sign_qtls_dtfm.loc[SNPID, :]

        # Will be time-consuming, as reload genotype (dosage) and information again
        gntp_obj = Genotype(self.wk_dir, self.gntp, self.gntp_ptrn)
        gntp_obj.load_gntp().load_gtif(sep=" ").encode(sign_qtls_snpids)

        gntp_dtfm = gntp_obj.gntp_dtfm
        gtif_dtfm = gntp_obj.get_dtfm(dtfm_sig="gtif")

        pntp_obj = Phenotype(self.wk_dir, self.pntp)
        pntp_dtfm = pntp_obj.load_files(dtfm_sig="pntp").get_dtfm(dtfm_sig="pntp")

        sign_qtls_dosage = gntp_dtfm.loc[sign_qtls_snpids, :]
        sign_qtls_dosage[SNPID] = sign_qtls_dosage.index

        sign_qtls_info = gntp_dtfm.loc[sign_qtls_snpids, :]
        sign_qtls_info[SNPID] = sign_qtls_info.index

        qtls_dosage = sign_qtls_dtfm.merge(sign_qtls_dosage, how="inner", on=SNPID)
        qtls_dosage_info = qtls_dosage.merge(sign_qtls_info, how="inner", on=SNPID)

        qtls_dosage_info_groupby_pntp_and_chrom = qtls_dosage_info.groupby([PNTP, CHROM])
        for group_name, group_idx in qtls_dosage_info_groupby_pntp_and_chrom.groups.items():
            pntp_name, chrom = group_name
            sub_dtfm_by_pntp_chrom = qtls_dosage_info.loc[group_idx, :]


        return self

    def dump_affected_genes(self, anno_path=None, chrom="1", pos=-1,
                            distance=2.5e5, top_snps_only=True,
                            extra_snps=None, opt_name=None, **kwargs):
        """Fetch genes around a loci by distance
        Args:
            anno_path (str): The path to annotation file (GFF/GTF)
            chrom (str): The chromosome.
            pos (int): The position of candidate QTL(s)
            distance (int): Maximum distance to a gene that will be fetched
            top_snps_only (bool): Whether only fetch around genes for top SNPs per chrom.
            extra_snps (list): A list of SNP ids fetch around genes
        Raises:
        Returns:
        To-dos:
            1. A directory to keep the constructed database, if it's not :memory:

        """
        if self.annos_path is None and anno_path is None:
            log_me.warning("No annotation file is supplied!!!")
            return self

        if self.annos_path is None:
            self.annos_path = anno_path

        if "dbfn" not in kwargs:
            _dbfn = ":memory:"
        else:
            _dbfn = kwargs["dbfn"]
            kwargs.pop("dbfn")

        gfdb = gfu.create_db(self.annos_path, dbfn=_dbfn, **kwargs)

        _target = ",".join(
            [
                "seqid", "source", "featuretype", "start", "end", "score",
                "strand", "frame", "attributes"
            ]
        )
        _conditions = chrom, pos - distance, pos + distance, "gene"
        _order_by = ",".join(["seqid", "start"])

        _query = "SELECT {} FROM features ".format(_target)
        _query += "WHERE seqid=={} AND start>={} AND end<={} feature=={} ".format(*_conditions)
        _query += "ORDER BY {}".format(_order_by)

        if opt_name is None:
            _flnm = "genesAroundQTL_{}_{}_{}_{}".format(snpid, chrom, start, end)
            opt_name = pjoin(self.wk_dir, "QTLInformation", _flnm)

        with open(opt_name, 'a') as opth:
            for sql_row in gfdb.execute(_query):
                chrom, source, feature, start, end, score, strand, _, attr = tuple(sql_row)
                id = attr.get("ID")
                name = attr.get("Name")
                if isinstance(name, list):  # Perhaps not
                    name = name[0]
                gene_info = "\t".join([id, name, chrom, source, feature, start, end, score, strand])
            opth.writelines(gene_info)

        return self

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
