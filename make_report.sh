#!/bin/bash

ipd=${1:?The first argument should be the dir including QTL mapping results}
eqtldb=${2:?The second argument should be eQTLGen database(file)}

# qtlifohd="snps,beta,FDR,cytokine,pvalue,statistic,SequenceName,Position,EffectAllele,AlternativeAllele"
qtlifohd="snps,beta,cytokine,pvalue,statistic,SequenceName,Position,EffectAllele,AlternativeAllele"
eqtlifohd="Pvalue	SNP	SNPChr	SNPPos	AssessedAllele	OtherAllele	Zscore	Gene	GeneSymbol	GeneChr	GenePos	NrCohorts	NrSamples	FDR"

echo "${qtlifohd},covariates,eQTL(ensembl_id(Pvalue|AssessedAllele|OtherAllele|Zscore|GeneSymbol|GeneChr|GenePos|FDR))"
for x in $(ls ${ipd}); do
	cvrt=${x##*_}
	qtlifod=${ipd}/${x}/reports/QtlInformation
	topsnpd=${ipd}/${x}/reports/GenotypeLevelPlot

	for topsnp_dsgplt in $(ls ${topsnpd}); do
		snpid=${topsnp_dsgplt/.pdf/}
		snpid=${snpid##*_}
		pntpf=${qtlifod}/${topsnp_dsgplt%_*}.csv
		grep -w ${snpid} ${pntpf} | cut -d',' -f1,2,4- | tr "\n" ","
		echo ${cvrt} | tr "\n" ","
		grep -w ${snpid} ${eqtldb} \
			| cut -f1,5-11,14 \
			| awk '{printf $5"("$1"|"$2"|"$3"|"$4"|"$6"|"$7"|"$8"|"$9");"} END{printf $3"("$1"|"$2"|"$4"|"$5"|"$6"|"$7"|"$8"|"$9")\n"}'
	done

done

