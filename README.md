qtlmapping
===


Introduction
---
A wrapper of `MatrixEQTL` and some post-mapping scripts. All the scripts are
written to analyze the latent HIV reservoir in 200 HIV project. However, it
should work for other dataset.


Dependencies
---
- `R (>= 3.5)`
- `MatrixEQTL`
- `data.table`
- `ggplot2`
- `optparse`


Usage
---

**CLI example**  
``` {bash}
Rscript ./qtlmapping_v0.2.0.R \
    --draw-pwcor \
    --trps-cvrt-dtfm \
    --trps-pntp-dtfm \
    --work-dir test_run \
    --pntp-file PHENOTYPE-FILE.tsv \
    --target-pntp TARGET-GENOTYPE-NAME \
    --cvrt-file COVARIATES-FILE.tsv \
    --gntp-dosage-file GENOTYPE-DOSAGE-FILE.gz \
    --gntp-info-file GENOTYPE-INFORMATION-FILE.gz
```

**Help**  
The help here could be out of date, as the scripts could be updated frequently
while the `README.md` isn't.

``` {bash}
Usage: qtlmapping_v0.2.0.R [options]
A QTL mapping script based on MatrixEQTL

Options:
        -h, --help
                Show this help message and exit

        -w WORK-DIR, --work-dir=WORK-DIR
                Output direcotry.

        --run-flag=RUN-FLAG
                Running flag which help to discrimnate different runs.

        -p PNTP-FILE, --pntp-file=PNTP-FILE
                Phenotype input file.

        --padding=PADDING
                The times of standard deviation as the boundary of outliers.

        --target-pntp=TARGET-PNTP
                The phenotypes will be used, if more than one, using comma as
                delimiter.

        --trps-pntp-dtfm
                Whether should transpose the data.frame of phenotypes to make
                the columns as sample ID.

        --pntp-idx-col=PNTP-IDX-COL
                The id column in phenotype file

        -c CVRT-FILE, --cvrt-file=CVRT-FILE
                Covariates file.

        --target-cvrt=TARGET-CVRT
                Covariates will be used, if more than one, using comma as
                delimiter.

        --trps-cvrt-dtfm
                Whether should transpose the data.frame of covariates.

        --cvrt-idx-col=CVRT-IDX-COL
                The id column in covariates file

        -t PWCOR-TRAIT, --pwcor-trait=PWCOR-TRAIT
                Traits will be correlated with in paire-wised way. The options
                will be ignored if --draw-pwcor isn't given.

        -d GNTP-DOSAGE-FILE, --gntp-dosage-file=GNTP-DOSAGE-FILE
                The genotype dosage file (could be compressed).

        --genotype-dosage-idx-col=GENOTYPE-DOSAGE-IDX-COL
                The id column in genotype file, usually its the name of column
                of SNP id.
                Default: id

        -i GNTP-INFO-FILE, --gntp-info-file=GNTP-INFO-FILE
                The genotype information file (could be compressed).

        --gntp-info-cols=GNTP-INFO-COLS
                The columns will be used. 
                Default: rsID,SequenceName,Position,EffectAllele,AlternativeAllele

        --maf-thrd=MAF-THRD
                Minor allele frequency.

        --pm-times=PM-TIMES
                How many times of permutations should be done. If it's less than
                1, no permutation will be performed but only 'raw' data will be
                used in the mapping.
                Default: 0

        --pm-seed=PM-SEED
                The random seed for permutation.
                Defautl: 31415

        --mhtn-fig-p-thrd=MHTN-FIG-P-THRD
                The threshold of p-value for Manhattan plot.
                Default: 0.05

        --draw-pwcor
                If the flag is given, the script will draw the pair-wise
                correlation plot for genotypes and covariates.
```


Notice
---


To-do list
---
- [ ] Make it into a R package
- [ ] Do permutation by multiple threads
- [ ] Add more details to the Usage.Help section


Licence
---

`爱(ài)干(gàn)嘛(má)干(gàn)嘛(má)。如果非得要就用MIT就行。`  
Do anything you want with it. If you insist, please use MIT one.
