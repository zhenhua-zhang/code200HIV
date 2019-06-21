# Two hundred HIV project


## Plasmarkers
Data source:
	- 200HIV_Plasmamarkers_20190528_YANG.csv

Preprocessing by Wouter:
	- None

Covariates:
	- Age
	- Gender
	- BMI (remove)
	- Smoking(packYears, remove)

Transformed subjects by `log2`
	- AAT
	- Adiponectin
	- CRP
	- IL-18
	- IL-18bpa
	- Leptin
	- Resistin
	- IFABP
	- sCD163
	- sCD14
	- IL-10
	- IL-1RA
	- TNF-A
	- d.dimer
	- CCL5
	- PF4
	- CXCL7
	- Beta_TG
	- TAT2
	- vWf
	- SerotoninPLT
	- IL6_LOD
	- BAFF_Result

Transformed subjects by `Inverse-rank`
	- IL-1A
	- IL1b_LLOD

Discarded subjects: 


## HIV reservoir
Data source:
	- 20190524_HIVreservoir_GENT.csv

Preprocessing by Wouter:
	- `log10` tranformed but still with an outlier

Variable:
	- log10(RNA)
	- log10(DNA)
	- log10(RNA/DNA)

Covariates:
	- Age
	- Gender
	- CD4_NADIR

Transformed subjects by `Inverse-rank`

Discarded subjects: 


## Innate Cytokine (24h)
Data source:
	- B_cytokineData200Hiv24h_boxCorData_log10_2019-03-20.csv

Preprocessing by Wouter:
	- `log10` transform and Box correction

Covariates:
	- Age
	- Gender
	- % monocytes
	- seasonnality/linear time variable
	- BMI (remove)
	- Smoking(packYears, remove)

Transformed subjects by `Inverse-rank`
	- C.gattii_PBMC_24h_IL.1b
	- C.gattii_PBMC_24h_IL.6
	- C.gattii_PBMC_24h_TNF.a
	- E.coli_PBMC_24h_TNF.a
	- Influenza_PBMC_24h_IL.1b
	- Influenza_PBMC_24h_IL.6
	- Influenza_PBMC_24h_TNF.a
	- LPS.1ng_PBMC_24h_TNF.a
	- M.tuberculosis_PBMC_24h_TNF.a
	- OxLDL.LPS_PBMC_24h_IL.1b
	- OxLDL.LPS_PBMC_24h_IL.6
	- OxLDL.LPS_PBMC_24h_TNF.a
	- P.gingivalis_PBMC_24h_IL.6
	- P.gingivalis_PBMC_24h_TNF.a
	- Poly.I.C_PBMC_24h_IL.1b
	- Poly.I.C_PBMC_24h_IL.6
	- Poly.I.C_PBMC_24h_TNF.a
	- RPMI.Serum_PBMC_24h_IL.1b
	- RPMI.Serum_PBMC_24h_IL.6
	- RPMI.Serum_PBMC_24h_TNF.a
	- S.pneumoniae_PBMC_24h_TNF.a
	- S.typhi_PBMC_24h_IL.1b
	- S.typhimurium_PBMC_24h_IL.1b

Discarded subjects: 
	- None


## Stimilated Cytokine (7d)
Data source:
	- C_cytokineData200Hiv7d_boxCorData_log10_2019-03-20.csv

Preprocessing by Wouter:
	- `log10` transform and Box correction

Covariates:
	- Age
	- Gender
	- seasonality/linear time variable
	- CD4P_T_cells_LMI1_grandparentPerc\_\_Index7 + CD8P_T_cells_LMI1_grandparentPerc\_\_Index7

	- BMI (remove)
	- Smoking (packyears, remove)

Transformed subjects by `Inverse-rank`
	- C.albicans.hyphae_PBMC_7d_IFNg
	- C.albicans.hyphae_PBMC_7d_IL.17
	- C.albicans.hyphae_PBMC_7d_IL.22
	- C.albicans.yeast_PBMC_7d_IFNg
	- C.albicans.yeast_PBMC_7d_IL.17
	- C.albicans.yeast_PBMC_7d_IL.22
	- C.gattii_PBMC_7d_IFNg
	- C.gattii_PBMC_7d_IL.17
	- C.gattii_PBMC_7d_IL.22
	- Imiquimod_PBMC_7d_IFNg
	- M.tuberculosis_PBMC_7d_IFNg
	- M.tuberculosis_PBMC_7d_IL.17
	- M.tuberculosis_PBMC_7d_IL.22
	- RPMI.Serum_PBMC_7d_IFNg
	- S.aureus_PBMC_7d_IFNg
	- S.aureus_PBMC_7d_IL.17
	- S.aureus_PBMC_7d_IL.22
	- S.enteritidis_PBMC_7d_IFNg
	- S.enteritidis_PBMC_7d_IL.17
	- S.enteritidis_PBMC_7d_IL.22
	- S.pneumoniae_PBMC_7d_IFNg
	- S.pneumoniae_PBMC_7d_IL.17
	- S.pneumoniae_PBMC_7d_IL.22
	- S.typhi_PBMC_7d_IFNg
	- S.typhi_PBMC_7d_IL.17
	- S.typhi_PBMC_7d_IL.22
	- S.typhimurium_PBMC_7d_IFNg
	- S.typhimurium_PBMC_7d_IL.17
	- S.typhimurium_PBMC_7d_IL.22


## Anti-inflammatory Cytokine
- Data source:
	- L_antiInflamCytokineData200Hiv_boxCorData_log10_2019-03-20.csv

- Preprocessing by Wouter:
	- `log10` transform and Box correction

- Covariates:
	- Age
	- Gender
	- seasonality/linear time variable
	- BMI (remove)
	- Smoking (packyears, remove)

- Transformed subjects by `Inverse-rank`
	- None

- Discarded subjects


## Chemokine Cytokine
- Data source: 
	- M_chemokineCytokineData200Hiv_boxCorData_log10_2019-03-20.csv

- Preprocessing by Wouter:
	- `log10` transform and Box correction

- Covariates:
	- Age
	- Gender
	- seasonality/linear time variable
	- BMI (remove)
	- Smoking (packyears, remove)

- Transformed subjects by `Inverse-rank`

- Discarded subjects



## Overview of Covariates
  - In the metadata
	- age
	- gender
	- BMI (remove)
	- packYears (smoking, remove)
	- smoking (binary)
	- numDaysFromJan2015 (linear time)
	- ssnlt_sin (seasonality sin part, sin(2*pi*numDaysFromJan2015/365))
	- ssnlt_cos (seasonality cos part, cos(2*pi*numDaysFromJan2015/365))
	- prcnt_mnct (percentage of monocytes, cell counts, 200HIV50FG_SysmexPBMC_curated20190227.csv)
	- D4P_TC_LMI1 (D4+ T cell, cell counts)
	- D8P_TC_LMI1 (D8+ T cell, cell counts)

  - Unknown source
	- HIV_DURATION


