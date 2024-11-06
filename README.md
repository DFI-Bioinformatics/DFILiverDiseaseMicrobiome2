# Data Repo for Fecal butyrate and deoxycholic acid concentrations correlate with mortality in hospitalized patients with liver disease


## Folder structure

1. **code**: contains all R scripts for generating heatmap, volcano plots and/or taxonomy barplots, etc. Scripts are organized by number. Resultant plots are labeled with same number as the script. 
2. **results**: contains all unmodified graphs in PDF format. Graphs published along with the paper have been slightly rearranged and beautified in Adobe Illustrator.
3. **data**: contains raw or derived data or template files
- `LD847_HD22.meta.quant.metabolomics.csv`: the most important file. It contains clinical variables, such as stool consistency, lactulose, Bifidobacterium group (10% is the cutoff), SBP, bacteremia status, or metagenomically derived toxin RPKM values, quant metabolomics readings from SCFA (unit: mM) and bile acid (unit: ug/ml) panels, etc.
- `LD847_HD22.metaphlan.csv`: taxonomy annotation from [metaphlan4](https://github.com/biobakery/MetaPhlAn)
- `LD847_HD22.qual.metabolomics.csv`: a long format of qualitative metabolomics values. There is not unit associated with them.
- `LD847.ratios.quant.metabolomics.csv`: primary to secondary bile acid ratios calculated in unit of mM.
- `taxumap_embedding.LD847.n375.csv`: [taxUmap](https://github.com/jsevo/taxumap) coordinates
- `In vitro and mouse data_for repo_230713.xlsx`: data for *in vitro* and mouse experiment. They are organized by tabs.
- `90day_survival_demographics_230518.csv`: data used for survival analysis and Cox Proportional Harzard model for initial sample of each lactulose exposed patient.
- `qPCR_Calculations_MattO_220708 GF.xlsx`: qPCR 16S copy number for Supp Fig 6, D.
- `MMF.16S.263_MattOdenwald.finalPhy_rdp.rds`: phyloseq object for Supp Fig 6, D. 
- `qual_compounds.csv`: a list of all qualitative compounds to use in the heatmap
- `qual_heatmap_lookup.csv`: a template for all qualitative compounds and their classes
- `mouse_consorita_phyloseq.rds`: phyloseq object for Supp Fig 7
- `mouse_consortia_metab_quant.csv`: quantitative metabolomic data for Supp Fig 7
