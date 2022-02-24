
### Analysis code for the publication
## Computational drug repurposing against SARS-CoV-2 reveals plasma membrane cholesterol depletion as key factor of antiviral drug activity

Code repository for [Computational drug repurposing against SARS-CoV-2 reveals plasma membrane cholesterol depletion as key factor of antiviral drug activity](https://www.biorxiv.org/content/10.1101/2021.09.10.459786v1) on bioRxiv.
### Abstract

Comparing SARS-CoV-2 infection-induced gene expression signatures to drug treatment-induced gene expression signatures is a promising bioinformatic tool to repurpose existing drugs against SARS-CoV-2. Based on the general hypothesis of signature based drug repurposing, drugs with inverse similarity to a disease signature can reverse disease phenotype, thus can be effective against the investigated disease. However, in case of viral infection diseases, like SARS-CoV-2, infected cells also activate adaptive, antiviral pathways, so the relationship between effective drug and disease signature can be more ambiguous.
To address this question, we analysed gene expression data from in vitro SARS-CoV-2 infected cell lines, and gene expression signatures of drugs showing anti-SARS-CoV-2 activity. Our extensive functional genomic analysis showed that both infection and treatment with in vitro effective drugs leads to activation of antiviral pathways like NFkB and JAK-STAT. Based on the similarity - and not inverse similarity - between drug and infection-induced gene expression signatures, we were able to predict the in vitro antiviral activity of drugs. We also identified SREBF1/2, key regulators of lipid metabolising enzymes, as the most activated transcription factors by several in vitro effective antiviral drugs. Using a fluorescently labeled cholesterol sensor, we showed that these drugs decrease the cholesterol levels of plasma-membrane. Supplementing drug-treated cells with cholesterol reversed the in vitro antiviral effect, suggesting the depleting plasma-membrane cholesterol plays a key role in virus inhibitory mechanism.
Our results can help to more effectively repurpose approved drugs against SARS-CoV-2, and also highlights key mechanisms behind their antiviral effect. 

<br/>
<p align="center">
<img style="" src="https://github.com/comp-sys-pharm/SARS-CoV-2-cholesterol/raw/main/figures/schematic/sars_chol_biorender_final.png" alt="Graphical abstract" width="50%"/></p>


<p align="center">

*Graphical abstract: Schematic figure of the hypothesis that antiviral drugs block virus entry into cells by cholesterol depletion from plasma membrane, and are leading to a compensatory increased SREBF1/2 activity.*

</p>
</br>

### Description of the analysis

- **Preparation of transcriptomic data of virus infection and drug treatment**
  -  1_1_Transcriptomics_data_preparation.ipynb
  -  1_2_Drug_signatures_data_preparation.ipynb
- **Functional genomic analysis and investigation of similarities of drug treatment and virus infection-induced signatures**
  -  2_Functional_analysis_of_virus_and_drug_signatures.ipynb
  -  3_Calculation_of_similarities.ipynb
  -  4_Prediction_with_machine_learning.ipynb

- **Investigation of effect on plasma membrane cholesterol of SREBF activating drugs**
  -  5_Effect_of_SREBF_activating_drugs.ipynb

- **Analysis of cholesterol sensor transfected cell lines**
  -  Cholesterol_sensor_analysis.md

- **CARNIVAL on A549 and CALU-3 cell lines after SARS-CoV-2 infections**
  -  Carnival_virus_signatures.md

- **Additional analysis with Whole genome Influenza A CRISPR screen and LINCS CRISPR signatures**
  - 6_Additional_analysis_CRISPR_KO_InfluenzaA_similary
