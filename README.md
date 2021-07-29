<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

### Code for the manuscript
# Patient stratification reveals the molecular basis of disease comoridities

medRxiv: <a href="https://https://doi.org/10.1101/2021.07.22.21260979">https://doi.org/10.1101/2021.07.22.21260979</a>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#manuscript">Manuscript</a>
    </li>
    <li>
      <a href="#web-application">Web Application</a>
    </li>
    <li><a href="#code">Code</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## Manuscript

Study of origin: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66117">GSE66117</a>

- Number of cases: 43
- Number of controls: 3
- Tissue: Blood

Raw counts where downloaded from <a href="http://www.ilincs.org/apps/grein/">GREIN</a>



<!-- GETTING STARTED -->
## Web Application

### Methods

First, samples with a percentage of aligned reads to the genome lower than 70% were removed, as well as studies with less than 3 cases (from now on patients) and control samples meeting the mentioned requirement. Secondly, and in order to perform the analysis at the disease level, gene counts and metadata for each disease were integrated (only studies with the disease, tissue and disease state information were considered). We performed quality controls using the edgeR pipeline (Robinson et al., 2009) and we applied within-sample normalization by considering the logarithm of the counts-per-million (log2CPM). Afterwards, we filtered out lowly expressed genes (those with less than 1 log2CPM in more than 20% of the samples) and we applied between-sample normalization using the trimmed mean of M values (TMM) method (Robinson and Oshlack, 2010). After performing batch effect identification, we used the limma pipeline (Ritchie et al., 2015) for differential expression analysis. Specifically, we built a model considering sample type (case vs. control) as our outcome of interest and adjusting for the study effect, as it is the most descriptive independent variable (tissue, platform and others depend on the study of origin). Genes with an FDR <=0.05 were considered significantly differentially expressed genes (sDEGs). Moreover, we used Combat (Johnson et al., 2007) and QR Decomposition (Ritchie et al., 2015) batch effect removal methods to check the clustering of the samples with t-Distributed Stochastic Neighbor embedding (tSNEs) (Van Der Maaten and Hinton, 2008).

### Results

File:  ChronicLymphocyticLeukemia_DEGs.txt 

Columns of interest:
- symbol: gene symbol 
- logFC: logFC of the mean expression in cases versus controls
- adj.P.Val: adjusted p-value

Genes with and adjusted p-value <0.05 can be considered significantly Differentially Expressed Genes (sDEGs).

Genes with a logFC > 1 are more expressed in the cases than in the controls (overexpressed) whereas genes with a logFC<1 are less expressed in the cases than in the controls (underexpressed). 

## Code

### Methods

To measure gene expression variability - the dispersion of expression of a given gene in a single condition- we used the Distance to Median (DM) (Newman et al., 2006). Since we wanted to capture the changes of gene expression variability associated with diseases (comparing patients vs. controls) and DM takes values in ℝ (including zero), we defined the differential variability value for each gene and disease as the ∆DM; that is, the difference between the DM of the patients and the DM of the controls. We used randomizations to establish the significance of the differential variability between patients and controls for each gene in a given disease. For each disease, we generated a ΔDM null model by sampling with replacement 100000 times patients and controls from the entire set of samples (maintaining the number of patients and controls respectively) and computing the ΔDM value for each gene. We computed the p-values by comparing the observed ΔDM of each gene with the null model. After correcting for multiple testing (FDR <= 0.05), we obtained the significantly Differentially Variable Genes (sDVGs).

### Results

File:  ChronicLymphocyticLeukemia_DVGs.txt 

Columns of interest:
- rn: gene symbol 
- dm_dif: delta DM (ΔDM)
- FDR_corrected: adjusted p-value
