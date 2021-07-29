<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

# Patient stratification reveals the molecular basis of disease comoridities
Beatriz Urda-García<a href="https://orcid.org/0000-0002-3845-5751">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" />
</a>, Jon Sánchez-Valle<a href="https://orcid.org/0000-0001-7959-6326">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" />
</a>, Rosalba Lepore<a href="https://orcid.org/0000-0002-9481-2557">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" />
</a>, Alfonso Valencia<a href="https://orcid.org/0000-0002-8937-6789">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" />
</a>

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
    <li><a href="#frequently-used-terms">Frequently used terms</a></li>
    <li><a href="#code">Code</a></li>
  </ol>
</details>



<!-- MANUSCRIPT INFORMATION -->
## Manuscript

Preprint available in medRxiv at <a href="https://https://doi.org/10.1101/2021.07.22.21260979">https://doi.org/10.1101/2021.07.22.21260979</a>

### Abstract
Epidemiological evidence shows that some diseases tend to co-occur; more exactly, certain groups of patients with a given disease are at a higher risk of developing a specific secondary condition. Despite the considerable interest, only a small number of connections between comorbidities and molecular processes have been identified.

Here we develop a new approach to generate a disease network that uses the accumulating RNA-seq data on human diseases to significantly match a large number of known comorbidities, providing plausible biological models for such co-occurrences. Furthermore, 64% of the known disease pairs can be explained by analysing groups of patients with similar expression profiles, highlighting the importance of patient stratification in the study of comorbidities.

These results solidly support the existence of molecular mechanisms behind many of the known comorbidities. All the information can be explored on a large scale and in detail at <a href="http://disease-perception.bsc.es/rgenexcom/">http://disease-perception.bsc.es/rgenexcom/</a>. 

<!-- WEB APPLICATION -->
## Web Application: rgenexcom

#### rgenexcom (RNA-seq Gene Expression Comorbidities)
Link: <a href="http://disease-perception.bsc.es/rgenexcom/">http://disease-perception.bsc.es/rgenexcom/</a>.


### Description
To facilitate the visualization and exploration of the generated networks, we implemented a web application that displays the Disease Similarity Network (DSN) and the Stratified Similarity Network (SSN) in a dynamic manner. The user can filter the networks by the type of interactions (positive or negative) and by selecting a minimum and maximum threshold for the edge’s weight. Community detection algorithms (greedy modularity optimization or random walks can be applied to the filtered network and interactions involving specific nodes can be filtered and highlighted. Furthermore, the molecular mechanisms behind diseases and disease interactions can be easily inspected and compared.

### Github
rgenexcom has its own Github repository, where all the code is accessible. 

Link: <a href="https://github.com/bsc-life/rgenexcom">https://github.com/bsc-life/rgenexcom</a>.

## Frequently used terms
- Disease Similarity Network (DSN): disease-disease network. Two diseases are connected if their gene expression profiles correlate significantly (positively or negatively)
- Meta-patients: groups of patients from a given disease with a similar gene expression profile. 
- Stratified Similarity Network (SSN): extension of the DSN network thar includes the meta-patients. Hence, the SSN contains three types of links: 
      (1) links between diseases
      (2) links between diseases and meta-patients
      (3) links between meta-patients.

<!-- CODE -->
## Code

### Data Preparation
First, uniformly processed gene counts were dowloaded from the <a href="http://www.ilincs.org/apps/grein/">GREIN platform</a> and Summarized Emperiment Objects (SE) from Bioconductor were constructed for each study. Finally, SE objects corresponding to the same disease where merged (merge_same_disease_se_objects.R).

### RNA-seq pipeline
Then, we applied an RNA-seq pipeline to each disease separately and in parallel (run_rnaseq_pipeline_for_disease.R). Then, we clustered diseases based on their significantly dysregulated pathways (molecular_insight_heatmap.R).

### DSN generation
1. First, we computed distances between diseases (Network_building/build_disease_level_network.py)
2. We used the generated distances to obtain the Disease Similarity Network (DSN) (generating_networks.R)

### DSN overlap
1. Obtain the SE object for each icd9 (generating_SE_objects_icd9_level.R)
2. Run the RNA-seq pipeline at the icd9 level (run_rnaseq_pipeline_for_disease.R)
3. Compute distances between diseases (build_ICD_level_network.R)
4. Use the obtained distances to generate the ICD9 level DSN network (generating_networks.R)
5. Compute the overlap with the epidemiological network from Hidalgo et al. (network_overlap_icd.R)

In the alternative approach, we computed the overlap of the DSN with the epidemiological network from Hidalgo et al. (network_overlap.R) directly from the DSN (by transforming the disease names into ICD9 codes)

### DSN analysis and comparison with other molecular networks
1. Topological analysis of the networks and how they compare to the epidemiological network (analyzing_networks.R)
2. Overlap of other molecular networks with the epidemiological network from Hidalgo et al. (Other_molecular_networks/other_networks_overlap.R)
3. Generating the barplot with the overlaps of the DSN and other molecular networks (barplot_overlap_molecular_networks.R)
4. Comparison of the PPI network with the DSN (Other_molecular_networks/comparison_with_ppi.R)
5. Comparison of the microarrays network with the DSN (Overlap_microarrays.R)

### Molecular mechanisms behind comorbidities
Inspecting the molecular mechanisms behind comorbidities (exploring_underlying_molecular_mechanisms.py and pathways_count_plots.R)

### Meta-patient definition and characterization. 
1. We used PAM and WARD algorithms to define meta-patients for each disease (groups of patients with a similar expression profile) (defining_meta_patients.R)
2. Then, we applied the RNA-seq pipeline for each meta-patient (DEanalysis_for_metapatients.R)

### SSN generation, analysis and overlap computation
1. First, we computed distances between diseases (Network_building/build_metapatient_disease_network.py)
2. We used the generated distances to obtain the Disease Similarity Network (DSN) (generating_networks.R)
3. We computed the overlap of the DSN with the epidemiological network from Hidalgo et al. (network_overlap_SSN.R)
4. Topological analysis of the DSN (analyzing_networks.R)

### Meta-patients increase the detection power
Randomizations were performed to check the significance of the increase in the detection power achieved with the definition of meta-patients (obtaining_random_meta_patients.R, DEA_random_metapatients.R, Network_building/build_metapatient_dis_network_randomization.py, Randomization/obtain_pvalues_increased_power_metapatients.R)














