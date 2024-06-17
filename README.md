<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

# Patient stratification reveals the molecular basis of disease co-occurrences
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
Epidemiological evidence shows that some diseases tend to co-occur; more exactly, certain groups of patients with a given disease are at a higher risk of developing a specific secondary condition. Here we develop a new approach to generate a disease network that uses the accumulating RNA-seq data on human diseases to significantly match an unprecedented proportion of known comorbidities, providing plausible biological models for such co-occurrences and effectively mirroring the underlying structure of complex disease relationships. Furthermore, 64% of the known disease pairs can be explained by analyzing groups of patients with similar expression profiles, highlighting the importance of patient stratification in the study of comorbidities. These results solidly support the existence of molecular mechanisms behind many of the known comorbidities, with most captured co-occurrences implicating the immune system.  Additionally, we identified new and potentially underdiagnosed comorbidities, providing molecular insights that could inform targeted therapeutic strategies. We provide a functional and comprehensive web application to explore diseases, disease co-occurrences and their underlying molecular processes at different resolution levels at <a href="http://disease-perception.bsc.es/rgenexcom/">http://disease-perception.bsc.es/rgenexcom/</a>. 


<!-- WEB APPLICATION -->
## Web Application: rgenexcom

#### rgenexcom (RNA-seq Gene Expression Comorbidities)
Link: <a href="http://disease-perception.bsc.es/rgenexcom/">http://disease-perception.bsc.es/rgenexcom/</a>.

![Captura de pantalla 2021-07-29 a las 17 20 21](https://user-images.githubusercontent.com/46993728/127519898-c808ac1c-34f7-4f2c-b266-3fdd1f5853cb.png)


### Description
To facilitate the visualization and exploration of the generated networks, we implemented a web application that displays the Disease Similarity Network (DSN) and the Stratified Similarity Network (SSN) in a dynamic manner. The user can filter the networks by the type of interactions (positive or negative) and by selecting a minimum and maximum threshold for the edge’s weight. Community detection algorithms (greedy modularity optimization or random walks can be applied to the filtered network and interactions involving specific nodes can be filtered and highlighted. Furthermore, the molecular mechanisms behind diseases and disease interactions can be easily inspected and compared.

### Github
rgenexcom has its own Github repository, where all the code is accessible. 

Link: <a href="https://github.com/bsc-life/rgenexcom">https://github.com/bsc-life/rgenexcom</a>.


## Frequently used terms
- **Disease Similarity Network (DSN):** disease-disease network. Two diseases are connected if their gene expression profiles correlate significantly (positively or negatively).
- **Meta-patients:** groups of patients from a given disease with a similar gene expression profile. 
- **Stratified Similarity Network (SSN):** extension of the DSN network that includes the meta-patients. Hence, the SSN contains three types of links: 
      (1) links between diseases
      (2) links between diseases and meta-patients
      (3) links between meta-patients.


<!-- CODE -->
## Code

### Data Preparation
First, uniformly processed gene counts were dowloaded from the <a href="http://www.ilincs.org/apps/grein/">GREIN platform</a> and Summarized Emperiment Objects (SE) from Bioconductor were constructed for each study. Finally, SE objects corresponding to the same disease were merged (<code>merge_same_disease_se_objects.R</code>).

### RNA-seq pipeline
Then, we applied an RNA-seq pipeline to each disease separately and in parallel (<code>run_rnaseq_pipeline_for_disease.R</code>). Then, we clustered diseases based on their significantly dysregulated pathways (<code>molecular_insight_heatmap.R</code>).

### DSN generation
1. First, we computed distances between diseases (<code>Network_building/build_disease_level_network.py</code>)
2. We used the generated distances to obtain the Disease Similarity Network (DSN) (<code>generating_networks.R</code>)

### DSN overlap
1. Obtain the SE object for each icd9 (<code>generating_SE_objects_icd9_level.R</code>)
2. Run the RNA-seq pipeline at the icd9 level (<code>run_rnaseq_pipeline_for_disease.R</code>)
3. Compute distances between diseases (<code>build_ICD_level_network.R</code>)
4. Use the obtained distances to generate the ICD9 level DSN network (<code>generating_networks.R</code>)
5. Compute the overlap with the epidemiological network from Hidalgo et al. (<code>network_overlap_icd.R</code>)

In the alternative approach, we computed the overlap of the DSN with the epidemiological network from Hidalgo et al. (<code>network_overlap.R</code>) directly from the DSN (by transforming the disease names into ICD9 codes)

### DSN analysis and comparison with other molecular networks
1. Topological analysis of the networks and how they compare to the epidemiological network (<code>analyzing_networks.R</code>)
2. Overlap of other molecular networks with the epidemiological network from Hidalgo et al. (<code>Other_molecular_networks/other_networks_overlap.R</code>)
3. Generating the barplot with the overlaps of the DSN and other molecular networks (<code>barplot_overlap_molecular_networks.R</code>)
4. Comparison of the PPI network with the DSN (<code>Other_molecular_networks/comparison_with_ppi.R</code>)
5. Comparison of the microarrays network with the DSN (<code>Overlap_microarrays.R</code>)

### Molecular mechanisms behind comorbidities
Inspecting the molecular mechanisms behind comorbidities (<code>exploring_underlying_molecular_mechanisms.py</code> and <code>pathways_count_plots.R</code>)

### Meta-patient definition and characterization. 
1. We used PAM and WARD algorithms to define meta-patients for each disease (groups of patients with a similar expression profile) (<code>defining_meta_patients.R</code>)
2. Then, we applied the RNA-seq pipeline for each meta-patient (<code>DEanalysis_for_metapatients.R</code>)

### SSN generation, analysis and overlap computation
1. First, we computed distances between diseases (<code>Network_building/build_metapatient_disease_network.py</code>)
2. We used the generated distances to obtain the Disease Similarity Network (DSN) (<code>generating_networks.R</code>)
3. We computed the overlap of the DSN with the epidemiological network from Hidalgo et al. (<code>network_overlap_SSN.R</code>)
4. Topological analysis of the SSN (<code>analyzing_networks.R</code>)

### Meta-patients increase the detection power
Randomizations were performed to check the significance of the increase in the detection power achieved with the definition of meta-patients (N=1000). 
Steps:
1. For each disease, create 1000 random meta-patients maintaining the number of size of the meta-patients for each disease (<code>obtaining_random_meta_patients.R</code>)
2. Run the RNA-seq pipeline for each meta-patient (<code>DEA_random_metapatients.R</code>)
3. Build 1000 networks with the random meta-patients (<code>Network_building/build_metapatient_dis_network_randomization.py</code>)
4. Use the generated networks to establish the significance of the increased in detection power observed with the original meta-patients (<code>Randomization/obtain_pvalues_increased_power_metapatients.R</code>)



