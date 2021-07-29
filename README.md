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
</a>,

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


<!-- CODE -->
## Code


