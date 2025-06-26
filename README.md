<!-- README.md is generated from README.Rmd. Please edit that file -->

# 🧬 ParallelForensicLR

<!-- badges: start -->
<!-- badges: end -->

## 🔍 Introduction

In the field of **Forensic Genetics**, mainly in *Disaster Victim
Identification (DVI)* and *Missing Persons Identification (MPI)*,
current open-source software solutions face significant limitations when
performing genetic screening, particularly with large-scale datasets.
These tools often struggle to handle the volume and complexity of
massive short tandem repeat (**STR**) and single nucleotide polymorphism
(**SNP**) comparisons essential for identifying individuals in complex
cases.

To address this, we developed **ParallelForensicLR** that leverages
automated ***parallelization technology***, a powerful feature that
enables the simultaneous use of multiple processors on a single machine.
By distributing genetic comparison tasks across multiple cores, this
approach significantly accelerates processing times, making high-volume
genetic screening not only feasible but efficient. The speed of results
will ultimately depend on the machine’s resources, but this parallelized
approach maximizes computational power to deliver faster, more effective
outcomes in genetic searches for missing persons.

## ⚙️ Requirements

This script relies on a `Familias3` (.fam) file containing essential
data such as the *Population Frequency Table* (for STR or SNP markers),
*Family Pedigrees*, and the *Persons Of Interest* (POIs) to be evaluated
for potential familial relationships.

- `Familias3` file with the _Population Frequency Table_ (STR/SNP) and the
  Ante-Mortem (AM) and Post-Mortem (PM) genetics profiles.
  To install the required software, go to this link: [Familias3](https://familias.no/)

- R version 4.2.0 or higher.
- R packages:
  - **Forensic genetics**:
    - pedtools
    - forrel
    - pedmut
  - **Paralelization**:
    - doParallel
  - **Data cleaning**:
    - tidyverse
  - **User interaction**:
    - svDialogs
  - **Output data**:
    - matrixStats

## ▶️ Usage Instructions

**Sample Data**: 
To perform the calculations, it is necessary to
    have a pre-configured `Familias3` file containing the _Frequency
    Table_, _Family Pedigrees_, and _POIs_.

**Run the Script**: 

1.    Make sure all necessary libraries are installed and loaded before execution
``` r
install.packages(c("pedtools", "forrel", "pedmut", "tidyverse", "doParallel", "svDialogs", "matrixStats")
``` 
2.  The script begins by automatically detects the `Familias3`
     file in "data" folder.
3.  Once is loaded, there is a precheck step where
    -   its been evaluated data-integrity check using `connectedPed()`,
        ensuring that every pedigree forms one connected component.
    -   standardises the input with `pedFormat()` and `pedNames()`, harmonising order
        and renaming _POIs_ and _Family Pedigrees_.
    -   verifies the genetic inheritance patterns within each _Family Pedigrees_
        to ensure no inconsistencies are present by runing `mendelian()`.
4.  The script will prompt the user twice:
       - **POI Selection**: Decide whether to analyse **all Persons of
         Interest** in the input file or only a specific subset.
    
       - **Pedigree Selection**: Choose whether to compare against **all
         available family pedigrees** or restrict the analysis to one or
         several selected families.
5.   Finally, set the main functions and parameters:
    - **Inconsistency Filter**: Enter the *maximum number of Mendelian
         inconsistencies* (e.g. `3`) permitted between a _POI_ and a _Family Pedigree_
         for the LR calculation. POIs exceeding this threshold for a given
         family are skipped, accelerating the analysis and reducing noise.
  ``` r
  computeExclusions(pm, am, maxIncomp = 3, ncores = detectCores() - 1)
  ```
  - **Parallelized Genetic Comparisons**:
  The `ParallelForensicLR()` function is the core engine of the workflow. 
  It takes the filtered list of _POI_ and _Family Pedigree_, splits the POIs into chunks,
  and distributes LR calculations across multiple CPU cores for maximum efficiency. 
  For each POI–Family pair that passes the inconsistency filter, it computes the likelihood ratio.
  ``` r
  ParallelForensicLR(pm, am, exclusions_list, ncores = detectCores() - 1, chunk_size = 500)
  ```
6.    Once calculations are complete, results are saved in an
        _csv file_, with POIs sorted by descending likelihood ratio (LR) for
        easier analysis. This output format allows for quick identification
        of the most likely familial matches based on LR values.
  
        The results are then collated into a single, ordered data frame containing:
      - **POI**: Identifier of the Person of Interest  
      - **Sex**: Sex of the POI, inferred via `getSexID()`  
      - **One column per family**: Formatted as `LR_total (nMarkers)`  

## 🧪 A toy MPI example

The `toy_MPI.fam` file included in this repository is a simplified
example representing a basic scenario for genetic screening. It is
structured similarly to real MPI/DVI cases but scaled down for
demonstration purposes. This file includes a **Population Frequency Table
from South America with 70 STR markers, 15 Family Pedigrees, and 15,000
POIs** to be evaluated. Users are encouraged to experiment with this
example to understand the workflow and adapt the format for larger or
more complex cases.

## 📊 Example Output

| POI       | Sex | FAM1          | FAM2 | FAM3 | FAM4 | FAM5          | …   |
|-----------|-----|---------------|------|------|------|---------------|-----|
| poi_13863 | F   | 0             | 0    | 0    | 0    | 2.13E+58 (70) | …   |
| poi_04493 | M   | 4.46E+52 (70) | 0    | 0    | 0    | 0             | …   |
| …         | …   | …             | …    | …    | …    | …             | …   |

## 🚧 Future Tasks

To enhance the usability and accessibility of this script, future
development will focus on integrating it into a **Shiny application**. This
Shiny app will provide an intuitive graphical user interface (GUI) that
simplifies data input, parameter selection, and result visualization. By
transitioning the script to a Shiny app, users will be able to interact
with the tool more easily without requiring knowledge of R code, making
it accessible to a broader audience in forensics genetic fields.
