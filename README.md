
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ParallelForensicLR <img src="man/figures/logo.png" align="right" height=100/>

<!-- badges: start -->
<!-- badges: end -->

## Introduction

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

## Requirements

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
    - writexl

## A toy MPI example

The `toy_MPI.fam` file included in this repository is a simplified
example representing a basic scenario for genetic screening. It is
structured similarly to real MPI/DVI cases but scaled down for
demonstration purposes. This file includes a **Population Frequency Table
from South America with 70 STR markers, 15 Family Pedigrees, and 15,000
POIs** to be evaluated. Users are encouraged to experiment with this
example to understand the workflow and adapt the format for larger or
more complex cases.

## Usage Instructions

1.  **Ensure Dependencies are Installed**: Before running the script,
    confirm all required R packages are installed. You may use the
    following command in your R console to install them:

``` r
install.packages(c("pedtools", "forrel", "pedmut", "tidyverse", "doParallel", "svDialogs", "matrixStats", "writexl")
```

2.  **Sample Data**: To perform the calculations, it is necessary to
    have a pre-configured `Familias3` file containing the _Frequency
    Table_, _Family Pedigrees_, and _POIs_.

3.  **Run the Script**: The script automatically detects the `Familias3`
    file in the current directory. Upon execution, all necessary
    libraries are loaded, and functions are initialized. The script will
    prompt the user twice:

    - Core Selection: Specify the number of CPU cores to use (the script
      will indicate the maximum available).
    - Pedigree Selection: Choose whether to compare against all
      available pedigrees or only specific ones.

The script includes a **Mendelian Consistency Check** function for the
*Family Pedigrees*, which verifies the genetic inheritance patterns
within each pedigree to ensure no inconsistencies are present.

Additionally, the **Mutation Rate** can also be set or removed directly
within the script if required for specific analyses.

After these selections, the parallelized genetic comparisons are run
according to the chosen parameters.

4.  **Output**: Once calculations are complete, results are saved in an
    Excel file, with POIs sorted by descending likelihood ratio (LR) for
    easier analysis. This output format allows for quick identification
    of the most likely familial matches based on LR values.

## Example Output

| POI       | Sex | FAM1     | FAM2 | FAM3 | FAM4 | FAM5     | …   |
|-----------|-----|----------|------|------|------|----------|-----|
| poi_13863 | F   | 0        | 0    | 0    | 0    | 2.13E+58 | …   |
| poi_04493 | M   | 4.46E+52 | 0    | 0    | 0    | 0        | …   |
| …         | …   | …        | …    | …    | …    | …        | …   |

## Limitations

This tool provides the likelihood ratio (LR) per pedigree and indicates
the sex of each POI, it currently lacks some of the advanced output
features available in `Familias3`. Specifically, the script does not
include details on the number of markers analyzed or any inconsistencies
detected within the _Family Pedigrees_ and _POIs_.

## Future Tasks

To enhance the usability and accessibility of this script, future
development will focus on integrating it into a **Shiny application**. This
Shiny app will provide an intuitive graphical user interface (GUI) that
simplifies data input, parameter selection, and result visualization. By
transitioning the script to a Shiny app, users will be able to interact
with the tool more easily without requiring knowledge of R code, making
it accessible to a broader audience in forensics genetic fields.
