
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Research Article

This repository contains the code and data used to carry out analyses
for this research article:

**Doré et al., 2023 - Mutualistic interactions shape global spatial
congruence and climatic niche evolution in Neotropical mimetic
butterflies**

***DOI will be provided once the paper is online***

Müllerian mimicry, where defended prey species converge on similar
warning signals, sharing the cost of educating predators, was first
described nearly 150 years ago. However, the implications of this type
of mutualism on large-scale biodiversity patterns and species niche
evolution remain little explored. Here, we show that Müllerian mimicry
is key to the maintenance of biodiversity in the c. 400 species of the
Neotropical butterfly tribe Ithomiini. Mimicry drives strong spatial
association among species and has channeled the convergence of their
climatic niches. While shared climatic niches may limit community
disassembly, highly adaptively-assembled communities tied by positive
interactions remain particularly vulnerable to extinction cascades,
where the loss of one species could result in the loss of its
mutualistic partners.

All content is available on
[GitHub](https://github.com/MaelDore/ithomiini_convergence) and on
[Zenodo](https://doi.org/10.5281/zenodo.6277769)
([![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6277769.svg)](https://doi.org/10.5281/zenodo.6277769)).

## Contents

-   [:file_folder: **figures_for_article**](figures_for_article/)
    directory contains high quality PDF and TIF versions of the
    Graphical Abstract and Figures displayed in the main text of the
    article.

-   [:file_folder: **functions**](functions/) directory contains all the
    homemade functions called in the scripts during the analyses.

-   [:file_folder: **graphs**](graphs/) directory contains all the
    graphs generated by the scripts during the analyses.

-   [:file_folder: **input_data**](input_data/) directory contains the
    raw and transformed data used in the analyses. Sub-folders include
    phylogeny from [Chazot et al.,
    2019](https://doi.org/10.1111/geb.12919), and species distribution
    maps from [Doré et al., 2021](https://doi.org/10.1111/ddi.13455).

-   [:file_folder: **maps**](models/) directory contains maps generated
    to evaluate community structure.

-   [:file_folder: **outputs**](outputs/) directory contains all files
    generated by the scripts that are not maps, graphs, or tables.

-   [:file_folder: **packages**](packages/) directory contains tar
    archives for packages used in the scripts.

-   [:file_folder: **renv**](renv/) directory contains the full library
    of packages used for this analyses to ensure reproductibility.

-   [:file_folder: **scripts**](scripts/) directory contains the scripts
    used to run the analyses

-   [:file_folder: **supplementaries**](supplementaries/) directory
    contains the outputs used to generate the Figures for the
    Supplementary Materials.

-   [:file_folder: **tables**](tables/) directory contains the outputs
    used to generate the Tables for the article and the Supplementary
    Materials.

## How to run it

This research has been developed using the statistical programming
language R. To run the analyses, you will need installed on your
computer the [R software](https://cloud.r-project.org/) itself and
optionally [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

You can download the entire project as a `.zip` from [this
URL](/archive/master.zip). After unzipping:

-   Open the `ithomiini_convergence.Rproj` file, found at the root of
    the project, in RStudio

-   Run sequentially the scripts found in the [:file_folder:
    **scripts**](scripts/) folder. It will rebuild the outputs, maps,
    graphs, figures and tables, including the final ones presented in
    the main text of the article.

## How to cite

Please cite this research article as:

> Doré, M., Willmott, K., Lavergne, S., Chazot, N., Freitas, A. V. L.,
> Fontaine, C., & Elias, M. (2023). Mutualistic interactions shape
> global spatial congruence and climatic niche evolution in Neotropical
> mimetic butterflies.

***Full reference will be provided once the article is published.***

## Associated online archives

The mimicry classification of ithomiine subspecies used in this study is
publicly available at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5497876.svg)](https://doi.org/10.5281/zenodo.5497876).

The occurrences data used in this study is publicly available at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4696055.svg)](https://doi.org/10.5281/zenodo.4696055).

The distribution maps used in the analyses are publicly available at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4673446.svg)](https://doi.org/10.5281/zenodo.4673446)
