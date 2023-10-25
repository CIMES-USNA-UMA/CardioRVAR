  <!-- badges: start -->
  [![R-CMD-check](https://github.com/CIMES-USNA-UMA/CardioRVAR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CIMES-USNA-UMA/CardioRVAR/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->


## Description

CardioRVAR is a tool developed by Alvaro Chao-Ecija (intern student at the 
Department of Physiology and the Autonomic Nervous System Unit at CIMES, University of
Malaga) and PhD. Marc Stefan Dawid-Milner (head of the unit and supervisor of the project). 
The functions provided by this R package allow the analysis of closed-loop 
cardiovascular interactions.

## Installation

To install the package, use the following code line in R (package remotes is required):

```ruby
remotes::install_github("CIMES-USNA-UMA/CardioRVAR", build_vignettes = TRUE)
```

## Shiny application

There is a complementary shiny application available for this package. To access it, you can either install
it and access it using the following commands (you may need to install package *shiny* separately):

```ruby
remotes::install_github("CIMES-USNA-UMA/CardioRVARapp", build_vignettes = TRUE)
CardioRVARapp::StartCardioRVARapp()
```

Or if you prefer use the following code line (you may need to install *shiny*, *ggplot2* and *CardioRVAR*):

```ruby
shiny::runGitHub("CardioRVARapp", "CIMES-USNA-UMA", subdir = "inst/app", launch.browser = TRUE)
```

## Tutorial

A tutorial for this package is available as a vignette. To access the tutorial vignette,
first ensure that packages *knitr* and *rmarkdown* are installed:

```ruby
install.packages("knitr")
install.packages("rmarkdown")
```

Then, install package *CardioRVAR* as described above (with the option to build vignettes turned on as
described). Then, you can use the following commands to access the tutorial:

```ruby
library(CardioRVAR)
vignette("CardioRVAR-Tutorial")

# Or alternatively:

vignette("CardioRVAR-Tutorial", package = "CardioRVAR")
```
This will show the vignette in the *help* section of *RStudio*. If you want to access
the HTML version, use the following command:

```ruby
browseVignettes("CardioRVAR")
```

An access for the vignette *CardioRVAR Tutorial* will be shown. To see the vignette, click on *HTML*.

## Help

To access information about a particular function, use *?* next to the function you wish to
check. For example, the following expression is used to access the "help" guide for function
*GetExpectedValues()*:

```ruby
?GetExpectedValues
```
To access a general description of the package, use the following command:

```ruby
?CardioRVAR
```

## Issues and requests

Please use the following link to create an issue or request:

https://github.com/CIMES-USNA-UMA/CardioRVAR/issues

## Contact information

Email: alvarochaoecija.rprojects@gmail.com

ORCID: https://orcid.org/0000-0002-2691-6936




