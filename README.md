# Analysis code for **Linking Social and Personal Preferences: Theory and Experiment**

Authors: 

- [William R. Zame](http://www.econ.ucla.edu/zame/)
- [Bertil Tungodden](https://sites.google.com/view/bertiltungodden/home)
- [Erik Ø. Sørensen](https://www.statsokonomen.no/erik-o-sorensen-cv/) [contact person for code and data, erik.sorensen@nhh.no]
- [Shachar Kariv](https://eml.berkeley.edu//~kariv/)
- [Alexander W. Cappelen](https://sites.google.com/view/alexander-w-cappelen/home)

**Abstract:** The attitudes of a Decision Maker toward riskless and risky choices—both personal
choices and social choices—enter virtually every realm of individual decision-making.
This paper asks when it is possible to link these attitudes. We provide a simple formalization of this question and necessary and sufficient conditions that such a link exists.
We also offer an experimental test of the theory in which subjects were confronted
with choices (involving monetary outcomes) in three domains: risky personal choices,
riskless social choices and risky social choices. Revealed preference tests show that
subject choices are generally consistent with utility maximization within each choice
domain but frequently involve at least some errors. We test for consistency across
choice domains using a novel nonparametric revealed preference test that accounts for
these errors.



## Overview

The master file this replication package will:

1. Install the required versions of the necessary `R` packages from CRAN.
2. Downloads the necessary datafile from Harvard Dataverse.
3. Create all the displays in the paper as separate files (documented below).
4. Create markdown documents for numbers referenced in the paper but not explicitly part of produced tables.

On a powerful desktop computer, the replicator should expect the code to run for about 6 hours.

## Data Availability

The experimental data used to support the findings of this study were collected by the authors. All data, 
with documentation, have been deposited in the public domain at Harvard Dataverse:

- Zame, William R.; Tungodden, Bertil; Sørensen, Erik Ø.; Kariv, Shachar; Cappelen, Alexander W., 2021, "Replication Data for: Linking Social and Personal Preferences: Theory and Experiment", https://doi.org/10.7910/DVN/LUF59R, Harvard Dataverse, V2.

We certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript, and the data are licensed under a Creative Commons/CC0 license. See [LICENSE_CC0.txt](LICENSE_CC0.txt) for details.

The data file is downloaded when the `targets` plan is first run.


## Computational requirements


### Software Requirements


- `R` (code was last run with version 4.1.2)
- `renv` (0.14.0)
    - system `curl` is a dependency of `renv`, should come with Windows 10 and Mac OSX, easily installable on linux.
- The libraries and versions specified by the `renv.lock` file will be installed into a project specific library when the master file is run.
- `rmarkdown` will require system `pandoc`, which is easiest installed with `Rstudio`.


### Memory and Runtime Requirements

Approximate time needed to reproduce all the analyses on a desktop machine (2021) would be between 5 and 10 hours.

The code was last run on a Ubuntu 20.04.3 system, with an AMD Ryzen 9 3950X 16-Core Processor and 32 GB memory. 
Memory requirements are limited, at around 200 MB per worker node. 
Disk use is very limited (at less than 0.5GB).

With a less powerful system, it would be good to adjust the following line in `main.R`:

```
tar_make_future(workers = 26)
```

The number of workers should not be larger than the number of threads the computer can comfortably run in parallel. 
There are 16 long running worker nodes. 
Parallelization is handled by Henrik Bengtsson's `future` library. 
The full potential efficiency gains from parallelization are not achieved, ordinary 
desktop computers will be fully occupied.



## Description of programs/code

- `main.R`: Top level script to run and generate all outputs.
- `_targets.R`: Specification of DAG for creating all outputs using Will Landau's `targets` package.
- `functions.R`: Function definitions called by the targets defined in `_targets.R`.
- `renv.lock`: Specification of necessary libraries (and version) for running the code.
- `.Rprofile`: Ensures that `renv` loads local libraries.


The graphical displays are produced in the `graphs/` directory (as pdf-files).
The other numbers and descriptions produced are part of the vignettes that are rendered in the `vignettes/` directory.

### License for Code

The code is licensed under a BSD-3-Clause license. See [LICENSE_BSD-3-Clause.txt](LICENSE_BSD-3-Clause.txt) for details.

## Instructions to Replicators

From the command line, the replicator should run the `main.R` script:

```
Rscript main.R
```

The first time it is run, this will install a local library of R packages from CRAN, 
download the data from Harvard Dataverse, and build all the targets specified in `_targets.R` 
(a specification of a directed acyclic graph of dependencies for generating all
necessary displays and descriptions). 

The `targets` system is smart about
caching intermediate results, so while running `main.R` takes a considerable amount of 
time for the first run, minor adjustments to the output routines in the vignettes are 
do not require the heavy computations to be re-run, and running `main.R` for the second
time is almost free of costs with respect to changes in the display layer. 


## List of tables and programs


The provided code reproduces all numbers provided in text in the paper. 
The table below indicates where and how the displays are produced.
Note that while the point at which displays are save is recorded in the table below,
there are dependencies (as described in `_targets.R`) such that it is always safest 
to run the `main.R` script to generate the displays. 


| Figure/Table #    | Program                  |Line Number | Output file                      | Note                            |
|-------------------|--------------------------|-------------|----------------------------------|---------------------------------|
| Table 1           |            |             |                ||
| Table 2           |  |           |                         ||
| Figure 1          | n.a. (no data)           |             |                                  | (theoretical illustration)         |
| Figure 2          | n.a. (no data)           |             |                                  | (theoretical illustration)          |
| Figure 3          |vignettes/aggregate_behavior.Rmd      |             | graphs/aggregate_choices.pdf       |       |
| Figure 4          |vignettes/individual_behavior.Rmd |                | graphs/logprice_scatters.pdf | | 
| Figure 5          |vignettes/testing_rationality.Rmd |                | graphs/empirical_cceis.pdf | |
| Figure 6          |vignettes/testing_rationality.Rmd |                | graphs/empirical_cceis_and_Bronars.pdf  | |
| Figure 7          |vignettes/testing_theory.Rmd      |                | graphs/prop3_permutations.pdf | | 
