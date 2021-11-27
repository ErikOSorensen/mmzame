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

On a powerful desktop computer, the replicator should expect the code to run for about 16 hours.

## Data Availability

The experimental data used to support the findings of this study were collected by the authors. All data, 
with documentation, have been deposited in the public domain at Harvard Dataverse:

- Zame, William R.; Tungodden, Bertil; Sørensen, Erik Ø.; Kariv, Shachar; Cappelen, Alexander W., 2021, "Replication Data for: Linking Social and Personal Preferences: Theory and Experiment", https://doi.org/10.7910/DVN/LUF59R, Harvard Dataverse, V2.

We certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript, and the data are licensed under a Creative Commons/CC0 license. See [LICENSE_CC0.txt](LICENSE_CC0.txt) for details.

The data file is downloaded when the `targets` plan is first run.


## Computational requirements


### Software Requirements


- R (code was last run with version 4.1.2)
  - `renv` (XX)
  - The libraries and versions specified by the `renv.lock` file will be installed into a project specific library when the master file is run.

### Memory and Runtime Requirements

Approximate time needed to reproduce the analyses on a (2021) desktop machine would be between 12 and 24 hours.

The code was last run on a Ubuntu 20.04.3 system, with an AMD Ryzen 9 3950X 16-Core Processor and 32 GB memory. 

With a less powerful system, it would be good to adjust the following line in `main.R`:

```
tar_make_future(workers = 26)
```

The number of workers should not be larger than the number of threads the computer can comfortably run in parallel.



## Description of programs/code

- `main.R`: Top level script to run and generate all outputs.
- `_targets.R`: Specification of DAG for creating all outputs using Will Landau's `targets` package.
- `renv.lock`: Specification of necessary libraries (and version) for running the code.
- `functions.R`: Function definitions called by the targets defined in `_targets.R`.

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

## List of tables and programs

> INSTRUCTIONS: Your programs should clearly identify the tables and figures as they appear in the manuscript, by number. Sometimes, this may be obvious, e.g. a program called "`table1.do`" generates a file called `table1.png`. Sometimes, mnemonics are used, and a mapping is necessary. In all circumstances, provide a list of tables and figures, identifying the program (and possibly the line number) where a figure is created.
>
> NOTE: If the public repository is incomplete, because not all data can be provided, as described in the data section, then the list of tables should clearly indicate which tables, figures, and in-text numbers can be reproduced with the public material provided.

The provided code reproduces:

- [ ] All numbers provided in text in the paper
- [ ] All tables and figures in the paper
- [ ] Selected tables and figures in the paper, as explained and justified below.


| Figure/Table #    | Program                  | Line Number | Output file                      | Note                            |
|-------------------|--------------------------|-------------|----------------------------------|---------------------------------|
| Table 1           | 02_analysis/table1.do    |             | summarystats.csv                 ||
| Table 2           | 02_analysis/table2and3.do| 15          | table2.csv                       ||
| Table 3           | 02_analysis/table2and3.do| 145         | table3.csv                       ||
| Figure 1          | n.a. (no data)           |             |                                  | Source: Herodus (2011)          |
| Figure 2          | 02_analysis/fig2.do      |             | figure2.png                      ||
| Figure 3          | 02_analysis/fig3.do      |             | figure-robustness.png            | Requires confidential data      |

