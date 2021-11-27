# Analysis code for **Linking Social and Personal Preferences: Theory and Experiment**

Authors: 

- [William R. Zame](http://www.econ.ucla.edu/zame/)
- [Bertil Tungodden](https://sites.google.com/view/bertiltungodden/home)
- [Erik Ø. Sørensen](https://www.statsokonomen.no/erik-o-sorensen-cv/)
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

> INSTRUCTIONS: The typical README in social science journals serves the purpose of guiding a reader through the available material and a route to replicating the results in the research paper. Start by providing a brief overview of the available material and a brief guide as to how to proceed from beginning to end.

Example: The code in this replication package constructs the analysis file from the three data sources (Ruggles et al, 2018; Inglehart et al, 2019; BEA, 2016) using Stata and Julia. Two master files run all of the code to generate the data for the 15 figures and 3 tables in the paper. The replicator should expect the code to run for about 14 hours.

## Data Availability and Provenance Statements

The experimental data used to support the findings of this study were collected by the authors. All data, 
with documentation, have been deposited in the public domain at  Dataverse:

- Zame, William R.; Tungodden, Bertil; Sørensen, Erik Ø.; Kariv, Shachar; Cappelen, Alexander W., 2021, "Replication Data for: Linking Social and Personal Preferences: Theory and Experiment", https://doi.org/10.7910/DVN/LUF59R, Harvard Dataverse, V2.

We certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript, and the data are licensed under a Creative Commons/CC0 license. See [LICENSE_CC0.txt](LICENSE_CC0.txt) for details.

- The data file is downloaded when the `targets` plan is first run.



## Computational requirements

> INSTRUCTIONS: In general, the specific computer code used to generate the results in the article will be within the repository that also contains this README. However, other computational requirements - shared libraries or code packages, required software, specific computing hardware - may be important, and is always useful, for the goal of replication. Some example text follows. 

> INSTRUCTIONS: We strongly suggest providing setup scripts that install/set up the environment. Sample scripts for [Stata](https://github.com/gslab-econ/template/blob/master/config/config_stata.do),  [R](https://github.com/labordynamicsinstitute/paper-template/blob/master/programs/global-libraries.R),  [Python](https://pip.readthedocs.io/en/1.1/requirements.html), [Julia](https://github.com/labordynamicsinstitute/paper-template/blob/master/programs/packages.jl) are easy to set up and implement.

### Software Requirements

> INSTRUCTIONS: List all of the software requirements, up to and including any operating system requirements, for the entire set of code. It is suggested to distribute most dependencies together with the replication package if allowed, in particular if sourced from unversioned code repositories, Github repos, and personal webpages. In all cases, list the version *you* used. 

- R (code was last run with version 4.1.2)
  - `renv` (XX)
  - The libraries and versions specified by the `renv.lock` file will be installed into a project specific library.

### Memory and Runtime Requirements

> INSTRUCTIONS: Memory and compute-time requirements may also be relevant or even critical. Some example text follows. It may be useful to break this out by Table/Figure/section of processing. For instance, some estimation routines might run for weeks, but data prep and creating figures might only take a few minutes.

#### Summary

Approximate time needed to reproduce the analyses on a standard (CURRENT YEAR) desktop machine:

- [ ] <10 minutes
- [ ] 10-60 minutes
- [X] 1-8 hours
- [ ] 8-24 hours
- [ ] 1-3 days
- [ ] 3-14 days
- [ ] > 14 days
- [ ] Not feasible to run on a desktop machine, as described below.

#### Details

The code was last run on a **4-core Intel-based laptop with MacOS version 10.14.4**. 

Portions of the code were last run on a **32-core Intel server with 1024 GB of RAM, 12 TB of fast local storage**. Computation took 734 hours. 

Portions of the code were last run on a **12-node AWS R3 cluster, consuming 20,000 core-hours**.  

> INSTRUCTIONS: Identifiying hardware and OS can be obtained through a variety of ways:
> Some of these details can be found as follows:
>
> - (Windows) by right-clicking on "This PC" in File Explorer and choosing "Properties"
> - (Mac) Apple-menu > "About this Mac"
> - (Linux) see code in [tools/linux-system-info.sh](https://github.com/AEADataEditor/replication-template/blob/master/tools/linux-system-info.sh)`


## Description of programs/code

> INSTRUCTIONS: Give a high-level overview of the program files and their purpose. Remove redundant/ obsolete files from the Replication archive.

### (Optional, but recommended) License for Code

The code is licensed under a BSD-3-Clause license. See [LICENSE_BSD-3-Clause.txt](LICENSE_BSD-3-Clause.txt) for details.

## Instructions to Replicators

> INSTRUCTIONS: The first two sections ensure that the data and software necessary to conduct the replication have been collected. This section then describes a human-readable instruction to conduct the replication. This may be simple, or may involve many complicated steps. It should be a simple list, no excess prose. Strict linear sequence. If more than 4-5 manual steps, please wrap a master program/Makefile around them, in logical sequences. Examples follow.




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

## References

> INSTRUCTIONS: As in any scientific manuscript, you should have proper references. For instance, in this sample README, we cited "Ruggles et al, 2019" and "DESE, 2019" in a Data Availability Statement. The reference should thus be listed here, in the style of your journal:
