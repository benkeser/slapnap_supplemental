
# `slapnap_supplement/multiple_sensitivty`

> Validation of slapnap for use with multiple antibodies

**Authors:** [David Benkeser](https://www.github.com/benkeser/), Brian
D. Williamson, Craig A. Magaret, Sohail Nizam, Peter B. Gilbert

-----

## Description

This directory contains all code needed to reproduce the `slapnap` analysis of multiple antibody regimens. The dependency structure of the code can be examined by studying the `Makefile`. See `make help` for more details on the components of the analysis.

The code requires the following installed software:
* `make`
* `docker`
* `R` (with `Rscript` executable in `PATH`)
* `R` packages: `here`, `SuperLearner`, `cvAUC`, `ggplot2`

With this software is installed, the analysis can be reproduced by executing `make auc_fig` or simply `make` at the command line. However, do note that this will simultaneously run 32 Docker containers to obtain results. Consequently, you may need to increase the available memory available to your Docker program and/or ensure that you are on a computer that has adequate resources for executing this task.

-----

## Details

The analysis is executed in several stages. First `slapnap_output` is made by running the `slapnap` container for each of 16 combination antibodies. For each antibody, two containers are run: one to trainer the learner based on the individual antibody data, the second to obtain a data set for the combination antibody. These results are saved in `bash_output/` in folders with the name of each antibody separated by underscores and a `_2` to denote the folder where data for the combination antibody is saved. 

Next, `R_output/all_nabs.txt` is created, which simply formats the antibody names to be used in later `R` scripts. 

Next, `R_output/rslt_df.RData` is created, by executing `R/combine_results.R`. Note that this execution occurs *inside* the slapnap container, to ensure that there are no compiler errors associated with using `slapnap`-trained learners to predict on features. 

Next, `R_output/n_sens_df.RData` is created, by executing `R/describe_combo_sens.R`. The file contains information on the number of resistant sequences for each measured combination antibody regimen. 

Finally, all results are pieced together by executing `R/make_auc_fig.R` (locally) with output being an image saved in `fig/auc_fig.tiff`.

-----

## Issues

If you encounter any bugs, please [file an issue](https://github.com/benkeser/slapnap_supplement/issues).

-----

## License

© 2020- David Benkeser

The contents of this repository are distributed under the MIT license:

    The MIT License (MIT)
    
    Copyright (c) 2020- David Benkeser
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
