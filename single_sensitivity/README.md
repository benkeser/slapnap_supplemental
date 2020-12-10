
# `slapnap_supplement/single_sensitivty`

> Validation of slapnap for use with single antibodies

**Authors:** [David Benkeser](https://www.github.com/benkeser/), [Brian
D. Williamson](https://www.github.com/bdwilliamson), Craig A. Magaret, Sohail Nizam, Peter B. Gilbert

-----

## Description

This directory contains all code needed to reproduce the `slapnap` analysis of single antibody regimens. The dependency structure of the code can be examined by studying the `Makefile`. See `make help` for more details on the components of the analysis.

The code requires the following installed software:
* `make`
* `docker`
* `R` (with `Rscript` executable in `PATH`)
* `R` packages: `here`, `SuperLearner`, `cvAUC`, `dplyr`, `tibble`, `ggplot2`, `cowplot`

With this software is installed, the analysis can be reproduced by executing `make auc_fig` or simply `make` at the command line. You may need to remove the existing output figure `fig/auc_fig.tiff` first.

Do note that this will simultaneously run 33 Docker containers to obtain results. Consequently, you may need to increase the available memory available to your Docker program and/or ensure that you are on a computer that has adequate resources for executing this task.

Also, the "make" command may fail, since this code takes some time to run and the Docker containers are running in the background. If this occurs, wait for the Docker containers to finish running, edit the Makefile (making sure you have a copy of the original somewhere) to remove dependence on the Docker containers, and re-run "make" from that point.

-----

## Details

The analysis is executed in several stages. First `slapnap_output` is made by running the `slapnap` container for each of 33 antibodies. The results are saved in `bash_output/` in folders with the name of each antibody.

Next, `R_output/all_nabs.txt` is created, which simply formats the antibody names to be used in later `R` scripts.

Next, `R_output/rslt_df.RData` is created, by executing `R/combine_results.R`. Note that this execution occurs *inside* the slapnap container, to ensure that there are no compiler errors associated with using `slapnap`-trained learners to predict on features.

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
