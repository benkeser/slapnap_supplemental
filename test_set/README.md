
# `slapnap_supplement/test_set`

> Validation of slapnap on test data external to CATNAP

**Authors:** [David Benkeser](https://www.github.com/benkeser/), [Brian
D. Williamson](https://www.github.com/bdwilliamson), Craig A. Magaret, Courtney Simmons, Sohail Nizam, Peter B. Gilbert

-----

## Description

This directory contains all code needed to reproduce the `slapnap` test-set analyses of 3BNC117 + 10-1074.

The code requires the following installed software:
* `docker`
* `R` (with `Rscript` executable in `PATH`)
* `R` packages: `here`, `SuperLearner`, `cvAUC`, `dplyr`, `tibble`, `ggplot2`, `cowplot`

With this software is installed, the analysis can be reproduced by:

1. Executing `bash chmod u+x bash/run_docker_containers.sh` from the command line;
2. Executing `bash bash/run_docker_containers.sh` and waiting for the results to finish (this may take a while);
3. Executing ;
4.

-----

## Details

The analysis is executed in several stages. First the `slapnap` container is run for the combination of antibodies 10-1074 and 3BNC117. The results are saved in `bash_output/` in a folder with the name of this antibody combination.

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
