May 18 2019

This directory contains the modified files from the `scanpy` package that is uesed in the paper.  We approximate the scanpy directory structure here, and include the scanpy License.

* The file `scanpy/tools/_rank_genes_groups.py` contains fixes to several bugs that were present in the version of scanpy that we examined (see the paper for the exact commit).  These bugs have since been fixed in the general version of scanpy.
* The file `scanpy/utils.py` contains some hacks to allow for integer valued class names in the differential expression calculations.  We used this when looking at two groups (e.g. for the simulated zheng data).
* The `LICENSE` file is scanpy's license file.  See the `scanpy` github page at  `github.com/theislab/scanpy` for the full package.

The file `scanpy-pval.ipynb` contains the code that we used to find the p-values for the Wilcoxon, t-test, and logistic regression methods on the Paul, Zeisel, ZhengFull, ZhengFilt, and ZhengSim data sets.
