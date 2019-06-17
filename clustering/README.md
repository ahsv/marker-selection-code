04 June 2019

This directory contains the data used for the clustering metrics.

* The file `clustering-metric.ipynb` contains all of the calculations, starting with the markers.  In the case of Spa and RankCorr, the markers are computed on-the-fly in the notebook.  The markers are all included in this repo; however, the paths to the markers in this notebook are probably incorrect.  This notebook produces the data files that are in this directory - all of the ari, ami, and fms (or aveScore) files.

* The file `clustering-metric-graphs.ipynb` contains the graphs that appear in the paper.  This notebook uses the data files that are in this directory and thus should run quite easily.
