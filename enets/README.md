Alexander Vargo
23 February 2019

This directory contains information about Elastic Nets.

* The notebook "paul-elastic-nets.ipynb" contains the code to fit elastic nets
  to the paul data set.
    * Each file of the form "paul15-nets-fold{}-marks.npz" contains the selected
      markers for a fold in 'arr_0'
    * Each file of the form "paul15-nets-fold{}-coefs.npz" contains the coefficients 
      for the fitted enet for a fold in 'arr_0'

* The notebook "zeisel-elastic-nets.ipynb" contains the code to fit elastic nets
  to the zeisel data set.  This is worth looking at - it shows the magic `precomputed=True`
  option that you need to use for faster convergence
    * Each file of the form "zeisel-nets-fold{}-marks.npz" contains the selected
      markers for a fold in 'arr_0'
    * Each file of the form "zeisel-nets-fold{}-coefs.npz" contains the coefficients 
      for the fitted enet for a fold in 'arr_0'
