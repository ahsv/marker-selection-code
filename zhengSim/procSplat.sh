#!/bin/bash

# process the cells to get the groups into the y vector
awk 'match($4, /([1-9])/, m) {print m[1]}' cellInfo.dat > y.dat

# proces the genes to get a list of differentially expressed genes and the DE
# factor.  Use `NR-2` since we write.table so the first record is the column
# names and awk starts indexing records at 1.
awk '$6 != 1 && NR > 1 {print NR-2, $6}' geneInfo.dat > deGenes.dat


# process the counts data: delete the first line and first column This should
# be fine since we are not currently planning to simulate huge data sets.  If
# the data sets get larger, you may want a more elegant solution here
awk 'NR > 1 {print}' counts.dat > tmp.dat
cut -d" " -f 2- tmp.dat > counts.dat
