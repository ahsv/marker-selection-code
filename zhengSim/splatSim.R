
# some parameters
nCells = 5000
#nGenes = 32738
paramFile = "bCell-params.rds"
# for the pre-filtered genes...
nGenes = 4999 

# need splatter
library(splatter)


# make simulation parameters - need to set the seed since it is not random :/
# also check to see if we have the simulated parameters already in the
# directory this should save A LOT of time once the parameters have already
# been generated for your data set.

seed <- sample(1:1000000,1)
write(seed, "seed.txt")

if (file.exists(paramFile)) {
    attach(paramFile)
    params <- setParams(params, seed = seed)
} else {
    params <- newSplatParams(seed=seed)

    # load the data
    mat <- scan('bCells.dat')
    mat <- matrix(mat, ncol=nGenes, byrow=TRUE)
    storage.mode(mat) <- "integer"
    mat <- t(mat)

    print(dim(mat))

    params <- splatEstimate(mat, params)
}
 
print(seed)

# here we are not using dropout, since it is stochasitc, not attenuated.
params <- setParams(params, batchCells = nCells, group.prob = c(0.5, 0.5), de.prob = c(0.1, 0), nGenes = 5000)
                    #, dropout.type="experiment")

# save the parameters if they are not found locally
if (!file.exists(paramFile)) {
    save(params, file=paramFile)
}


# simulate
sim <- splatSimulate(params, method="groups")

# write output to files
write.table(rowData(sim), file="geneInfo.dat", sep=" ")
write.table(colData(sim), file="cellInfo.dat", sep=" ")
write.table(counts(sim), file="counts.dat", sep=" ")
