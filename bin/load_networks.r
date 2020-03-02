#########
# To load full network from binary: 
netfile="agg.rerank.Rdata"
label="agg.rerank"
load(netfile)
network = diag(length(genes.t))
bottom = row(network) > col(network)
colnames(network) = genes.t
rownames(network) = genes.t
network[bottom] = temp
network = network + t(network)
diag(network) = 1


#########
# To load full network from HDF5
require(rhdf5)
genes <- rhdf5::h5read("generic.genes.h5", "genes" )
network <- rhdf5::h5read("generic.net.h5", "net" )

# To access individual gene(s) from network
require(rhdf5)
genes <- rhdf5::h5read("generic.genes.h5", "genes" )
index = which(genes[,2] == "KCNH2")
coexp <- rhdf5::h5read("generic.net.h5", "net", index = list(index, NULL))
colnames(network) = genes[,1] # 1 = entrezIDs, 2 = symbols, 3 = ensemblIDs
 rownames(network) = genes[,1]

m = match(c("KCNH2", "CACNA1C", "SNC5A"), genes[,2])
index = m[!is.na(m)]
subnetwork <- rhdf5::h5read("generic.net.h5", "net", index = list(index, NULL))

