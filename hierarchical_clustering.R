suppressMessages(library(cba))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library(RSvgDevice))

#dat = read.table("/Users/Fiziev/projects/chrom_compare/data/tonis/concatenated.with_replicates.12_states/comparisons/all_vs_all_distances.csv", sep="\t", row.names=1, header=TRUE)
dat = read.table("/Users/Fiziev/projects/chrom_compare/data/roadmap/all_vs_all/all_vs_all_distances.csv", sep="\t", row.names=1, header=TRUE)


metric = "correlation"
clust_method <- "complete"

#dissimilarity <- 1 - cor(dat)
#rd <- dist(dat, "correlation")
dat <- dat - min(dat)
rd <- as.dist(dat)
#rd <- dist((1 - cor(dat))/2)

rc <- hclust(rd,method=clust_method)
ro <- order.optimal(rd, rc$merge)
rc$merge <- ro$merge
rc$order <- ro$order

#cd <- dist(t(dat), "correlation")
cd <- as.dist(dat)
#cd <- dist((1 - cor(t(dat)))/2)

cc <- hclust(cd,method=clust_method)
co <- order.optimal(cd, cc$merge)
cc$merge <- co$merge
cc$order <- co$order

rowv = as.dendrogram(rc)
colv = as.dendrogram(cc)

heatmap.2(data.matrix(dat),
          Rowv=rowv,
          keysize=1, 
          Colv=colv,
          scale="none",
          trace = "none", 
          margins = c(10, 10),
          density.info = "none",
          dendrogram = "col",
          col=brewer.pal(9,"OrRd"))
