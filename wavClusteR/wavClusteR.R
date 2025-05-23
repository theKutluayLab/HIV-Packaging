library(wavClusteR)

bam_path <- "/path/to/bam_files" #no slash
out_path <- "/path/to/wavClusteR_output" #no slash
exp <- "Experiment name"

dir.create(file.path(out_path, "plots", exp), recursive = T)

Bam <- readSortedBam(filename = file.path(bam_path, paste0(exp, ".sorted.bam")))

countTable <- countTable[ !elementMetadata( countTable )[, 'substitutions'] %in% c( 'AN', 'CN', 'GN', 'TN', 'NA', 'NC', 'NG', 'NT' ) ]
emd <- elementMetadata( countTable )
subst <- emd[, 'substitutions']
count <- emd[, 'count']
countPos <- table( subst )
n <- length( countPos )
substNames <- names( countPos )
highlight<-"TC"
posHL <- which( substNames == highlight )
percentagePos <- round( countPos[ posHL ] / sum( countPos ) * 100, 2 )
col <- rep( 'gray60', n )
col[ posHL ] <- 'skyblue2'
countReads <- sapply( split( count, subst ), sum )
percentageReads <- round( countReads[ posHL ] / sum( countReads ) * 100, 2 )

# count all substitutions
countTable <- getAllSub(sortedBam = Bam, minCov = 10)
head(countTable)

png(filename = file.path(out_path, "plots", exp, paste0(exp, "_sub_dist.png")), width = 2000, height = 2000, res = 300)
plotSubstitutions(countTable = countTable, highlight = "TC")
dev.off()


# countReads plot (number of reads with substitution)
png(filename = file.path(out_path, "plots", exp, paste0(exp, "_sub_nbr_reads.png")), 
    width = 2000, height = 2000, res = 300)

barplot(countReads, col=col, border="black", 
        main="",
        xlab="Substitutions", 
        ylab="Number of reads with substitution", 
        ylim=c(0, max(countReads)*1.3))

highlighted_value <- countReads[highlight]

dev.off()

# estimate a model
model <- fitMixtureModel(countTable = countTable, substitution = "TC")
str(model)

png(filename = file.path(out_path, "plots", exp, paste0(exp, "_model.png")), width = 3000, height = 2000, res = 300)
support <- getExpInterval(model = model, bayes = T)
dev.off()

png(filename = file.path(out_path, "plots", exp, paste0(exp, "_sub_dist_model.png")), width = 3000, height = 2000, res = 300)
plotSubstitutions(countTable = countTable, model = model, highlight = "TC")
dev.off()


# get high confidence PAR-CLIP induced transitions
highConfSub <- getHighConfSub(countTable = countTable, support = support, substitution = "TC")
head(highConfSub)

# identify protein binding sites (clusters)
coverage <- coverage(Bam)
coverage$chrX

# hard thresholding
clusters <- getClusters(highConfSub = highConfSub, coverage = coverage, sortedBam = Bam, threshold = 1)
# note ways to estimate threshold choice

# filterClusters
require(BSgenome.Hsapiens.UCSC.hg19)
wavclusters <- filterClusters(clusters = clusters,
                              highConfSub = highConfSub,
                              coverage = coverage,
                              model = model,
                              genome = Hsapiens,
                              refBase = "T",
                              minWidth = 12)

exportClusters(clusters = wavclusters,
               filename = file.path(out_path, "data", paste0(exp, "_wavClusters.bed")),
               trackname = "wavClusters",
               description = "wavClusters")




### PLOTS ###

require(GenomicFeatures)
txDB <- makeTxDbFromUCSC(genome = "hg19", tablename = "ensGene")
png(filename = file.path(out_path, "plots", exp, paste0(exp, "_anno_clusters.png")), width = 2000, height = 2000, res = 300)
annotateClusters(clusters = wavclusters, txDB = txDB,
                 plot = T, verbose = T)
dev.off()


png(filename = file.path(out_path, "plots", exp, paste0(exp, "_metagene.png")), width = 2000, height = 2000, res = 300)
metagene <- getMetaGene(clusters = wavclusters, txDB = txDB,
                        upstream = 1e3, downstream = 1e3,
                        nBins = 40, nBinsUD = 10,
                        minLength = 1, plot = TRUE, verbose = TRUE )
dev.off()


png(filename = file.path(out_path, "plots", exp, paste0(exp, "_metaTSS.png")), width = 2000, height = 2000, res = 300)
metaTSS <- getMetaTSS(sortedBam = Bam, txDB = txDB,
                      upstream = 1e3, downstream = 1e3, nBins = 40,
                      unique = FALSE, plot = TRUE, verbose = TRUE )
dev.off()


png(filename = file.path(out_path, "plots", exp, paste0(exp, "_size_dist.png")), width = 2000, height = 2000, res = 300)
plotSizeDistribution(clusters = wavclusters, showCov = TRUE, col = "skyblue2")
dev.off()


png(filename = file.path(out_path, "plots", exp, paste0(exp, "_stats.png")), width = 2000, height = 2000, res = 300)
plotStatistics(clusters = wavclusters, corMethod = "spearman", lower = panel.smooth)
dev.off()