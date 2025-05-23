library(tidyverse)
library(RColorBrewer)
library(pals)
library(ggthemes)
library(ggrepel)
library(gtable)
require(grid)
library(gridExtra)

# script for viral density map

setwd('/path/to/directory/') #with slash
title <- 'experiment name'

# read in the annotaion file
anno = read_tsv('/path/to/HIV_NL4-3.gtf')
anno = data.frame(anno)

cds = anno[anno$Name != "RefSeq", 1:8]
cds = cds[cds$Strand == "+",]
cds$Frame = as.numeric(cds$Frame)

cds$Feature[which(cds$Feature=='five_prime_UTR')] = "5'UTR"
cds$Feature[which(cds$Feature=='three_prime_UTR')] = "3'UTR"
cds$Name[which(cds$Name == 'RefSeq')] = cds$Feature[which(cds$Name=='RefSeq')]

# change the range of ORF1A
# Note: there is an overlapping region between ORF1A and ORF1B
# the overlapping region will be labeled ORF1B for the purpose of plotting
cds$End[2] = cds$Start[3]-1

#read in the counts_repaired file
filelist = list.files(path=".",pattern='sense_counts_repaired')


if (!dir.exists("plots")) {
  dir.create("plots")
}

for (i in 1:length(filelist)) {
filename = unlist(strsplit(filelist[i],'_'))
filename <- sub("\\.txt$", "", filelist[i]) #used for .txt files
savepath = paste('plots',filename[1],sep='/')
print(filename)
print(savepath)

sense = read.table(filelist[i], header=TRUE, sep='\t')
sense = data.frame(sense$bp, rowSums(sense[2:5]))
colnames(sense) = c('bp','count')


# plotting data
plotting.data = data.frame(bp = 1:tail(cds$End,1))

plotting.data = merge(plotting.data, sense, all=TRUE)
colnames(plotting.data)[2] = 'sense_count'

plotting.data$sense_count[is.na(plotting.data$sense_count)] <- 0
plotting.data$gene = ''

# add 3'UTR and 5'UTR
cds <- rbind(data.frame(Chromosome = "NL4-3", Name = "5'UTR", Feature = "CDS", Start = 1, End = 335, Score = ".", Strand = "+", Frame = "0"), cds)
cds <- rbind(cds, data.frame(Chromosome = "NL4-3", Name = "3'UTR", Feature = "CDS", Start = 8714, End = length(plotting.data$bp), Score = ".", Strand = "+", Frame = "0"))

cds1<-data.frame(Name=cds$Name)
cds1<-rbind(cds1,data.frame(Name="na"))
cds1<-rbind(cds1,data.frame(Name="na1"))
cds1<-rbind(cds1,data.frame(Name="na2"))

for (i in 1:length(cds$Name)) {
  plotting.data$gene[cds$Start[i]:cds$End[i]] = cds$Name[i]
}

plotting.data$gene[1:335] <- "5'UTR"
plotting.data$gene[8714:length(plotting.data$bp)]<-"3'UTR"


plotting.data$gene1 <- plotting.data$gene
plotting.data$gene1[plotting.data$gene1 == ""] <- "na"
plotting.data$gene1[8342] <- "na1" #change and use this line if needed only for aesthetics!!!!! Check na_data to see which data should be named na1 or na2. Add more nas if necessary.
plotting.data$gene1[5592:5607] <- "na2" #change and use this line if needed only for aesthetics!!!!!

na_rows <- which(plotting.data$gene1 == "na")
na_data <- plotting.data[na_rows, ]

y_max = max(sense$count) 
y_lim = max(log10(plotting.data$sense_count))
print(max(log10(plotting.data$sense_count)))

### PLOT TYPE I

#mixPalette = c('#FFFFFF', sample(as.vector(tableau20(14))))
#colPalette = mixPalette
colPalette = c("#FFFFFF", "#AEC7E87F", "#83c9817f", 
               "#C49C947F", "#C5B0D57F", "#ffd2787f", "#d6272762",
               "#ff7e0e6a", "#9467BD7F", "#E377C27F", 
               "#F7B6D27F", "#2e40b366","#2ca0986f","#a02c6864",
               "#aba5a8a2", "#aba5a8a2","#aba5a8a2","#aba5a8a2")


breaks <- c(0, 336, seq(1000, 8714, 1000), 8714, length(plotting.data$bp))
labels <- c(0, 336, seq(1000, 8714, 1000), 8714, length(plotting.data$bp))

### plotting sense-counts
# sense count with log-norm
ps <- (ggplot(data = plotting.data, aes(x = bp, y = log10(sense_count))) 
      + geom_area(aes(fill = gene1))
      + scale_fill_manual(breaks = c('',cds1$Name),
                          values = colPalette,
                          labels = c('',cds1$Name))
      + theme_classic()
      + xlim(0,length(plotting.data$bp))
      + ylim(0, y_lim)
      + labs(y = 'log10(counts)')
      + theme(legend.position = 'none',
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=10),
              axis.text.x = element_text(size=8),
              axis.text.y = element_text(size=8))
      + scale_x_continuous(breaks = breaks, labels = labels)
      )


# sense count WITHOUT log-norm
psl <- (ggplot(data = plotting.data, aes(x = bp, y = sense_count)) 
       + geom_area(aes(fill = gene1))
       + scale_fill_manual(breaks = c('',cds1$Name),
                           values = colPalette,
                           labels = c('',cds1$Name))
       + theme_classic()
       + xlim(336,length(plotting.data$bp))
       + ylim(0,y_max)
       + labs(y = 'counts')
       + theme(legend.position = 'none',
               axis.title.x = element_blank(),
               axis.title.y = element_text(size=10),
               axis.text.x = element_text(size=8),
               axis.text.y = element_text(size=8))
      + scale_x_continuous(breaks = breaks, labels = labels)

       )

### plotting genome CDS regions
gme <- (ggplot(data = plotting.data, aes(x = bp))
      + theme_tufte()
      + xlim(0,length(plotting.data$bp))
      + ylim(0,1)
      + geom_rect(aes(NULL,NULL,xmin=bp-1,xmax=bp,fill=gene),
                  ymin=0, ymax=0.1)

      + scale_fill_manual(breaks = (c('',cds$Name)),
                          values = colPalette,
                          labels = (c('',cds$Name)))
      + labs(title = title)
      + theme(plot.title = element_text(family='sans', size=13, hjust=0.5),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = c(0.5,0.3),
              legend.direction = 'horizontal',
              legend.title = element_blank(),
              legend.key.size = unit(0.5,'cm'),
              legend.text = element_text(family='sans', size=10))
      + guides(fill = guide_legend(nrow = 1, order = 2))
      )


ps_histogram <- (ggplot(data = plotting.data, aes(x = bp, y = log10(sense_count))) 
      + geom_bar(stat="identity", width=1, aes(fill=gene1))
      + scale_fill_manual(breaks = c('',cds1$Name),
                          values = colPalette,
                          labels = c('',cds1$Name))
      + theme_classic()
      + xlim(0,length(plotting.data$bp))
      + ylim(0, y_lim)
      + labs(y = 'log10(counts)')
      + theme(legend.position = 'none',
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=10),
              axis.text.x = element_text(size=8),
              axis.text.y = element_text(size=8))
      + scale_x_continuous(breaks = breaks, labels = labels)
      )

psl_histogram <- (ggplot(data = plotting.data, aes(x = bp, y = sense_count)) 
            + geom_bar(stat="identity", width=1, aes(fill=gene1))
            + scale_fill_manual(breaks = c('',cds1$Name),
                                values = colPalette,
                                labels = c('',cds1$Name))
            + theme_classic()
            + xlim(336,length(plotting.data$bp))
            + ylim(0,y_max)
            + labs(y = 'counts')
            + theme(legend.position = 'none',
               axis.title.x = element_blank(),
               axis.title.y = element_text(size=10),
               axis.text.x = element_text(size=8),
               axis.text.y = element_text(size=8))
            + scale_x_continuous(breaks = breaks, labels = labels)

)

# combining plots
g1 <- ggplotGrob(gme)
g2 <- ggplotGrob(ps_histogram)
g3 <- ggplotGrob(psl_histogram)


# log-normed plot [Genome region bar + log-norm sense]
g <- rbind(g1, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g3$widths)

png(paste(savepath,'_vgenome.png',sep=''), 12, 3, units='in', res=300)
grid.newpage()
grid.draw(g)
dev.off()

# full plot [Genome region bar + raw sense + log-norm sense]
g <- rbind(g1, g3, g2, size = "first")
g$widths <- unit.pmax(g1$widths, g3$widths, g2$widths)

png(paste(savepath,'_vgenome_full.png',sep=''), 12, 3, units='in', res=300)
grid.newpage()
grid.draw(g)
dev.off()
