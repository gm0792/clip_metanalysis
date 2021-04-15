# ICLIP ANALYSIS 
# 2nd March 2021
setwd("/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/")


#DISTRIBUTION OF XLINKS OVER GENOME FRATURE-----------------------
# Running this code requires a bed file with 'crosslink sites' output by BEDtools (*.xl.bed.gz)
# Path to bed file with curated binding sites
#files previously unzipped in terminal 

library(dplyr)    
library(ggplot2) 
library(DESeq2)
library(data.table)
library(dplyr)
library(tidyverse)
library(tidyr)
library(reshape2)
library(BiocManager)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(plotly)
library(ggplot2)
library(tximport)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(readr)
library(scales) 



#FUNCTIONS----------------------

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}



#COLOUR PALETTES----------------------
mypal.dedup<-c("input_reads"= "#5FB1D4",  "output_reads"= "#39373B") 
mypal.dedup.ratio<-c("ratio"= "#5FB1D4")
mypal.region<-c("intron" = "#454241","CDS" ="#F0421C", "intergenic" ="#DDD3D1", "ncRNA" = "#3DB75D", "UTR5" = "#3DB7E6", "UTR3"= "#D644C7")
mypal.model<-c("FTLD-TDP_human_brain" = "#1A7387","SH-SY5Y_neuroblastoma"="#32A8BD","293Flp _exp_GFP-TDP43" ="#E66405", "embrionic_stem_cells" ="#319786")


metadata$model

# DEDUPLICATION------------------------------------

#import dedup.log Tollervey et al data--------------------

tollervey.dedup.li = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/tollervey/data/run/results/dedup", pattern = ".log$", full.names = TRUE) #store file paths for each file in a list
tollervey.fi.li.names = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/tollervey/data/run/results/dedup", pattern = ".log$", full.names = FALSE) %>%
  gsub(".log","",.)  #extracts files names in a seprate list
tollervey.dedup.list = list() #create empty list
for (i in 1:length(tollervey.dedup.li)){
  temp = read_delim(tollervey.dedup.li[[i]],"\t", escape_double = FALSE, trim_ws = TRUE)
  tollervey.dedup.list[[i]] = temp
  
}

names(tollervey.dedup.list) <-tollervey.fi.li.names 

#import dedup.log G.Rot et al data--------------------
#path to Nobby's run 
#/camp/lab/luscomben/home/users/chakraa2/projects/giulia/results/dedup
grot.dedup.li = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/grot/prova/run-imp/results/dedup", pattern = ".log$", full.names = TRUE) #store file paths for each file in a list
grot.fi.li.names = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/grot/prova/run-imp/results/dedup", pattern = ".log$", full.names = FALSE) %>%
  gsub(".log","",.)  #extracts files names in a seprate list
grot.dedup.list = list() #create empty list
for (i in 1:length(grot.dedup.li)){
  temp = read_delim(grot.dedup.li[[i]],"\t", escape_double = FALSE, trim_ws = TRUE)
  grot.dedup.list[[i]] = temp
  
}

names(grot.dedup.list) <-grot.fi.li.names #assign samples name to list's elements

#build dedup reast in/out/ratio dataset Tollervey et al----------
tollervey.dedup.input.li = list() #list input reads
tollervey.dedup.output.li = list() #list output reads
for (i in 1:length(dedup.list)){
  tollervey.dedup.input.li[i]<-as.integer(sub(".*?Input Reads:.*?(\\d+).*", "\\1", tollervey.dedup.list[[i]]))
  tollervey.dedup.output.li[i]<-as.integer(sub(".*?Number of reads out:.*?(\\d+).*", "\\1", tollervey.dedup.list[[i]]))
}

names(tollervey.dedup.input.li) <-tollervey.fi.li.names #assign samples name to list's elements
names(tollervey.dedup.output.li) <-tollervey.fi.li.names #assign samples name to list's elements
tollervey.dedup.input.df<-do.call(rbind, tollervey.dedup.input.li) %>% as.data.table( .,keep.rownames=TRUE) #transform to df -input reads df 
colnames(tollervey.dedup.input.df)<- c("sample","input_reads")
tollervey.dedup.output.df<-do.call(rbind, tollervey.dedup.output.li) %>% as.data.table( .,keep.rownames=TRUE) #transform to df -output reads df 
colnames(tollervey.dedup.output.df)<- c("sample","output_reads")

tollervey.dedup.df<-left_join(tollervey.dedup.input.df,tollervey.dedup.output.df) %>% as.data.frame() #join df by sample names

tollervey.sample.order<-c("tollervey_ESC", "tollervey_SHSY5Y1a","tollervey_SHSY5Y1b","tollervey_SHSY5Y2","tollervey_SHSY5Y3",
                "tollervey_SHSY5Y_cyt","tollervey_SHSY5Y_nucl","tollervey_brain1","tollervey_brain2","tollervey_brain3",
                "tollervey_brain4","tollervey_brain5","tollervey_brain6.high","tollervey_brain7.low")

#setting df levels to reorder samples
#input reads and output reads 
tollervey.dedup.df<-tollervey.dedup.df[match(stollervey.sample.order, tollervey.dedup.df$sample),]
tollervey.dedup.df$sample <- factor(tollervey.dedup.df$sample, levels = tollervey.sample.order)
tollervey.dedup.tidy<-melt(tollervey.dedup.df, variable.name = "read", value.name = "count")
tollervey.dedup.tidy$sample <- factor(tollervey.dedup.tidy$sample, levels = tollervey.sample.order)

#ratio df 
tollervey.dedup.ratio.df<-left_join(tollervey.dedup.input.df,tollervey.dedup.output.df) %>% mutate(ratio= input_reads/output_reads) %>% dplyr::select(sample,ratio) %>% as.data.table( .,keep.rownames=TRUE)
tollervey.dedup.ratio.df<-tollervey.dedup.ratio.df[match(tollervey.sample.order, tollervey.dedup.ratio.df$sample),]
tollervey.dedup.ratio.df$sample <- factor(tollervey.dedup.ratio.df$sample,levels = stollervey.sample.order)
tollervey.dedup.radio.tidy<-melt(tollervey.dedup.ratio.df, variable.name = "read", value.name = "count")


#build dedup reast in/out/ratio dataset G.Rot et al----------
grot.dedup.input.li = list() #list input reads
grot.dedup.output.li = list() #list output reads
for (i in 1:length(dedup.list)){
  grot.dedup.input.li[i]<-as.integer(sub(".*?Input Reads:.*?(\\d+).*", "\\1", grot.dedup.list[[i]]))
  grot.dedup.output.li[i]<-as.integer(sub(".*?Number of reads out:.*?(\\d+).*", "\\1", grot.dedup.list[[i]]))
}

names(grot.dedup.input.li) <-grot.fi.li.names #assign samples name to list's elements
names(grot.dedup.output.li) <-grot.fi.li.names #assign samples name to list's elements
grot.dedup.input.df<-do.call(rbind, grot.dedup.input.li) %>% as.data.table( .,keep.rownames=TRUE) #transform to df -input reads df 
colnames(grot.dedup.input.df)<- c("sample","input_reads")
grot.dedup.output.df<-do.call(rbind, grot.dedup.output.li) %>% as.data.table( .,keep.rownames=TRUE) #transform to df -output reads df 
colnames(grot.dedup.output.df)<- c("sample","output_reads")

grot.dedup.df<-left_join(grot.dedup.input.df,grot.dedup.output.df) %>% as.data.frame() #join df by sample names

grot.sample.order<-c("grot_293fl_1", "grot_293fl_2")

#setting df levels to reorder samples
#input reads and output reads 
grot.dedup.df<-grot.dedup.df[match(grot.sample.order, grot.dedup.df$sample),]
grot.dedup.df$sample <- factor(grot.dedup.df$sample, levels = grot.sample.order)
grot.dedup.tidy<-melt(grot.dedup.df, variable.name = "read", value.name = "count")
grot.dedup.tidy$sample <- factor(grot.dedup.tidy$sample, levels = grot.sample.order)

#ratio df 
grot.dedup.ratio.df<-left_join(grot.dedup.input.df,grot.dedup.output.df) %>% mutate(ratio= input_reads/output_reads) %>% as.data.frame( .,keep.rownames=TRUE) %>% dplyr::select(sample,ratio)

grot.dedup.ratio.df<-grot.dedup.ratio.df[match(grot.sample.order, grot.dedup.ratio.df$sample),]
grot.dedup.ratio.df$sample <- factor(grot.dedup.ratio.df$sample,levels = grot.sample.order)
grot.dedup.radio.tidy<-melt(grot.dedup.ratio.df, variable.name = "read", value.name = "count")


#merging tidy df from different datasets-----------------------------

#input reads and output reads df
df.dedup.in.out<-dplyr::bind_rows(tollervey.dedup.tidy,grot.dedup.tidy)

#reads ratio df
df.dedup.ratio<-dplyr::bind_rows(tollervey.dedup.radio.tidy,grot.dedup.radio.tidy)


#Dedup Plots-----------------------------
#input reads and output reads
dedup.in.out.plot<-ggplot(df.dedup.in.out) +
  geom_bar(stat = "identity", aes(x=sample, y=count, fill=read), position = "dodge") +
  scale_color_manual(values=mypal.dedup) +
  scale_fill_manual(values=alpha(c(mypal.dedup))) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  coord_flip()
dedup.in.out.plot

#reads ratio
dedup.ratio.plot<-ggplot(df.dedup.ratio) +
  geom_bar(stat = "identity", aes(x=sample, y=count, fill=read), position = "dodge") +
  scale_color_manual(values=mypal.dedup.ratio) +
  scale_fill_manual(values=alpha(c(mypal.dedup.ratio))) +
  coord_flip()

#INSERT THRESHOLD! 
threshold.dedup = 10
dedup.ratio.plot<-dedup.ratio.plot + geom_hline(yintercept = threshold.dedup , linetype = "dashed", color ="black") 
dedup.ratio.plot





#GENE LENGTH AND GC CONTENT TABLE----------------------------
GC_lengths_regions <- read_csv("/Users/manferg/Documents/ref/GC_lengths_regions.csv")
colnames(GC_lengths_regions)<-c("gene_id","length","GC")
GC_lengths<-as.data.frame(GC_lengths_regions)

#grep("ENSG00000251562.8",GC_lengths$gene_id)

#BED INTERSECTED FILES------------------------------------
#Tollervey - import intersected files genialis----------------------
tollervey.fi.li = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/tollervey/data/run/results/xlinks/intersected-genialis", pattern = ".xl.bed$", full.names = TRUE) #store file paths for each file in a list
tollervey.fi.li.names = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/tollervey/data/run/results/xlinks/intersected-genialis", pattern = "xl.bed$", full.names = FALSE) %>%
  gsub(".xl.bed","",.) %>%
  gsub("intersect_","",.)#extracts files names in a seprate list
tollervey.intersect.bed.li = list() #create empty list
for (i in 1:length(tollervey.fi.li)){
  temp = read.table(tollervey.fi.li[[i]])
  tollervey.intersect.bed.li[[i]] = temp
  tollervey.intersect.bed.li[[i]]<-dplyr::select(tollervey.intersect.bed.li[[i]], -V6,-V8,-V10,-V12) #remove variable columns
  colnames(tollervey.intersect.bed.li[[i]]) = c("seqname","start","end","score","region","gene_id","gene_name","biotype") #rename columns
  tollervey.intersect.bed.li[[i]]$gene_id <- sub(".", "NA", tollervey.intersect.bed.li[[i]]$gene_id)
  tollervey.intersect.bed.li[[i]]$gene_id <- sub("NANSG","ENSG", tollervey.intersect.bed.li[[i]]$gene_id)
    # Add `sample` column to the data frame
  }
names(tollervey.intersect.bed.li) <-tollervey.fi.li.names #adding sample names to list elements


for (i in 1:length(tollervey.intersect.bed.li)){
tollervey.intersect.bed.li[[i]]<- as.data.frame(tollervey.intersect.bed.li[[i]])
sample<-as.character(names(tollervey.intersect.bed.li[i])) # Create a new vector with sample names
tollervey.intersect.bed.li[[i]]$sample <- sample
tollervey.intersect.bed.li[[i]]<-left_join(tollervey.intersect.bed.li[[i]],GC_lengths,by="gene_id")
}



#G.Rot - import intersected files genialis----------------------
grot.fi.li = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/grot/prova/run-imp/results/xlinks/intersected/", pattern = ".xl.bed$", full.names = TRUE) #store file paths for each file in a list
grot.fi.li.names = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/grot/prova/run-imp/results/xlinks/intersected/", pattern = "xl.bed$", full.names = FALSE) %>%
  gsub(".xl.bed","",.) %>%
  gsub("intersect_","",.)#extracts files names in a seprate list
grot.intersect.bed.li = list() #create empty list
for (i in 1:length(grot.fi.li)){
  temp = read.table(grot.fi.li[[i]])
  grot.intersect.bed.li[[i]] = temp
  grot.intersect.bed.li[[i]]<-dplyr::select(grot.intersect.bed.li[[i]], -V6,-V8,-V10,-V12) #remove variable columns
  colnames(grot.intersect.bed.li[[i]]) = c("seqname","start","end","score","region","gene_id","gene_name","biotype") #rename columns
  grot.intersect.bed.li[[i]]$gene_id <- sub(".", "NA", grot.intersect.bed.li[[i]]$gene_id)
  grot.intersect.bed.li[[i]]$gene_id <- sub("NANSG","ENSG", grot.intersect.bed.li[[i]]$gene_id)
}
names(grot.intersect.bed.li) <-grot.fi.li.names #adding sample names to list elements


for (i in 1:length(grot.intersect.bed.li)){
  grot.intersect.bed.li[[i]]<- as.data.frame(grot.intersect.bed.li[[i]]) # Add `sample` column to the data frame
  sample<-as.character(names(grot.intersect.bed.li[i])) # Create a new vector with sample names
  grot.intersect.bed.li[[i]]$sample <- sample
  grot.intersect.bed.li[[i]]<-left_join(grot.intersect.bed.li[[i]],GC_lengths,by="gene_id")
}

#main df on BED intersected for Tollervey et al + Grot et al ---------------------------

#main df on BED intersected for Tollervey--------------------
tollervey.intresected.df<-do.call(rbind, tollervey.intersect.bed.li)
tollervey.intresected.df.chromosomes<-tollervey.intresected.df[grep("chr", tollervey.intresected.df$seqname),] #select only chromosomes / discard scaffolds 


#main df on BED intersected for Grot---------------------
grot.intresected.df<-do.call(rbind, grot.intersect.bed.li)
grot.intresected.df
grot.intresected.df.chromosomes<-grot.intresected.df[grep("chr", grot.intresected.df$seqname),] #select only chromosomes / discard scaffolds 




#main df----------------

#intresected.df<-dplyr::bind_rows(tollervey.intresected.df,grot.intresected.df)
#write.csv(intresected.df,file="/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/intresected.df.csv")
#intresected.chr.df<-dplyr::bind_rows(tollervey.intresected.df.chromosomes,grot.intresected.df.chromosomes)
#write.csv(intresected.chr.df,file="/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/intresected.chr.df.csv")

intresected.chr.df<-read_csv("/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/intresected.chr.df.csv")
intresected.df<-read_csv("/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/intresected.df.csv")

intresected.df<-intresected.df %>% as.data.frame() %>% dplyr::select(-X1)
intresected.chr.df<-intresected.chr.df %>% as.data.frame() %>% dplyr::select(-X1)


#METADATA----------------------


metadata<- read_csv("/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/meta_metadata.csv") %>% as.data.frame %>%  na.omit %>% column_to_rownames(var="meta_id")
rownames(metadata)

reorder_idx <- (unique(intresected.df$sample))
metadata<- metadata[reorder_idx,]
metadata$meta_id <- rownames(metadata)


intresected.df<-left_join(intresected.df, metadata, by=c("sample" = "meta_id")) 
intresected.chr.df<-left_join(intresected.chr.df, metadata, by=c("sample" = "meta_id"))


intresected.df.fil<- dplyr::select(intresected.df,-species,-technology, -study_id, -barcode) #filter metada columns to exclude
intresected.chr.df.fil<-dplyr::select(intresected.chr.df,-species,-technology, -study_id, -barcode)

#main df list----------------
main.li = split(intresected.df.fil,intresected.df.fil$sample) #transform df into list 
main.chr.li = split(intresected.chr.df.fil,intresected.chr.df.fil$sample) #transform df into list 














#number of xlink events per region--------------------
library(janitor)

#region dupes list------- 
# main.dupe.region.li =list()
# for (i in 1:length(main.li)){
# main.dupe.region.li[[i]]<-main.li[[i]] %>% get_dupes(region) #dupe counts (or number of xlinks) for each region
# }

#counts and perc (or number of xlinks) for each region------------------
# main.region.counts.li =list()
# for (i in 1:length(main.li)){
#   main.region.counts.li[[i]]<-main.li[[i]] %>%  group_by(region) %>% summarize(n=n()) %>%  mutate(perc = n * 100/ sum(n)) %>% arrange( .,desc(n)) 
# sample<-as.character(names(main.li[i])) # Create a new vector with sample names
# main.region.counts.li[[i]]$sample <- sample
# }

# main.region.counts.df<-do.call(rbind,main.region.counts.li) #convert back to df to plot

#CHR - counts (or number of xlinks) for each region------------------
chr.region.counts.li =list()
for (i in 1:length(main.chr.li)){
  chr.region.counts.li[[i]]<-main.chr.li[[i]] %>% group_by(region) %>% summarize(n=n()) %>%  mutate(perc = n * 100/ sum(n)) %>% arrange( .,desc(n)) 
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  chr.region.counts.li[[i]]$sample <- sample
}
chr.region.counts.df<-do.call(rbind,chr.region.counts.li) %>% as.data.frame()#convert back to df to plot

chr.region.counts.df<-left_join(chr.region.counts.df,metadata,by=c("sample" = "meta_id"))


#plot of xlink events per region------------------

# xlink.region<-ggplot(main.region.counts.df, aes(x = reorder(region, -n), y = n, fill=region)) + 
#   geom_bar(stat = "identity") +
#   #scale_color_manual(values=mypal.region) +
#   #scale_fill_manual(values=alpha(c(mypal.region))) +
#   facet_wrap(~sample,scales = "free") + #each facet has a different scale
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ggtitle("Cross-links events distribution") +
#   ylab("xlink events counts") +
#   theme_classic() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())

xlink.region.chr<-ggplot(chr.region.counts.df, aes(x = reorder(region, -n), y = n, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  facet_wrap(~sample,scales = "free") + #each facet has a different scale
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("xlink events counts") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(xlink.region.chr, filename = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/r-plots/xlink.region.chr.png", height = 6, width = 6)

xlink.region.chr.perc<-ggplot(chr.region.counts.df, aes(x = reorder(region, -perc), y = perc, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  facet_wrap(~sample,scales = "free") + #each facet has a different scale
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution -% ") +
  ylab("xlink events %") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(xlink.region.chr.perc, filename = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/r-plots/xlink.region.chr.perc.png", height = 6, width = 6)

xlink.region.chr.perc.model<-ggplot(chr.region.counts.df,aes(x = region, y = perc, fill=model)) + 
  geom_bar(stat = "identity",position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values=mypal.model) +
  scale_fill_manual(values=alpha(c(mypal.model))) +
  ggtitle("Cross-links events distribution -% ") +
  ylab("xlink events %") 


ggsave(xlink.region.chr.perc.model, filename = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/clip_metanalysis/r-plots/xlink.region.chr.perc.model.png", height = 6, width = 6)

#the number of scaffolds events is limited therefore there is no massive differrence in the number of cross-link in intergenic regions 
#between chr only and main dataset 

#========XLINK EVENTS PER GENE===========#--------

#----MAIN--------
#gene_name dupes list------- MAIN
main.dupe.gene_name.li =list()
for (i in 1:length(main.li)){
  main.dupe.gene_name.li[[i]]<-main.li[[i]] %>% get_dupes(gene_name) #dupe counts (or number of xlinks) for each gene
}

#counts (or number of xlinks) for each gene------------------MAIN
main.xlink.events.gene.li =list()
top.main.xlink.events.gene.li=list()
for (i in 1:length(main.li)){
  main.xlink.events.gene.li[[i]]<-main.li[[i]] %>% group_by(gene_name) %>% summarize(n=n()) %>% arrange( .,desc(n)) #counts (or number of xlinks) for each gene 
  sample<-as.character(names(main.li[i])) # Create a new vector with sample names
  main.xlink.events.gene.li[[i]]$sample <- sample
  top.main.xlink.events.gene.li[[i]]<-main.xlink.events.gene.li[[i]][2:10,] #exclude "None" and list top 20 genes
}

main.gene_name.counts.df<-do.call(rbind,main.xlink.events.gene.li) #convert back to df to plot
top.main.xlink.events.gene.df<-do.call(rbind,top.main.xlink.events.gene.li)


#plot of xlink events per gene------------MAIN
xlink.genes<-ggplot(top.main.xlink.events.gene.df, aes(x = reorder(gene_name, n), y = n, fill=gene_name)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  coord_flip()

#counts (or number of xlinks) for each gene/region------------------MAIN
main.xlink.events.gene.region.li =list()
top.main.xlink.events.gene.region.li=list()
for (i in 1:length(main.li)){
  main.xlink.events.gene.region.li[[i]]<-main.li[[i]] %>% group_by(gene_name,region) %>% summarize(n=n()) %>% arrange( .,desc(n)) #counts (or number of xlinks) for each gene 
  sample<-as.character(names(main.li[i])) # Create a new vector with sample names
  main.xlink.events.gene.region.li[[i]]$sample <- sample
  top.main.xlink.events.gene.region.li[[i]]<-main.xlink.events.gene.region.li[[i]][2:10,] #exclude "None" and list top 20 genes
}

main.gene_name.region.counts.df<-do.call(rbind,main.xlink.events.gene.region.li) #convert back to df to plot
top.main.xlink.events.gene.region.df<-do.call(rbind,top.main.xlink.events.gene.region.li)

#plot of xlink events per gene/region-----------MAIN

xlink.genes.region<-ggplot(top.main.xlink.events.gene.region.df, aes(x = reorder(gene_name, n), y = n, fill=region)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  coord_flip()


#----CHR--------

#gene_name dupes list------- CHR
main.dupe.gene_name.chr.li =list()
for (i in 1:length(main.chr.li)){
  main.dupe.gene_name.chr.li[[i]]<-main.chr.li[[i]] %>% get_dupes(gene_name) #dupe counts (or number of xlinks) for each gene
}

#counts (or number of xlinks) for each gene------------------CHR
xlink.events.gene.chr.li =list()
top.xlink.events.gene.chr.li=list()
for (i in 1:length(main.chr.li)){
  xlink.events.gene.chr.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name) %>% summarize(n=n()) %>% arrange( .,desc(n)) #counts (or number of xlinks) for each gene 
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  xlink.events.gene.chr.li[[i]]$sample <- sample
  top.xlink.events.gene.chr.li[[i]]<-xlink.events.gene.chr.li[[i]][2:10,] #exclude "None" and list top 20 genes
}

xlink.events.gene.chr.df<-do.call(rbind,xlink.events.gene.chr.li) #convert back to df to plot
top.xlink.events.gene.chr.df<-do.call(rbind,top.xlink.events.gene.chr.li)


#plot of xlink events per gene------------CHR
xlink.genes.chr<-ggplot(top.xlink.events.gene.chr.df, aes(x = reorder_within(gene_name,n,sample), y = n, fill=gene_name)) + 
  geom_bar(stat = "identity") +
  scale_x_reordered() +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#counts (or number of xlinks) for each gene/region------------------CHR
xlink.events.gene.region.chr.li =list()
top.xlink.events.gene.region.chr.li=list()
for (i in 1:length(main.chr.li)){
  xlink.events.gene.region.chr.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,region) %>% summarize(n=n()) %>% arrange( .,desc(n)) #counts (or number of xlinks) for each gene.region. 
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  xlink.events.gene.region.chr.li[[i]]$sample <- sample
  top.xlink.events.gene.region.chr.li[[i]]<-xlink.events.gene.region.chr.li[[i]][2:10,]  #exclude "None" and list top 20 gene.region.s
  top.xlink.events.gene.region.chr.li[[i]]$gene_name <- factor(top.xlink.events.gene.region.chr.li[[i]]$gene_name , levels = top.xlink.events.gene.region.chr.li[[i]]$gene_name[order(top.xlink.events.gene.region.chr.li[[i]]$n)])

}

xlink.events.gene.region.chr.df<-do.call(rbind,xlink.events.gene.region.chr.li) #convert back to df to plot
top.xlink.events.gene.region.chr.df<-do.call(rbind,top.xlink.events.gene.region.chr.li)

#test
p<-intresected.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% filter(region == "intron")
dim(p)

view(top.xlink.events.gene.region.chr.df)
view(xlink.events.gene.region.chr.df)

#plot of xlink events per gene/region-----------CHR

xlink.genes.chr.region<-ggplot(top.xlink.events.gene.region.chr.df,aes(reorder_within(gene_name,n,sample), n, fill=region)) + 
  facet_wrap(~sample,scales = "free") +
  scale_x_reordered() +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  geom_bar(stat = "identity") +
  xlab("genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  theme_classic() +
  coord_flip()


#===gc=====#---------------


xlink.events.gene.region.chr.GC.li =list()
top.xlink.events.gene.region.chr.GC.li=list()
for (i in 1:length(main.chr.li)){
  xlink.events.gene.region.chr.GC.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,region,GC) %>% summarize(n=n()) %>% arrange( .,desc(n)) #counts (or number of xlinks) for each gene.region. 
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  xlink.events.gene.region.chr.GC.li[[i]]$sample <- sample
  top.xlink.events.gene.region.chr.GC.li[[i]]<-xlink.events.gene.region.chr.GC.li[[i]][1:10,]  #exclude "None" and list top 20 gene.region.s
  top.xlink.events.gene.region.chr.GC.li[[i]]$gene_name <- factor(top.xlink.events.gene.region.chr.GC.li[[i]]$gene_name , levels = top.xlink.events.gene.region.chr.GC.li[[i]]$gene_name[order(top.xlink.events.gene.region.chr.GC.li[[i]]$n)])
  
}




xlink.events.gene.region.chr.GC.df<-do.call(rbind,xlink.events.gene.region.chr.GC.li) #convert back to df to plot
top.xlink.events.gene.region.chr.GC.df<-do.call(rbind,top.xlink.events.gene.region.chr.GC.li)

#test
#p<-intresected.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% filter(region == "intron")
#dim(p)

#view(top.xlink.events.gene.region.chr.df)
#view(xlink.events.gene.region.chr.df)

#plot of xlink events per gene/region-----------CHR----GC

xlink.genes.chr.GC.region<-ggplot(top.xlink.events.gene.region.chr.GC.df,aes(reorder_within(gene_name,n,sample), n, fill=GC)) + 
  facet_wrap(~sample,scales = "free") +
  scale_fill_gradient(low = "#f7f7f7", high = "#99000d") +
  scale_x_reordered() +
  geom_bar(stat = "identity") +
  xlab("genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  theme_classic() +
  coord_flip()



#===gene length normalisation=====#---------------

xlink.events.gene.region.chr.norm.li =list()
top.xlink.events.gene.region.chr.norm.li=list()
xlink.chr.normalised.li=list()
for (i in 1:length(main.chr.li)){
  xlink.events.gene.region.chr.norm.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,region,GC,lenght) %>% summarize(n=n()) %>% arrange( .,desc(n)) #counts (or number of xlinks) for each gene.region. 
  xlink.chr.normalised.li[[i]]<-xlink.events.gene.region.chr.norm.li[[i]] %>% mutate( .,xlink.enrichment=n/lenght) %>% arrange( .,desc(xlink.enrichment))
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  xlink.chr.normalised.li[[i]]$sample <- sample
  top.xlink.events.gene.region.chr.norm.li[[i]]<-xlink.chr.normalised.li[[i]][1:10,]  #exclude "None" and list top 20 gene.region.s
  top.xlink.events.gene.region.chr.norm.li[[i]]$gene_name <- factor(top.xlink.events.gene.region.chr.norm.li[[i]]$gene_name , levels = top.xlink.events.gene.region.chr.norm.li[[i]]$gene_name[order(top.xlink.events.gene.region.chr.norm.li[[i]]$xlink.enrichment)])
  
}




xlink.chr.normalised.df<-do.call(rbind,xlink.chr.normalised.li) #convert back to df to plot
top.xlink.events.gene.region.chr.norm.df<-do.call(rbind,top.xlink.events.gene.region.chr.norm.li)

#test
#p<-intresected.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% filter(region == "intron")
#dim(p)

#view(top.xlink.events.gene.region.chr.df)
#view(xlink.events.gene.region.chr.df)

#plot of xlink events per gene/region-----------CHR

xlink.genes.chr.norm.region<-ggplot(top.xlink.events.gene.region.chr.norm.df,aes(reorder_within(gene_name,xlink.enrichment,sample), xlink.enrichment, fill=GC)) + 
  facet_wrap(~sample,scales = "free") +
  #scale_color_manual(values=mypal.region) +
  #scale_fill_manual(values=alpha(c(mypal.region))) +
  scale_fill_gradient(low = "#f7f7f7", high = "#99000d") +
  scale_x_reordered() +
  geom_bar(stat = "identity") +
  xlab("genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  theme_classic() +
  coord_flip()



















#===========================SCORE PER GENE=================#



score.events.gene.region.chr.norm.li =list()
top.score.events.gene.region.chr.norm.li=list()
score.chr.normalised.li=list()
for (i in 1:length(main.chr.li)){
  score.events.gene.region.chr.norm.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,region,GC,lenght,score) %>% summarize(n=n()) %>% arrange( .,desc(score))#counts (or number of scores) for each gene.region. 
  
}

for (i in 1:length(score.events.gene.region.chr.norm.li)){
  score.chr.normalised.li[[i]]<-score.events.gene.region.chr.norm.li[[i]]%>% mutate(score.enrichment=score/lenght) %>% arrange( .,desc(score.enrichment))
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  score.chr.normalised.li[[i]]$sample <- sample
  top.score.events.gene.region.chr.norm.li[[i]]<-score.chr.normalised.li[[i]][1:10,]
  top.score.events.gene.region.chr.norm.li[[i]]$gene_name <- factor(top.score.events.gene.region.chr.norm.li[[i]]$gene_name , levels = top.score.events.gene.region.chr.norm.li[[i]]$gene_name[order(top.score.events.gene.region.chr.norm.li[[i]]$score.enrichment)])
  
}




score.chr.normalised.df<-do.call(rbind,score.chr.normalised.li) #convert back to df to plot
top.score.events.gene.region.chr.norm.df<-do.call(rbind,top.score.events.gene.region.chr.norm.li)

#test
#p<-intresected.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% filter(region == "intron")
#dim(p)

#view(top.score.events.gene.region.chr.df)
#view(score.events.gene.region.chr.df)

#plot of score events per gene/region-----------CHR

score.genes.chr.norm.region<-ggplot(top.score.events.gene.region.chr.norm.df,aes(reorder_within(gene_name,score.enrichment,sample), score.enrichment, fill=region)) + 
  facet_wrap(~sample,scales = "free") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  #scale_fill_gradient(low = "#f7f7f7", high = "#99000d") +
  scale_x_reordered() +
  geom_bar(stat = "identity") +
  xlab("genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  theme_classic() +
  coord_flip()















#score count per region--------------------
intresected.score.region.df<- intresected.df %>% dplyr::select(score,region,sample) #smaller df for plot
score<-ggplot(intresected.score.region.df, aes(x = reorder(region, -score), y = score, fill=region)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~sample,scales = "free") + #each facet has a different scale
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links score distribution") +
  ylab("total score counts") 
#counts (or number of xlinks) for each gene------------------MAIN
main.score.events.gene.li =list()
top.main.score.events.gene.li=list()
for (i in 1:length(main.li)){
  main.score.events.gene.li[[i]]<-main.li[[i]] %>% group_by(gene_name,score,region) %>% summarize(n=n()) %>% arrange( .,desc(score)) #counts (or number of scores) for each gene 
  sample<-as.character(names(main.li[i])) # Create a new vector with sample names
  main.score.events.gene.li[[i]]$sample <- sample
  top.main.score.events.gene.li[[i]]<-main.score.events.gene.li[[i]][2:10,] #exclude "None" and list top 20 genes
}

main.gene_name.counts.df<-do.call(rbind,main.score.events.gene.li) #convert back to df to plot
top.main.score.events.gene.df<-do.call(rbind,top.main.score.events.gene.li)


#gene_name dupes list------- CHR
main.dupe.gene_name.chr.li =list()
for (i in 1:length(main.chr.li)){
  main.dupe.gene_name.chr.li[[i]]<-main.chr.li[[i]] %>% get_dupes(gene_name) #dupe counts (or number of scores) for each gene
}

#counts (or number of scores) for each gene------------------CHR
score.events.gene.chr.li =list()
top.score.events.gene.chr.li=list()
for (i in 1:length(main.chr.li)){
  score.events.gene.chr.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,score,region) %>% summarize(n=n()) %>% arrange( .,desc(score)) #counts (or number of scores) for each gene 
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  score.events.gene.chr.li[[i]]$sample <- sample
  top.score.events.gene.chr.li[[i]]<-score.events.gene.chr.li[[i]][2:10,] #exclude "None" and list top 20 genes
}

score.events.gene.chr.df<-do.call(rbind,score.events.gene.chr.li) #convert back to df to plot
top.score.events.gene.chr.df<-do.call(rbind,top.score.events.gene.chr.li)




#plot of score events per gene
score.genes<-ggplot(top.main.score.events.gene.df, aes(x = reorder(gene_name, score), y = score, fill=gene_name)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  coord_flip()



score.genes.chr<-ggplot(top.score.events.gene.chr.df, aes(x = reorder(gene_name, score), y = score, fill=gene_name)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  coord_flip()














#=============================================================================
#genome level xlinks counts/region per chromosomes-------------------------------------- 

#counts (or number of xlinks) for each chromosome------------------
main.xlink.genome.li =list()
for (i in 1:length(main.li)){
  main.xlink.genome.li[[i]]<-main.li[[i]] %>%  group_by(seqname,region) %>% summarize(n=n()) %>% arrange( .,desc(seqname))
  sample<-as.character(names(main.li[i])) # Create a new vector with sample names
  main.xlink.genome.li[[i]]$sample <- sample
}

main.xlink.genome.df<-do.call(rbind,main.xlink.genome.li) #convert back to df to plot

#CHR - counts (or number of xlinks) for each chromosome------------------
chr.xlink.genome.li =list()
for (i in 1:length(main.chr.li)){
  chr.xlink.genome.li[[i]]<-main.chr.li[[i]] %>% group_by(seqname,region) %>% summarize(n=n()) %>% arrange( .,desc(seqname))
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  chr.xlink.genome.li[[i]]$sample <- sample
}

chr.xlink.genome.df<-do.call(rbind,chr.xlink.genome.li) #convert back to df to plot


# dupe count list for chromosomes--------------
chr.xlink.duped.genome.li =list()
for (i in 1:length(main.chr.li)){
  chr.xlink.duped.genome.li[[i]]<-main.chr.li[[i]] %>% get_dupes(seqname) #dupe counts (or number of xlinks) for each gene
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  chr.xlink.duped.genome.li[[i]]$sample <- sample
}

chr.xlink.duped.genome.df<-do.call(chr.xlink.duped.genome.li)

#xlinks counts at genome level plots--------------------------------------
xlink.events.xlink.genome.main<-ggplot(main.xlink.genome.df, aes(x = reorder(seqname,n), y = n, fill=region)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links xlink distribution across the genome") +
  ylab("events counts") + 
  heme_classic()+
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  coord_flip()


xlink.events.xlink.genome.chr<-ggplot(chr.xlink.genome.df, aes(x = reorder(seqname,n), y = n, fill=region)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links xlink distribution across the genome") +
  ylab("events counts") + 
  theme_classic()+
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  coord_flip()

#genome level SCORE counts/region per chromosomes-------------------------------------- 

#score for each chromosome------------------
main.score.genome.li =list()
for (i in 1:length(main.li)){
  main.score.genome.li[[i]]<-main.li[[i]] %>%  group_by(score,seqname,region) %>% summarize(n=n()) %>% arrange( .,desc(score))
  sample<-as.character(names(main.li[i])) # Create a new vector with sample names
  main.score.genome.li[[i]]$sample <- sample
}

main.score.genome.df<-do.call(rbind,main.score.genome.li) #convert back to df to plot

#CHR - counts (or number of scores) for each chromosome------------------
chr.score.genome.li =list()
for (i in 1:length(main.chr.li)){
  chr.score.genome.li[[i]]<-main.chr.li[[i]] %>% group_by(score,seqname,region) %>% summarize(n=n()) %>% arrange( .,desc(score))
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  chr.score.genome.li[[i]]$sample <- sample
}

chr.score.genome.df<-do.call(rbind,chr.score.genome.li) #convert back to df to plot

#score counts at genome level plots--------------------------------------


xlink.events.score.genome.main<-ggplot(main.score.genome.df, aes(x = reorder(seqname,score), y = score, fill=region)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links score distribution across the genome") +
  ylab("events counts") + 
  theme_classic()+
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  coord_flip()


xlink.events.score.genome.chr<-ggplot(chr.score.genome.df, aes(x = reorder(seqname,score), y = score, fill=region)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links score distribution across the genome") +
  ylab("events counts") + 
  coord_flip()

#=============================================================================













# import summary tables----------------------
#old summary tables in "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/tollervey/data/run/results/xlinks/summary"
#New summary tables with regions.gtf
tab.li = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/tollervey/data/run/results/xlinks/summary-genialis", pattern = "_summary.tab$", full.names = TRUE) #store file paths for each file in a list
fi.li.names = list.files(path = "/Volumes/lab-luscomben/home/users/manferg/projects/nf/clip/tollervey/data/run/results/xlinks/summary-genialis", pattern = "_summary.tab$", full.names = FALSE) %>%
  gsub(".xl.bed_summary.tab","",.) %>%
  gsub("sorted_","",.) #extracts files names in a seprate list
summary.li = list() #create empty list
for (i in 1:length(tab.li)){
  temp = read_delim(tab.li[[i]],"\t", escape_double = FALSE, trim_ws = TRUE)
  summary.li[[i]] = temp
}
names(summary.li) <-fi.li.names #adding sample names to list elements

# add column with sample names
for (i in 1:length(summary.li)){
  summary.li[[i]]<- as.data.frame(summary.li[[i]])
  sample<- as.character(names(summary.li[i])) # Create a new vector with sample names
  summary.li[[i]]$sample <- sample  # Add `sample` column to the data frame
}


df<-do.call(rbind, summary.li)


#tidy dataset--------------------
df.tidy<-melt(data = df, id.vars = c("type", "length", "length %", "sites #", "sites %", "sites enrichment", "events #","events %", "events enrichment"), measure.vars = c("sample"))
df.tidy<-dplyr::select(df.tidy, -variable)
colnames(df.tidy) <- c("type", "length", "length.p", "sites.n", "sites.p", "sites.enrichment", "events.n","events.p","events.enrichment", "sample_name")
head(df.tidy)


#regions x counts-------------------------
xlink.events<-ggplot(df.tidy, aes(x = reorder(type, -events.n), y = events.n, fill=sample_name)) + 
  geom_bar(stat = "identity") +
  xlab("regions") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links event distribution") +
  ylab("counts")


xlink.events.log<-ggplot(df.tidy, aes(x = reorder(type, -events.n), y = events.n, fill=sample_name)) + 
  geom_bar(stat = "identity") +
  scale_y_log10() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links event distribution") +
  xlab("regions") +
  ylab("counts [log10]")


#samples x counts------------------

xlink.event<-ggplot(df.tidy, aes(x = sample_name, y = events.n, fill=type)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-link distribution") +
  coord_flip() +
  ylab("events/xlinks number")

xlink.event.log<-ggplot(df.tidy, aes(x = sample_name, y = events.n, fill=type)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-link distribution") +
  scale_y_log10() +
  coord_flip() +
  xlab("labels") + 
  ylab("counts [log10]")


#samples x events perc

xlink.event.perc<-ggplot(df.tidy, aes(x = sample_name, y = events.p, fill=type)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-link distribution") +
  coord_flip() +
  ylab("events/xlinks %")

#regions x events perc
regions.perc<-ggplot(df.tidy, aes(x = reorder(type, -events.p), y = events.p)) + 
  geom_col(position = "stack", aes(fill= sample_name)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Events number") +
  coord_flip() +
  xlab("labels") + 
  ylab("counts [log10]")

#samples x enriched 

events.enrich<-ggplot(df.tidy, aes(x = sample_name, y = events.enrichment, fill=type)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-link distribution") +
  coord_flip() +
  ylab("events/xlinks number")

sites.enrich<-ggplot(df.tidy, aes(x = sample_name, y = sites.enrichment, fill=type)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-link distribution") +
  coord_flip() +
  xlab("labels") + 
  ylab("counts")





# import intersected files new (gtf.regions intersect)-------------------------------


# filter xlink bed files-------------------------------
view(bed.li[[1]])

x<-bed.li[[1]][seqnames(bed.li[[1]]== "^chr")]
view(x)
x<-bed.li[[1]][seqnames(bed.li[[1]]) == "chr1"]
x<-bed.li[[1]] %>% select(grepl('^chr', seqnames))


### Object cleaning
xl.li = list()
for (i in 1:length(bed.li)){
  temp = select 
    xl.li[[i]] = temp
}


x<-bed.li[[1]]

# Import annotation file (GTF)---------------------------

annotation.file = "/Volumes/lab-luscomben/home/users/manferg/ref/gencode.v36.annotation.gtf"
anno = import(annotation.file, format = "GTF") #import GTF as GRange object


head(anno)
#view(anno)
class(anno) #GRanges 


##read gtf file as data frame
#gtf <- readGFF("/Volumes/lab-luscomben/home/users/manferg/ref/gencode.v36.annotation.gtf")
#save gtf as RDS
#saveRDS(gtf, file = "gencode.v36.annotation.gtf.rds")
#restore the object
#gtf.rds<-readRDS(file = "gencode.v36.annotation.gtf.rds")


### Filter feature level annotation--------------------------------------
anno = anno[anno$level != 3] #filter in only gene support level 1 and 2 

view(anno)
### Filter transcript level annotation
anno$transcript_support_level[is.na(anno$transcript_support_level)] = 0 #assign level 0 to NA 
anno$transcript_support_level[anno$transcript_support_level == "NA"] = 10 #assign level 10 to NA in italics
anno = anno[anno$transcript_support_level == 0
            | anno$transcript_support_level == 1
            | anno$transcript_support_level == 2
            | anno$transcript_support_level == 3 ]

### Create txdb databse from filtered annotations
#The GenomicFeatures package uses TxDb objects to store transcript metadata

anno.db = makeTxDbFromGRanges(anno)

class(anno.db) #TxDb
str(anno.db)


### 
### Genomic targets overview
### 

### Select annotated genes
gns = genes(anno.db)
view(gns)
view(anno)
idx = match(gns$gene_id, anno$gene_id) #match genes 
elementMetadata(gns) = cbind(elementMetadata(gns), elementMetadata(anno)[idx,])# retrieve metadata to genes

### Gene perspective
df = data.frame(gene_type = subsetByOverlaps(gns, x)$gene_type) %>% #subset ranges from annotation that overlaps with ranges in xlink and extract gene_types
  table %>% 
  as.data.frame

ggplot(df, aes(x = reorder(df[,1], -df[,2]), y = df[,2])) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Gene target fractions") + 
  scale_y_log10() + 
  xlab("labels") + 
  ylab("counts [log10]")


### Binding site perspective
bs.tg = subsetByOverlaps(bs.filter,gns)
bs.tg = bs.tg[countOverlaps(bs.tg,targets) == 1]
overlaps = findOverlaps(bs.tg, targets) %>% as.data.frame
df = data.frame(gene_type = targets$gene_type[overlaps$subjectHits]) %>%
  table %>%
  as.data.frame

ggplot(df, aes(x = reorder(df[,1], -df[,2]), y = df[,2])) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Binding site fractions") +
  scale_y_log10() +
  xlab("labels") + 
  ylab("counts [log10]")