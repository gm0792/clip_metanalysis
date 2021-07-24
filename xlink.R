# ICLIP ANALYSIS 
# 2nd March 2021
setwd("/Users/manferg/clip_metanalysis/")


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
my.pal.lib<-c("output_reads"= "#39373B")
mypal.region<-c("intron" = "#454241","CDS" ="#F0421C", "intergenic" ="#DDD3D1", "ncRNA" = "#3DB75D", "UTR5" = "#3DB7E6", "UTR3"= "#D644C7")
mypal.model<-c("FTLD-TDP_human_brain" = "#6E6EDB","Healthy_human_brain" = "#1A7497","SH-SY5Y"="#6DD6F2","293Flp" ="#E86E68", "hES" ="#41EB81")


#SAMPLE ORDER------------------------

reorder_sample_idx <- c("grot_293fl_1", "grot_293fl_2", "tollervey_brain1","tollervey_brain2","tollervey_brain6.high","tollervey_brain7.low","tollervey_brain3",
                        "tollervey_brain4","tollervey_brain5","tollervey_ESC","tollervey_SHSY5Y1a","tollervey_SHSY5Y1b","tollervey_SHSY5Y2","tollervey_SHSY5Y3","tollervey_SHSY5Y_cyt","tollervey_SHSY5Y_nucl")





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

tollervey.sample.order<-c("tollervey_brain1","tollervey_brain2","tollervey_brain6.high","tollervey_brain7.low","tollervey_brain3",
                          "tollervey_brain4","tollervey_brain5","tollervey_ESC","tollervey_SHSY5Y1a","tollervey_SHSY5Y1b",
                          "tollervey_SHSY5Y2","tollervey_SHSY5Y3","tollervey_SHSY5Y_cyt","tollervey_SHSY5Y_nucl")


#setting df levels to reorder samples
#input reads and output reads 
tollervey.dedup.df<-tollervey.dedup.df[match(tollervey.sample.order, tollervey.dedup.df$sample),]
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
  coord_flip() + 
  theme_classic()

ggsave(dedup.in.out.plot, filename = "/Users/manferg/clip_metanalysis/r-plots/dedup.in.out.plot.png", height = 6, width = 6)

#reads ratio
dedup.ratio.plot<-ggplot(df.dedup.ratio) +
  geom_bar(stat = "identity", aes(x=sample, y=count, fill=read), position = "dodge") +
  scale_color_manual(values=mypal.dedup.ratio) +
  scale_fill_manual(values=alpha(c(mypal.dedup.ratio))) +
  coord_flip() +
theme_classic()
ggsave(dedup.ratio.plot, filename = "/Users/manferg/clip_metanalysis/r-plots/dedup.ratio.plot.png", height = 6, width = 6)


#INSERT THRESHOLD! 
threshold.dedup = 10
dedup.ratio.plot<-dedup.ratio.plot + geom_hline(yintercept = threshold.dedup , linetype = "dashed", color ="black") 
dedup.ratio.plot

#Library size plot-----------
#output reads (unique reads are considered as library size)

df.dedup.in<-df.dedup.in.out[df.dedup.in.out$read == "input_reads",]
df.dedup.out.lib<-df.dedup.in.out[df.dedup.in.out$read == "output_reads",] #library sizes

library.size.plot<-ggplot(df.dedup.out.lib)+
  geom_bar(aes(x= sample,y=count,fill=read), stat ='identity') +
  scale_color_manual(values=my.pal.lib) +
  scale_fill_manual(values=my.pal.lib) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic() +
  coord_flip() +
  ggtitle("Library size") 
ggsave(library.size.plot, filename = "/Users/manferg/clip_metanalysis/r-plots/library.size.plot.png", height = 6, width = 6)



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

intresected.chr.df<-read_csv("/Users/manferg/clip_metanalysis/intresected.chr.df.csv")
intresected.df<-read_csv("/Users/manferg/clip_metanalysis/intresected.df.csv")

intresected.df<-intresected.df %>% as.data.frame() %>% dplyr::select(-X1)
intresected.chr.df<-intresected.chr.df %>% as.data.frame() %>% dplyr::select(-X1)


#METADATA----------------------

write.csv()

metadata<- read_csv("/Users/manferg/clip_metanalysis/meta_metadata.csv") %>% as.data.frame %>%  na.omit %>% column_to_rownames(var="meta_id")
rownames(metadata)
colnames(metadata)






metadata<- metadata[reorder_sample_idx ,]
metadata$meta_id <- rownames(metadata)


intresected.d<-left_join(intresected.df, metadata, by=c("sample" = "meta_id")) 
intresected.chr.df<-left_join(intresected.chr.df, metadata, by=c("sample" = "meta_id"))



intresected.df.fil<- dplyr::select(intresected.d,-species,-technology, -study_id, -barcode) #filter metada columns to exclude
intresected.chr.df.fil<-dplyr::select(intresected.chr.df,-species,-technology, -study_id, -barcode)

head(intresected.chr.df.fil)


#main df list----------------
main.li = split(intresected.df.fil,intresected.d$sample) #transform df into list 
main.chr.li = split(intresected.chr.df.fil,intresected.chr.df$sample) #transform df into list 


main.chr.li<-main.chr.li[reorder_sample_idx] #reprder liste elements 
names(main.chr.li) #check order




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

#reorder levels as samples order
chr.region.counts.df$sample <- factor(chr.region.counts.df$sample , levels=unique(chr.region.counts.df$sample ))



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

xlink.region.chr<-ggplot(chr.region.counts.df, aes(x = reorder_within(region, -n, sample), y = n, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  facet_wrap(~sample + model,scales = "free") + #each facet has a different scale
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("xlink events counts") + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
         axis.ticks.x=element_blank()) #+
  # theme(strip.background =element_rect(fill= c("black")))+
  # theme(strip.text = element_text(colour = 'white'))
  #       
ggsave(xlink.region.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.region.chr.png", height = 8, width = 9)

xlink.region.chr.perc<-ggplot(chr.region.counts.df, aes(x = reorder(region, -perc), y = perc, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  facet_wrap(~sample + model,scales = "free") + #each facet has a different scale
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution -% ") +
  ylab("xlink events %") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(xlink.region.chr.perc, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.region.chr.perc.png", height = 8, width = 9)

xlink.region.chr.perc.model<-ggplot(chr.region.counts.df,aes(x = region, y = perc, fill=model)) + 
  geom_bar(stat = "identity",position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values=mypal.model) +
  scale_fill_manual(values=alpha(c(mypal.model))) +
  ggtitle("Cross-links events distribution -% ") +
  ylab("xlink events %") + 
  theme_classic()


ggsave(xlink.region.chr.perc.model, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.region.chr.perc.model.png", height = 6, width = 6)

#the number of scaffolds events is limited therefore there is no massive differrence in the number of cross-link in intergenic regions 
#between chr only and main dataset 

#========XLINK EVENTS PER GENE===========#--------

#----MAIN--------
#gene_name dupes list------- MAIN
# main.dupe.gene_name.li =list()
# for (i in 1:length(main.li)){
#   main.dupe.gene_name.li[[i]]<-main.li[[i]] %>% get_dupes(gene_name) #dupe counts (or number of xlinks) for each gene
# }

#counts (or number of xlinks) for each gene------------------MAIN
# main.xlink.events.gene.li =list()
# top.main.xlink.events.gene.li=list()
# for (i in 1:length(main.li)){
#   main.xlink.events.gene.li[[i]]<-main.li[[i]] %>% group_by(gene_name) %>% summarize(n=n()) %>% arrange( .,desc(n)) #counts (or number of xlinks) for each gene 
#   sample<-as.character(names(main.li[i])) # Create a new vector with sample names
#   main.xlink.events.gene.li[[i]]$sample <- sample
#   top.main.xlink.events.gene.li[[i]]<-main.xlink.events.gene.li[[i]][1:10,] #exclude "None" and list top 20 genes
# }
# 
# main.gene_name.counts.df<-do.call(rbind,main.xlink.events.gene.li) #convert back to df to plot
# top.main.xlink.events.gene.df<-do.call(rbind,top.main.xlink.events.gene.li)


#plot of xlink events per gene------------MAIN
# xlink.genes<-ggplot(top.main.xlink.events.gene.df, aes(x = reorder(gene_name, n), y = n, fill=gene_name)) + 
#   geom_bar(stat = "identity") +
#   xlab("genes") +
#   facet_wrap(~sample,scales = "free") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ggtitle("Cross-links events distribution") +
#   ylab("events counts") + 
#   coord_flip()

#counts (or number of xlinks) for each gene/region------------------MAIN
# main.xlink.events.gene.region.li =list()
# top.main.xlink.events.gene.region.li=list()
# for (i in 1:length(main.li)){
#   main.xlink.events.gene.region.li[[i]]<-main.li[[i]] %>% group_by(gene_name,region) %>% summarize(n=n()) %>% arrange( .,desc(n)) #counts (or number of xlinks) for each gene 
#   sample<-as.character(names(main.li[i])) # Create a new vector with sample names
#   main.xlink.events.gene.region.li[[i]]$sample <- sample
#   top.main.xlink.events.gene.region.li[[i]]<-main.xlink.events.gene.region.li[[i]][2:10,] #exclude "None" and list top 20 genes
# }
# 
# main.gene_name.region.counts.df<-do.call(rbind,main.xlink.events.gene.region.li) #convert back to df to plot
# top.main.xlink.events.gene.region.df<-do.call(rbind,top.main.xlink.events.gene.region.li)
# 
# #plot of xlink events per gene/region-----------MAIN
# 
# xlink.genes.region<-ggplot(top.main.xlink.events.gene.region.df, aes(x = reorder(gene_name, n), y = n, fill=region)) + 
#   geom_bar(stat = "identity") +
#   xlab("genes") +
#   facet_wrap(~sample,scales = "free") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ggtitle("Cross-links events distribution") +
#   ylab("events counts") + 
#   coord_flip()


#----CHR--------

#gene_name dupes list------- CHR
# main.dupe.gene_name.chr.li =list()
# for (i in 1:length(main.chr.li)){
#   main.dupe.gene_name.chr.li[[i]]<-main.chr.li[[i]] %>% get_dupes(gene_name) #dupe counts (or number of xlinks) for each gene
# }




#counts (or number of xlinks) for each gene/region------------------CHR
xlink.events.gene.region.chr.li =list()
top.xlink.events.gene.region.chr.li=list()
for (i in 1:length(main.chr.li)){
  xlink.events.gene.region.chr.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,region,GC,length) %>% summarize(n=n()) %>% arrange( .,desc(n)) %>% as.data.frame()#counts (or number of xlinks) for each gene.region. 
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  xlink.events.gene.region.chr.li[[i]]$sample <- sample
  top.xlink.events.gene.region.chr.li[[i]]<-xlink.events.gene.region.chr.li[[i]][2:10,]  #exclude "None" and list top 20 gene.region.s
  top.xlink.events.gene.region.chr.li[[i]]$gene_name <- factor(top.xlink.events.gene.region.chr.li[[i]]$gene_name , levels = top.xlink.events.gene.region.chr.li[[i]]$gene_name[order(top.xlink.events.gene.region.chr.li[[i]]$n)])

}

xlink.events.chr.df<-do.call(rbind,xlink.events.gene.region.chr.li) #convert back to df to plot



#reorder levels as samples order
xlink.events.chr.df$sample <- factor(xlink.events.chr.df$sample , levels=unique(xlink.events.chr.df$sample ))


top.xlink.events.chr.df<-do.call(rbind,top.xlink.events.gene.region.chr.li)
top.xlink.events.chr.df<-left_join(top.xlink.events.chr.df,metadata,by=c("sample" = "meta_id"))
#reorder levels as samples order
top.xlink.events.chr.df$sample <- factor(top.xlink.events.chr.df$sample , levels=unique(top.xlink.events.chr.df$sample ))

#data exploration xlink df----------

dim(xlink.events.chr.df) #390292

length(unique(xlink.events.chr.df$gene_name)) #35248 unique gene_name

xlink.events.unique.gene<-xlink.events.chr.df %>% group_by(gene_name,sample) %>% summarize(n=sum(n)) %>% arrange( .,desc(n)) %>% as.data.frame() #total xlinks per gene (sum of all region)
dim(xlink.events.unique.gene)

#p<-xlink.events.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% summarize(n=sum(n))
#p<-xlink.events.unique.gene %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1")




# normalize raw x-links per Library size ## NO NEED TO NORMALISE X-LINK DATA PER LIBRARY SIZE?

# head(df.dedup.out.lib)
# head(xlink.events.chr.df)
# xlink.events.chr.lb.norm.df<-left_join(xlink.events.unique.gene, df.dedup.out.lib, by = "sample")
# xlink.events.chr.lb.norm.df<-xlink.events.chr.lb.norm.df %>% mutate(ls.factor=count/1e+06) %>% mutate(n.norm = n/ls.factor)
# 
# tidy<-xlink.events.chr.lb.norm.df %>% dplyr::select(sample, gene_name,n.norm) %>% arrange( .,desc(n.norm)) %>% filter( .,gene_name != "None")
# 
# top_genes<-tidy %>% dplyr::select(gene_name) %>% unique()
# 
# spread.df<-spread(tidy,sample,n.norm)







#spread.df.na.rm<-spread.df %>% na.omit()
#dim(spread.df.na.rm)

#HEATH-MAP---------THIS HEATHMAP IS PLOTTED ON LIBRARY SIZE NORM X-LINKS
spread.df[is.na(spread.df)] = 0 #transform NA values into 0 

head(spread.df)

mat<-as.matrix(spread.df[,2:17])#pheatmap only takes the numeric matrix object as input. So, we need to transfer the numeric part of the data frame to a matrix by removing the first 5 columns of categorical data.
rownames(mat)<-as.character(spread.df$gene_name)
mat_scale = scale(mat)

dim(mat_scale)

library(pheatmap)
mn_select <- top_genes[1:30,]
mat <- mat[,mn_reorder_idx]
mat <-mat_scale[mn_select,] 
mn_reorder_idx <- match(reorder_sample_idx, colnames(mat_scale), nomatch = 0) #REORDER SAMPLES ORDER 


#heathmap annotation-------
ann_df = data.frame(metadata$meta_id, metadata$model, row.names = TRUE)
colnames(ann_df)<-c("model")
rownames(ann_df) <- colnames(mat_scale) # name matching

anno_colors = list(mypal.model)


#heath_map
#xl_matrix_col <- pheatmap(mat, annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, main = paste("Top 20 xlinked genes"),fontsize = 10)
xl_matrix_row <- pheatmap(mat,scale = "row", annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE, main = paste("whole dataset"),fontsize = 10)
ggsave(xl_matrix_row, filename = "/Users/manferg/clip_metanalysis/r-plots/xl_matrix_row.png", height = 10, width = 10)


xl_matrix_row_30 <- pheatmap(mat,scale = "row", annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, main = paste("Top 30 xlinked genes"),fontsize = 10)
ggsave(xl_matrix_row_30, filename = "/Users/manferg/clip_metanalysis/r-plots/xl_matrix_row_30.png", height = 10, width = 10)






# #test
# p<-intresected.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% filter(region == "intron")
# dim(p)

# view(top.xlink.events.gene.region.chr.df)
# view(xlink.events.gene.region.chr.df)


#plot of xlink events per gene------------CHR
#plot of xlink events per gene
xlink.genes.chr<-ggplot(top.xlink.events.chr.df, aes(x = reorder_within(gene_name,n,sample), y = n, fill=gene_name)) + 
  geom_bar(stat = "identity") +
  scale_x_reordered() +
  xlab("genes") +
  facet_wrap(~sample + model,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution") +
  ylab("events counts") + 
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(xlink.genes.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.genes.chr.png", height = 6, width = 6)

#plot of xlink events per gene/region
xlink.genes.region.chr<-ggplot(top.xlink.events.chr.df, aes(x = reorder_within(gene_name,n,sample), y = n, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_x_reordered() +
  xlab("genes") +
  facet_wrap(~sample + model,scales = "free") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution -genomic regions") +
  ylab("events counts") + 
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(xlink.genes.region.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.genes.region.chr.png", height = 8, width = 11)



#plot of xlink events per gene/GC
xlink.genes.GC.chr<-ggplot(top.xlink.events.chr.df, aes(x = reorder_within(gene_name,n,sample), y = n, fill=GC)) + 
  geom_bar(stat = "identity") +
  scale_x_reordered() +
  xlab("genes") +
  scale_fill_gradient(low = "#f7f7f7", high = "#99000d") +
  facet_wrap(~sample + model,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution-GC content") +
  ylab("events counts") + 
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(xlink.genes.GC.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.genes.GC.chr.png", height = 8, width = 11)





#plot of xlink events per gene/length
xlink.genes.length.chr<-ggplot(top.xlink.events.chr.df, aes(x = reorder_within(gene_name,n,sample), y = n, fill=length)) + 
  geom_bar(stat = "identity") +
  scale_x_reordered() +
  xlab("genes") +
  scale_fill_gradient(low = "#f7f7f7", high = "#3B6A8D") +
  facet_wrap(~sample + model,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution -gene length") +
  ylab("events counts") + 
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(xlink.genes.length.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.genes.length.chr.png", height = 8, width = 11)


#========SCATTER PLOTS AND FREQUENCY PLOTS=====================================================================

# library(scattermore)
# 
# 
# 
# 
# 
# 
# 
# 
# xlink.events.chr.tidy<-melt(data = xlink.events.chr.df)
# dim(xlink.events.chr.tidy)
# 
# x<-ggplot(xlink.events.chr.tidy) +
#   geom_scattermore(aes(x = gene_name, y = value, color = sample))
# ggsave(x, filename = "/Users/manferg/clip_metanalysis/r-plots/x.png", height = 6, width = 6)
# 
# 
# +
#   scale_y_log10() +
#   facet_wrap(~sample,scales = "free") +
#   xlab("Genes") +
#   ylab("Normalized Counts") +
#   ggtitle("scatterplot") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(plot.title=element_text(hjust=0.5))
# 
# 
# 
# total.xlink.counts<- ggplot(xlink.events.chr.tidy) + 
#   geom_bar(aes(x= sample,y=value,fill=region), stat ='identity') +
#   theme_classic() +
#   coord_flip() +
#   scale_color_manual(values=mypal.region) +
#   scale_fill_manual(values=alpha(c(mypal.region))) +
#   ggtitle("Total Cross-link counts") 
# ggsave(total.xlink.counts, filename = "/Users/manferg/clip_metanalysis/r-plots/total.xlink.counts.png", height = 6, width = 6)
# 


#===gene length normalisation=====#---------------

xlink.events.gene.region.chr.norm.li =list()
top.xlink.events.gene.region.chr.norm.li=list()
xlink.chr.normalised.li=list()
for (i in 1:length(main.chr.li)){
  xlink.events.gene.region.chr.norm.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,region,GC,length) %>% summarize(n=n()) %>% arrange( .,desc(n)) %>% as.data.frame() #counts (or number of xlinks) for each gene.region. 
  xlink.chr.normalised.li[[i]]<-xlink.events.gene.region.chr.norm.li[[i]] %>% mutate( .,xlink.enrichment=n/length) %>% arrange( .,desc(xlink.enrichment)) %>% as.data.frame()
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  xlink.chr.normalised.li[[i]]$sample <- sample
  top.xlink.events.gene.region.chr.norm.li[[i]]<-xlink.chr.normalised.li[[i]][1:10,]  #exclude "None" and list top 20 gene.region.s
  top.xlink.events.gene.region.chr.norm.li[[i]]$gene_name <- factor(top.xlink.events.gene.region.chr.norm.li[[i]]$gene_name , levels = top.xlink.events.gene.region.chr.norm.li[[i]]$gene_name[order(top.xlink.events.gene.region.chr.norm.li[[i]]$xlink.enrichment)])
  
}


xlink.chr.normalised.df<-do.call(rbind,xlink.chr.normalised.li) #convert back to df to plot
top.xlink.events.chr.norm.df<-do.call(rbind,top.xlink.events.gene.region.chr.norm.li)
top.xlink.events.chr.norm.df<-left_join(top.xlink.events.chr.norm.df,metadata,by=c("sample" = "meta_id"))
#reorder levels as samples order
top.xlink.events.chr.norm.df$sample <- factor(top.xlink.events.chr.norm.df$sample , levels=unique(top.xlink.events.chr.norm.df$sample ))


#test
#p<-intresected.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% filter(region == "intron")
#dim(p)

#view(top.xlink.events.gene.region.chr.df)
#view(xlink.events.gene.region.chr.df)


#plot of xlink events per gene/region----LENGTH NORM-----CHR
xlink.genes.region.norm.chr<-ggplot(top.xlink.events.chr.norm.df,aes(reorder_within(gene_name,xlink.enrichment,sample), xlink.enrichment, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_x_reordered() +
  xlab("genes") +
  facet_wrap(~sample + model,scales = "free") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution normalised per gene length -genomic regions") +
  ylab("events counts") + 
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(xlink.genes.region.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.genes.region.norm.chr.png",height = 8, width = 11)


#plot of xlink events per gene/GC----LENGTH NORM-----CHR

xlink.genes.GC.norm.chr<-ggplot(top.xlink.events.chr.norm.df,aes(reorder_within(gene_name,xlink.enrichment,sample), xlink.enrichment, fill=GC)) + 
  facet_wrap(~sample + model,scales = "free") +
  #scale_color_manual(values=mypal.region) +
  #scale_fill_manual(values=alpha(c(mypal.region))) +
  scale_fill_gradient(low = "#f7f7f7", high = "#99000d") +
  scale_x_reordered() +
  geom_bar(stat = "identity") +
  xlab("genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution normalised per gene length-GCcontent") +
  ylab("events counts") + 
  theme_classic() +
  coord_flip()
ggsave(xlink.genes.GC.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.genes.GC.norm.chr.png", height = 8, width = 11)


xlink.genes.length.norm.chr<-ggplot(top.xlink.events.chr.norm.df,aes(reorder_within(gene_name,xlink.enrichment,sample), xlink.enrichment, fill=length)) + 
  facet_wrap(~sample + model,scales = "free") +
  #scale_color_manual(values=mypal.region) +
  #scale_fill_manual(values=alpha(c(mypal.region))) +
  scale_fill_gradient(low = "#f7f7f7", high = "#3B6A8D") +
  scale_x_reordered() +
  geom_bar(stat = "identity") +
  xlab("genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution normalised per gene length-GCcontent") +
  ylab("events counts") + 
  theme_classic() +
  coord_flip()
ggsave(xlink.genes.length.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.genes.length.norm.chr.png", height = 8, width = 11)



# normalize raw x-links per Library size and gene length

head(df.dedup.out.lib)
head(xlink.chr.normalised.df)


xlink.enrich.unique.gene<-xlink.chr.normalised.df %>% group_by(gene_name,sample) %>% summarize(n=sum(xlink.enrichment)) %>% arrange( .,desc(n)) %>% as.data.frame() #total xlinks per gene (sum of all region)
dim(xlink.enrich.unique.gene)

#p<-xlink.events.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% summarize(n=sum(n))
#p<-xlink.events.unique.gene %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1")




# normalize raw x-links per Library size 

head(df.dedup.out.lib)
head(xlink.enrich.unique.gene)
xlink.enrich.chr.lb.norm.df<-left_join(xlink.enrich.unique.gene, df.dedup.out.lib, by = "sample")
xlink.enrich.chr.lb.norm.df<-xlink.enrich.chr.lb.norm.df %>% mutate(ls.factor=count/1e+06) %>% mutate(enrich.norm = n/ls.factor)

tidy.enrich.norm<-xlink.enrich.chr.lb.norm.df %>% dplyr::select(sample, gene_name,enrich.norm) %>% arrange( .,desc(enrich.norm)) %>% filter( .,gene_name != "None")

top_genes.enrich.norm<-tidy.enrich.norm %>% dplyr::select(gene_name) %>% unique()

spread.enrich.norm<-spread(tidy.enrich.norm,sample,enrich.norm)


dim(top_genes.enrich.norm) #35247
#spread.df.na.rm<-spread.df %>% na.omit()
#dim(spread.df.na.rm)

#HEATH-MAP---------
spread.enrich.norm[is.na(spread.enrich.norm)] = 0 #transform NA values into 0 



mat<-as.matrix(spread.enrich.norm[,2:17])#pheatmap only takes the numeric matrix object as input. So, we need to transfer the numeric part of the data frame to a matrix by removing the first 5 columns of categorical data.
rownames(mat)<-as.character(spread.enrich.norm$gene_name)
mat_scale = scale(mat)

dim(mat_scale)

library(pheatmap)
mn_select <- top_genes.enrich.norm[1:30,]
mat <- mat[,mn_reorder_idx]
mat <-mat_scale[mn_select,] 
mn_reorder_idx <- match(reorder_sample_idx, colnames(mat_scale), nomatch = 0) #REORDER SAMPLES ORDER 


#heathmap annotation-------
ann_df = data.frame(metadata$meta_id, metadata$model, row.names = TRUE)
colnames(ann_df)<-c("model")
rownames(ann_df) <- colnames(mat_scale) # name matching

anno_colors = list(mypal.model)




#heath_map
#xl_matrix_col <- pheatmap(mat, annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, main = paste("Top 20 xlinked genes"),fontsize = 10)
xl_matrix_row_enrich_whole <- pheatmap(mat,scale = "row", annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE, main = paste("whole dataset normalised"),fontsize = 10)
ggsave(xl_matrix_row_enrich_whole, filename = "/Users/manferg/clip_metanalysis/r-plots/xl_matrix_row_enrich_whole.png", height = 10, width = 10)

xl_matrix_row_enrich_30 <- pheatmap(mat,scale = "row", annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, main = paste("top 30 xlinkd genes - normalised"),fontsize = 10)
ggsave(xl_matrix_row_enrich_30, filename = "/Users/manferg/clip_metanalysis/r-plots/xl_matrix_row_enrich_30.png", height = 10, width = 10)
























#===========================SCORE PER GENE=================#


score.events.gene.region.chr.li =list()
top.score.events.gene.region.chr.li=list()
for (i in 1:length(main.chr.li)){
  score.events.gene.region.chr.li[[i]]<-main.chr.li[[i]] %>% group_by(score,region) %>% summarize(score=sum(score)) %>% mutate(perc = score * 100/ sum(score)) %>% arrange( .,desc(score)) %>% as.data.frame()#counts (or number of scores) for each gene.region. 
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  score.events.gene.region.chr.li[[i]]$sample <- sample
  
}
score.events.chr.df<-do.call(rbind,score.events.gene.region.chr.li )%>% as.data.frame()#convert back to df to plot #convert back to df to plot
score.events.chr.df<-left_join(score.events.chr.df,metadata,by=c("sample" = "meta_id"))


#reorder levels as samples order
score.events.chr.df$sample <- factor(score.events.chr.df$sample , levels=unique(score.events.chr.df$sample ))


#score per regions plot--------CHR---
score.region.chr<-ggplot(score.events.chr.df, aes(x = reorder(region, -score), y = score, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  facet_wrap(~sample + model,scales = "free") + #each facet has a different scale
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Score distribution") +
  ylab("score events counts") + 
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(score.region.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/score.region.chr.png", height = 8, width =9)

score.region.chr.perc<-ggplot(score.events.chr.df, aes(x = reorder(region, -perc), y = perc, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  facet_wrap(~sample + model,scales = "free") + #each facet has a different scale
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links score distribution -% ") +
  ylab("score events %") + 
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


ggsave(score.region.chr.perc, filename = "/Users/manferg/clip_metanalysis/r-plots/score.region.chr.perc.png",  height = 8, width =9)

score.region.chr.perc.model<-ggplot(score.events.chr.df,aes(x = region, y = perc, fill=model)) + 
  geom_bar(stat = "identity",position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values=mypal.model) +
  scale_fill_manual(values=alpha(c(mypal.model))) +
  ggtitle("Cross-links score distribution -% ") +
  ylab("score events %") +
theme_classic()


ggsave(score.region.chr.perc.model, filename = "/Users/manferg/clip_metanalysis/r-plots/score.region.chr.perc.model.png", height = 6, width = 6)


#score distribution across gene regions-------------------


score.gene.total.chr.li =list()

for (i in 1:length(main.chr.li)){
  score.gene.total.chr.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name) %>% summarize(score=sum(score)) %>% arrange( .,desc(score)) %>% as.data.frame()#counts (or number of scores) for each gene.region. 
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  score.gene.total.chr.li[[i]]$sample <- sample
  
}
score.gene.total.df<-do.call(rbind,score.gene.total.chr.li )%>% as.data.frame()#convert back to df to plot #convert back to df to plot
score.gene.total.df<-left_join(score.gene.total.df,metadata,by=c("sample" = "meta_id"))
#reorder levels as samples order
score.gene.total.df$sample <- factor(score.gene.total.df$sample , levels=unique(score.gene.total.df$sample ))


score.gene.region.chr.li =list()
for (i in 1:length(main.chr.li)){
  score.gene.region.chr.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,region) %>% summarize(score=sum(score)) %>% arrange( .,desc(score)) %>% as.data.frame()
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  score.gene.region.chr.li[[i]]$sample <- sample
  
}

score.gene.region.chr.df<-do.call(rbind,score.gene.region.chr.li )%>% as.data.frame()
score.gene.region.chr.df<-left_join(score.gene.region.chr.df,metadata,by=c("sample" = "meta_id"))
#reorder levels as samples order
score.gene.region.chr.df$sample <- factor(score.gene.region.chr.df$sample , levels=unique(score.gene.region.chr.df$sample ))




score_sumary.li =list()
for (i in 1:length(score.gene.total.chr.li)){
  
  score_sumary.li[[i]]<-left_join(score.gene.total.chr.li[[i]],score.gene.region.chr.li[[i]], by="gene_name") 
  
}

score_summary.df<-do.call(rbind,score_sumary.li )%>% dplyr::select(-sample.y) %>% as.data.frame() 
colnames(score_summary.df)<-c("gene_name","score_total","sample","region","score_region")



score_summary.tidy<-reshape2::melt(score_summary.df,variable.name = "score_name",value.name="score")

top.genes.score<-score.gene.total.df[2:10,]
score_summary.tidy<-score_summary.tidy[score_summary.tidy$gene_name %in% top.genes.score$gene_name,]



#plot of xlink events per gene/score/region #yes!
score.regions.per.gene<-ggplot(score_summary.tidy, aes(x = reorder(gene_name, score), y = score, fill=region)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  facet_wrap(~sample,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links score distribution") +
  ylab("score counts") + 
  theme_classic() +
  coord_flip()

ggsave(score.regions.per.gene, filename = "/Users/manferg/clip_metanalysis/r-plots/score.regions.per.gene.png", height = 10, width = 10)



#===gene length normalisation=====#-----SCORE----------


score.events.gene.region.chr.norm.li =list()
top.score.events.gene.region.chr.norm.li=list()
score.chr.normalised.li=list()
for (i in 1:length(main.chr.li)){
  score.events.gene.region.chr.norm.li[[i]]<-main.chr.li[[i]] %>% group_by(gene_name,region,GC,length) %>% summarize(score=sum(score)) %>% arrange( .,desc(score)) %>% as.data.frame() #counts (or number of scores) for each gene.region. 
  score.chr.normalised.li[[i]]<-score.events.gene.region.chr.norm.li[[i]] %>% mutate( .,score.enrichment=score/length) %>% arrange( .,desc(score.enrichment)) %>% as.data.frame()
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  score.chr.normalised.li[[i]]$sample <- sample
  top.score.events.gene.region.chr.norm.li[[i]]<-score.chr.normalised.li[[i]][1:10,]  #exclude "None" and list top 20 gene.region.s
  top.score.events.gene.region.chr.norm.li[[i]]$gene_name <- factor(top.score.events.gene.region.chr.norm.li[[i]]$gene_name , levels = top.score.events.gene.region.chr.norm.li[[i]]$gene_name[order(top.score.events.gene.region.chr.norm.li[[i]]$score.enrichment)])
  
}


score.chr.normalised.df<-do.call(rbind,score.chr.normalised.li) #convert back to df to plot
top.score.events.chr.norm.df<-do.call(rbind,top.score.events.gene.region.chr.norm.li)
top.score.events.chr.norm.df<-left_join(top.score.events.chr.norm.df,metadata,by=c("sample" = "meta_id"))
#reorder levels as samples order
top.score.events.chr.norm.df$sample <- factor(top.score.events.chr.norm.df$sample , levels=unique(top.score.events.chr.norm.df$sample ))



#test
#p<-intresected.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% filter(region == "intron")
#dim(p)

#view(top.score.events.gene.region.chr.df)
#view(score.events.gene.region.chr.df)

#plot of score events per gene/region-----------CHR

#plot of score events per gene/region----LENGTH NORM-----CHR
score.genes.region.norm.chr<-ggplot(top.score.events.chr.norm.df,aes(reorder_within(gene_name,score.enrichment,sample), score.enrichment, fill=region)) + 
  geom_bar(stat = "identity") +
  scale_x_reordered() +
  xlab("genes") +
  facet_wrap(~sample + model,scales = "free") +
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution normalised per gene length -genomic regions") +
  ylab("events counts") + 
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(score.genes.region.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/score.genes.region.norm.chr.png", height = 8, width =11)


#plot of score events per gene/GC----LENGTH NORM-----CHR

score.genes.GC.norm.chr<-ggplot(top.score.events.chr.norm.df,aes(reorder_within(gene_name,score.enrichment,sample), score.enrichment, fill=GC)) + 
  facet_wrap(~sample + model,scales = "free") +
  #scale_color_manual(values=mypal.region) +
  #scale_fill_manual(values=alpha(c(mypal.region))) +
  scale_fill_gradient(low = "#f7f7f7", high = "#99000d") +
  scale_x_reordered() +
  geom_bar(stat = "identity") +
  xlab("genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution normalised per gene length-GCcontent") +
  ylab("events counts") + 
  theme_classic() +
  coord_flip()
ggsave(score.genes.GC.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/score.genes.GC.norm.chr.png",height = 8, width =11)


score.genes.length.norm.chr<-ggplot(top.score.events.chr.norm.df,aes(reorder_within(gene_name,score.enrichment,sample), score.enrichment, fill=length)) + 
  facet_wrap(~sample + model,scales = "free") +
  #scale_color_manual(values=mypal.region) +
  #scale_fill_manual(values=alpha(c(mypal.region))) +
  scale_fill_gradient(low = "#f7f7f7", high = "#3B6A8D") +
  scale_x_reordered() +
  geom_bar(stat = "identity") +
  xlab("genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links events distribution normalised per gene length-GCcontent") +
  ylab("events counts") + 
  theme_classic() +
  coord_flip()
ggsave(score.genes.length.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/score.genes.length.norm.chr.png", height = 8, width =11)



#HEATH-MAP-----score length normalised----
score.enrich.unique.gene<-score.chr.normalised.df %>% group_by(gene_name,sample) %>% summarize(n=sum(score.enrichment)) %>% arrange( .,desc(n)) %>% as.data.frame() #total xlinks per gene (sum of all region)


score.enrich.chr.lb.norm.df<-left_join(score.enrich.unique.gene, df.dedup.out.lib, by = "sample")
score.enrich.chr.lb.norm.df<-score.enrich.chr.lb.norm.df %>% mutate(ls.factor=count/1e+06) %>% mutate(enrich.norm = n/ls.factor)

tidy.enrich.norm<-score.enrich.chr.lb.norm.df %>% dplyr::select(sample, gene_name,enrich.norm) %>% arrange( .,desc(enrich.norm)) %>% filter( .,gene_name != "None")

top_genes.enrich.norm<-tidy.enrich.norm %>% dplyr::select(gene_name) %>% unique()

spread.enrich.norm<-spread(tidy.enrich.norm,sample,enrich.norm)

dim(top_genes.enrich.norm)

#spread.df.na.rm<-spread.df %>% na.omit()
#dim(spread.df.na.rm)


spread.enrich.norm[is.na(spread.enrich.norm)] = 0 #transform NA values into 0 



mat<-as.matrix(spread.enrich.norm[,2:17])#pheatmap only takes the numeric matrix object as input. So, we need to transfer the numeric part of the data frame to a matrix by removing the first 5 columns of categorical data.
rownames(mat)<-as.character(spread.df$gene_name)
mat_scale = scale(mat)

dim(mat_scale)

library(pheatmap)
mn_select <- top_genes.enrich.norm[1:30,]
 mat <- mat[,mn_reorder_idx]
mat <-mat_scale[mn_select,]
mn_reorder_idx <- match(reorder_sample_idx, colnames(mat_scale), nomatch = 0) #REORDER SAMPLES ORDER 

#heathmap annotation-------
ann_df = data.frame(metadata$meta_id, metadata$model, row.names = TRUE)
colnames(ann_df)<-c("model")
rownames(ann_df) <- colnames(mat_scale) # name matching

anno_colors = list(mypal.model)


#heath_map
#xl_matrix_col <- pheatmap(mat, annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, main = paste("Top 20 xlinked genes"),fontsize = 10)
score_matrix_row_enrich_whole <- pheatmap(mat,scale = "row", annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE, main = paste("whole dataset normalised"),fontsize = 10)
ggsave(score_matrix_row_enrich_whole, filename = "/Users/manferg/clip_metanalysis/r-plots/score_matrix_row_enrich_whole.png", height = 10, width = 10)

score_matrix_row_enrich_30 <- pheatmap(mat,scale = "row", annotation = ann_df,annotation_colors = anno_colors,cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=TRUE, main = paste("top 30 xlinkd genes - normalised"),fontsize = 10)
ggsave(score_matrix_row_enrich_30, filename = "/Users/manferg/clip_metanalysis/r-plots/score_matrix_row_enrich_30.png", height = 10, width = 10)



#=====NORM ON GENE LENGHT AND LIBRARY SIZE========#



# #test
# #p<-intresected.chr.df %>% filter(sample == "grot_293fl_1") %>% filter(gene_name == "GSE1") %>% filter(region == "intron")
# #dim(p)
# 
# #view(top.score.events.gene.region.chr.df)
# #view(score.events.gene.region.chr.df)
# 
# #plot of score events per gene/region-----------CHR
# 
# #plot of score events per gene/region----LENGTH NORM-----CHR
# score.genes.region.norm.chr<-ggplot(top.score.events.chr.norm.df,aes(reorder_within(gene_name,score.enrichment,sample), score.enrichment, fill=region)) + 
#   geom_bar(stat = "identity") +
#   scale_x_reordered() +
#   xlab("genes") +
#   facet_wrap(~sample + model,scales = "free") +
#   scale_color_manual(values=mypal.region) +
#   scale_fill_manual(values=alpha(c(mypal.region))) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ggtitle("Cross-links events distribution normalised per gene length -genomic regions") +
#   ylab("events counts") + 
#   coord_flip() +
#   theme_classic() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# ggsave(score.genes.region.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/score.genes.region.norm.chr.png", height = 6, width = 6)
# 
# 
# #plot of score events per gene/GC----LENGTH NORM-----CHR
# 
# score.genes.GC.norm.chr<-ggplot(top.score.events.chr.norm.df,aes(reorder_within(gene_name,score.enrichment,sample), score.enrichment, fill=GC)) + 
#   facet_wrap(~sample + model,scales = "free") +
#   #scale_color_manual(values=mypal.region) +
#   #scale_fill_manual(values=alpha(c(mypal.region))) +
#   scale_fill_gradient(low = "#f7f7f7", high = "#99000d") +
#   scale_x_reordered() +
#   geom_bar(stat = "identity") +
#   xlab("genes") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ggtitle("Cross-links events distribution normalised per gene length-GCcontent") +
#   ylab("events counts") + 
#   theme_classic() +
#   coord_flip()
# ggsave(score.genes.GC.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/score.genes.GC.norm.chr.png", height = 6, width = 6)
# 
# 
# score.genes.length.norm.chr<-ggplot(top.score.events.chr.norm.df,aes(reorder_within(gene_name,score.enrichment,sample), score.enrichment, fill=length)) + 
#   facet_wrap(~sample + model,scales = "free") +
#   #scale_color_manual(values=mypal.region) +
#   #scale_fill_manual(values=alpha(c(mypal.region))) +
#   scale_fill_gradient(low = "#f7f7f7", high = "#3B6A8D") +
#   scale_x_reordered() +
#   geom_bar(stat = "identity") +
#   xlab("genes") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ggtitle("Cross-links events distribution normalised per gene length-GCcontent") +
#   ylab("events counts") + 
#   theme_classic() +
#   coord_flip()
# ggsave(score.genes.length.norm.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/score.genes.length.norm.chr.png", height = 6, width = 6)
# 
# 
# 










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
main.xlink.genome.df<-left_join(main.xlink.genome.df,metadata,by=c("sample" = "meta_id"))
#reorder levels as samples order
main.xlink.genome.df$sample <- factor(main.xlink.genome.df$sample , levels=unique(main.xlink.genome.df$sample ))



#CHR - counts (or number of xlinks) for each chromosome------------------
chr.xlink.genome.li =list()
for (i in 1:length(main.chr.li)){
  chr.xlink.genome.li[[i]]<-main.chr.li[[i]] %>% group_by(seqname,region) %>% summarize(n=n()) %>% arrange( .,desc(seqname))
  sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
  chr.xlink.genome.li[[i]]$sample <- sample
}

chr.xlink.genome.df<-do.call(rbind,chr.xlink.genome.li) #convert back to df to plot
chr.xlink.genome.df<-left_join(chr.xlink.genome.df,metadata,by=c("sample" = "meta_id"))
#reorder levels as samples order
chr.xlink.genome.df$sample <- factor(chr.xlink.genome.df$sample , levels=unique(chr.xlink.genome.df$sample ))



# dupe count list for chromosomes--------------
# chr.xlink.duped.genome.li =list()
# for (i in 1:length(main.chr.li)){
#   chr.xlink.duped.genome.li[[i]]<-main.chr.li[[i]] %>% get_dupes(seqname) #dupe counts (or number of xlinks) for each gene
#   sample<-as.character(names(main.chr.li[i])) # Create a new vector with sample names
#   chr.xlink.duped.genome.li[[i]]$sample <- sample
# }

#chr.xlink.duped.genome.df<-do.call(chr.xlink.duped.genome.li)





#xlinks counts at genome level plots--------------------------------------
xlink.events.xlink.genome.main<-ggplot(main.xlink.genome.df, aes(x = reorder(seqname,n), y = n, fill=region)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample + model,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links xlink distribution across the genome") +
  ylab("events counts") + 
  theme_classic()+
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  coord_flip()
ggsave(xlink.events.xlink.genome.main, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.events.xlink.genome.main.png", height = 11, width = 11)



xlink.events.xlink.genome.chr<-ggplot(chr.xlink.genome.df, aes(x = reorder(seqname,n), y = n, fill=region)) + 
  geom_bar(stat = "identity") +
  xlab("genes") +
  facet_wrap(~sample + model,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cross-links xlink distribution across the genome") +
  ylab("events counts") + 
  theme_classic()+
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  coord_flip()
ggsave(xlink.events.xlink.genome.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.events.xlink.genome.chr.png", height = 14, width = 14)

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
  chr.score.genome.li[[i]]<-main.chr.li[[i]] %>% group_by(score,seqname,region) %>% summarize(n=sum(score)) %>% arrange( .,desc(score))
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
  coord_flip()+
  scale_color_manual(values=mypal.region) +
  scale_fill_manual(values=alpha(c(mypal.region))) +
  theme_classic()
ggsave(xlink.events.score.genome.chr, filename = "/Users/manferg/clip_metanalysis/r-plots/xlink.events.score.genome.chr.png", height = 14, width = 14)












### Annotate our heatmap (optional)
annotation <- data.frame(sampletype=xlink.events.chr.df[sample,'sampletype'], 
                         row.names=rownames(xlink.events.chr.df))

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap
library(pheatmap)
p<-pheatmap(xlink.events.chr.df, color = heat_colors, cluster_rows = T, show_rownames=F,, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)
















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
