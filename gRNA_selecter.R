.libPaths("D:\\Programming\\R_working\\Rpackages")
setwd("D:\\Programming\\R_working\\CENH3")

library(data.table)
#library(ggbio)
library(ggthemes)
library(RColorBrewer)
library(ComplexHeatmap)
library(VennDiagram)
library(ggplot2)
library(ggsci)
library(GenomicRanges)

#file from java -jar FlashFry-assembly-1.10.jar discover
score=fread('variable_gRNAs.scored')
#file after intersection of off-target coordinates and sunflower genes
intersect_bed = fread('intersect_outputoffTarget_vs_HA412_annotaion.bed')

# cenh3 position in sunflower genome
target_coordinate_cenh3_chr = '07|chr|HA412.bronze.20140814'
start = 71337443	
end = 71340378

#remove gRNAs of CENH3 sites
on_target = which(intersect_bed$V1 == target_coordinate_cenh3_chr & intersect_bed$V2 > start & intersect_bed$V3 < end)
off_intersect_bed = intersect_bed[-on_target]

# count number of targets in genes
count = off_intersect_bed[V8 == "gene",.(n_gene_off_target = length(V1)), by = V4]

#merge score and gene off targets
mer_score_off_target = merge(score, by.x = 'target', count, by.y = 'V4')

write.table(mer_score_off_target[,-c(3,4,5,6,11,12,13,15,16)], "gRNA_CENH3_selected.tab", 
            sep = "\t", quote = F, row.names = F, dec = ',')

intersect_bed[V4 == 'CACAGGCGCTACGTGAGATTAGG']

