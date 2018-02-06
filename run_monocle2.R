#!/usr/bin/Rscript
##source("http://bioconductor.org/biocLite.R")
##biocLite()
##biocLite("monocle",lib.loc="/data/murphy/shared_libs/R/")

library (Seurat)
library (monocle)
library (plyr)
library (gplots)
library (RColorBrewer)


load ("/data/murphy/home/mplass/planaria/DGEs/tsne_multimapped/reg.Robj")
mycol    <- read.csv ("/data/murphy/home/mplass/planaria/DGEs/tsne_multimapped/paper+reg/clusters_and_colors2.txt",sep="\t",header=F)
mycol2 <- as.vector(mycol[,2])
names(mycol2) <- mycol[,1]

### run monocle for the full dataset
df.full <- reg@raw.data[rownames(reg@data),colnames(reg@data)]

df.pheno.full <- as.data.frame(reg@data.info)
df.pheno.full <- new ("AnnotatedDataFrame",data=df.pheno.full)
cds.full <- newCellDataSet(as.matrix(df.full),phenoData = df.pheno.full)

cds.full <- estimateSizeFactors(cds.full)
cds.full <- estimateDispersions(cds.full)

pData(cds.full)$cell_type <- pData(cds.full)$final_Id

options(repr.plot.width=4, repr.plot.height=3)
plot_pc_variance_explained(cds.full) + geom_vline(xintercept = c(10,20,25,30,35,40,45,50)) 
ggsave ("variance_expolained_PCs.full.pdf",width=5,height=5) 


valid_cds.full2 <- reduceDimension(cds.full,norm_method = 'log', verbose = F, max_components = 50) 
valid_cds.full2 <- orderCells(valid_cds.full2, reverse = T) 

save (cds.full,valid_cds.full2,mycol2, file="monocle2.test1.Robj") 


colsamp <-sample(unname(mycol2[-(1:12)]), 30)
colsamp[15] <- "forestgreen"
colsamp[23] <- "seagreen"

####################################################################################################
## we asisign the color of the most frequent cluster of the do a plot of monocle in which we color de monocle 

monocle_col <- c("#9ECAE1","NA","royalblue","dodgerblue", "hotpink", "grey80","NA","NA",
	       "#2171B5","#6BAED6","NA","indianred 3","#993404",
	       	"grey80","grey80","tomato2","grey80","grey80","grey80",
		"grey80","firebrick","mediumpurple3","#FFF7BC","#FEE391",
		"#FFF7BC","hotpink","mediumorchid1","deeppink","grey80","forestgreen")

names(monocle_col) <-1:30

####################################################################################################
## we plot the monocle output and the seurat output
pData(valid_cds.full2)$cell_type <- reg@data.info$final_Id

p <- plot_complex_cell_trajectory(valid_cds.full2, color_by ='cell_type',
  show_branch_points = T, cell_size = 1.5, cell_link_size =0.5, root_states = 29) +
    theme (legend.position="top", legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = mycol2)

p

ggsave ("monocle_tree_seurat.pdf", width=12,height=10)

p <- plot_complex_cell_trajectory(valid_cds.full2, color_by = 'State',
     show_branch_points = T,
                                   cell_size = 2, cell_link_size = 0.3,root_states =29) +
       theme (legend.position="top", legend.title=element_blank()) +
       guides(colour = guide_legend(override.aes = list(size = 5))) +
       scale_color_manual(values = monocle_col)



p
ggsave ("monocle_tree_monocle.pdf", width=12,height=10)

####################################################################################################
### paint a tSNE wwith monocle 2 states
load ("../DGEs/tsne_multimapped/reg.Robj")
monocle_groups <- (pData(valid_cds.full2))$State

reg@data.info$monocle <-monocle_groups

TSNEPlot(reg, group.by="monocle",do.label=T,colors.use=monocle_col)
ggsave ("tsne_reg_monocle.pdf", width=12,height=8)


####################################################################################################
### we compare the classification of cells in clusters and the Seurat clusters

### the order of out clusters
myident <- c("neoblast 1","neoblast 2","neoblast 3","neoblast 4","neoblast 5","neoblast 6","neoblast 7","neoblast 8","neoblast 9",
	           "neoblast 10","neoblast 11","neoblast 12","neoblast 13","epidermal neoblasts","epidermis DVb neoblast",
   	         "parenchyma progenitors","parenchyma 1","parenchyma 2","glia","parenchyma 3","parenchyma 4","pigment",
             "gut progenitors","phagocytes","goblet cells","gut unknown",
             "muscle progenitors","muscle body","muscle pharynx",
             "neural progenitors","ChAT neurons 1","ChAT neurons 2","spp-11+ neurons","npp-18+ neurons","cav-1+ neurons 1","cav-1+ neurons 2","neural 1","neural 2",
             "secretory 1","secretory 2","secretory 3","secretory 4",
             "early epidermal progenitors","activated early epidermal progenitors","late epidermal progenitors 1","late epidermal progenitors 2","epidermis",
             "epidermis DVb",
             "pharynx progenitors","pharynx",
             "protonephridia") 




mytable <- table(pData(valid_cds.full2)[, c('State', 'cell_type')]) 
mytab_sum <- apply (mytable,1,sum) 
mynewtable <- mytable*100/mytab_sum 

mymat <- matrix (data <- as.numeric(mynewtable),nrow=30,ncol=51,byrow=F)
colnames(mymat) <- dimnames(mynewtable)$cell_type
rownames(mymat) <- dimnames(mynewtable)$State
mymat[is.nan(mymat)] <-0

mymat <- mymat[,myident]

title <- "heatmap_percent_monocle_clusters_overlap.pdf"
pdf (title,width=10,height=10)
heatmap.2 (t(mymat),Rowv = FALSE,Colv=TRUE,scale="none",col=colorRampPalette(c("white","gold","red"),space = "Lab"),
	  sepcolor="black",colsep=0:31,rowsep=0:52,xlab="monocle Sates", ylab= "Seurat Clusters",
	  trace="none",dendrogram="none",margins = c(5, 15), srtCol=90,
	  density.info="none",key=T,keysize=1,key.title=NA, key.ylab=NA,key.ytickfun=NA,key.xlab="% overlap")

dev.off()


mytable <- table(pData(valid_cds.full2)[, c('State', 'cell_type')])
mytab_sum <- apply (mytable,2,sum)
mynewtable <- t(t(mytable)*100/mytab_sum)

mymat <- matrix (data <- as.numeric(mynewtable),nrow=30,ncol=51,byrow=F)
colnames(mymat) <- dimnames(mynewtable)$cell_type
rownames(mymat) <- dimnames(mynewtable)$State
mymat <- mymat[,myident]


title <- "heatmap_percent_seurat_clusters_overlap.pdf"
pdf (title,width=10,height=10)

heatmap.2 (t(mymat),Rowv = FALSE,Colv=TRUE,scale="none",col=colorRampPalette(c("white","gold","red"),space = "Lab"),
	  sepcolor="black",colsep=0:31,rowsep=0:52,xlab="monocle Sates", ylab= "Seurat Clusters",
	  trace="none",dendrogram="none",margins = c(5, 15), srtCol=90,
	  density.info="none",key=T,keysize=1,key.title=NA, key.ylab=NA,key.ytickfun=NA,key.xlab="% overlap")

dev.off()
