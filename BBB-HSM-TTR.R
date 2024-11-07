#加载R包
library(Seurat)
library(glmGamPoi)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#设置路径
setwd("xx")
folders<-list.files(".")
folders

###################################################
#创建Seurat对象
library(Seurat)
scelist<-lapply(folders, function(folder){
  CreateSeuratObject(counts = Read10X(folder),
                     project = folder)
})
scelist
#汇总样本得到6个对象的集合
Tdata<-merge(scelist[[1]],
             y=c(scelist[[2]],scelist[[3]],
                 scelist[[4]],scelist[[5]],scelist[[6]]),
             add.cell.ids = folders,
             project = "BBB")

#统计各样本细胞数
table(Tdata$orig.ident)

#稀疏矩阵#GetAssayData(Tdata)[1:20,1:25] 
Tdata<-JoinLayers(Tdata)
GetAssayData(object = Tdata, assay = "RNA", slot = "counts")[1:20,1:25]
save(Tdata,file = "Tdata.Rdata")
load("Tdata.Rdata")

min(Tdata$nCount_RNA)
min(Tdata$nFeature_RNA)

#########################################################
#过滤数据
Tdata$ratio<-Tdata$nFeature_RNA/Tdata$nCount_RNA
filtered_Tdata <- subset(x = Tdata, 
                         subset= (ratio>0.2) &
                           nFeature_RNA >= 250)
save(filtered_Tdata,file = "filtered_Tdata.Rdata")
load("filtered_Tdata.Rdata")
#统计过滤后各样本细胞数
min(filtered_Tdata$nCount_RNA)
min(filtered_Tdata$nFeature_RNA)
table(filtered_Tdata$orig.ident)
#绘图nFeature_RNA，nCount_RNA
VlnPlot(filtered_Tdata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot2 <- FeatureScatter(filtered_Tdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

#########################################################
#十分耗时十分占内存,放后台跑还可以
filtered_Tdata<-SCTransform(filtered_Tdata,method="glmGamPoi",verbose=FALSE)
save(filtered_Tdata,file = "filtered_Tdata_SCTransform.Rdata")
load("filtered_Tdata_SCTransform.Rdata")

#########################################################
#聚类分群
filtered_Tdata <- NormalizeData(filtered_Tdata)
filtered_Tdata <- ScaleData(filtered_Tdata)
filtered_Tdata<-RunPCA(filtered_Tdata,verbose = T)
filtered_Tdata<-RunUMAP(filtered_Tdata,dims = 1:30,verbose = T)
filtered_Tdata<-FindNeighbors(filtered_Tdata,dims = 1:30,verbose = T)
filtered_Tdata<-FindClusters(filtered_Tdata,verbose = T,resolution = 0.2)
DimPlot(filtered_Tdata,label = TRUE)
filtered_Tdata<-RunTSNE(filtered_Tdata,dims = 1:30)
DimPlot(filtered_Tdata,reduction = "tsne",label = TRUE)

#########################################################
#验证批次效应
table(filtered_Tdata@meta.data$orig.ident)
#定义condition
filtered_Tdata@meta.data$condition<-c(rep("T0",24887),
                                      rep("T1",25246),
                                      rep("T3",23237))
library(RColorBrewer)
DimPlot(filtered_Tdata,reduction = "umap",label = FALSE,group.by = "condition",
        cols = c("#E41A1C","#FFFF33","#3182BD"),
        pt.size = 0.00001)
DimPlot(filtered_Tdata,reduction = "tsne",label = FALSE,group.by = "condition",
        cols = c("#E41A1C","#FFFF33","#3182BD"),
        pt.size = 0.00001)

DimPlot(filtered_Tdata,reduction = "umap",label = T,split.by = "condition")

#########################################################

brewer.pal(3,"Set1")

#########################################################
#展示marker gene
FeaturePlot(filtered_Tdata, features = c("Ttr"),
            reduction="umap",
            cols = c("grey", "red"))
VlnPlot(filtered_Tdata, "Ttr")
FeaturePlot(filtered_Tdata, features = c("Cx3cr1","Aqp4","Cldn5","Vtn","Acta2","Enpp2"),
            reduction="umap",
            cols = c("grey", "red"))

#####################################################################
FeaturePlot(filtered_Tdata, features = c("Selplg","Tmem119","Ctss","Cx3cr1","P2ry12","Hexb"),
            reduction="umap",
            cols = c("grey", "red"))
FeaturePlot(filtered_Tdata, features = c("Dclk1","Aqp4","Gfap"),
            reduction="umap",
            cols = c("grey", "red"))
FeaturePlot(filtered_Tdata, features = c("Cldn5","Ly6a","Pglyrp1","Pecam1"),
            reduction="umap",
            cols = c("grey", "red"))
FeaturePlot(filtered_Tdata, features = c("Acta2","Sncg"),
            reduction="umap",
            cols = c("grey", "red"))
FeaturePlot(filtered_Tdata, features = c("Enpp2","Folr1"),
            reduction="umap",
            cols = c("grey", "red"))
FeaturePlot(filtered_Tdata, features = c("Pdgfra"),
            reduction="umap",
            cols = c("grey", "red"))
FeaturePlot(filtered_Tdata, features = c("Tubb3","Rbfox1","Meg3"),
            reduction="umap",
            cols = c("grey", "red"))

markergene<-c("Selplg","Tmem119","Ctss","Cx3cr1","P2ry12","Hexb",
              "Dclk1","Aqp4","Gfap",
              "Cldn5","Ly6a","Pglyrp1","Pecam1",
              "Vtn","Pdgfrb","Cspg4",
              "Acta2","Sncg",
              "Enpp2","Folr1",
              "Pdgfra",
              "Tubb3","Rbfox1","Meg3")

#dotplot图可视化
DotPlot(filtered_Tdata, features = markergene)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

#microglia   0,9,10
#astrocyte 1,3,7,17
#endothelium  2,8,12,14,11
#pericyte  4
#SMC  5
#choroid plexus   6,13
#OPC 16
#neuron 15,18

new.cluster.ids <- c("microglia",# 0
                     "astrocyte",# 1
                     "endothelium",# 2
                     "astrocyte",# 3
                     "pericyte",# 4
                     "SMC",# 5
                     "choroid plexus",# 6
                     "astrocyte",# 7
                     "endothelium",# 8
                     "microglia",# 9
                     "microglia",# 10
                     "endothelium",# 11
                     "endothelium",# 12
                     "choroid plexus",# 13
                     "endothelium",# 14
                     "neuron",# 15
                     "OPC",# 16
                     "astrocyte",# 17
                     "neuron")# 18

names(new.cluster.ids) <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
filtered_Tdata <- RenameIdents(filtered_Tdata, new.cluster.ids)

DimPlot(filtered_Tdata, reduction = "umap", label = TRUE,
        cols = brewer.pal(9,"Set1"))

saveRDS(filtered_Tdata,"filtered_Tdata_anno.rds")

#读取细胞注释后数据
filtered_Tdata <- readRDS("filtered_Tdata_anno.rds")

#################################################################
#统计各分组细胞数
table(filtered_Tdata$condition)

#统计各细胞类型占比
prop.table(table((Idents(filtered_Tdata))))

#统计各细胞类型数（condition）
table(Idents(filtered_Tdata),filtered_Tdata$condition)

#统计各细胞类型占比（condition）
Cellratio<-prop.table(table(Idents(filtered_Tdata),filtered_Tdata$condition),margin = 2)
Cellratio
Cellratio<-as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  scale_fill_brewer(palette = "Set1")+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
#细胞分亚群前过滤
meta<- filtered_Tdata@meta.data
a<- Idents(filtered_Tdata)
filtered_Tdata$cell_type <- a
filtered_Tdata[["percent.mt"]] <- PercentageFeatureSet(filtered_Tdata, pattern = "^mt-")
meta<- filtered_Tdata@meta.data

######################################################################
#小胶细胞分亚群
microglia<-subset(filtered_Tdata,cell_type=="microglia"&percent.mt<20)
microglia<-SCTransform(microglia,method="glmGamPoi",verbose=FALSE)
microglia<-RunPCA(microglia,verbose = T)
microglia<-RunUMAP(microglia,dims = 1:50,verbose = T)
microglia<-FindNeighbors(microglia,dims = 1:50,verbose = T)
microglia<-FindClusters(microglia,verbose = T,resolution = 0.2)
DimPlot(microglia,reduction = "umap",label = T)
saveRDS(microglia,"microglia.rds")
microglia <- readRDS("xx/microglia.rds")
FeaturePlot(microglia, features = c("Ttr"),
            reduction="umap",
            cols = c("grey", "red"))
VlnPlot(microglia, "Ttr")

FeaturePlot(microglia, features = c("Flt1"),
            reduction="umap",
            cols = c("grey", "red"))
VlnPlot(microglia, "Flt1")

DotPlot(microglia, features = "Ttr")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

microglia@meta.data$condition<-c(rep("T0",8617),
                                 rep("T1",8509),
                                 rep("T3",8183))


FeaturePlot(microglia,"Ttr",split.by = "condition")
VlnPlot(microglia,"Ttr",split.by = "condition")

#####################################################################
new.cluster.ids <- c("LSM",# 0
                     "LSM",# 1
                     "LSM",# 2
                     "HSM",# 3
                     "LSM",# 4
                     "LSM",# 5
                     "LSM",# 6
                     "LSM",# 7
                     "LSM",# 8
                     "LSM",# 9
                     "LSM",# 10
                     "LSM",# 11
                     "LSM",# 12
                     "LSM",# 13
                     "LSM",# 14
                     "LSM")# 15

names(new.cluster.ids) <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
microglia <- RenameIdents(microglia, new.cluster.ids)

DimPlot(microglia, reduction = "umap", label = TRUE,
        cols = brewer.pal(9,"Set1"))
FeaturePlot(microglia,"Ttr",split.by = "condition",
            cols = c("grey", "red"))
VlnPlot(microglia,"Ttr",split.by = "condition")
VlnPlot(microglia,"Ttr")

saveRDS(microglia,"microglia_anno.rds")
#########################################################################

#读取小胶细胞注释亚群后数据
microglia <- readRDS("xx/microglia_anno.rds")

table(Idents(microglia))

prop.table(table((Idents(microglia))))

table(Idents(microglia),microglia$condition)

Cellratio<-prop.table(table(Idents(microglia),microglia$condition),margin = 2)
Cellratio
Cellratio<-as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  scale_fill_brewer(palette = "Set1")+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

microglia.markers<-FindAllMarkers(object = microglia,only.pos = TRUE,min.pct = 0.25,thresh.use=0.25,
                                  logfc.threshold = 0.25)
print(microglia.markers)
write.csv(microglia.markers,file = "microglia_markers.csv")
microglia.markers <- read.csv("microglia_markers.csv", header = TRUE)
top10<-microglia.markers%>%group_by(cluster)%>%top_n(10,avg_log2FC)

DoHeatmap(object = microglia, features = top10$gene,
          group.colors = c("red", "#336699"))+ #设置组别颜色
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))#设置热图颜色

#################################################################
## barcode_ttr ##
meta <- microglia@meta.data
head(meta)
table(microglia$seurat_clusters)
DimPlot(microglia, reduction = "umap", label = TRUE) + NoLegend()
VlnPlot(microglia, features = c("Ttr"), group.by = "seurat_clusters")
meta <- meta %>% mutate(ttr_expr=case_when(seurat_clusters%in%c(0,1,2,4,5,6,7,8,9,10,11,12,13,14,15)~"low",
                                           TRUE~"high"))
barcode_HSM <- row.names(meta %>% filter(ttr_expr=="high"))
barcode_LSM <- row.names(meta %>% filter(ttr_expr=="low"))

meta<- filtered_Tdata@meta.data
celltype<- Idents(filtered_Tdata)
endo_barcode <-celltype[which(celltype=="endothelium")] 
meta<- meta %>% mutate(expr_com=case_when(rownames(meta)%in%barcode_HSM~"HSM",
                                          rownames(meta)%in%barcode_LSM~"LSM",
                                          rownames(meta)%in%names(endo_barcode)~"endothelium",
                                          TRUE~"no"))
filtered_Tdata$expr_com <- meta$expr_com
DimPlot(filtered_Tdata)
filter_data<- subset(filtered_Tdata,expr_com%in%c("HSM","LSM","endothelium"))
Idents(filter_data) <- "expr_com"
DimPlot(filter_data)
rm(filtered_Tdata)

### t3时期分析 #####
meta<- filter_data@meta.data
table(meta$orig.ident)
meta<- meta %>% mutate(day=case_when(grepl("T3",orig.ident)~"T3",
                                     TRUE~"OTHER"))
filter_data$day <- meta$day
filter_data <- subset(filter_data,day=="T3")
Idents(filter_data) <- "expr_com"
DimPlot(filter_data)
### 受体配体分析 ####
# 安装CellChat包#
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.19")  # 安装最新版本的 Bioconductor
BiocManager::install(c("Biobase"))
BiocManager::install(c("BiocNeighbors"))
BiocManager::install(c("BiocGenerics"))
BiocManager::install(c("ComplexHeatmap"))
library(Biobase)
library(BiocNeighbors)
library(BiocGenerics)
library(ComplexHeatmap)
devtools::install_github("jinworks/CellChat")
library(CellChat)
seurat_object <- filter_data
data.use <- GetAssayData(seurat_object, layer = "data", assay = "SCT")
labels <- Idents(seurat_object)
meta <- data.frame(labels = labels, row.names = names(labels))

cellChat <- createCellChat(object = seurat_object$SCT$data, group.by = "labels", assay = "SCT")
str(filter_data)
class(seurat_object)
slotNames(seurat_object)
seurat_object@assays

## 准备数据集
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
cellChat@DB <- CellChatDB
cellchat<- cellChat
##预步骤
cellchat <- subsetData(cellChat)
future::plan("multisession", workers = 4)
result <- presto::wilcoxauc(data.use, labels, groups_use = level.use)
print(result)
getAnywhere(identifyOverExpressedGenes)
# 调用 identifyOverExpressedGenes
cellchat <- identifyOverExpressedGenes(cellchat)
devtools::install_github('immunogenomics/presto')
cellchat <- identifyOverExpressedInteractions(cellchat)
### 构建网络 ##
options(future.globals.maxSize=4000000000) 
cellchat <- computeCommunProb(cellchat, type = "triMean")
saveRDS(cellchat,"T3_ttr_endo_process.rds")
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"t3_ttr_endo_communicate.csv")
### 绘图 ###
cellchat <- aggregateNet(cellchat)
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
pdf("t3_ttr_endo_bubble.pdf",height = 15,width = 6)
netVisual_bubble(cellchat, sources.use = c(1:3), targets.use = c(1:3), remove.isolate = FALSE)
dev.off()
pdf("t3_ttr_endo_circle.pdf",height = 15,width = 15)
netVisual_chord_gene(cellchat, sources.use = c(1:3), targets.use = c(1:3))
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat)

######################################################################
#提取小胶细胞，组间差异分析，分群
table(Idents(filtered_Tdata))
microglia<-filtered_Tdata[,Idents(filtered_Tdata)%in%c("microglia")]
table(Idents(microglia))
microglia<-RunPCA(microglia,verbose = T)
microglia<-RunUMAP(microglia,dims = 1:30,verbose = T)
microglia<-FindNeighbors(microglia,dims = 1:30,verbose = T)
microglia<-FindClusters(microglia,verbose = T,resolution = 0.01)
DimPlot(microglia,label = TRUE)#+NoLegend()
DimPlot(microglia, reduction = "umap", label = TRUE,
        cols = brewer.pal(9,"Set1"),
        split.by = "condition") 
microglia<-RunTSNE(microglia,dims = 1:30)
DimPlot(microglia, reduction = "tsne", label = TRUE,
        cols = brewer.pal(9,"Set1"),
        split.by = "condition") 

microglia.markers<-FindAllMarkers(object = microglia,only.pos = TRUE,min.pct = 0.25,thresh.use=0.25,
                                  logfc.threshold = 0.25)



library(dplyr)
#每个亚群挑10个marker基因进行展示
top10<-microglia.markers%>%group_by(cluster)%>%top_n(10,avg_log2FC)
DoHeatmap(object=microglia,features = top10$gene)

######################################################################

table(microglia$condition)

microgliadeg_T0_T1<-FindMarkers(microglia,group.by="condition",
                                ident.1="T0",ident.2="T1",
                                verbose=T,test.use="wilcox",min.pct=0.25)
write.csv(microgliadeg_T0_T1,file = "microgliadeg_T0_T1.csv")
microgliadeg_T1_T3<-FindMarkers(microglia,group.by="condition",
                                ident.1="T1",ident.2="T3",
                                verbose=T,test.use="wilcox",min.pct=0.25)
write.csv(microgliadeg_T1_T3,file = "microgliadeg_T1_T3.csv")