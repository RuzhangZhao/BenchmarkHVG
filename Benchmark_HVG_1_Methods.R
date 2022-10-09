# Author: Ruzhang Zhao
# Dataset loading can be found  under 
## https://github.com/RuzhangZhao/scRNAseqProcess/blob/main/DatatsetLoading.R
library(parallel)
library(scran)
library(MASS)
library(DuoClustering2018)
library(caret)
library(Seurat)
library(mclust)
library(FNN)
library(caret)
library(Rfast)
library(Matrix)
library(SeuratObject)
library(ggplot2)
library(data.table)
library(HDF5Array)
library(pdist)

# saved as paste0(save_path,dataset_name,"_pcalist_base.rds")
hvg_pca<-function(rna_mat,dataset_name,
                  rna_cell_label=NULL,# for plot
                  nfeatures = 2000,
                  save_path="benchmark/",
                  save_figure=FALSE){
    save_path = paste0(save_path,dataset_name,"/")

    seurat.obj0<-CreateSeuratObject(rna_mat)
    rna_mat_PFlog1pPF<-t(t(seurat.obj0@assays$RNA@counts)/colSums(seurat.obj0@assays$RNA@counts))*mean(colSums(seurat.obj0@assays$RNA@counts))
    rna_mat_PFlog1pPF<-log1p(rna_mat_PFlog1pPF)
    rna_mat_PFlog1pPF<-t(t(rna_mat_PFlog1pPF)/colSums(rna_mat_PFlog1pPF))*mean(colSums(rna_mat_PFlog1pPF))

    seurat.obj.pca<-list()
    var.seurat.obj<-list()

    #################################################
    #1.seurat selection: log(var(count))~log(mean(count))
    message("method1")
    seurat.obj1<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj1$cell_type<-rna_cell_label
    seurat.obj1<-NormalizeData(seurat.obj1,verbose = F)
    seurat.obj1<-FindVariableFeatures(seurat.obj1,nfeatures=nfeatures,verbose = F)
    seurat.obj1<-ScaleData(seurat.obj1,verbose = F)
    seurat.obj1<-RunPCA(seurat.obj1,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj1<-RunUMAP(seurat.obj1,dims=1:30,verbose = F)
        DimPlot(seurat.obj1,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_1seurat.png"))}
    seurat.obj.pca[[1]]<-seurat.obj1@reductions$pca@cell.embeddings
    var.seurat.obj[[1]]<-VariableFeatures(seurat.obj1)
    rm(seurat.obj1)

    #################################################
    #2.normalized selection: log(var(normalized))~log(mean(normalized))
    message("method2")
    seurat.obj2<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj2$cell_type<-rna_cell_label
    seurat.obj2<-NormalizeData(seurat.obj2,verbose = F)
    seurat.obj2@assays$RNA@counts<-seurat.obj2@assays$RNA@data
    seurat.obj2@assays$RNA@counts@x<-exp(seurat.obj2@assays$RNA@counts@x)-1
    seurat.obj2<-FindVariableFeatures(seurat.obj2,nfeatures=nfeatures,verbose = F)
    seurat.obj2<-ScaleData(seurat.obj2,verbose = F)
    seurat.obj2<-RunPCA(seurat.obj2,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj2<-RunUMAP(seurat.obj2,dims=1:30,verbose = F)
        DimPlot(seurat.obj2,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_2normalized.png"))}

    seurat.obj.pca[[2]]<-seurat.obj2@reductions$pca@cell.embeddings
    var.seurat.obj[[2]]<-VariableFeatures(seurat.obj2)
    rm(seurat.obj2)
    #################################################
    #3.log-normalized selection: log(var(log1p(normalized)))~log(mean(log1p(normalized)))
    message("method3")

    seurat.obj3<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj3$cell_type<-rna_cell_label
    seurat.obj3<-NormalizeData(seurat.obj3,verbose = F)
    seurat.obj3@assays$RNA@counts<-seurat.obj3@assays$RNA@data
    seurat.obj3<-FindVariableFeatures(seurat.obj3,nfeatures=nfeatures,verbose = F)
    seurat.obj3<-ScaleData(seurat.obj3,verbose = F)
    seurat.obj3<-RunPCA(seurat.obj3,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj3<-RunUMAP(seurat.obj3,dims=1:30,verbose = F)
        DimPlot(seurat.obj3,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_3lognormalized.png"))}


    seurat.obj.pca[[3]]<-seurat.obj3@reductions$pca@cell.embeddings
    var.seurat.obj[[3]]<-VariableFeatures(seurat.obj3)
    rm(seurat.obj3)
    #################################################
    #4.scran selection: var(log1p(normalized))~ mean(log1p(normalized))
    message("method4")
    ## scran
    library(scran)
    sce <- SingleCellExperiment(list(counts=rna_mat))
    sce <- logNormCounts(sce)
    dec <- modelGeneVar(sce)
    top.hvgs2 <- getTopHVGs(dec, n=2000)

    seurat.obj4<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj4$cell_type<-rna_cell_label
    seurat.obj4<-NormalizeData(seurat.obj4,verbose = F)
    VariableFeatures(seurat.obj4)<-top.hvgs2
    seurat.obj4<-ScaleData(seurat.obj4,verbose = F)
    seurat.obj4<-RunPCA(seurat.obj4,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj4<-RunUMAP(seurat.obj4,dims=1:30,verbose = F)
        DimPlot(seurat.obj4,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_4scran.png"))}

    seurat.obj.pca[[4]]<-seurat.obj4@reductions$pca@cell.embeddings
    var.seurat.obj[[4]]<-VariableFeatures(seurat.obj4)

    rm(seurat.obj4)
    rm(sce)
    rm(dec)
    #################################################
    #5.seurat
    message("method5")
    ## seurat mvp
    seurat.obj5<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj5$cell_type<-rna_cell_label
    seurat.obj5<-NormalizeData(seurat.obj5,verbose = F)
    seurat.obj5<-FindVariableFeatures(seurat.obj5,nfeatures=nfeatures,verbose = F,selection.method = "mvp")
    seurat.obj5<-ScaleData(seurat.obj5,verbose = F)
    seurat.obj5<-RunPCA(seurat.obj5,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj5<-RunUMAP(seurat.obj5,dims=1:30,verbose = F)
        DimPlot(seurat.obj5,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_5seurat_mvp.png"))}

    seurat.obj.pca[[5]]<-seurat.obj5@reductions$pca@cell.embeddings
    var.seurat.obj[[5]]<-VariableFeatures(seurat.obj5)

    rm(seurat.obj5)
    #################################################
    #6.cell ranger selection (dispersion based, scanpy):
    message("method6")
    library(reticulate)
    sc = import("scanpy")
    seurat.obj6<-CreateSeuratObject(rna_mat,verbose = F)
    scanpy_rna=sc$AnnData( t((seurat.obj6@assays$RNA@counts)))
    scanpy_rna$obs_names = as.character(1:ncol(rna_mat))
    scanpy_rna$var_names = rownames(rna_mat)
    #sc$pp$filter_cells(scanpy_rna, min_genes=200)
    sc$pp$filter_genes(scanpy_rna, min_cells=1)
    sc$pp$normalize_total(scanpy_rna, target_sum=1e4)
    sc$pp$log1p(scanpy_rna)
    sc$pp$highly_variable_genes(scanpy_rna,n_top_genes=as.integer(2000),flavor='cell_ranger')
    #[‘seurat’, ‘cell_ranger’, ‘seurat_v3’]
    top_scanpy_cell_ranger<-scanpy_rna$var_names[scanpy_rna$var$highly_variable]
    top_scanpy_cell_ranger<-sapply(1:length(top_scanpy_cell_ranger), function(i){
        top_scanpy_cell_ranger[i-1]
    })

    if(save_figure) seurat.obj6$cell_type<-rna_cell_label
    seurat.obj6<-NormalizeData(seurat.obj6,verbose = F)
    VariableFeatures(seurat.obj6)<-top_scanpy_cell_ranger
    seurat.obj6<-ScaleData(seurat.obj6,verbose = F)
    seurat.obj6<-RunPCA(seurat.obj6,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj6<-RunUMAP(seurat.obj6,dims=1:30,verbose = F)
        DimPlot(seurat.obj6,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_6scanpy_cell_ranger.png"))}

    seurat.obj.pca[[6]]<-seurat.obj6@reductions$pca@cell.embeddings
    var.seurat.obj[[6]]<-VariableFeatures(seurat.obj6)
    rm(seurat.obj6)
    rm(scanpy_rna)
    #################################################
    #7.log-PFlog1pPF selection : log(var(PFlog1pPF))~log(mean(PFlog1pPF))
    message("method7")

    seurat.obj7<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj7$cell_type<-rna_cell_label
    tmp<-seurat.obj7@assays$RNA@counts
    seurat.obj7@assays$RNA@counts<-rna_mat_PFlog1pPF
    seurat.obj7<-FindVariableFeatures(seurat.obj7,nfeatures=nfeatures,verbose = F)
    seurat.obj7@assays$RNA@counts<-tmp
    seurat.obj7<-NormalizeData(seurat.obj7,verbose = F)
    seurat.obj7<-ScaleData(seurat.obj7,verbose = F)
    seurat.obj7<-RunPCA(seurat.obj7,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj7<-RunUMAP(seurat.obj7,dims=1:30,verbose = F)
        DimPlot(seurat.obj7,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_7logPFlog1pPF.png"))}

    seurat.obj.pca[[7]]<-seurat.obj7@reductions$pca@cell.embeddings
    var.seurat.obj[[7]]<-VariableFeatures(seurat.obj7)
    rm(seurat.obj7)

    #################################################
    #8.PFlog1pPF selection: var(PFlog1pPF)~mean(PFlog1pPF)
    message("method8")
    library(scran)
    sce <- SingleCellExperiment(list(counts=rna_mat_PFlog1pPF))
    #sce <- logNormCounts(sce)
    sce@assays@data$logcounts<-sce@assays@data$counts
    dec <- modelGeneVar(sce)
    top.hvgs <- getTopHVGs(dec, n=2000)
    rm(sce)
    seurat.obj8<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj8$cell_type<-rna_cell_label
    #res_here<-hypervar_clip(rna_mat_PFlog1pPF)
    #VariableFeatures(seurat.obj8)<-res_here$data$feature[1:nfeatures]
    seurat.obj8<-NormalizeData(seurat.obj8,verbose = F)
    VariableFeatures(seurat.obj8)<-top.hvgs
    seurat.obj8<-ScaleData(seurat.obj8,verbose = F)
    seurat.obj8<-RunPCA(seurat.obj8,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj8<-RunUMAP(seurat.obj8,dims=1:30,verbose = F)
        DimPlot(seurat.obj8,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_8PFlog1pPF.png"))}

    seurat.obj.pca[[8]]<-seurat.obj8@reductions$pca@cell.embeddings
    var.seurat.obj[[8]]<-VariableFeatures(seurat.obj8)

    rm(seurat.obj8)

    ######## 9. Most original var count ~ mean count
    message("method9")
    library(scran)
    sce <- SingleCellExperiment(list(counts=rna_mat))
    #sce <- logNormCounts(sce)
    sce@assays@data$logcounts<-sce@assays@data$counts
    dec <- modelGeneVar(sce)
    top.hvgs <- getTopHVGs(dec, n=2000)
    rm(sce)
    seurat.obj9<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj9$cell_type<-rna_cell_label
    seurat.obj9<-NormalizeData(seurat.obj9,verbose = F)
    #res_here<-hypervar_clip(seurat.obj9@assays$RNA@counts)
    #VariableFeatures(seurat.obj9)<-res_here$data$feature[1:nfeatures]
    VariableFeatures(seurat.obj9)<-top.hvgs
    seurat.obj9<-ScaleData(seurat.obj9,verbose = F)
    seurat.obj9<-RunPCA(seurat.obj9,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj9<-RunUMAP(seurat.obj9,dims=1:30,verbose = F)
        DimPlot(seurat.obj9,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_9var_mean.png"))}

    seurat.obj.pca[[9]]<-seurat.obj9@reductions$pca@cell.embeddings
    var.seurat.obj[[9]]<-VariableFeatures(seurat.obj9)
    rm(seurat.obj9)
    #################################################
    #10.original normalized selection: var(normalized)~mean(normalized)
    message("method10")
    seurat.obj10<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj10$cell_type<-rna_cell_label
    seurat.obj10<-NormalizeData(seurat.obj10,verbose = F)
    rna_mat_norm<-seurat.obj10@assays$RNA@data
    rna_mat_norm@x<-exp(rna_mat_norm@x)-1
    library(scran)
    sce <- SingleCellExperiment(list(counts=rna_mat_norm))
    #sce <- logNormCounts(sce)
    sce@assays@data$logcounts<-sce@assays@data$counts
    dec <- modelGeneVar(sce)
    top.hvgs <- getTopHVGs(dec, n=2000)
    rm(sce)
    VariableFeatures(seurat.obj10)<-top.hvgs
    seurat.obj10<-ScaleData(seurat.obj10,verbose = F)
    seurat.obj10<-RunPCA(seurat.obj10,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj10<-RunUMAP(seurat.obj10,dims=1:30,verbose = F)
        DimPlot(seurat.obj10,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_10var_norm.png"))}

    seurat.obj.pca[[10]]<-seurat.obj10@reductions$pca@cell.embeddings
    var.seurat.obj[[10]]<-VariableFeatures(seurat.obj10)
    rm(seurat.obj10)
    #################################################
    #11.Seurat disp
    message("method11")
    seurat.obj11<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj11$cell_type<-rna_cell_label
    seurat.obj11<-NormalizeData(seurat.obj11,verbose = F)
    seurat.obj11<-FindVariableFeatures(seurat.obj11,nfeatures=nfeatures,verbose = F,selection.method = "disp")
    seurat.obj11<-ScaleData(seurat.obj11,verbose = F)
    seurat.obj11<-RunPCA(seurat.obj11,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj11<-RunUMAP(seurat.obj11,dims=1:30,verbose = F)
        DimPlot(seurat.obj11,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_11seurat_dispersion.png"))}

    seurat.obj.pca[[11]]<-seurat.obj11@reductions$pca@cell.embeddings
    var.seurat.obj[[11]]<-VariableFeatures(seurat.obj11)

    rm(seurat.obj11)

    #################################################
    #12.Seurat disp PFlog1pPF
    message("method12")
    seurat.obj12<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj12$cell_type<-rna_cell_label
    seurat.obj12<-NormalizeData(seurat.obj12,verbose = F)
    seurat.obj12@assays$RNA@data<-rna_mat_PFlog1pPF
    seurat.obj12<-FindVariableFeatures(seurat.obj12,nfeatures=nfeatures,verbose = F,selection.method = "disp")
    seurat.obj12<-NormalizeData(seurat.obj12,verbose = F)
    seurat.obj12<-ScaleData(seurat.obj12,verbose = F)
    seurat.obj12<-RunPCA(seurat.obj12,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj12<-RunUMAP(seurat.obj12,dims=1:30,verbose = F)
        DimPlot(seurat.obj12,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_12seurat_dispPF.png"))}

    seurat.obj.pca[[12]]<-seurat.obj12@reductions$pca@cell.embeddings
    var.seurat.obj[[12]]<-VariableFeatures(seurat.obj12)

    rm(seurat.obj12)

    savis_nth<- function(x, k) {
        if(sum(is.na(x))>0){
            x[is.na(x)]<-min(x[!is.na(x)])-0.1
        }
        ## might have problem when k is too large for nan case
        p <- length(x) - k
        if(p < 0){
            stop("savis_nth: input k too larger")
        }else if(p == 0){
            res<-1:length(x)
        }else{
            xp <- base::sort(x, partial=p)[p]
            res<-which(x > xp)
        }
        res
    }

    #################################################
    #13.highly expressed genes
    message("method13")
    seurat.obj13<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj13$cell_type<-rna_cell_label
    row_means_count<-rowMeans(seurat.obj13@assays$RNA@counts)
    hvg_rank<-savis_nth(row_means_count,nfeatures)
    VariableFeatures(seurat.obj13)<-rownames(seurat.obj13@assays$RNA@counts)[hvg_rank]
    seurat.obj13<-NormalizeData(seurat.obj13,verbose = F)
    seurat.obj13<-ScaleData(seurat.obj13,verbose = F)
    seurat.obj13<-RunPCA(seurat.obj13,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj13<-RunUMAP(seurat.obj13,dims=1:30,verbose = F)
        DimPlot(seurat.obj13,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_13mean_max_ct.png"))}

    seurat.obj.pca[[13]]<-seurat.obj13@reductions$pca@cell.embeddings
    var.seurat.obj[[13]]<-VariableFeatures(seurat.obj13)

    rm(seurat.obj13)


    #################################################
    #14.highly expressed genes normalized count
    message("method14")
    seurat.obj14<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj14$cell_type<-rna_cell_label
    seurat.obj14<-NormalizeData(seurat.obj14,verbose = F)
    seurat.obj14@assays$RNA@data@x<-exp(seurat.obj14@assays$RNA@data@x)-1
    row_means_count<-rowMeans(seurat.obj14@assays$RNA@data)
    hvg_rank2<-savis_nth(row_means_count,nfeatures)
    VariableFeatures(seurat.obj14)<-rownames(seurat.obj14@assays$RNA@counts)[hvg_rank2]
    seurat.obj14<-ScaleData(seurat.obj14,verbose = F)
    seurat.obj14<-RunPCA(seurat.obj14,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj14<-RunUMAP(seurat.obj14,dims=1:30,verbose = F)
        DimPlot(seurat.obj14,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_14mean_max_nc.png"))}

    seurat.obj.pca[[14]]<-seurat.obj14@reductions$pca@cell.embeddings
    var.seurat.obj[[14]]<-VariableFeatures(seurat.obj14)

    rm(seurat.obj14)


    #################################################
    #15.highly expressed genes lognormalized count
    message("method15")
    seurat.obj15<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj15$cell_type<-rna_cell_label
    seurat.obj15<-NormalizeData(seurat.obj15,verbose = F)
    row_means_count<-rowMeans(seurat.obj15@assays$RNA@data)
    hvg_rank3<-savis_nth(row_means_count,nfeatures)
    VariableFeatures(seurat.obj15)<-rownames(seurat.obj15@assays$RNA@counts)[hvg_rank3]
    seurat.obj15<-ScaleData(seurat.obj15,verbose = F)
    seurat.obj15<-RunPCA(seurat.obj15,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj15<-RunUMAP(seurat.obj15,dims=1:30,verbose = F)
        DimPlot(seurat.obj15,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_15mean_max_lognc.png"))}

    seurat.obj.pca[[15]]<-seurat.obj15@reductions$pca@cell.embeddings
    var.seurat.obj[[15]]<-VariableFeatures(seurat.obj15)

    rm(seurat.obj15)


    #################################################
    #16.scran selection: possion: var(log1p(normalized))~ mean(log1p(normalized))
    message("method16")
    ## scran
    library(scran)
    sce <- SingleCellExperiment(list(counts=rna_mat))
    sce <- logNormCounts(sce)
    dec<- modelGeneVarByPoisson(sce)
    #dec <- modelGeneVar(sce)
    #dec <- modelGeneCV2(sce)
    top.hvgs3 <- getTopHVGs(dec, n=2000)
    seurat.obj16<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj16$cell_type<-rna_cell_label
    seurat.obj16<-NormalizeData(seurat.obj16,verbose = F)
    VariableFeatures(seurat.obj16)<-top.hvgs3
    seurat.obj16<-ScaleData(seurat.obj16,verbose = F)
    seurat.obj16<-RunPCA(seurat.obj16,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj16<-RunUMAP(seurat.obj16,dims=1:30,verbose = F)
        DimPlot(seurat.obj16,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_16scran_possion.png"))}

    seurat.obj.pca[[16]]<-seurat.obj16@reductions$pca@cell.embeddings
    var.seurat.obj[[16]]<-VariableFeatures(seurat.obj16)

    rm(seurat.obj16)
    rm(sce)
    rm(dec)


    #################################################
    #17.Seurat SCT
    message("method17")
    seurat.obj17<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj17$cell_type<-rna_cell_label
    seurat.obj17<-SCTransform(seurat.obj17,variable.features.n=2000,verbose = F)
    seurat.obj17<-RunPCA(seurat.obj17,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj17<-RunUMAP(seurat.obj17,dims=1:30,verbose = F)
        DimPlot(seurat.obj17,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_17seurat_SCT.png"))}

    seurat.obj.pca[[17]]<-seurat.obj17@reductions$pca@cell.embeddings
    var.seurat.obj[[17]]<-VariableFeatures(seurat.obj17)

    #18.random
    message("method18")

    hvg<-sample(rownames(rna_mat)[which(rowSums(rna_mat)>0)],nfeatures)

    seurat.obj18<-CreateSeuratObject(rna_mat,verbose = F)
    if(save_figure) seurat.obj18$cell_type<-rna_cell_label
    seurat.obj18<-NormalizeData(seurat.obj18,verbose = F)
    VariableFeatures(seurat.obj18)<-hvg
    seurat.obj18<-ScaleData(seurat.obj18,verbose = F)
    seurat.obj18<-RunPCA(seurat.obj18,npcs=30,verbose=F)
    if(save_figure){
        seurat.obj18<-RunUMAP(seurat.obj18,dims=1:30,verbose = F)
        DimPlot(seurat.obj18,group.by = "cell_type")
        ggsave(paste0(save_path,dataset_name, "_18random.png"))}

    seurat.obj.pca[[18]]<-seurat.obj18@reductions$pca@cell.embeddings
    var.seurat.obj[[18]]<-VariableFeatures(seurat.obj18)

    newList<-list("seurat.obj.pca"=seurat.obj.pca,
                  "var.seurat.obj"=var.seurat.obj)

    saveRDS(newList,paste0(save_path,dataset_name,"_pcalist_base.rds"))

    newList<-list("seurat.obj.pca"=seurat.obj.pca,
                  "var.seurat.obj"=var.seurat.obj)
}


# saved as paste0(save_path,dataset_name,"_pcalist_mixture.rds")
mixture_hvg_pca<-function(rna_mat,dataset_name,
                          nfeatures = 2000,
                          save_path="benchmark/"){

    save_path = paste0(save_path,dataset_name,"/")
    seurat.obj.pca<-list()
    var.seurat.obj<-list()

    seurat.obj0<-CreateSeuratObject(rna_mat,verbose = F)
    rna_mat<-seurat.obj0@assays$RNA@counts
    rm(seurat.obj0)
    rna_mat_PFlog1pPF<-t(t(rna_mat)/colSums(rna_mat))*mean(colSums(rna_mat))
    rna_mat_PFlog1pPF<-log1p(rna_mat_PFlog1pPF)
    rna_mat_PFlog1pPF<-t(t(rna_mat_PFlog1pPF)/colSums(rna_mat_PFlog1pPF))*mean(colSums(rna_mat_PFlog1pPF))


    library(scran)
    sce <- SingleCellExperiment(list(counts=rna_mat))
    sce <- logNormCounts(sce)
    dec <- modelGeneVar(sce)
    dec.var <- dec@listData$bio
    dec.keep <- !is.na(dec.var) & dec.var > 0
    v_scran<-dec.var
    v_scran[!dec.keep]<-0
    rm(sce)
    rm(dec)

    seurat.obj1<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj1<-NormalizeData(seurat.obj1,verbose = F)
    seurat.obj1<-FindVariableFeatures(seurat.obj1,nfeatures=nfeatures,verbose = F)
    v_seuratv3<-seurat.obj1@assays$RNA@meta.features$vst.variance.standardized

    v_maxct<-rowMeans(seurat.obj1@assays$RNA@data)
    rm(seurat.obj1)
    seurat.obj11<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj11<-NormalizeData(seurat.obj11,verbose = F)
    seurat.obj11<-FindVariableFeatures(seurat.obj11,nfeatures=nfeatures,verbose = F,selection.method = "disp")
    v_disp<-seurat.obj11@assays$RNA@meta.features$mvp.dispersion
    v_disp[is.na(v_disp)]<-0
    v_disp[v_disp<0]<-0
    rm(seurat.obj11)

    seurat.obj10<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj10<-NormalizeData(seurat.obj10,verbose = F)
    rna_mat_norm<-seurat.obj10@assays$RNA@data
    rna_mat_norm@x<-exp(rna_mat_norm@x)-1
    library(scran)
    sce <- SingleCellExperiment(list(counts=rna_mat_norm))
    #sce <- logNormCounts(sce)
    sce@assays@data$logcounts<-sce@assays@data$counts
    dec <- modelGeneVar(sce)
    top.hvgs <- getTopHVGs(dec, n=2000)
    dec.var <- dec@listData$bio
    dec.keep <- !is.na(dec.var) & dec.var > 0
    v_scran_nc<-dec.var
    v_scran_nc[!dec.keep]<-0
    rm(sce)
    rm(dec)

    sce <- SingleCellExperiment(list(counts=rna_mat))
    sce <- logNormCounts(sce)
    dec<- modelGeneVarByPoisson(sce)
    dec.var <- dec@listData$bio
    dec.keep <- !is.na(dec.var) & dec.var > 0
    v_scran_pos<-dec.var
    v_scran_pos[!dec.keep]<-0
    rm(sce)
    rm(dec)

    mixture_rank<-function(input_lst){
        input_lst_order<-list()
        for(i in 1:length(input_lst)){
            input_lst_order[[i]]<-order(order(input_lst[[i]],decreasing = T))
        }

        apply(matrix(unlist(input_lst_order),
                     ncol=length(input_lst_order),byrow = FALSE),1,FUN = min)
    }
    mixture_mat<-cbind(v_scran_nc,v_scran,v_scran_pos,v_disp,v_maxct,v_seuratv3)
    mixture_index_list<-list(c(1,2),c(1,3),c(1,4),c(1,5),c(2,3),c(2,4),c(2,5),c(3,4),c(3,5),c(4,5),
                          c(1,2,3),c(1,2,4),c(1,2,5),c(1,3,4),c(1,3,5),c(1,4,5),
                          c(2,3,4),c(2,3,5),c(2,4,5),c(3,4,5),
                          c(1,2,3,4),c(1,2,3,5),c(1,2,4,5),c(1,3,4,5),c(2,3,4,5),
                          c(1,2,3,4,5),c(1,2,3,4,5,6))

    hvg<-list()

    for(kk in 1:length(mixture_index_list)){
        cur_mixture_index<-mixture_index_list[[kk]]
        cur_mixture_list<-list()
        cur_mixture_list_count<-1
        for( i in cur_mixture_index){
            cur_mixture_list[[cur_mixture_list_count]]<-mixture_mat[,i]
            cur_mixture_list_count<-cur_mixture_list_count+1
        }
        hvg[[kk]]<-order(mixture_rank(cur_mixture_list))[1:nfeatures]
    }

    print("Method_Mixture")
    rna_mat_norm<-NormalizeData(rna_mat,verbose = FALSE)

    for( i in 1:length(hvg)){
        message(mixture_index_list[i])
        cur_index<-i
        rna_mat_norm_hvg<-rna_mat_norm[hvg[[i]],]
        rna_mat_scale<-ScaleData(rna_mat_norm_hvg,verbose = F)
        suppressWarnings(rna_mat_pca <- RunPCA(
            object = rna_mat_scale,
            features = rownames(rna_mat_scale),
            npcs = 30,
            verbose = FALSE)@cell.embeddings)
        seurat.obj.pca[[cur_index]]<-rna_mat_pca
        var.seurat.obj[[cur_index]]<-rownames(rna_mat)[hvg[[i]]]
    }


    newList<-list("seurat.obj.pca"=seurat.obj.pca,
                  "var.seurat.obj"=var.seurat.obj)

    saveRDS(newList,paste0(save_path,dataset_name,"_pcalist_mixture.rds"))
    print("Finish Methods Mixture!")
    return(newList)
}


