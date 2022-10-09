# Author: Ruzhang Zhao
# Evaluate with cell sorting
# Criterion 1: ARI
ARILouvain<-function(
        embedding,
        cell_label,
        maxit = 25
){
    N_label<-length(unique(cell_label))
    snn_<- FindNeighbors(object = embedding,
                         nn.method = "rann",
                         verbose = F)$snn
    cluster_label <- FindClusters(snn_,
                                  resolution = 0,
                                  verbose = F)[[1]]
    N_cluster0<-length(unique(cluster_label))
    if(N_cluster0>=N_label){
        print("Resolution 0 still larger than true")
        return(c(N_cluster0,
                 adjustedRandIndex(cluster_label,cell_label)))
    }

    snn_<- FindNeighbors(object = embedding,
                         nn.method = "rann",
                         verbose = F)$snn
    cluster_label <- FindClusters(snn_,
                                  resolution = 1,
                                  verbose = F)[[1]]
    N_cluster1<-length(unique(cluster_label))
    if(N_cluster1 <= N_label ){
        return(c(N_cluster1,
                 adjustedRandIndex(cluster_label,cell_label)))
    }else if (N_cluster1 > N_label){
        a<-0
        b<-1
    }
    keepsign<-TRUE
    iter_<-0
    while (keepsign) {
        c<-(a+b)/2
        snn_<- FindNeighbors(object = embedding,
                             nn.method = "rann",
                             verbose = F)$snn
        cluster_label <- FindClusters(snn_,
                                      resolution = c,
                                      verbose = F)[[1]]
        N_clusterc<-length(unique(cluster_label))
        if( N_clusterc == N_label){
            keepsign<-FALSE
            return(c(N_clusterc,
                     adjustedRandIndex(cluster_label,cell_label)))
        }else if (N_clusterc < N_label){
            a<-c
        }else{
            b<-c
        }
        iter_<-iter_+1
        if(iter_ > maxit){
            keepsign<-FALSE
        }
    }
    return(c(N_clusterc,
             adjustedRandIndex(cluster_label,cell_label)))
}

# Criterion 2: Variance Ratio
withinVSbetween<-function(embedding,cell_label){
    cell_label <- as.character(cell_label)
    # sorted unique cluster label
    label_index<-sort(unique(as.character(cell_label)))
    if(sum(label_index == "") >0){
        label_index<-label_index[-which(label_index == "")]
    }
    # number of cluster label
    N_label<-length(label_index)
    within<-sapply(1:N_label, function(i){
        index_i<-which(cell_label == label_index[i])
        (length(index_i)-1)*Rfast::colVars(embedding[index_i,])
    })

    all<-(sum(nrow(embedding)-1)*colVars(embedding))
    prop1<-sum(rowSums(within))/(sum(all)- sum(rowSums(within)))
    list("variance"=prop1)
}

# Criterion 3: Nearest Neighbor Accuracy
knn_accuracy<-function(embedding,cell_label){
    N_cell<-nrow(embedding)
    pre_label<-sapply(c(1:N_cell), function(cell){
        train_label<-cell_label[-cell]
        a<-FNN::knn(train = embedding[-cell,],
                    test = embedding[cell,],
                    cl = train_label,k = 3,prob = T)
        a[1]
    })
    sum(pre_label == cell_label)/length(cell_label)
}

## Comprehensive Evaluation for cell sorting.
evaluate_hvg_sorting<-function(label,
                               dataset_name,nfeatures = 2000,
                               save_path="benchmark/",
                               save = TRUE,evaluate_type="base"
){
    save_path = paste0(save_path,dataset_name,"/")
    if(evaluate_type=="base"){
        newList<-readRDS(paste0(save_path,dataset_name,"_pcalist_base.rds"))
    }else if(evaluate_type=="mixture"){
        newList<-readRDS(paste0(save_path,dataset_name,"_pcalist_mixture.rds"))
    }else{
        newList<-readRDS(paste0(save_path,dataset_name,"_pcalist.rds"))
    }

    seurat.obj.pca<-newList$seurat.obj.pca
    var.seurat.obj<-newList$var.seurat.obj


    Num_method<-length(seurat.obj.pca)
    len_between<-matrix(0,nrow = Num_method,ncol = Num_method)
    diag(len_between)<-nfeatures
    for(i in 1:(Num_method-1) ){
        for(j in (i+1):Num_method){
            len_between[j,i]<-length(intersect(var.seurat.obj[[i]],var.seurat.obj[[j]]))
            len_between[i,j]<-len_between[j,i]
        }
    }
    #################################################
    #Within CellType Variance / Between CellType Variance (Ratio)
    ## Part 1 within and between variance
    message("var_ratio")
    start_t<-Sys.time()
    withinVSbetween<-function(embedding,cell_label){
        cell_label <- as.character(cell_label)
        # sorted unique cluster label
        label_index<-sort(unique(as.character(cell_label)))
        if(sum(label_index == "") >0){
            label_index<-label_index[-which(label_index == "")]
        }
        # number of cluster label
        N_label<-length(label_index)
        within<-sapply(1:N_label, function(i){
            index_i<-which(cell_label == label_index[i])
            (length(index_i)-1)*Rfast::colVars(embedding[index_i,])
        })

        all<-(sum(nrow(embedding)-1)*colVars(embedding))
        prop1<-sum(rowSums(within))/(sum(all)- sum(rowSums(within)))
        list("variance"=prop1)
    }

    variance_ratio<-list()
    library(caret)
    Nosample<-FALSE
    if(length(label)<10000){
        Nosample<-TRUE
    }

    if(!Nosample){
        set.seed(1)
        index_sample_pca<-createDataPartition(label,p = min(1,10000/length(label)))$Resample1
        set.seed(2)
        index_sample_pca1<-createDataPartition(label,p = min(1,10000/length(label)))$Resample1
    }

    for(i in 1:Num_method){
        if(Nosample){
            withinbetween_ratio<- withinVSbetween(seurat.obj.pca[[i]],label)
            variance_ratio[['CP10k']][[i]]<-(withinbetween_ratio$variance)

        }else{
            withinbetween_ratio<- withinVSbetween(seurat.obj.pca[[i]][index_sample_pca,],label[index_sample_pca])
            withinbetween_ratio1<- withinVSbetween(seurat.obj.pca[[i]][index_sample_pca1,],label[index_sample_pca1])
            variance_ratio[['CP10k']][[i]]<-(withinbetween_ratio$variance+withinbetween_ratio1$variance)/2
        }
    }

    if(evaluate_type=="base"){
        dataset_name<-paste0(dataset_name,"_base")
    }else if(evaluate_type=="mixture"){
        dataset_name<-paste0(dataset_name,"_mixture")
    }


    end_t<-Sys.time()
    print(paste0("VAR Time"))
    print(end_t-start_t)

    if(save){
        newList<-list(
            "len_between"=len_between,
            "var_ratio"=variance_ratio)
        saveRDS(newList,paste0(save_path,dataset_name,"CP10kPFfinal6.rds"))
    }
    #################################################
    #ARI with Louvain clustering  same number
    message("ari_louvain")
    start_t<-Sys.time()
    ARILouvain<-function(
        embedding,
        cell_label,
        maxit = 25
    ){
        N_label<-length(unique(cell_label))
        snn_<- FindNeighbors(object = embedding,
                             nn.method = "rann",
                             verbose = F)$snn
        cluster_label <- FindClusters(snn_,
                                      resolution = 0,
                                      verbose = F)[[1]]
        N_cluster0<-length(unique(cluster_label))
        if(N_cluster0>=N_label){
            print("Resolution 0 still larger than true")
            return(c(N_cluster0,
                     adjustedRandIndex(cluster_label,cell_label)))
        }

        snn_<- FindNeighbors(object = embedding,
                             nn.method = "rann",
                             verbose = F)$snn
        cluster_label <- FindClusters(snn_,
                                      resolution = 1,
                                      verbose = F)[[1]]
        N_cluster1<-length(unique(cluster_label))
        if(N_cluster1 <= N_label ){
            return(c(N_cluster1,
                     adjustedRandIndex(cluster_label,cell_label)))
        }else if (N_cluster1 > N_label){
            a<-0
            b<-1
        }
        keepsign<-TRUE
        iter_<-0
        while (keepsign) {
            c<-(a+b)/2
            snn_<- FindNeighbors(object = embedding,
                                 nn.method = "rann",
                                 verbose = F)$snn
            cluster_label <- FindClusters(snn_,
                                          resolution = c,
                                          verbose = F)[[1]]
            N_clusterc<-length(unique(cluster_label))
            if( N_clusterc == N_label){
                keepsign<-FALSE
                return(c(N_clusterc,
                         adjustedRandIndex(cluster_label,cell_label)))
            }else if (N_clusterc < N_label){
                a<-c
            }else{
                b<-c
            }
            iter_<-iter_+1
            if(iter_ > maxit){
                keepsign<-FALSE
            }
        }
        return(c(N_clusterc,
                 adjustedRandIndex(cluster_label,cell_label)))
    }

    ari_list<-list()
    for(i in 1:Num_method){
        print(i)
        if(Nosample){
            ari_list[['CP10k']][[i]]<-(ARILouvain(seurat.obj.pca[[i]],label)[2] )
        }else{
            ari_list[['CP10k']][[i]]<-(ARILouvain(seurat.obj.pca[[i]][index_sample_pca,],label[index_sample_pca])[2] +ARILouvain(seurat.obj.pca[[i]][index_sample_pca1,],label[index_sample_pca1])[2] )/2
        }

    }

    end_t<-Sys.time()
    print(paste0("ARI Time"))
    print(end_t-start_t)
    if(save){
        newList<-list(
            "len_between"=len_between,
            "var_ratio"=variance_ratio,
            "ari"=ari_list)

        saveRDS(newList,paste0(save_path,dataset_name,"CP10kPFfinal6.rds"))
    }
    #################################################
    # 3NN accuracy
    message("3nn")
    start_t<-Sys.time()

    set.seed(3)
    if(length(label)>1e4){
        index_sample_pca<-createDataPartition(label,p = 1e4/length(label))$Resample1
    }else{
        index_sample_pca<-1:length(label)
    }

    knn_accuracy<-function(embedding,cell_label){
        N_cell<-nrow(embedding)
        pre_label<-sapply(c(1:N_cell), function(cell){
            train_label<-cell_label[-cell]
            a<-FNN::knn(train = embedding[-cell,],
                        test = embedding[cell,],
                        cl = train_label,k = 3,prob = T)
            a[1]
        })
        sum(pre_label == cell_label)/length(cell_label)
    }

    nn_list<-list()
    nn_listPF<-list()
    for(i in 1:Num_method){
        nn_list[['CP10k']][[i]]<-knn_accuracy(seurat.obj.pca[[i]][index_sample_pca,],label[index_sample_pca])

    }
    end_t<-Sys.time()
    print(paste0("KNN Time"))
    print(end_t-start_t)
    if(save){
        newList<-list(
            "len_between"=len_between,
            "var_ratio"=variance_ratio,
            "ari"=ari_list,
            "3nn"=nn_list)

        saveRDS(newList,paste0(save_path,dataset_name,"CP10kPFfinal6.rds"))
    }

    #ARI with Kmeans same number
    message("ari_kmeans")
    start_t<-Sys.time()
    ARIKmeans<-function(
        embedding,
        cell_label
    ){
        N_label<-length(unique(cell_label))
        clu_res<-kmeans(embedding,centers = N_label)
        adjustedRandIndex(clu_res$cluster,cell_label)
    }

    ari_kmeans_list<-list()
    ari_kmeans_listPF<-list()
    for(i in 1:Num_method){

        ari_kmeans_list[['CP10k']][[i]]<-ARIKmeans(seurat.obj.pca[[i]][index_sample_pca,],label[index_sample_pca])

    }
    end_t<-Sys.time()
    print(paste0("ARI Kmeans Time"))
    print(end_t-start_t)
    newList<-list(
        "len_between"=len_between,
        "var_ratio"=variance_ratio,
        "ari"=ari_list,
        "3nn"=nn_list,
        "ari_kmeans"=ari_kmeans_list)
    if(save){
        saveRDS(newList,paste0(save_path,dataset_name,"CP10kPFfinal6.rds"))
    }
    return(newList)
}

# Evaluate with CITEseq & MultiomeATAC
# Criterion 1: Variance Ratio
withinVSbetween_variance_distance<-function(embedding,surface_markers,all_var,dmat,resolution = 0.8){

    snn_<- FindNeighbors(object = embedding,
                         nn.method = "rann",
                         verbose = F)$snn
    cluster_label <- FindClusters(snn_,
                                  resolution = resolution,
                                  verbose = F)[[1]]

    cluster_label <- as.numeric(as.character(cluster_label))
    # sorted unique cluster label
    label_index<-sort(as.numeric(
        unique(as.character(cluster_label))))
    # number of cluster label
    N_label<-length(label_index)
    print(paste0("Current Label Number is ",N_label))
    within<-sapply(1:N_label, function(i){
        index_i<-which(cluster_label == label_index[i])
        (length(index_i)-1)*Rfast::rowVars(surface_markers[,index_i])
    })

    var_r<-rowSums(within)/(all_var- rowSums(within))


    return(list("var"=var_r))
}

# Criterion 2: Nearest Neighbor Mean Square Error
knn_regression<-function(embedding,surface_markers){
    pre_surface_markers<-sapply(1:nrow(embedding), function(cell){
        surface_markers1<-surface_markers[,-cell]
        a<-FNN::knn(train=embedding[-cell,],
                    test= embedding[cell,],
                    cl = rep(1,nrow(embedding)-1),
                    k=3)
        rowMeans(surface_markers1[,attr(a,"nn.index")[1,]])
    })

    sum((pre_surface_markers - surface_markers)^2)
}

# Criterion 3: Nearest Neighbor Distance Ratio
knn_ratio<-function(embedding,surface_markers,scale_surface_markers){
    N_cell<-nrow(embedding)
    N_suf<-nrow(surface_markers)
    surface_markers<-t(surface_markers)
    scale_surface_markers<-t(scale_surface_markers)
    pre_label<-sapply(c(1:N_cell), function(cell){
        train_label<-rep(1,nrow(embedding)-1)
        a<-FNN::knn(train = embedding[-cell,],
                    test = embedding[cell,],
                    cl = train_label,k = 100,prob = T)
        nn_index<-c(1:N_cell)[-cell][attr(a,"nn.index")[1,]]
        a<-sapply(1:ncol(surface_markers), function(i){
            mean(abs(surface_markers[cell,i]-surface_markers[nn_index,i]))/
                max(mean(abs(surface_markers[cell,i]-surface_markers[-c(nn_index,cell),i])),1e-8)
        })
        a1<-mean(pdist::pdist(scale_surface_markers[cell,],scale_surface_markers[nn_index,])@dist)/
            mean(pdist::pdist(scale_surface_markers[cell,],scale_surface_markers[-c(nn_index,cell),])@dist)
        c(a1,a)
    })
    rowMeans(pre_label)
}


## Comprehensive Evaluation for CITEseq & MultiomeATAC
## For MultiomeATAC, pro = LSI and set scale_protein = FALSE
evaluate_hvg_CITEseq<-function(pro,
                               dataset_name,nfeatures = 2000,
                               save_path="benchmark/",
                               save = TRUE,scale_protein=TRUE,
                               evaluate_type="base"){

    save_path = paste0(save_path,dataset_name,"/")
    if(evaluate_type=="base"){
        newList<-readRDS(paste0(save_path,dataset_name,"_pcalist_base.rds"))
    }else if(evaluate_type=="mixture"){
        newList<-readRDS(paste0(save_path,dataset_name,"_pcalist_mixture.rds"))
    }else{
        newList<-readRDS(paste0(save_path,dataset_name,"_pcalist.rds"))
    }
    if(evaluate_type=="base"){
        dataset_name<-paste0(dataset_name,"_base")
    }else if(evaluate_type=="mixture"){
        dataset_name<-paste0(dataset_name,"_mixture")
    }
    seurat.obj.pca<-newList$seurat.obj.pca
    var.seurat.obj<-newList$var.seurat.obj

    Num_method<-length(seurat.obj.pca)

    len_between<-matrix(0,nrow = Num_method,ncol = Num_method)
    diag(len_between)<-nfeatures
    for(i in 1: (Num_method-1)){
        for(j in (i+1):Num_method){
            len_between[j,i]<-length(intersect(var.seurat.obj[[i]],var.seurat.obj[[j]]))
            len_between[i,j]<-len_between[j,i]
        }
    }
    #################################################
    #Within CellType Variance / Between CellType Variance (Ratio)
    ## Part 1 within and between variance
    message("knn_ratio")
    #Within CellType Distance / Between CellType Distance (Ratio2)
    #message("distance_ratio")
    Nosample<-FALSE
    if(nrow(seurat.obj.pca[[1]])>10000){
        set.seed(1)
        index_sample_pca<-sample(1:ncol(pro),size = 10000)
        set.seed(2)
        index_sample_pca1<-sample(1:ncol(pro),size = 10000)
        set.seed(3)
        index_sample_pca2<-sample(1:ncol(pro),size = 10000)
    }else{
        index_sample_pca<-1:ncol(pro)
        Nosample<-TRUE
    }
    if(scale_protein){
        scale_pro<-CreateSeuratObject(pro)
        scale_pro <- NormalizeData(scale_pro, normalization.method = "CLR", margin = 2,verbose=F)
        scale_pro <- ScaleData(scale_pro)
        scale_pro <- scale_pro@assays$RNA@scale.data
    }else{
        scale_pro<-pro
    }

    knn_ratio<-function(embedding,surface_markers,scale_surface_markers){
        N_cell<-nrow(embedding)
        N_suf<-nrow(surface_markers)
        surface_markers<-t(surface_markers)
        scale_surface_markers<-t(scale_surface_markers)
        pre_label<-sapply(c(1:N_cell), function(cell){
            train_label<-rep(1,nrow(embedding)-1)
            a<-FNN::knn(train = embedding[-cell,],
                        test = embedding[cell,],
                        cl = train_label,k = 100,prob = T)
            nn_index<-c(1:N_cell)[-cell][attr(a,"nn.index")[1,]]
            a<-sapply(1:ncol(surface_markers), function(i){
                mean(abs(surface_markers[cell,i]-surface_markers[nn_index,i]))/
                    max(mean(abs(surface_markers[cell,i]-surface_markers[-c(nn_index,cell),i])),1e-8)
            })
            a1<-mean(pdist::pdist(scale_surface_markers[cell,],scale_surface_markers[nn_index,])@dist)/
                mean(pdist::pdist(scale_surface_markers[cell,],scale_surface_markers[-c(nn_index,cell),])@dist)
            c(a1,a)
        })
        rowMeans(pre_label)
    }

    knnratio_ratio<-list()
    start_t<-Sys.time()

    for(i in 1:Num_method){

        withinbetween_ratio<- knn_ratio(seurat.obj.pca[[i]][index_sample_pca,],pro[,index_sample_pca],scale_pro[,index_sample_pca])
        knnratio_ratio[['CP10k']][[i]]<-list("rep0"=as.matrix(withinbetween_ratio))

    }
    end_t<-Sys.time()
    print(paste0("KNNRATIO Time"))
    print(end_t-start_t)
    if(save){
        newList<-list(
            "len_between"=len_between,
            "knn_ratio"=knnratio_ratio)
        saveRDS(newList,paste0(save_path,dataset_name,"CP10kPFfinal6.rds"))
    }

    withinVSbetween_variance_distance<-function(embedding,surface_markers,all_var,dmat,resolution = 0.8){

        snn_<- FindNeighbors(object = embedding,
                             nn.method = "rann",
                             verbose = F)$snn
        cluster_label <- FindClusters(snn_,
                                      resolution = resolution,
                                      verbose = F)[[1]]

        cluster_label <- as.numeric(as.character(cluster_label))
        # sorted unique cluster label
        label_index<-sort(as.numeric(
            unique(as.character(cluster_label))))
        # number of cluster label
        N_label<-length(label_index)
        print(paste0("Current Label Number is ",N_label))
        within<-sapply(1:N_label, function(i){
            index_i<-which(cluster_label == label_index[i])
            (length(index_i)-1)*Rfast::rowVars(surface_markers[,index_i])
        })

        var_r<-rowSums(within)/(all_var- rowSums(within))


        return(list("var"=var_r))
    }
    message("var_ratio")
    start_t<-Sys.time()
    all_var<-(ncol(scale_pro[,index_sample_pca])-1)*rowVars(scale_pro[,index_sample_pca])
    dmat<-dist(t(scale_pro[,index_sample_pca]))
    if(!Nosample){
        all_var1<-(ncol(scale_pro[,index_sample_pca1])-1)*rowVars(scale_pro[,index_sample_pca1])
        dmat1<-dist(t(scale_pro[,index_sample_pca1]))
        all_var2<-(ncol(scale_pro[,index_sample_pca2])-1)*rowVars(scale_pro[,index_sample_pca2])
        dmat2<-dist(t(scale_pro[,index_sample_pca2]))
    }


    variance_ratio<-list()

    for(i in 1:Num_method){
        if(Nosample){
            withinbetween_ratio<- withinVSbetween_variance_distance(seurat.obj.pca[[i]][index_sample_pca,],scale_pro[,index_sample_pca],all_var,dmat)
            variance_ratio[['CP10k']][[i]]<-
                list("rep0"=withinbetween_ratio$var)
        }else{
            withinbetween_ratio<- withinVSbetween_variance_distance(seurat.obj.pca[[i]][index_sample_pca,],scale_pro[,index_sample_pca],all_var,dmat)
            withinbetween_ratio1<- withinVSbetween_variance_distance(seurat.obj.pca[[i]][index_sample_pca1,],scale_pro[,index_sample_pca1],all_var1,dmat1)
            withinbetween_ratio2<- withinVSbetween_variance_distance(seurat.obj.pca[[i]][index_sample_pca2,],scale_pro[,index_sample_pca2],all_var2,dmat2)

            variance_ratio[['CP10k']][[i]]<-
                list("rep0"=withinbetween_ratio$var,
                     "rep1"=withinbetween_ratio1$var,
                     "rep2"=withinbetween_ratio2$var)
            sil_ratio[['CP10k']][[i]]<-list("rep0"=withinbetween_ratio$sil,
                                            "rep1"=withinbetween_ratio1$sil,
                                            "rep2"=withinbetween_ratio2$sil)
        }
    }
    end_t<-Sys.time()
    print("VAR Time")
    print(end_t-start_t)
    if(save){
        newList<-list(
            "len_between"=len_between,
            "knn_ratio"=knnratio_ratio,
            "var_ratio"=variance_ratio)

        saveRDS(newList,paste0(save_path,dataset_name,"CP10kPFfinal6.rds"))
    }


    #################################################
    # 3NN accuracy
    message("3nn")

    knn_regression<-function(embedding,surface_markers){
        pre_surface_markers<-sapply(1:nrow(embedding), function(cell){
            surface_markers1<-surface_markers[,-cell]
            a<-FNN::knn(train=embedding[-cell,],
                        test= embedding[cell,],
                        cl = rep(1,nrow(embedding)-1),
                        k=3)
            rowMeans(surface_markers1[,attr(a,"nn.index")[1,]])
        })

        sum((pre_surface_markers - surface_markers)^2)
    }

    nn_list<-list()
    nn_listPF<-list()
    start_t<-Sys.time()
    for(i in 1:Num_method){
        nn_list[['CP10k']][[i]]<-knn_regression(seurat.obj.pca[[i]][index_sample_pca,],scale_pro[,index_sample_pca])
    }
    end_t<-Sys.time()
    print(paste0("KNN Time"))
    print(end_t-start_t)

    newList<-list(
        "len_between"=len_between,
        "knn_ratio"=knnratio_ratio,
        "var_ratio"=variance_ratio,
        "3nn"=nn_list)
    if(save){
        saveRDS(newList,paste0(save_path,dataset_name,"CP10kPFfinal6.rds"))
    }
    return(newList)
}

