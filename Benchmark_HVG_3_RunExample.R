# Author: Ruzhang Zhao
## Example: CITEseq
## bmcite
library(SeuratData)
bm <- LoadData(ds = "bmcite")
expr<-bm@assays$RNA@counts
pro<-bm@assays$ADT@data
# baseline methods
.<-hvg_pca(expr,dataset_name="bmcite")
# mixture methods
.<-mixture_hvg_pca(expr,dataset_name="bmcite")
rm(expr)
.<-evaluate_hvg_CITEseq(pro=pro,dataset_name="bmcite",evaluate_type="base")
.<-evaluate_hvg_CITEseq(pro=pro,dataset_name="bmcite",evaluate_type="mixture")

## Example: MultiomeATAC
## pbmc3k_multi
## Assume we have processed pbmc3k_multi
## and get pbmc3k_multi_rna_mat and pbmc3k_multi_lsi
# baseline methods
.<-hvg_pca(pbmc3k_multi_rna_mat,dataset_name="pbmc3k_multi")
# mixture methods
.<-mixture_hvg_pca(pbmc3k_multi_rna_mat,dataset_name="pbmc3k_multi")
rm(pbmc3k_multi_rna_mat)
.<-evaluate_hvg_CITEseq(pro=pbmc3k_multi_lsi,dataset_name="pbmc3k_multi",evaluate_type="base",scale_protein=FALSE)
.<-evaluate_hvg_CITEseq(pro=pbmc3k_multi_lsi,dataset_name="pbmc3k_multi",evaluate_type="mixture",scale_protein=FALSE)





## Example: Cell Sorting
## duo4_pbmc
library(DuoClustering2018)
sce_names<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
sn<-sce_names[1]
sce<-do.call(paste0("sce_full_",sn),list())
expr<-sce@assays$data@listData$counts
cell_label<-sce$phenoid
# baseline methods
.<-hvg_pca(expr,dataset_name="duo4_pbmc")
# If we need plots
#.<-hvg_pca(expr,dataset_name="duo4_pbmc",rna_cell_label=cell_label,save_figure=TRUE)
# mixture methods
.<-mixture_hvg_pca(expr,dataset_name="duo4_pbmc")
rm(expr)
.<-evaluate_hvg_sorting(label=cell_label,dataset_name="duo4_pbmc",evaluate_type="base")
.<-evaluate_hvg_sorting(label=cell_label,dataset_name="duo4_pbmc",evaluate_type="mixture")
