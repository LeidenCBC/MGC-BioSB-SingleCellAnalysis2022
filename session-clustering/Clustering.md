Clustering
================

Created by: Ahmed Mahfouz

Edited by: Mohammed Charrout, Lieke Michielsen

# Overview

In this tutorial we will look at different approaches to cluster
scRNA-seq datasets in order to characterize the different subgroups of
cells. Using unsupervised clustering, we will try to identify groups of
cells based on the similarities of the transcriptomes without any prior
knowledge of the labels.

Load required packages:

``` r
suppressMessages(require(Seurat))
```

## Datasets

Again, we will continue with dataset the you have preprocessed and
visualized in the previous practicals. Let’s start by loading the data
again.

During the previous practical, we have already selected highly variable
genes. This step is also to decide which genes to use when clustering
the cells. Single cell RNA-seq can profile a huge number of genes in a
lot of cells. But most of the genes are not expressed enough to provide
a meaningful signal and are often driven by technical noise. Including
them could potentially add some unwanted signal that would blur the
biological variation. Moreover gene filtering can also speed up the
computational time for downstream analysis.

``` r
pbmc <- readRDS('../session-dimensionalityreduction/pbmc3k.rds')
```

## Clustering

### Hierarchical clustering

``` r
# Get scaled counts from the Seurat object
scaled_pbmc <- pbmc@assays$RNA@scale.data

# Calculate Distances (default: Euclidean distance)
distance_euclidean <- dist(t(scaled_pbmc))

#Perform hierarchical clustering using ward linkage
ward_hclust_euclidean <- hclust(distance_euclidean,method = "ward.D2")
plot(ward_hclust_euclidean, main = "dist = eucledian, Ward linkage", labels=FALSE)
```

![](Clustering_files/figure-gfm/hierarchical_eucledian_ward-1.png)<!-- -->

Now cut the dendrogram to generate 10 clusters and plot the cluster
labels and the previously given celltype labels on the t-SNE plot. For
now, we just pick 10, but you can of course vary this number to see how
it influences your results.

``` r
#Cutting the cluster tree to make 10 groups
cluster_hclust <- cutree(ward_hclust_euclidean,k = 10)
pbmc@meta.data$cluster_hclust <- factor(cluster_hclust)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "cluster_hclust")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/hierarchical_eucledian_ward_pcaplot-1.png)<!-- -->

Now let’s try a different distance measure. A commonly used distance
measure is 1 - correlation.

``` r
# Calculate Distances (1 - correlation)
C <- cor(scaled_pbmc)

# Run clustering based on the correlations, where the distance will 
# be 1-correlation, e.g. higher distance with lower correlation.
distance_corr <- as.dist(1-C) 
    
#Perform hierarchical clustering using ward linkage
ward_hclust_corr <- hclust(distance_corr,method="ward.D2")
plot(ward_hclust_corr, main = "dist = 1-corr, Ward linkage", labels=FALSE)
```

![](Clustering_files/figure-gfm/hierarchical_corr_ward-1.png)<!-- -->

Again, let’s cut the dendrogram to generate 10 clusters and plot the
cluster labels on the t-SNE plot.

``` r
#Cutting the cluster tree to make 10 groups
cluster_hclust <- cutree(ward_hclust_corr,k = 10)
pbmc@meta.data$cluster_hclust <- factor(cluster_hclust)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "cluster_hclust")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/hierarchical_corr_ward_pcaplot-1.png)<!-- -->

Instead of changing the distance metric, we can change the linkage
method. Instead of using Ward’s method, let’s use complete linkage.

``` r
#Perform hierarchical clustering using complete linkage & euclidean distance
comp_hclust_eucledian <- hclust(distance_euclidean,method = "complete")
plot(comp_hclust_eucledian, main = "dist = euclidean, complete linkage", labels=FALSE)
```

![](Clustering_files/figure-gfm/hierarchical_eucledian_complete-1.png)<!-- -->

Once more, let’s cut the dendrogram to generate 10 clusters and plot the
cluster labels on the t-SNE plot.

``` r
#Cutting the cluster tree to make 10 groups
cluster_hclust <- cutree(comp_hclust_eucledian,k = 10)
pbmc@meta.data$cluster_hclust <- factor(cluster_hclust)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "cluster_hclust")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/hierarchical_eucledian_complete_pcaplot-1.png)<!-- -->
As you can see, these linkage methods and distances cluster the data
differently. If you want, there are even more distance measures and
linkage methods to play around with.

### K-means

Next, we will try the k-means algorithm on the scaled data.

``` r
pbmc_kmeans <- kmeans(x = t(scaled_pbmc), centers = 10)
pbmc@meta.data$cluster_kmeans <- factor(pbmc_kmeans$cluster)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "cluster_kmeans")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/kmeans-1.png)<!-- -->

### Graph based clustering

The clustering algorithm of Seurat itself is based on graph based
clustering. The output of the clustering, will be saved automatically in
the metadata as ‘seurat_clusters’. As explained in the lecture, the
resolution parameter is related to the number of clusters. You can play
around with this parameters to see how it influences the results.

``` r
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.25, verbose = FALSE)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "seurat_clusters")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/graph_clust-1.png)<!-- -->

## Visualizing marker genes and annotating the cells

Once, you are satisfied with the clusters, these can be annotated by
visualizing known marker genes or by looking at differentially expressed
genes. In a later practical, you will learn how to select these, for now
we will just focus on known marker genes. A commonly used approach is
that the data is annotated in a hierarchical fashion. First the data is
annotated at a low resolution (e.g. only 2-3 cell types) and afterwards
each cluster is subsetted from the data, clustered and annotated again.
This process can continue until you’re satisfied with the resolution.

``` r
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.01, verbose = FALSE)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "seurat_clusters")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/graph_clust_lowres-1.png)<!-- -->

So now that we have clustered the data at a low resolution, we can
visualize some marker genes: CD19 (B cells), CD3D (T cells), CD14
(Monocytes), NKG7 (NK cells).

``` r
FeaturePlot(pbmc, reduction='tsne', features=c('CD19', 'CD3D', 'CD14', 'NKG7'))
```

![](Clustering_files/figure-gfm/featplot-1.png)<!-- -->

For a new, more complex dataset, you will probably need to visualize
more genes before you can label a cluster. For now, we will assume that
cluster 0 are NK and T cells, cluster 1 are Monocytes and cluster 2 are
B cells. In the code below, you will assign these labels to your
cluster.

``` r
new.cluster.ids <- c("NK and T cells", "Monocytes", "B cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE) + NoLegend()
```

![](Clustering_files/figure-gfm/clusnames-1.png)<!-- -->

If you want to cluster the cells at a higher resolution, you could for
instance subset the data now and repeat these steps. For now, we will
just save the object for the next practicals.

``` r
saveRDS(pbmc, file = "pbmc3k.rds")
```

### Session info

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] sp_1.5-0           SeuratObject_4.1.2 Seurat_4.2.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6         
    ##   [4] ellipsis_0.3.2        ggridges_0.5.4        rstudioapi_0.14      
    ##   [7] spatstat.data_2.2-0   farver_2.1.1          leiden_0.4.3         
    ##  [10] listenv_0.8.0         ggrepel_0.9.1         fansi_1.0.3          
    ##  [13] codetools_0.2-18      splines_4.2.1         knitr_1.40           
    ##  [16] polyclip_1.10-0       jsonlite_1.8.2        ica_1.0-3            
    ##  [19] cluster_2.1.3         png_0.1-7             rgeos_0.5-9          
    ##  [22] uwot_0.1.14           shiny_1.7.2           sctransform_0.3.5    
    ##  [25] spatstat.sparse_2.1-1 compiler_4.2.1        httr_1.4.4           
    ##  [28] assertthat_0.2.1      Matrix_1.5-0          fastmap_1.1.0        
    ##  [31] lazyeval_0.2.2        cli_3.4.1             later_1.3.0          
    ##  [34] htmltools_0.5.3       tools_4.2.1           igraph_1.3.5         
    ##  [37] gtable_0.3.1          glue_1.6.2            RANN_2.6.1           
    ##  [40] reshape2_1.4.4        dplyr_1.0.10          Rcpp_1.0.9           
    ##  [43] scattermore_0.8       vctrs_0.4.2           nlme_3.1-157         
    ##  [46] progressr_0.11.0      lmtest_0.9-40         spatstat.random_2.2-0
    ##  [49] xfun_0.33             stringr_1.4.1         globals_0.16.1       
    ##  [52] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.2      
    ##  [55] irlba_2.3.5.1         goftest_1.2-3         future_1.28.0        
    ##  [58] MASS_7.3-57           zoo_1.8-11            scales_1.2.1         
    ##  [61] spatstat.core_2.4-4   promises_1.2.0.1      spatstat.utils_2.3-1 
    ##  [64] parallel_4.2.1        RColorBrewer_1.1-3    yaml_2.3.5           
    ##  [67] reticulate_1.26       pbapply_1.5-0         gridExtra_2.3        
    ##  [70] ggplot2_3.3.6         rpart_4.1.16          stringi_1.7.8        
    ##  [73] highr_0.9             rlang_1.0.6           pkgconfig_2.0.3      
    ##  [76] matrixStats_0.62.0    evaluate_0.16         lattice_0.20-45      
    ##  [79] ROCR_1.0-11           purrr_0.3.5           tensor_1.5           
    ##  [82] labeling_0.4.2        patchwork_1.1.2       htmlwidgets_1.5.4    
    ##  [85] cowplot_1.1.1         tidyselect_1.1.2      parallelly_1.32.1    
    ##  [88] RcppAnnoy_0.0.19      plyr_1.8.7            magrittr_2.0.3       
    ##  [91] R6_2.5.1              generics_0.1.3        DBI_1.1.3            
    ##  [94] mgcv_1.8-40           pillar_1.8.1          fitdistrplus_1.1-8   
    ##  [97] survival_3.3-1        abind_1.4-5           tibble_3.1.8         
    ## [100] future.apply_1.9.1    KernSmooth_2.23-20    utf8_1.2.2           
    ## [103] spatstat.geom_2.4-0   plotly_4.10.0         rmarkdown_2.16       
    ## [106] grid_4.2.1            data.table_1.14.2     digest_0.6.29        
    ## [109] xtable_1.8-4          tidyr_1.2.1           httpuv_1.6.6         
    ## [112] munsell_0.5.0         viridisLite_0.4.1
