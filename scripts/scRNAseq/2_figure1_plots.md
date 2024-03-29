2_figure1_plots
================
2022-07-10

## Overview

Figure 1: Single cell transcriptional landscape of CD4+ and CD8+ T cells
during spontaneous and anti-PD-1-induced T1D in NOD mice. The purpose of
this figure is to show the single cell transcriptional landscape of CD4+
and CD8+ T cells during spontaneous and anti-PD-1-induced T1D in NOD
mice across the tissues examined: pancreatic lymph node (pLN), pancreas,
and blood.

``` r
require(dplyr, quietly = T) #simple manipulation of data.frames
require(Seurat, quietly = T) #scRNA-seq analysis
require(reshape, quietly = T) #melting and casting data.frames
require(ggplot2, quietly = T) #plotting
require(rstatix, quietly = T) #adding p values to plots
require(ggprism, quietly = T) #adding p values to plots
require(ggsci, quietly = T) #IGV color palette
require(patchwork, quietly = T) #adding p values to plots
require(magrittr, quietly = T) #adding p values to plots
require(gridExtra, quietly = T) #display multiple plots
require(ggrepel, quietly = T) #labels for plots
require(tibble, quietly = T) #rownames_to_column function
require(useful, quietly = T) #corner function to peek at data.frames
require(viridis, quietly = T) #viridis color palettes
require(readxl, quietly = T) #reading Excel files
require(ggseqlogo, quietly = T) #for TCR amino acid analysis of NRP-v7+ T cells 
source("scripts/utils.R") #custom color palettes and other functions
set.seed(1) #set a seed for reproducibility
```

Read in the integrated seurat object containing CD4+ and CD8+ T cells
from all tissues:

``` r
seurat_integrated <- readRDS("objects/seurat_integrated_clustered.rds")
pan8 <- readRDS("objects/seurat_cd8_panonly_clustered.rds")
pan4 <- readRDS("objects/seurat_cd4_panonly_clustered.rds")
```

We examined the incidence of acinar inflammation in the three treatment
groups. Here we perform a Fisher’s Exact test for differences in acinar
inflammation by creating three contingency tables:

``` r
igg <- data.frame(acinar = c(1,6), not_acinar = c(19,39))
rownames(igg) <- c("IgG","not_IgG")

spt <- data.frame(acinar = c(0,7), not_acinar = c(20,38))
rownames(spt) <- c("Spt","not_Spt")

pd1 <- data.frame(acinar = c(6,1), not_acinar = c(19,39))
rownames(pd1) <- c("PD1","not_PD1")

res <- data.frame(treatment = c("IgG", "Spt", "PD1"), 
            fisher.p.value = p.adjust(c(fisher.test(igg)$p.value, fisher.test(spt)$p.value, fisher.test(pd1)$p.value)),
            chisq.p.value = p.adjust(c(chisq.test(igg)$p.value, chisq.test(spt)$p.value, chisq.test(pd1)$p.value)))
```

    ## Warning in chisq.test(igg): Chi-squared approximation may be incorrect

    ## Warning in chisq.test(spt): Chi-squared approximation may be incorrect

    ## Warning in chisq.test(pd1): Chi-squared approximation may be incorrect

``` r
res
```

    ##   treatment fisher.p.value chisq.p.value
    ## 1       IgG     0.42258839     0.5708221
    ## 2       Spt     0.17969304     0.3032732
    ## 3       PD1     0.03259754     0.0628019

All of the T cells (both CD4 and CD8) will be clustered on a single
UMAP. Additional UMAPs showing distribution of tissues, treatment
groups, cell types, and clonotype expansion are also plotted to show the
distribution of the samples. To label the initial UMAP with the Seurat
clusters, the centroids of each cluster must be measured.

``` r
#define the centroids of each cluster for labeling
integrated_cls <- seurat_integrated@meta.data %>% 
                    dplyr::select(allcells_UMAP_1, 
                                  allcells_UMAP_2, 
                                  seurat_clusters) %>% 
                    group_by(seurat_clusters) %>% 
                    summarise(U1 = mean(allcells_UMAP_1),
                              U2 = mean(allcells_UMAP_2)) %>% 
                    data.frame()
integrated_cls
```

    ##        seurat_clusters        U1         U2
    ## 1               CD8_cm -1.712744  4.8496332
    ## 2              CD8_mem -3.776977  2.5400027
    ## 3           CD8_effmem -3.636370  1.1572626
    ## 4             CD8_slec -7.657589 -0.5356249
    ## 5              CD8_eff -5.382211  1.1655254
    ## 6             CD8_pexh -5.144333  1.6715835
    ## 7             CD8_texh -5.275911 -0.3912251
    ## 8              Tcon_cm  3.596482  2.8960248
    ## 9             Tcon_mem  1.679384  1.1377887
    ## 10         Tcon_effmem -1.034969  0.8894641
    ## 11            Tcon_eff -1.670464 -0.8911204
    ## 12           Tcon_prog  1.004587 -1.2863303
    ## 13           Tcon_Th21  2.554923 -1.7345814
    ## 14            Tcon_Tfh  1.605868 -2.7830132
    ## 15        Treg_resting  4.592014 -2.2855875
    ## 16            Treg_eff  1.732288 -4.2231346
    ## 17 Acinar_contaminated -1.247957  0.5114292
    ## 18  Recently_activated  4.733982 -0.9761269
    ## 19  Interferon_sensing -2.486465 -2.4838935
    ## 20       Proliferating -3.620606 -3.9032418

``` r
#plot the UMAP from the integrated Seurat containing the cluster labels
p_allcells_cluster <- seurat_integrated@meta.data %>% 
                        ggplot(aes(x = allcells_UMAP_1, 
                                   y = allcells_UMAP_2, 
                                   color = seurat_clusters)) + 
                        geom_point(shape = 16, 
                                   alpha = 0.3) +
                        theme_minimal() +
                        theme(axis.title = element_blank(),
                              axis.text = element_blank(),
                              text = element_text(size = 15),
                              panel.grid = element_blank()) +
                        labs(color= NULL,
                             title = NULL) +
                        geom_label_repel(data = integrated_cls, #these centroids will change depending on Seurat used
                                         mapping = aes(x = U1,
                                                       y = U2, 
                                                       label = seurat_clusters), 
                                         size = 5, 
                                         min.segment.length = Inf, 
                                         show.legend = FALSE) +
                        guides(color = guide_legend(override.aes = list(alpha = 1, #custom legend key
                                                                      size = 5, 
                                                                      shape = 19))) + 
                        scale_color_manual(values = custom_colors("clusters"))

#plot the UMAP of the different treatment groups
p_treatment <- seurat_integrated@meta.data %>% 
                ggplot(aes(x = allcells_UMAP_1, 
                           y = allcells_UMAP_2, 
                           color = treatment)) + 
                geom_point(shape = 16, 
                           alpha = 0.3,
                           size = 0.5) +
                theme_minimal() +
                theme(axis.title = element_blank(),
                      axis.text = element_blank(),
                      panel.grid = element_blank(),
                      text = element_text(size = 15),
                      legend.position = "none") +
                labs(color= NULL,
                     title = NULL) +
                guides(color = guide_legend(override.aes = list(alpha = 1, 
                                                              size = 5, 
                                                              shape = 19))) + 
                scale_color_manual(values = custom_colors("treatments"))

#plot the UMAP of the different tissues
p_tissues <- seurat_integrated@meta.data %>% 
              ggplot(aes(x = allcells_UMAP_1, 
                         y = allcells_UMAP_2, 
                         color = tissue)) + 
              geom_point(shape = 16, 
                         alpha = 0.3,
                           size = 0.5) +
              theme_minimal() +
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    panel.grid = element_blank(),
                    text = element_text(size = 15),
                    legend.position = "none") +
              labs(color= NULL,
                   title = NULL) +
              guides(color = guide_legend(override.aes = list(alpha = 1, 
                                                            size = 5, 
                                                            shape = 19))) + 
              scale_color_manual(values = custom_colors("tissues"))

#plot the UMAP of NRP-v7+ CD8 T cells
p_tetramer <- seurat_integrated@meta.data %>% 
                ggplot(aes(x = allcells_UMAP_1, 
                           y = allcells_UMAP_2, 
                           color = tetramer_specificity)) + 
                geom_point(shape = 16, 
                           alpha = 0.2) +
                theme_minimal() +
                theme(axis.title = element_blank(),
                      axis.text = element_blank(),
                      panel.grid = element_blank(),
                      text = element_text(size = 15),
                      legend.position = "none") +
                labs(color = NULL,
                     title = "All T cells: NRP-v7 tetramer positive") +
                guides(color = guide_legend(override.aes = list(alpha = 1, 
                                                              size = 5, 
                                                              shape = 19))) + 
                scale_color_manual(values = custom_colors("tetramer"))

#plot the UMAP of the different celltypes
p_celltypes <- seurat_integrated@meta.data %>% 
                ggplot(aes(x = allcells_UMAP_1, 
                           y = allcells_UMAP_2, 
                           color = celltype)) + 
                geom_point(shape = 16, 
                           alpha = 0.3,
                           size = 0.5) +
                theme_minimal() +
                theme(axis.title = element_blank(),
                      axis.text = element_blank(),
                      panel.grid = element_blank(),
                      text = element_text(size = 15),
                      legend.position = "none") +
                labs(color= NULL,
                     title = NULL) +
                guides(color = guide_legend(override.aes = list(alpha = 1, 
                                                              size = 5, 
                                                              shape = 19))) + 
                scale_color_manual(values = custom_colors("celltypes"))

#plot the UMAP of clonotype expansion
p_expansion <- seurat_integrated@meta.data %>% 
                ggplot(aes(x = allcells_UMAP_1, 
                           y = allcells_UMAP_2, 
                           color = log10(clonotype_count))) + 
                geom_point(shape = 16, 
                           alpha = 0.5) +
                theme_minimal() +
                theme(axis.title = element_blank(),
                      axis.text = element_blank(),
                      panel.grid = element_blank(),
                      text = element_text(size = 15),
                      legend.position = "none",
                      legend.key.height = unit(0.25,"cm")) +
                labs(color= NULL,
                     title = NULL) +
                scale_color_gradient(na.value ="grey90", low = "black", high = "cyan", breaks = c(0.5,2)) +
                guides(color = guide_colorbar(ticks = F)) #removes ticks from colorbar

p_coinh <- seurat_integrated@meta.data %>% 
            ggplot(aes(x = allcells_UMAP_1, 
                       y = allcells_UMAP_2, 
                       color = coinh_count,
                       alpha = coinh_count)) + 
            geom_point(shape = 16) +
            theme_minimal() +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  panel.grid = element_blank(),
                  legend.position = "none",
                  text = element_text(size = 15)) +
            labs(color= NULL,
                 title = NULL) +
            guides(color = guide_legend(override.aes = list(alpha = 1, 
                                                          size = 5, 
                                                          shape = 19))) +
            scale_color_gradient(low = "bisque", high = "darkslategray")

p_prog <- seurat_integrated@meta.data %>%
            ggplot(aes(x = allcells_UMAP_1,
                       y = allcells_UMAP_2,
                       color = Miller2019_progenitor_exh1)) +
            geom_point(shape = 16) +
            theme_minimal() +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  panel.grid = element_blank(),
                  text = element_text(size = 15),
                  legend.position = "none",
                  legend.key.height = unit(0.25,"cm")) +
            labs(color= NULL,
                 title = NULL) +
            scale_color_gradient(low = "cornsilk", high = "#CE3D32FF")

p_term <- seurat_integrated@meta.data %>% 
            ggplot(aes(x = allcells_UMAP_1, 
                       y = allcells_UMAP_2, 
                       color = Miller2019_terminally_exh1)) + 
            geom_point(shape = 16) +
            theme_minimal() +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  panel.grid = element_blank(),
                  text = element_text(size = 15),
                  legend.position = "none",
                  legend.key.height = unit(0.25,"cm")) +
            labs(color= NULL,
                 title = NULL) +
            scale_color_gradient(low = "cornsilk", high = "#1A0099FF") +
            guides(color = guide_colorbar(ticks = F)) #removes ticks from colorbar

#show TOX
data.frame(t(as.matrix(seurat_integrated@assays$RNA["Tox",]))) -> new_metadata
seurat_integrated <- AddMetaData(seurat_integrated, new_metadata)

p_tox <- seurat_integrated@meta.data %>% 
    ggplot(aes(x = allcells_UMAP_1, 
               y = allcells_UMAP_2, 
               color = Tox,
               alpha = Tox)) + 
    geom_point(shape = 16) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 15),
          legend.position = "none",
          legend.key.height = unit(0.25,"cm")) +
    labs(color= NULL,
         title = NULL) +
    scale_color_gradient(low = "papayawhip", high = "sienna4") +
    guides(color = guide_colorbar(ticks = F)) #removes ticks from colorbar


grid.arrange(p_celltypes, p_tox, p_coinh, p_term, p_expansion, p_tissues, p_treatment,   ncol = 2)
```

![](2_figure1_plots_files/figure-gfm/Improved%20cluster%20UMAP-1.png)<!-- -->

Plot the UMAP visualization of some critical markers of T cells used to
identify cell subsets (CD4, CD8a, FoxP3):

``` r
data.frame(t(as.matrix(seurat_integrated@assays$RNA[c("Cd4","Cd8a","Foxp3"),]))) -> new_metadata
seurat_integrated <- AddMetaData(seurat_integrated, new_metadata)

p_cd4 <- seurat_integrated@meta.data %>% 
          ggplot(aes(x = allcells_UMAP_1, 
                     y = allcells_UMAP_2, 
                     color = Cd4)) + 
          geom_point(shape = 16, 
                     alpha = 0.6) +
          theme_minimal() +
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                panel.grid = element_blank(),
                text = element_text(size = 15),
                legend.position = c(0.2, 0.1),
                legend.key.height = unit(0.25,"cm")) +
          labs(color= NULL,
               title = NULL) +
          scale_color_gradient(low = "cornsilk", high = "#011f4b") +
          guides(color = guide_colorbar(ticks = F)) #removes ticks from colorbar

p_cd8a <- seurat_integrated@meta.data %>% 
            ggplot(aes(x = allcells_UMAP_1, 
                       y = allcells_UMAP_2, 
                       color = Cd8a)) + 
            geom_point(shape = 16, 
                       alpha = 0.6) +
            theme_minimal() +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  panel.grid = element_blank(),
                  text = element_text(size = 15),
                  legend.position = c(0.2, 0.1),
                  legend.key.height = unit(0.25,"cm")) +
            labs(color= NULL,
                 title = NULL) +
            scale_color_gradient(low = "cornsilk", high = "deeppink2") +
            guides(color = guide_colorbar(ticks = F)) #removes ticks from colorbar

p_foxp3 <- seurat_integrated@meta.data %>% 
            ggplot(aes(x = allcells_UMAP_1, 
                       y = allcells_UMAP_2, 
                       color = Foxp3)) + 
            geom_point(shape = 16, 
                       alpha = 0.6) +
            theme_minimal() +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  panel.grid = element_blank(),
                  text = element_text(size = 15),
                  legend.position = c(0.2, 0.1),
                  legend.key.height = unit(0.25,"cm")) +
            labs(color= NULL,
                 title = NULL) +
            scale_color_gradient(low = "cornsilk", high = "skyblue") +
            guides(color = guide_colorbar(ticks = F)) #removes ticks from colorbar


grid.arrange(p_cd4, p_cd8a, p_foxp3, ncol = 3)
```

![](2_figure1_plots_files/figure-gfm/UMAP%20of%20markers-1.png)<!-- -->

Next, create a stacked bar chart showing the distribution of tissues in
each Seurat cluster:

``` r
p_tissue_per_seurat_cluster <- seurat_integrated@meta.data %>% 
                                ggplot(aes(x = factor(seurat_clusters, 
                                                    levels = rev(levels(seurat_integrated$seurat_clusters))), #orders clusters properly
                                                    fill = tissue)) + 
                                geom_bar(stat = "count", 
                                         position = "fill") + 
                                scale_x_discrete(expand = c(0, 0)) +
                                scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                                   expand = c(0, 0)) +
                                theme_minimal() +
                                theme(panel.grid = element_blank(),
                                      text = element_text(size = 15),
                                      axis.text = element_blank(),
                                      legend.position = "none") +
                                labs(title = NULL, 
                                     y = NULL, 
                                     x = NULL,
                                     color = "Tissue") +
                                scale_fill_manual(values = custom_colors("tissues"))

p_celltype_per_seurat_cluster <- seurat_integrated@meta.data %>% 
                                  ggplot(aes(x = factor(seurat_clusters, 
                                                      levels = rev(levels(seurat_integrated$seurat_clusters))), #orders clusters properly
                                                      fill = celltype)) + 
                                  geom_bar(stat = "count", 
                                           position = "fill") + 
                                  scale_x_discrete(expand = c(0, 0)) +
                                  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                                     expand = c(0, 0)) +
                                  theme_minimal() +
                                  theme(panel.grid = element_blank(),
                                        text = element_text(size = 15),
                                        axis.text = element_blank(),
                                        legend.position = "none") +
                                  labs(title = NULL, 
                                       y = NULL, 
                                       x = NULL,
                                       color = "Tissue") +
                                  scale_fill_manual(values = custom_colors("celltypes"))

p_treatment_per_seurat_cluster <- seurat_integrated@meta.data %>% 
                                    ggplot(aes(x = factor(seurat_clusters, 
                                                        levels = rev(levels(seurat_integrated$seurat_clusters))), #orders clusters properly
                                                        fill = treatment)) + 
                                    geom_bar(stat = "count", 
                                             position = "fill") + 
                                    scale_x_discrete(expand = c(0, 0)) +
                                    scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                                       expand = c(0, 0)) +
                                    theme_minimal() +
                                    theme(panel.grid = element_blank(),
                                          text = element_text(size = 15),
                                          axis.text = element_blank(),
                                          legend.position = "none") +
                                    labs(title = NULL, 
                                         y = NULL, 
                                         x = NULL,
                                         color = "Tissue") +
                                    scale_fill_manual(values = custom_colors("treatments"))

grid.arrange(p_tissue_per_seurat_cluster, p_celltype_per_seurat_cluster, p_treatment_per_seurat_cluster, nrow = 3)
```

![](2_figure1_plots_files/figure-gfm/stacked%20bar%20chart%20of%20tissues%20and%20cell%20type%20per%20seurat%20cluster-1.png)<!-- -->

Plot a bar chart of NRP-v7+ T cells per cluster

``` r
p_tetramer_per_cluster <- seurat_integrated@meta.data %>%
                            filter(!grepl("Tcon|Treg", seurat_clusters)) %>%
                            mutate(tetramer_tissue = case_when(tetramer_specificity == "NRPv7" & tissue == "pancreas" ~ "pancreas NRPv7",
                                     tetramer_specificity == "NRPv7" & tissue != "pancreas" ~ "periphery NRPv7",
                                     TRUE ~ paste0(tetramer_specificity))) %>%
                            ggplot(aes(x = factor(seurat_clusters, 
                                                levels = rev(levels(seurat_integrated$seurat_clusters))), #orders clusters properly
                                                fill = tetramer_tissue)) + 
                            geom_bar(stat = "count", 
                                     position = "fill") + 
                            scale_x_discrete(expand = c(0, 0)) +
                            scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                               expand = c(0, 0)) +
                            theme_minimal() +
                            theme(panel.grid = element_blank(),
                                  text = element_text(size = 15),
                                  legend.position = "right") +
                            labs(title = "CD8 T cells", 
                                 y = "Percent of Seurat Cluster", 
                                 x = NULL,
                                 fill = "Tetramer\nspecificity") +
                            scale_fill_manual(values = c("unknown" = "grey90", "pancreas NRPv7" = "Firebrick", "periphery NRPv7" = "hotpink")) +
                            coord_flip()
p_tetramer_per_cluster
```

![](2_figure1_plots_files/figure-gfm/plot%20NRP-v7%20staining-1.png)<!-- -->

``` r
## TOP 50 CLONOTYPES
## needs to be fixe to factor each group separately
cumulative_sums <- pan8@meta.data %>% filter(!is.na(clonotype_id)) %>% 
                      group_by(clonotype_id, tetramer_specificity) %>% 
                      tally() %>% 
                      group_by(tetramer_specificity) %>% 
                      mutate(prcnt = n*100/sum(n)) %>% 
                      group_by(tetramer_specificity) %>% 
                      arrange(tetramer_specificity, desc(prcnt)) %>% 
                      mutate(cs = cumsum(prcnt))

cumulative_sums$clonotype_id <- factor(cumulative_sums$clonotype_id, 
                                       levels = unique(cumulative_sums$clonotype_id))

p_top50_clonos <- cumulative_sums %>% 
                    group_by(tetramer_specificity) %>% 
                    slice_head(n = 50) %>% 
                    ggplot(aes(x = clonotype_id, 
                               y = cs, 
                               color = tetramer_specificity)) + 
                    geom_point() + 
                    geom_segment(aes(xend = clonotype_id, 
                                     yend = cs, 
                                     x = clonotype_id, 
                                     y = 0, 
                                     color = tetramer_specificity)) + 
                    scale_color_manual(values = custom_colors("tetramer")) +
                    scale_y_continuous(expand = c(0, 0)) +
                    theme_minimal() +
                    theme(panel.grid = element_blank(),
                          panel.spacing.x = unit(0, "lines"),
                          text = element_text(size = 15),
                          axis.text.x = element_blank(),
                          panel.border = element_rect(colour = "black", 
                                                      fill = NA, 
                                                      size = 1),
                          legend.position = "none") +
                    labs(x = "Top 50 clonotypes", 
                         y = "Cumulative % CD8 TCRs",
                         title = "CD8+ T cells in pancreas") +
                    facet_wrap(~ tetramer_specificity, scales = "free_x")
p_top50_clonos
```

![](2_figure1_plots_files/figure-gfm/plot%20NRP-v7%20staining-2.png)<!-- -->

To plot a heatmap of genes of interest across the clusters in Morpheus,
the following manually annotated list of genes of interest, **goi**,
will be examined. The average expression for each cluster will be
determined and this data will be used to plot in Morpheus (Broad
institute). The average expression for every gene will be determined to
plot the pearson correlations for the Seurat clusters across tissues as
well.

``` r
## average gene expression of each tissue for each cluster/tissue - for Pearson correlation
avgs <- as.data.frame(t(AverageExpression(seurat_integrated, 
                                          assay = "RNA", 
                                          group.by = c("tissue", "seurat_clusters"))$RNA)) #get average expression for each gene (column) in each tissue/cluster pair (row)
avgs <- avgs %>% rownames_to_column("tissue_and_cluster") #the generated rownames include the tissue and cluster
avgs <- avgs %>% 
          mutate(tissue = sapply(tissue_and_cluster, function(x) strsplit(x, "_")[[1]][1]), #parse the tissue out to a separate column
                 seurat_clusters = sapply(tissue_and_cluster, function(x) paste(strsplit(x, "_")[[1]][-1], collapse = "_"))) %>%  #parse the cluster name to a separate column
          relocate(tissue_and_cluster, tissue, seurat_clusters) #move columns to left of data frame
#write.table(avgs, "output/seurat_cluster_avg_expression_all_genes_cluster_tissue_for_morpheus.tsv", sep = "\t", quote = F, row.names = F) #write the data to table

## average expression for each cluster for heatmap
gene_annot <- read.table("output/morpheus/figure1/gene_cluster_heatmap/genes_for_heatmap.tsv", sep = "\t", header = T) #custom curated list of annotated genes of interest
avgs <- as.data.frame(AverageExpression(seurat_integrated, 
                                        features = gene_annot$gene, 
                                        assay = "RNA", 
                                        group.by = "seurat_clusters")$RNA) #get average expression for each gene of interest
avgs <- avgs %>% rownames_to_column("gene")
avgs <- merge(gene_annot, avgs, by = "gene")
#write.table(avgs, "output/morpheus/figure1/gene_cluster_heatmap/seurat_cluster_avg_expression_curated_genes_for_morpheus.tsv", sep = "\t", quote = F, row.names = F) #write the data to table
#write.table(seurat_integrated@meta.data %>% group_by(seurat_clusters) %>% summarize(coinh_count = mean(coinh_count)), "output/morpheus/figure1/gene_cluster_heatmap/seurat_cluster_avg_coinh_count_for_morpheus.tsv", row.names = F, quote = F, sep = "\t")

#get the markers that can be used to define each cluster
# seurat_markers <- FindAllMarkers(seurat_integrated, 
#                                  min.pct = 0.25, 
#                                  only.pos = T, 
#                                  logfc.threshold = 0.25) %>%
#                      group_by(cluster) %>%
#                      filter(p_val_adj < 0.05) %>%
#                      group_by(cluster) %>%
#                      arrange(desc(avg_log2FC), .by_group = T) %>%
#                      top_n(n = 50, 
#                            wt = avg_log2FC)
#saveRDS(seurat_markers, "objects/integrated_seurat_FindAllMarkers.rds")

# avgs <- as.data.frame(AverageExpression(seurat_integrated, 
#                                         assay = "RNA", 
#                                         group.by = "seurat_clusters")$RNA) #get average expression for each gene (column) in each cluster (row)
# avgs <- avgs %>% rownames_to_column("gene")

## write a table containing the top 20 DE genes for each cluster
# top20_genes <- seurat_markers %>% group_by(cluster) %>%
#                       arrange(desc(avg_log2FC), .by_group = T) %>%
#                       top_n(n = 20, wt = avg_log2FC) %>%
#                       dplyr::select(gene, cluster)

#write.table(top20_genes, "output/integrated_seurat_FindAllMarkers_top20.tsv", quote = F, sep = "\t", row.names = F)
```
