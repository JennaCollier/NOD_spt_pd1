geomx_analysis
================
2022-06-29

## GeoMX analysis

In this analysis, we looked at various proteins within beta islets of
NOD mice age 16-22 weeks of age. The three treatment groups examined
include healthy **isotype antibody-treated (IgG)**; **spontaneously
diabetic (Spt)**; and **PD-1 blockade-treated diabetic (PD1)** mice.
Both IgG and PD1 mice were treated with two doses of 200ug of antibody
(anti-PD-1 clone 29F.1A12 or isotype IgG2a control) on days 0,2 before
sacrifice on day 4. Pancreata were excised and immediately preserved in
10% formalin for 24h, stored in 70% ethanol, embedded in paraffin, and
sectioned for H&E staining or staining with anti-INS or anti-CD45
antibodies. Both CD45 and INS staining was used to differentiate various
regions within the islets. Islets were scored as healthy with no immune
infiltrate (score 0), per-insulitis (score 1-2), or overt insulitis
(score 3-4).

``` r
require(dplyr, quietly = T) #simple manipulation of data.frames
require(reshape, quietly = T) #melting and casting data.frames
require(ggplot2, quietly = T) #plotting
require(readxl, quietly = T) #reading excel spreadsheets
require(rstatix, quietly = T) #to add significance to p values
require(ggprism, quietly = T) #adding p values to plots
require(patchwork, quietly = T) #adding p values to plots
require(magrittr, quietly = T) #adding p values to plots
require(gridExtra, quietly = T) #display multiple plots
require(vegan , quietly = T) #principal component analysis
require(nlme, quietly = T) #mixed linear model analysis
require(ggrepel, quietly = T) #labels for pca plots
require(tibble, quietly = T) #rownames_to_column function
```

``` r
treatment_colors <- c("IgG" = "black", 
                      "Spt" = rgb(253,128,8, maxColorValue = 255), 
                      "PD1" = rgb(128,0,128, maxColorValue = 255))
```

Read in the data files. One data file contains the code-friendly protein
names, description, and associated mouse genes. Two data files for CD45+
and INS+ regions of interest describing the normalized (using the GeoMX
software) quantification of protein expression. The data for the CD45
and ins regions are merged, as they are already distinguished based on
staining for CD45 and INS.

``` r
#read in CD45 data that has been normalized
cd45 <- as_tibble(t(read_excel(file.path("data/geoMX/norm_geoMean_CD45_trimmed.xlsx"))))
```

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
    ## Using compatibility `.name_repair`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

``` r
colnames(cd45) <- cd45[1,] #set column names manually (doesn't work in read_excel)
cd45 <- cd45[-1,] #remove row with column names
cd45$islet_score <- gsub("\\-", ".", cd45$islet_score) #change dashes in Peri-insulitis to periods (Peri.insulitis) in data
colnames(cd45) <- gsub("\\-", ".", colnames(cd45)) #change dashes in protein colnames in data

#read in CD45 data that has been normalized
ins <- as_tibble(t(read_excel(file.path("data/geoMX/norm_geoMean_INS_trimmed.xlsx"))))
colnames(ins) <- ins[1,] #set column names manually (doesn't work in read_excel)
ins <- ins[-1,] #remove row with column names
ins$islet_score <- gsub("\\-", ".", ins$islet_score) #change dashes in Peri-insulitis to periods (Peri.insulitis) in data
colnames(ins) <- gsub("\\-", ".", colnames(ins)) #change dashes in protein colnames in data

#merge data
geomx <- rbind(cd45,ins)
geomx <- geomx[,-c(1:8)] #remove excess column information
head(geomx)
```

    ## # A tibble: 6 x 74
    ##   treatment mouse_id islet_score label PD.L1   CD45  Rb_IgG CTLA4 Rt_IgG2b GZMB 
    ##   <chr>     <chr>    <chr>       <chr> <chr>   <chr> <chr>  <chr> <chr>    <chr>
    ## 1 PD1       437      Insulitis   CD45  63.127~ 2433~ 175.1~ 87.9~ 16.3109~ 144.~
    ## 2 Spt       535      Insulitis   CD45  140.98~ 1638~ 36.01~ 156.~ 6.93852~ 726.~
    ## 3 Spt       669      Insulitis   CD45  105.76~ 1073~ 43.27~ 135.~ 2.62767~ 92.9~
    ## 4 Spt       669      Insulitis   CD45  260.89~ 1806~ 162.8~ 229.~ 32.8831~ 193.~
    ## 5 PD1       101      Insulitis   CD45  152.26~ 1099~ 35.58~ 146.~ 10.3706~ 896.~
    ## 6 PD1       440      Insulitis   CD45  96.717~ 7766~ 49.04~ 114.~ 7.81749~ 131.~
    ## # ... with 64 more variables: Histone.H3 <chr>, CD8a <chr>, CD19 <chr>,
    ## #   CD3e <chr>, CD11c <chr>, F480 <chr>, SMA <chr>, CD4 <chr>, PanCk <chr>,
    ## #   CD11b <chr>, GAPDH <chr>, Fibronectin <chr>, Ki.67 <chr>, PD.1 <chr>,
    ## #   S6 <chr>, MHCII <chr>, Rt_IgG2a <chr>, VISTA <chr>, B7.H3 <chr>,
    ## #   Tim.3 <chr>, OX40L <chr>, LAG3 <chr>, GITR <chr>, CD44 <chr>, CD86 <chr>,
    ## #   CD40L <chr>, CD40 <chr>, CD27 <chr>, CD12 <chr>, ICOS <chr>, CD28 <chr>,
    ## #   CD34 <chr>, `Ly6G&Ly6C` <chr>, FOXP3 <chr>, CD14 <chr>, CD163 <chr>, ...

``` r
#read in gene annotations
gene_annot <- read_excel(file.path("data/geoMX/gene_annotations.xlsx"), col_names = T)
gene_annot$protein <- gsub("\\-", ".", gene_annot$protein)
gene_annot$target_group <- gsub("All Targets;","",gene_annot$target_group) #remove extraneous text from description
head(gene_annot)
```

    ## # A tibble: 6 x 6
    ##   target_group                        info  code_class description protein gene 
    ##   <chr>                               <chr> <chr>      <chr>       <chr>   <chr>
    ## 1 Myeloid Activation;Checkpoint       chec~ Endogenous PD-L1       PD.L1   Cd274
    ## 2 Total Immune                        immu~ Endogenous CD45        CD45    Cd45 
    ## 3 Background                          cont~ Negative   Rb IgG      Rb_IgG  <NA> 
    ## 4 T cell Activation;T cells;Th cells~ chec~ Endogenous CTLA4       CTLA4   Ctla4
    ## 5 Background                          cont~ Negative   Rt IgG2b    Rt_IgG~ <NA> 
    ## 6 T cell Activation;Cytotoxicity      T ce~ Endogenous GZMB        GZMB    Gzmb

Tally the number of samples. Note that there is only one normal islet in
the CD45 group that was analyzed.

``` r
geomx %>% group_by(treatment,islet_score,label) %>% tally() %>% cast(label+islet_score ~ treatment, value = "n", fill = 0)
```

    ##   label    islet_score IgG PD1 Spt
    ## 1  CD45      Insulitis  30  16  24
    ## 2  CD45         Normal   1   0   0
    ## 3  CD45 Peri.insulitis  22   4   8
    ## 4   INS      Insulitis  16   0  10
    ## 5   INS         Normal   7   0   0
    ## 6   INS Peri.insulitis  19   0   6

Examine protein expression broadly across all of the samples to confirm
that expression levels generally make sense:

``` r
#melt data to long format for plotting
colnames(geomx)
```

    ##  [1] "treatment"        "mouse_id"         "islet_score"      "label"           
    ##  [5] "PD.L1"            "CD45"             "Rb_IgG"           "CTLA4"           
    ##  [9] "Rt_IgG2b"         "GZMB"             "Histone.H3"       "CD8a"            
    ## [13] "CD19"             "CD3e"             "CD11c"            "F480"            
    ## [17] "SMA"              "CD4"              "PanCk"            "CD11b"           
    ## [21] "GAPDH"            "Fibronectin"      "Ki.67"            "PD.1"            
    ## [25] "S6"               "MHCII"            "Rt_IgG2a"         "VISTA"           
    ## [29] "B7.H3"            "Tim.3"            "OX40L"            "LAG3"            
    ## [33] "GITR"             "CD44"             "CD86"             "CD40L"           
    ## [37] "CD40"             "CD27"             "CD12"             "ICOS"            
    ## [41] "CD28"             "CD34"             "Ly6G&Ly6C"        "FOXP3"           
    ## [45] "CD14"             "CD163"            "CD31"             "BatF3"           
    ## [49] "S100B"            "Epcam"            "ER"               "Pmel17"          
    ## [53] "Her2"             "AR"               "AhR"              "IFNGR"           
    ## [57] "GFP"              "BIM"              "p53"              "BCLXL"           
    ## [61] "PARP"             "p21"              "Perforin"         "Cleaved_Caspase3"
    ## [65] "BAD"              "gamma.H2AX"       "pS6"              "pGSK3A"          
    ## [69] "pPRAS40"          "pAMPK.alpha"      "pAKT1"            "Pan.AKT"         
    ## [73] "PLCG1"            "MET"

``` r
melted <- data.frame(geomx) %>% melt(id.vars = c("treatment","mouse_id","islet_score","label"))
colnames(melted)[which(colnames(melted) == "variable")] <- "protein"
colnames(melted)[which(colnames(melted) == "value")] <- "expression"
head(melted)
```

    ##   treatment mouse_id islet_score label protein         expression
    ## 1       PD1      437   Insulitis  CD45   PD.L1 63.127837555183959
    ## 2       Spt      535   Insulitis  CD45   PD.L1 140.98315394086532
    ## 3       Spt      669   Insulitis  CD45   PD.L1 105.76616280048388
    ## 4       Spt      669   Insulitis  CD45   PD.L1 260.89662556491243
    ## 5       PD1      101   Insulitis  CD45   PD.L1 152.26059349306291
    ## 6       PD1      440   Insulitis  CD45   PD.L1 96.717817435207635

``` r
melted$expression <- as.numeric(melted$expression)
melted$protein <- gsub("\\-", ".", melted$protein) #replace periods "." with dash "-" in protein names eg: PD.L1 to PD-L1
#set factor levels based on average expression
melted$protein <- factor(melted$protein, levels = melted %>% 
                           group_by(protein) %>% 
                           summarize(avg = mean(expression)) %>% 
                           arrange(desc(avg)) %>% 
                           pull(var = protein))
melted$treatment <- factor(melted$treatment, levels = c("IgG","Spt","PD1"))
melted$islet_score <- factor(melted$islet_score)

melted <- merge(melted, gene_annot[,c("protein","info")], by = c("protein")) #add protein annotation info

#calculate the p values for plotting using t tests
p_vals <- melted %>%
  rstatix::group_by(protein) %>%
  rstatix::t_test(expression ~ label) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "protein", dodge = 0.8) %>% # important for positioning!
  filter(p.adj <= 0.05)
p_vals$y.position <- log2(p_vals$y.position) #log2 transformed data will be plotted, so djust pvalue position

#plot the data as violin plots
melted %>% 
  ggplot(aes(x=protein, y=log2(expression), color=label, fill=label)) + 
  geom_violin(alpha=0.3, adjust=2, scale = "width", draw_quantiles = c(0.5), position = position_dodge(0.8)) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank()) +
  labs(title = "GeoMX protein expression") +
  add_pvalue(p_vals, 
             color = "black",
               xmin = "xmin", 
               xmax = "xmax",
               label = "p.adj.signif",
               tip.length = 0,
             show.legend = F,
             inherit.aes = F)
```

<img src="geomx_analysis_files/figure-gfm/protein expression by region of interest-1.png" width="900px" />

Expression of CD45 and other lymphocyte-associated proteins (CD19,
CD11c, F4/80, CD3e, etc) are higher in the CD45 regions which would be
expected. The controls look fairly consistent, so they can be removed to
focus on the proteins of interest from the melted data frame.

look at differences between islets with peri-insulitis and overt
insulitis. The single normal islet was dropped for this analysis.

``` r
p1_vals <- melted %>%
  filter(label == "CD45") %>%
  filter(islet_score != "Normal") %>%
  rstatix::group_by(protein) %>%
  rstatix::t_test(expression ~ islet_score) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "protein", dodge = 0.8) # important for positioning!
p1_vals$y.position <- log2(p1_vals$y.position)
p1_vals
```

    ## # A tibble: 69 x 16
    ##    protein    .y.        group1  group2    n1    n2 statistic    df      p p.adj
    ##    <fct>      <chr>      <chr>   <chr>  <int> <int>     <dbl> <dbl>  <dbl> <dbl>
    ##  1 S6         expression Insuli~ Peri.~    70    34   -2.76    56.1 0.0078 0.179
    ##  2 Histone.H3 expression Insuli~ Peri.~    70    34   -0.0964  48.1 0.924  0.996
    ##  3 CD45       expression Insuli~ Peri.~    70    34    2.14    80.0 0.0354 0.244
    ##  4 SMA        expression Insuli~ Peri.~    70    34    1.66    80.9 0.0998 0.574
    ##  5 Pan.AKT    expression Insuli~ Peri.~    70    34    0.847   67.9 0.4    0.852
    ##  6 GAPDH      expression Insuli~ Peri.~    70    34    1.23    56.5 0.225  0.707
    ##  7 IFNGR      expression Insuli~ Peri.~    70    34    0.614   61.3 0.542  0.92 
    ##  8 BAD        expression Insuli~ Peri.~    70    34   -0.600   64.2 0.551  0.92 
    ##  9 CD44       expression Insuli~ Peri.~    70    34   -0.445   80.2 0.657  0.940
    ## 10 CD4        expression Insuli~ Peri.~    70    34    2.27    87.7 0.0255 0.220
    ## # ... with 59 more rows, and 6 more variables: p.adj.signif <chr>,
    ## #   y.position <dbl>, groups <named list>, x <dbl>, xmin <dbl>, xmax <dbl>

``` r
# p1 <- melted %>%
#     filter(label == "CD45") %>%
#     filter(islet_score != "Normal") %>% 
#     ggplot(aes(x=protein,y=log2(expression),color=islet_score,fill=islet_score)) + 
#     geom_violin(alpha=0.3, adjust=2, scale = "width", draw_quantiles = c(0.5), position = position_dodge(0.8)) + 
#     theme_minimal() + 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.title.x = element_blank()) +
#     labs(title = "GeoMX protein expression (CD45+ regions)") +
#     add_pvalue(p1_vals, 
#                color = "black",
#                xmin = "xmin", 
#                xmax = "xmax",
#                label = "p.adj.signif",
#                tip.length = 0,
#                show.legend = F,
#                inherit.aes = F)

p2_vals <- melted %>%
  filter(label == "INS") %>%
  filter(islet_score != "Normal") %>%
  rstatix::group_by(protein) %>%
  rstatix::t_test(expression ~ islet_score) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "protein", dodge = 0.8) %>% # important for positioning!
  filter(p.adj <= 0.05)
p2_vals$y.position <- log2(p2_vals$y.position)
p2_vals
```

    ## # A tibble: 0 x 16
    ## # ... with 16 variables: protein <fct>, .y. <chr>, group1 <chr>, group2 <chr>,
    ## #   n1 <int>, n2 <int>, statistic <dbl>, df <dbl>, p <dbl>, p.adj <dbl>,
    ## #   p.adj.signif <chr>, y.position <dbl>, groups <named list>, x <dbl>,
    ## #   xmin <dbl>, xmax <dbl>

``` r
# p2 <- melted %>%
#     filter(label == "INS") %>%
#     filter(islet_score != "Normal") %>% 
#     ggplot(aes(x=protein,y=log2(expression),color=islet_score,fill=islet_score)) + 
#     geom_violin(alpha=0.3,adjust=2,scale = "width",draw_quantiles = c(0.5),position = position_dodge(0.8)) + 
#     theme_minimal() + 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.title.x = element_blank()) +
#     labs(title = "GeoMX protein expression (INS+ regions)") +
#     add_pvalue(p2_vals, 
#                color = "black",
#                xmin = "xmin", 
#                xmax = "xmax",
#                label = "p.adj.signif",
#                tip.length = 0,
#                show.legend = F,
#                inherit.aes = F)
# 
# grid.arrange(p1, p2, nrow=2)
```

The two differences in CD45+ regions is that Granzyme B and Epcam is
increased in peri-insulitis. This suggests that T cells are making more
granzyme B in peri-insulitis lesions during initial destruction of the
beta-islets than at later stages of destruction (perhaps they are
becoming exhausted or dysfunctional as beta-islets become increasingly
destroyed). Epcam is also increased - this is a gene involved in
epithelial cell adhesions which suggests that the cell-to-cell adhesions
are reduced as the beta-islet degrades.

``` r
df_for_pca <- geomx %>% filter(label=="CD45")
df_for_pca[,-c(1:4)] <- sapply(df_for_pca[,-c(1:4)],as.numeric)
df_for_pca[,-c(1:4)] <- sapply(df_for_pca[,-c(1:4)],log2)

prin_comp <- rda(df_for_pca[,-c(1:4)])
pca_scores <- scores(prin_comp)

# Define shapes for plot
shapes = c(16, 17, 15) 
shapes <- shapes[as.numeric(factor(df_for_pca$islet_score))]

#plot using base R
plot(pca_scores$sites[,1],
     pca_scores$sites[,2],
     pch = shapes,
     col = paste(treatment_colors[df_for_pca$treatment]),
     xlim = c(-3,3),
     ylim = c(-3,3),
     cex=3)
ordiellipse(prin_comp, factor(df_for_pca$treatment, levels = c("IgG","Spt","PD1")), conf=0.90, col = treatment_colors, lwd = 2)
```

![](geomx_analysis_files/figure-gfm/PCA%20analysis%20CD45-1.png)<!-- -->

``` r
#keep data to later plot in ggplot
df_all <- cbind(df_for_pca, pca_scores$sites)
pca_weights <- data.frame(pca_scores$species, label = "CD45") %>% rownames_to_column("protein")

#collect variance explained for each PC for labeling
variance_explained_labels <- c(paste0("PC1 (", 
                                      round(summary(prin_comp)$cont$importance[2, "PC1"]*100, 2),
                                      "% of variance)"),
                               paste0("PC2 (", 
                                      round(summary(prin_comp)$cont$importance[2, "PC2"]*100, 2),
                                      "% of variance)"))
```

``` r
df_for_pca <- geomx %>% filter(label=="INS")
df_for_pca[,-c(1:4)] <- sapply(df_for_pca[,-c(1:4)],as.numeric)
df_for_pca[,-c(1:4)] <- sapply(df_for_pca[,-c(1:4)],log2)

prin_comp <- rda(df_for_pca[,-c(1:4)])
pca_scores <- scores(prin_comp)

# Define shapes for plot
shapes = c(16, 17, 15) 
shapes <- shapes[as.numeric(factor(df_for_pca$islet_score))]

plot(pca_scores$sites[,1],
     pca_scores$sites[,2],
     pch = shapes,
     col = paste(treatment_colors[df_for_pca$treatment]),
     xlim = c(-4,3),
     ylim = c(-2,4),
     cex=3)
ordiellipse(prin_comp, factor(df_for_pca$treatment, levels = c("IgG","Spt","PD1")), conf=0.90, col = treatment_colors, lwd = 2)
```

![](geomx_analysis_files/figure-gfm/PCA%20analysis%20INS-1.png)<!-- -->

``` r
#gather data for PCA plot in ggplot
df_for_pca <- cbind(df_for_pca, pca_scores$sites)
df_all <- rbind(df_all, df_for_pca)
pca_weights <- rbind(pca_weights, data.frame(pca_scores$species, label = "INS") %>% rownames_to_column("protein"))

#collect variance explained for each PC (as a percentage of 100, rounded to two decimals) for labeling axes in plots
variance_explained_labels <- c(variance_explained_labels,
                               paste0("PC1 (", 
                                      round(summary(prin_comp)$cont$importance[2, "PC1"]*100, 2),
                                      "% of variance)"),
                               paste0("PC2 (", 
                                      round(summary(prin_comp)$cont$importance[2, "PC2"]*100, 2),
                                      "% of variance)"))
names(variance_explained_labels) <- c("PC1_CD45","PC2_CD45","PC1_INS","PC2_INS")
```

These plots can be generated using ggplot2 which looks nicer and allows
for increased customization. Both the PCA plots will be generated
separately to allow for the axes to be labeled with the customized
variance labels (rather than facetting the plots).

``` r
pca_cd45 <- df_all %>%
    filter(label == "CD45") %>%
    ggplot(aes(x = PC1, y= PC2, color = treatment, shape = islet_score, group = treatment)) +
    geom_point(size = 4, alpha = 0.6) +
    stat_ellipse(level = 0.9) + #90% confidence intervals
    scale_color_manual(values = treatment_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 15),
          legend.position = "none") + #this plot will be on the left so it doesn't need a legend
    labs(title = "CD45+ regions of islets",
         x = variance_explained_labels["PC1_CD45"],
         y = variance_explained_labels["PC2_CD45"])
    
pca_ins <- df_all %>%
    filter(label == "INS") %>%
    ggplot(aes(x = PC1, y= PC2, color = treatment, shape = islet_score, group = treatment)) +
    geom_point(size = 4, alpha = 0.6) +
    stat_ellipse(level = 0.9) + #90% confidence intervals
    scale_color_manual(values = treatment_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 15)) +
    labs(shape = "Islet score", 
         color = "Treatment",
         title = "INS+ regions of islets",
         x = variance_explained_labels["PC1_INS"],
         y = variance_explained_labels["PC2_INS"])

grid.arrange(pca_cd45, pca_ins, nrow = 1, widths = c(1.5, 2))
```

![](geomx_analysis_files/figure-gfm/ggplot%20of%20PCA%20plots-1.png)<!-- -->

``` r
weights_cd45 <- pca_weights %>%
    filter(label == "CD45") %>%
    ggplot(aes(x = PC1, y= PC2)) +
    facet_wrap(~ label, scales = "free") +
    scale_color_manual(values = treatment_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 15)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(4, "mm")),
                 color = "black",
                 alpha = 0.5) +
    geom_text_repel(aes(label = protein), 
                    color = "black",
                    max.overlaps = 4) +
    labs(x = variance_explained_labels["PC1_CD45"],
         y = variance_explained_labels["PC2_CD45"])

weights_ins <- pca_weights %>%
    filter(label == "INS") %>%
    ggplot(aes(x = PC1, y= PC2)) +
    facet_wrap(~ label, scales = "free") +
    scale_color_manual(values = treatment_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 15)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(4, "mm")),
                 color = "black",
                 alpha = 0.5) +
    geom_text_repel(aes(label = protein), 
                    color = "black",
                    max.overlaps = 4) +
    labs(x = variance_explained_labels["PC1_INS"],
         y = variance_explained_labels["PC2_INS"])

grid.arrange(weights_cd45, weights_ins, nrow = 1)
```

    ## Warning: ggrepel: 58 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 56 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](geomx_analysis_files/figure-gfm/ggplot%20of%20PCA%20weights-1.png)<!-- -->

Next, the differences between treatment groups will be analyzed by
generating volcano plots. **NOTE** need to fix this to use mixed linear
models instead

``` r
avgs <- melted %>% 
  group_by(treatment,islet_score,label,info,protein) %>% 
  summarize(avg = mean(expression), .groups = "keep")

p_vals1 <- melted %>% filter(label == "CD45" & islet_score == "Peri.insulitis") %>%
    group_by(protein) %>% 
    wilcox_test(expression~treatment) %>% 
    adjust_pvalue() %>%
    add_significance("p.adj") %>% 
    mutate(comparison = paste0(group1,"_",group2), label = "CD45", islet_score = "Peri.insulitis")
p_vals2 <- melted %>% filter(label == "INS" & islet_score == "Peri.insulitis") %>%
    group_by(protein) %>% 
    wilcox_test(expression~treatment) %>% 
    adjust_pvalue() %>%
    add_significance("p.adj") %>% 
    mutate(comparison = paste0(group1,"_",group2), label = "INS", islet_score = "Peri.insulitis")
p_vals3 <- melted %>% filter(label == "CD45" & islet_score == "Insulitis") %>%
    group_by(protein) %>% 
    wilcox_test(expression~treatment) %>% 
    adjust_pvalue() %>%
    add_significance("p.adj") %>% 
    mutate(comparison = paste0(group1,"_",group2), label = "CD45", islet_score = "Insulitis")
p_vals4 <- melted %>% filter(label == "INS" & islet_score == "Insulitis") %>%
    group_by(protein) %>% 
    wilcox_test(expression~treatment) %>% 
    adjust_pvalue() %>%
    add_significance("p.adj") %>% 
    mutate(comparison = paste0(group1,"_",group2), label = "INS", islet_score = "Insulitis")

p_vals <- rbind(p_vals1,p_vals2,p_vals3,p_vals4)
head(p_vals)
```

    ## # A tibble: 6 x 13
    ##   protein  .y.   group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
    ##   <fct>    <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
    ## 1 S6       expr~ IgG    Spt       22     8        45 4.5 e-2 1      ns          
    ## 2 S6       expr~ IgG    PD1       22     4        35 5.6 e-1 1      ns          
    ## 3 S6       expr~ Spt    PD1        8     4        22 3.68e-1 1      ns          
    ## 4 Histone~ expr~ IgG    Spt       22     8       164 8.82e-5 0.0183 *           
    ## 5 Histone~ expr~ IgG    PD1       22     4        84 2   e-3 0.408  ns          
    ## 6 Histone~ expr~ Spt    PD1        8     4        19 6.83e-1 1      ns          
    ## # ... with 3 more variables: comparison <chr>, label <chr>, islet_score <chr>

``` r
avgs_tbl <- merge(avgs, p_vals[,c(1,9:13)], by=c("protein","islet_score","label"))
avgs_tbl <- avgs_tbl %>% 
    cast(protein+islet_score+label+comparison+p.adj+p.adj.signif~treatment,value="avg") %>%
    mutate(log2FC = ifelse(comparison == "IgG_Spt", log2(Spt/IgG), ifelse(comparison == "IgG_PD1", log2(PD1/IgG), log2(PD1/Spt))))
head(avgs_tbl)
```

    ##   protein    islet_score label comparison      p.adj p.adj.signif      IgG
    ## 1      S6      Insulitis  CD45    IgG_PD1 1.00000000           ns 3557.392
    ## 2      S6      Insulitis  CD45    IgG_Spt 0.00099086          *** 3557.392
    ## 3      S6      Insulitis  CD45    Spt_PD1 1.00000000           ns 3557.392
    ## 4      S6      Insulitis   INS    IgG_Spt 1.00000000           ns 5149.459
    ## 5      S6 Peri.insulitis  CD45    IgG_PD1 1.00000000           ns 6260.097
    ## 6      S6 Peri.insulitis  CD45    IgG_Spt 1.00000000           ns 6260.097
    ##        Spt      PD1      log2FC
    ## 1 7090.272 4791.471  0.42964862
    ## 2 7090.272 4791.471  0.99502087
    ## 3 7090.272 4791.471 -0.56537226
    ## 4 4998.275       NA -0.04299055
    ## 5 9755.671 6369.982  0.02510415
    ## 6 9755.671 6369.982  0.64005608

``` r
avgs_tbl %>%
    ggplot(aes(x=log2FC,y=-log10(p.adj))) + 
    theme_minimal() + 
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color="grey10", fill = NA)) +
    facet_wrap(~label+islet_score+comparison) +
    geom_hline(yintercept = -log10(0.05), color = "grey90") + 
    geom_vline(xintercept = c(-0.5,0.5), color = "grey90") +
    geom_point() + 
    geom_text_repel(data = subset(avgs_tbl, p.adj < 0.05 & abs(log2FC) > 0.5), aes(label=protein)) +
    scale_x_continuous(limits = c(-4,4))
```

![](geomx_analysis_files/figure-gfm/Protein%20expression%20between%20treatment%20groups-1.png)<!-- -->
