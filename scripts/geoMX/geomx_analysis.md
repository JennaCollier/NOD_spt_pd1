geomx_analysis
================
2022-06-29

## Introduction

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
require(useful, quietly = T) #corner function to peek at data.frames
source("scripts/utils.R") #custom color palettes
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
# cd45$islet_score <- gsub("\\-", ".", cd45$islet_score) #change dashes in Peri-insulitis to periods (Peri.insulitis) in data if not already
# colnames(cd45) <- gsub("\\-", ".", colnames(cd45)) #change dashes in protein colnames in data if not already

#read in CD45 data that has been normalized
ins <- as_tibble(t(read_excel(file.path("data/geoMX/norm_geoMean_INS_trimmed.xlsx"))))
colnames(ins) <- ins[1,] #set column names manually (doesn't work in read_excel)
ins <- ins[-1,] #remove row with column names
#ins$islet_score <- gsub("\\-", ".", ins$islet_score) #change dashes in Peri-insulitis to periods (Peri.insulitis) in data if not already
#colnames(ins) <- gsub("\\-", ".", colnames(ins)) #change dashes in protein colnames in data if not already

#merge the CD45 and INS data
geomx <- rbind(cd45,ins)
colnames(geomx)
```

    ##  [1] "ROI"                                              
    ##  [2] "ROI_label"                                        
    ##  [3] "segment"                                          
    ##  [4] "segment_label"                                    
    ##  [5] "ROI_X_Coordinate"                                 
    ##  [6] "ROI_Y_Coordinate"                                 
    ##  [7] "segment_tags"                                     
    ##  [8] "AOI"                                              
    ##  [9] "AOI_surface_area"                                 
    ## [10] "AOI_nuclei_count"                                 
    ## [11] "fov_count"                                        
    ## [12] "fov_counted"                                      
    ## [13] "bindingDensity"                                   
    ## [14] "scan_name"                                        
    ## [15] "slide"                                            
    ## [16] "scan_plexset"                                     
    ## [17] "scan_id"                                          
    ## [18] "roi_id"                                           
    ## [19] "LOT_Mouse_Immune_Cell_Profiling_Protein_Core"     
    ## [20] "LOT_Mouse_IO_Drug_Target_Protein_Module"          
    ## [21] "LOT_Mouse_Immune_Activation_Status_Protein_Module"
    ## [22] "LOT_Mouse_Immune_Cell_Typing_Protein_Module"      
    ## [23] "LOT_Mouse_Pan_Tumor__Protein_Module"              
    ## [24] "LOT_Mouse_nC_Cell_Death_Protein"                  
    ## [25] "LOT_Mouse_nC_PI3K_AKT_Signaling_Protein"          
    ## [26] "treatment"                                        
    ## [27] "mouse"                                            
    ## [28] "islet_score"                                      
    ## [29] "label"                                            
    ## [30] "PD.L1"                                            
    ## [31] "CD45"                                             
    ## [32] "Rb_IgG"                                           
    ## [33] "CTLA4"                                            
    ## [34] "Rt_IgG2b"                                         
    ## [35] "GZMB"                                             
    ## [36] "Histone.H3"                                       
    ## [37] "CD8a"                                             
    ## [38] "CD19"                                             
    ## [39] "CD3e"                                             
    ## [40] "CD11c"                                            
    ## [41] "F480"                                             
    ## [42] "SMA"                                              
    ## [43] "CD4"                                              
    ## [44] "PanCk"                                            
    ## [45] "CD11b"                                            
    ## [46] "GAPDH"                                            
    ## [47] "Fibronectin"                                      
    ## [48] "Ki.67"                                            
    ## [49] "PD.1"                                             
    ## [50] "S6"                                               
    ## [51] "MHCII"                                            
    ## [52] "Rt_IgG2a"                                         
    ## [53] "VISTA"                                            
    ## [54] "B7.H3"                                            
    ## [55] "Tim3"                                             
    ## [56] "OX40L"                                            
    ## [57] "LAG3"                                             
    ## [58] "GITR"                                             
    ## [59] "CD44"                                             
    ## [60] "CD86"                                             
    ## [61] "CD40L"                                            
    ## [62] "CD40"                                             
    ## [63] "CD27"                                             
    ## [64] "CD127"                                            
    ## [65] "ICOS"                                             
    ## [66] "CD28"                                             
    ## [67] "CD34"                                             
    ## [68] "Ly6G&Ly6C"                                        
    ## [69] "FOXP3"                                            
    ## [70] "CD14"                                             
    ## [71] "CD163"                                            
    ## [72] "CD31"                                             
    ## [73] "BatF3"                                            
    ## [74] "S100B"                                            
    ## [75] "Epcam"                                            
    ## [76] "ER"                                               
    ## [77] "Pmel17"                                           
    ## [78] "Her2"                                             
    ## [79] "AR"                                               
    ## [80] "AhR"                                              
    ## [81] "IFNGR"                                            
    ## [82] "GFP"                                              
    ## [83] "BIM"                                              
    ## [84] "p53"                                              
    ## [85] "BCLXL"                                            
    ## [86] "PARP"                                             
    ## [87] "p21"                                              
    ## [88] "Perforin"                                         
    ## [89] "Cleaved_Caspase3"                                 
    ## [90] "BAD"                                              
    ## [91] "gamma.H2AX"                                       
    ## [92] "pS6"                                              
    ## [93] "pGSK3A"                                           
    ## [94] "pPRAS40"                                          
    ## [95] "pAMPK.alpha"                                      
    ## [96] "pAKT1"                                            
    ## [97] "Pan.AKT"                                          
    ## [98] "PLCG1"                                            
    ## [99] "MET"

``` r
geomx <- geomx[,c("scan_id","ROI",colnames(geomx)[26:99])] #remove excess columns
corner(geomx)
```

    ## # A tibble: 5 x 5
    ##   scan_id     ROI                                  treatment mouse islet_score
    ##   <chr>       <chr>                                <chr>     <chr> <chr>      
    ## 1 437,669,535 194aaf9d-eb8a-4754-a0f8-08d6722f5200 PD1       437   Insulitis  
    ## 2 437,669,535 00f36e18-66f6-4920-9dc9-7b058271765d Spt       535   Insulitis  
    ## 3 437,669,535 5b6fff74-162f-4be9-8f74-c4b09dc6c71b Spt       669   Insulitis  
    ## 4 437,669,535 4d2d541a-2f9c-49f8-8f14-ba8e0c48a16a Spt       669   Insulitis  
    ## 5 101,440,647 b0a1b2ef-30a4-4d49-b1f6-0ad1434e693e PD1       101   Insulitis

``` r
#read in gene annotations file that includes additional information for proteins detected
gene_annot <- read_excel(file.path("data/geoMX/gene_annotations.xlsx"), col_names = T)
#gene_annot$protein <- gsub("\\-", ".", gene_annot$protein) #change dashes in protein colnames in data if not already
gene_annot$target_group <- gsub("All Targets;","",gene_annot$target_group) #remove extraneous text from description
head(gene_annot)
```

    ## # A tibble: 6 x 5
    ##   target_group                              code_class description protein gene 
    ##   <chr>                                     <chr>      <chr>       <chr>   <chr>
    ## 1 Myeloid Activation;Checkpoint             Endogenous PD-L1       PD.L1   Cd274
    ## 2 Total Immune                              Endogenous CD45        CD45    Cd45 
    ## 3 Background                                Negative   Rb IgG      Rb_IgG  <NA> 
    ## 4 T cell Activation;T cells;Th cells;Check~ Endogenous CTLA4       CTLA4   Ctla4
    ## 5 Background                                Negative   Rt IgG2b    Rt_IgG~ <NA> 
    ## 6 T cell Activation;Cytotoxicity            Endogenous GZMB        GZMB    Gzmb

Tally the number of samples to get an idea of what the data looks like.
Note that there is only one normal islet in the CD45 group that was
analyzed, so that will be removed before any statistical analyses as it
doesn’t make sense - normal islets should have no CD45+ infiltrate. Note
that there aren’t many samples for the PD1 treatment group.

``` r
geomx %>% group_by(treatment,islet_score,label) %>% tally() %>% cast(label+islet_score ~ treatment, value = "n", fill = 0)
```

    ##   label    islet_score IgG PD1 Spt
    ## 1  CD45      Insulitis  30  16  24
    ## 2  CD45         Normal   1   0   0
    ## 3  CD45 Peri.insulitis  22   4   8
    ## 4   INS      Insulitis  16   0  10
    ## 5   INS         Normal   7   0   0
    ## 6   INS Peri-insulitis  19   0   6

``` r
geomx <- geomx %>% filter(islet_score != "Normal" & label == "CD45" | label == "INS")
```

## Differences by score and region classification

Examine protein expression broadly across all of the samples based on
their region of classification (CD45+ or INS+ region) to confirm that
expression levels generally make sense. CD45+ regions should have higher
expression of immune-cell associated proteins (CD3e, CD45, CD8a, PD-1,
etc.) compared to INS+ regions. Multiple t tests are used with a
Benjamini & Hochberg “BH” correction.

``` r
#melt data to long format for plotting
colnames(geomx)
```

    ##  [1] "scan_id"          "ROI"              "treatment"        "mouse"           
    ##  [5] "islet_score"      "label"            "PD.L1"            "CD45"            
    ##  [9] "Rb_IgG"           "CTLA4"            "Rt_IgG2b"         "GZMB"            
    ## [13] "Histone.H3"       "CD8a"             "CD19"             "CD3e"            
    ## [17] "CD11c"            "F480"             "SMA"              "CD4"             
    ## [21] "PanCk"            "CD11b"            "GAPDH"            "Fibronectin"     
    ## [25] "Ki.67"            "PD.1"             "S6"               "MHCII"           
    ## [29] "Rt_IgG2a"         "VISTA"            "B7.H3"            "Tim3"            
    ## [33] "OX40L"            "LAG3"             "GITR"             "CD44"            
    ## [37] "CD86"             "CD40L"            "CD40"             "CD27"            
    ## [41] "CD127"            "ICOS"             "CD28"             "CD34"            
    ## [45] "Ly6G&Ly6C"        "FOXP3"            "CD14"             "CD163"           
    ## [49] "CD31"             "BatF3"            "S100B"            "Epcam"           
    ## [53] "ER"               "Pmel17"           "Her2"             "AR"              
    ## [57] "AhR"              "IFNGR"            "GFP"              "BIM"             
    ## [61] "p53"              "BCLXL"            "PARP"             "p21"             
    ## [65] "Perforin"         "Cleaved_Caspase3" "BAD"              "gamma.H2AX"      
    ## [69] "pS6"              "pGSK3A"           "pPRAS40"          "pAMPK.alpha"     
    ## [73] "pAKT1"            "Pan.AKT"          "PLCG1"            "MET"

``` r
melted <- data.frame(geomx) %>% melt(id.vars = c("treatment","mouse","islet_score","label","scan_id","ROI"))
colnames(melted)[which(colnames(melted) == "variable")] <- "protein"
colnames(melted)[which(colnames(melted) == "value")] <- "expression"
head(melted)
```

    ##   treatment mouse islet_score label     scan_id
    ## 1       PD1   437   Insulitis  CD45 437,669,535
    ## 2       Spt   535   Insulitis  CD45 437,669,535
    ## 3       Spt   669   Insulitis  CD45 437,669,535
    ## 4       Spt   669   Insulitis  CD45 437,669,535
    ## 5       PD1   101   Insulitis  CD45 101,440,647
    ## 6       PD1   440   Insulitis  CD45 101,440,647
    ##                                    ROI protein         expression
    ## 1 194aaf9d-eb8a-4754-a0f8-08d6722f5200   PD.L1 63.127837555183959
    ## 2 00f36e18-66f6-4920-9dc9-7b058271765d   PD.L1 140.98315394086532
    ## 3 5b6fff74-162f-4be9-8f74-c4b09dc6c71b   PD.L1 105.76616280048388
    ## 4 4d2d541a-2f9c-49f8-8f14-ba8e0c48a16a   PD.L1 260.89662556491243
    ## 5 b0a1b2ef-30a4-4d49-b1f6-0ad1434e693e   PD.L1 152.26059349306291
    ## 6 dd5988b6-8ed2-418d-87e1-8d3012d88cd0   PD.L1 96.717817435207635

``` r
#correct the class of some of the variables
melted$expression <- as.numeric(melted$expression)
melted$protein <- factor(melted$protein, levels = melted %>% 
                           group_by(protein) %>% 
                           summarize(avg = mean(expression)) %>% 
                           arrange(desc(avg)) %>% 
                           pull(var = protein))
melted$treatment <- factor(melted$treatment, levels = c("IgG","Spt","PD1"))
melted$islet_score <- factor(melted$islet_score)

melted <- merge(melted, gene_annot[,c("protein","gene","code_class")], by = c("protein")) #add protein annotation info

#calculate the p values for plotting using t tests
p_vals <- melted %>%
  rstatix::group_by(protein) %>%
  rstatix::t_test(expression ~ label) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "protein", dodge = 0.8) # important for positioning!
p_vals$y.position <- log2(p_vals$y.position) #log2 transformed data will be plotted, so djust pvalue position

#plot the expression data as violin plots grouped by region classification
melted %>% 
  ggplot(aes(x = protein, 
             y = log2(expression), #log2 expression values are used
             color = label, 
             fill = label)) + 
  geom_violin(alpha = 0.3, 
              adjust = 2, 
              scale = "width", 
              draw_quantiles = c(0.5), #draw the median as a line
              position = position_dodge(0.8)) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank()) +
  labs(title = "GeoMX protein expression: CD45+/INS+ regions") +
  add_pvalue(p_vals %>% filter(p.adj < 0.05),  #only plot significant values
             color = "black",
               xmin = "xmin", 
               xmax = "xmax",
               label = "p.adj.signif",
               tip.length = 0,
             show.legend = F,
             inherit.aes = F)
```

<img src="geomx_analysis_files/figure-gfm/differences in protein expression by region classification (CD45/INS)-1.png" width="900px" />

Expression of CD45 and other lymphocyte-associated proteins (CD19,
CD11c, F4/80, CD3e, etc) are higher in the CD45 regions as expected.
Look at differences between islets with peri-insulitis and overt
insulitis. Multiple t tests are used with a Benjamini & Hochberg “BH”
correction.

``` r
p1_vals <- melted %>%
  filter(label == "CD45") %>%
  rstatix::group_by(protein) %>%
  rstatix::t_test(expression ~ islet_score) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "protein", dodge = 0.8) #important for positioning!
p1_vals$y.position <- log2(p1_vals$y.position) #log2 transformed expression values are used, so this adjustment is needed
p1_vals <- merge(p1_vals, gene_annot[,c("protein","code_class")], by = "protein")

p1 <- melted %>%
    filter(label == "CD45") %>%
    ggplot(aes(x = protein, 
               y = log2(expression), 
               color = islet_score, 
               fill = islet_score)) +
    geom_violin(alpha=0.3, 
                adjust=2, 
                scale = "width", 
                draw_quantiles = c(0.5), 
                position = position_dodge(0.8)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1),
          axis.title.x = element_blank()) +
    scale_color_manual(values = custom_colors("islets")) +
    scale_fill_manual(values = custom_colors("islets")) +
    { if (sum(p1_vals$p.adj < 0.05) > 0) #check if there are any significant p values, otherwise, don't plot any NS
      add_pvalue(p1_vals %>% filter(p.adj < 0.05),  #only plot significant values
               color = "black",
               xmin = "xmin",
               xmax = "xmax",
               label = "p.adj.signif",
               tip.length = 0,
               show.legend = F,
               inherit.aes = F) } +
    labs(title = "GeoMX protein expression (CD45+ regions)")


p2_vals <- melted %>%
  filter(label == "INS") %>%
  rstatix::group_by(protein) %>%
  rstatix::t_test(expression ~ islet_score) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "protein", dodge = 0.8) # important for positioning!
p2_vals$y.position <- log2(p2_vals$y.position) #log2 transformed expression values are used, so this adjustment is needed
p2_vals <- merge(p2_vals, gene_annot[,c("protein","code_class")], by = "protein")

p2_vals %>% filter(p.adj < 0.05) %>% arrange(p.adj) %>% head(15) #show some of the top DE proteins
```

    ##    protein        .y.    group1         group2 n1 n2 statistic       df
    ## 1    GAPDH expression Insulitis         Normal 26  7  5.510810 25.82472
    ## 2    IFNGR expression Insulitis         Normal 26  7  5.196953 29.47353
    ## 3     CD40 expression Insulitis         Normal 26  7  4.116087 28.23845
    ## 4     CD45 expression    Normal Peri-insulitis  7 25 -4.149129 29.96856
    ## 5      p21 expression Insulitis         Normal 26  7  4.300210 26.21041
    ## 6    VISTA expression Insulitis         Normal 26  7  4.038324 29.35055
    ## 7    BCLXL expression Insulitis         Normal 26  7  3.982615 28.95669
    ## 8    CD11c expression Insulitis         Normal 26  7  3.847352 28.88931
    ## 9    IFNGR expression    Normal Peri-insulitis  7 25 -3.892207 29.96216
    ## 10     MET expression Insulitis         Normal 26  7  3.857970 30.57686
    ## 11     CD4 expression Insulitis         Normal 26  7  3.740306 29.18747
    ## 12    CD8a expression Insulitis         Normal 26  7  3.801756 25.85006
    ## 13     MET expression    Normal Peri-insulitis  7 25 -3.760412 28.53623
    ## 14    CD3e expression Insulitis         Normal 26  7  3.583154 28.62778
    ## 15    LAG3 expression Insulitis         Normal 26  7  3.644435 26.01626
    ##           p      p.adj p.adj.signif y.position                 groups  x
    ## 1  8.97e-06 0.00144900           **   11.52768      Insulitis, Normal  6
    ## 2  1.40e-05 0.00144900           **   12.03331      Insulitis, Normal  7
    ## 3  3.04e-04 0.01221300            *   10.89488      Insulitis, Normal 44
    ## 4  2.53e-04 0.01221300            *   13.45210 Normal, Peri-insulitis  3
    ## 5  2.10e-04 0.01221300            *   10.94884      Insulitis, Normal 53
    ## 6  3.54e-04 0.01221300            *   10.89635      Insulitis, Normal 13
    ## 7  4.20e-04 0.01242000            *   11.18571      Insulitis, Normal 20
    ## 8  6.07e-04 0.01256490            *   11.04978      Insulitis, Normal 12
    ## 9  5.14e-04 0.01256490            *   13.20345 Normal, Peri-insulitis  7
    ## 10 5.50e-04 0.01256490            *   10.88190      Insulitis, Normal 40
    ## 11 8.00e-04 0.01273846            *   10.94928      Insulitis, Normal 10
    ## 12 7.88e-04 0.01273846            *   10.93289      Insulitis, Normal 14
    ## 13 7.78e-04 0.01273846            *   12.79930 Normal, Peri-insulitis 40
    ## 14 1.00e-03 0.01293750            *   10.93084      Insulitis, Normal 22
    ## 15 1.00e-03 0.01293750            *   10.90278      Insulitis, Normal 24
    ##         xmin      xmax code_class
    ## 1   5.733333  6.000000    Control
    ## 2   6.733333  7.000000 Endogenous
    ## 3  43.733333 44.000000 Endogenous
    ## 4   3.000000  3.266667 Endogenous
    ## 5  52.733333 53.000000 Endogenous
    ## 6  12.733333 13.000000 Endogenous
    ## 7  19.733333 20.000000 Endogenous
    ## 8  11.733333 12.000000 Endogenous
    ## 9   7.000000  7.266667 Endogenous
    ## 10 39.733333 40.000000 Endogenous
    ## 11  9.733333 10.000000 Endogenous
    ## 12 13.733333 14.000000 Endogenous
    ## 13 40.000000 40.266667 Endogenous
    ## 14 21.733333 22.000000 Endogenous
    ## 15 23.733333 24.000000 Endogenous

``` r
p2 <- melted %>%
    filter(label == "INS") %>%
    ggplot(aes(x = protein, 
               y = log2(expression), 
               color = islet_score, 
               fill = islet_score)) +
    geom_violin(alpha=0.3, 
                adjust=2, 
                scale = "width", 
                draw_quantiles = c(0.5), 
                position = position_dodge(0.8)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1),
          axis.title.x = element_blank()) +
    scale_color_manual(values = custom_colors("islets")) +
    scale_fill_manual(values = custom_colors("islets")) +
    { if (sum(p2_vals$p.adj < 0.05) > 0) #check if there are any significant p values, otherwise, don't plot any NS
      add_pvalue(p2_vals %>% filter(p.adj < 0.05), #only plot significant values
               color = "black",
               xmin = "xmin",
               xmax = "xmax",
               label = "p.adj.signif",
               tip.length = 0,
               show.legend = F,
               inherit.aes = F) } +
    labs(title = "GeoMX protein expression (INS+ regions)")

grid.arrange(p1, p2, nrow=2)
```

<img src="geomx_analysis_files/figure-gfm/differences in protein expression by islet score (Insulitis/Peri.Insulitis)-1.png" width="900px" />

## PCA analysis

There are no significant differences between CD45+ regions within islets
with insulitis vs. peri-insulitis. The only identified differences were
between INS+ normal islets (in IgG mice) and those with peri-insulitis
or insulitis. Most DE proteins in the INS+ regions are immune
cell-associated (CD40, CD45, VISTA, CD11c) which indicates that there is
some contamination of these regions with immune infiltrate. However, the
immune cells in this region are likely to be associated with the
“leading edge” of the infiltrate that are closest to the INS+ cells.

``` r
df_for_pca <- geomx %>% filter(label == "CD45")
df_for_pca[,-c(1:6)] <- sapply(df_for_pca[,-c(1:6)], as.numeric)
df_for_pca[,-c(1:6)] <- sapply(df_for_pca[,-c(1:6)], log2)

prin_comp <- rda(df_for_pca[,-c(1:6)])
pca_scores <- scores(prin_comp)

# Define shapes for plot
shapes = c(16, 17, 15) 
shapes <- shapes[as.numeric(factor(df_for_pca$islet_score))]

#plot using base R
plot(pca_scores$sites[,1],
     pca_scores$sites[,2],
     pch = shapes,
     col = paste(custom_colors("treatments")[df_for_pca$treatment]),
     xlim = c(-3,3),
     ylim = c(-3,3),
     cex=3)
ordiellipse(prin_comp, factor(df_for_pca$treatment, levels = c("IgG","Spt","PD1")), 
            conf=0.90, 
            col = custom_colors("treatments"), 
            lwd = 2)
```

![](geomx_analysis_files/figure-gfm/PCA%20analysis%20CD45%20regions-1.png)<!-- -->

``` r
#keep data to later plot in ggplot
df_all <- cbind(df_for_pca, pca_scores$sites)
pca_weights <- data.frame(pca_scores$species, label = "CD45") %>% 
                          rownames_to_column("protein")

#collect variance explained for each PC for labeling
variance_explained_labels <- c(paste0("PC1 (", 
                                      round(summary(prin_comp)$cont$importance[2, "PC1"]*100, 2),
                                      "% of variance)"),
                               paste0("PC2 (", 
                                      round(summary(prin_comp)$cont$importance[2, "PC2"]*100, 2),
                                      "% of variance)"))
```

``` r
df_for_pca <- geomx %>% filter(label == "INS")
df_for_pca[,-c(1:6)] <- sapply(df_for_pca[,-c(1:6)], as.numeric)
df_for_pca[,-c(1:6)] <- sapply(df_for_pca[,-c(1:6)], log2)

prin_comp <- rda(df_for_pca[,-c(1:6)])
pca_scores <- scores(prin_comp)

# Define shapes for plot
shapes = c(16, 17, 15) 
shapes <- shapes[as.numeric(factor(df_for_pca$islet_score))]

plot(pca_scores$sites[,1],
     pca_scores$sites[,2],
     pch = shapes,
     col = paste(custom_colors("treatments")[df_for_pca$treatment]),
     xlim = c(-4,3),
     ylim = c(-2,4),
     cex=3)
ordiellipse(prin_comp, 
            factor(df_for_pca$treatment, levels = c("IgG","Spt","PD1")), 
            conf=0.90, 
            col = custom_colors("treatments"), 
            lwd = 2)
```

![](geomx_analysis_files/figure-gfm/PCA%20analysis%20INS-1.png)<!-- -->

``` r
#gather data for PCA plot in ggplot
df_for_pca <- cbind(df_for_pca, pca_scores$sites)
df_all <- rbind(df_all, df_for_pca)
pca_weights <- rbind(pca_weights, 
                     data.frame(pca_scores$species, label = "INS") %>% 
                     rownames_to_column("protein"))

#collect variance explained for each PC (as a percentage of 100, rounded to two decimals) for labeling axes in plots
variance_explained_labels <- c(variance_explained_labels,
                               paste0("PC1 (", 
                                      round(summary(prin_comp)$cont$importance[2, "PC1"]*100, 2),
                                      "% of variance)"),
                               paste0("PC2 (", 
                                      round(summary(prin_comp)$cont$importance[2, "PC2"]*100, 2),
                                      "% of variance)"))
names(variance_explained_labels) <- c("PC1_CD45","PC2_CD45","PC1_INS","PC2_INS") #name the vector of axis titles
```

These plots can be generated using ggplot2 which looks nicer and allows
for increased customization. Both the PCA plots will be generated
separately to allow for the axes to be labeled with the customized
variance labels (rather than facetting the plots).

``` r
pca_cd45 <- df_all %>%
    filter(label == "CD45") %>%
    ggplot(aes(x = PC1, 
               y= PC2, 
               color = treatment, 
               shape = islet_score, 
               group = treatment)) +
    geom_point(size = 4, 
               alpha = 0.6) +
    stat_ellipse(level = 0.9) + #90% confidence intervals
    scale_color_manual(values = custom_colors("treatments")) +
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
    ggplot(aes(x = PC1, y= PC2, 
               color = treatment, 
               shape = islet_score, 
               group = treatment)) +
    geom_point(size = 4, 
               alpha = 0.6) +
    stat_ellipse(level = 0.9) + #90% confidence intervals
    scale_color_manual(values = custom_colors("treatments")) +
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
proteins_of_interest <- c("LAG3","Perforin","GZMB","CTLA4","Ki.67","ICOS","CD127","IFNGR","Tim3","BAD","GITR","PD.1","Icos","BCLXL","Cleaved_Caspase3","CD4","CD8a","VISTA")

pca_labels <- pca_weights %>% 
  filter(protein %in% proteins_of_interest)

weights_cd45 <- pca_weights %>%
    filter(label == "CD45") %>%
    ggplot(aes(x = PC1, y = PC2)) +
    facet_wrap(~ label, scales = "free") +
    scale_color_manual(values = custom_colors("treatments")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 15)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(4, "mm")),
                 color = "black",
                 alpha = 0.5) +
    geom_text_repel(data = pca_labels %>% filter(label=="CD45"), aes(label = protein), 
                    color = "black",
                    max.overlaps = 4) +
    labs(x = variance_explained_labels["PC1_CD45"],
         y = variance_explained_labels["PC2_CD45"])

weights_ins <- pca_weights %>%
    filter(label == "INS") %>%
    ggplot(aes(x = PC1, y= PC2)) +
    facet_wrap(~ label, scales = "free") +
    scale_color_manual(values = custom_colors("treatments")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 15)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(4, "mm")),
                 color = "black",
                 alpha = 0.5) +
    geom_text_repel(data = pca_labels %>% filter(label=="INS"), aes(label = protein), 
                    color = "black",
                    max.overlaps = 4) +
    labs(x = variance_explained_labels["PC1_INS"],
         y = variance_explained_labels["PC2_INS"])

grid.arrange(weights_cd45, weights_ins, nrow = 1)
```

    ## Warning: ggrepel: 5 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](geomx_analysis_files/figure-gfm/ggplot%20of%20PCA%20weights-1.png)<!-- -->

# Treatment-associated differences

Next, the differences between treatment groups will be analyzed by
generating volcano plots. Here is an example of the comparison between
PD1 and Spt. There are no significant differences.

``` r
meltednrml <- melted %>% 
                filter(label == "CD45") %>% #examine the CD45 gated regions
                filter(!protein %in% c("Rt_IgG2b","Rb_IgG","Rt_IgG2a","CD45","INS","S6","GFP","MHCII","GAPDH","Histone.H3","PanCk")) %>% #remove housekeeping genes and irrelevant genes not in the dataset (GFP is not in this mouse and MHCII is different in NOD)
                dplyr::select(-islet_score) 
meltednrml$treatment <- relevel(meltednrml$treatment, ref = "PD1") #use PD1 as the reference group
df <- data.frame()

for (current_protein in as.character(unique(meltednrml$protein))) {
  tmp <- meltednrml %>% 
            filter(protein == current_protein)
  tmp <- data.frame(summary(lme(expression ~ treatment, #mixed linear model to look at correlation of expression to treatment group
                                random = list(~1|mouse, #control for the mouse examined
                                              ~1|scan_id), #also control for each scan
                                data = tmp))$tTable[-1,]) %>% 
          rownames_to_column("comparison") %>% #move the treatment group being compared to a column
          mutate(protein = current_protein, ref = "PD1", comparison = gsub("treatment", "", comparison)) #clean up
  df <- rbind(df, tmp) #add the temporary df with one protein to the df with data for all the proteins
}

avgs <- meltednrml %>% 
          group_by(treatment, protein) %>% 
          summarize(avg = mean(expression), .groups = "keep") #calculate the mean expression for each protein within each treatment

tbl <- merge(df, avgs %>% cast(protein ~ treatment, value = "avg"), by = c("protein")) #merge the df with the casted averages (wide format)
tbl <- tbl %>% mutate(log2FC = ifelse(comparison == "Spt", log2(PD1/Spt), log2(PD1/Spt))) #calculate the log fold change based on the desired comparison
pd1_v_spt <- tbl %>% 
              filter(comparison == "Spt") %>% #only select the relevant comparison
              dplyr::select(-IgG) #remove column containing average values that aren't relevant
pd1_v_spt$p.adj <- p.adjust(pd1_v_spt$p.value, method = "BH") #adjust the p values using BH

log2FC_cutoff <- 1.25 #define a cutoff log2fc
pd1_v_spt$log2FC_signif <- case_when(pd1_v_spt$log2FC > log2FC_cutoff ~ "PD1",
                              pd1_v_spt$log2FC < -log2FC_cutoff ~ "Spt",
                              TRUE ~ "none")

## PLOT THE DATA
proteins_of_interest <- c("CTLA4","Ki.67","CD4","Tim3","LAG3","PD.1","CD8a")
pd1_v_spt %>%
  ggplot(aes(x = log2FC, 
             y = -log10(p.adj),
             color = log2FC_signif)) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "Average Log2 Fold Change (PD1 / Spt)",
       y = expression(-log[10]~(Adjusted~p~value)),
       title="CD45+ regions") +
  #geom_hline(yintercept = -log10(0.05), color = "grey90") + 
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), color = "grey90") +
  geom_point(size = 2) + 
  scale_x_continuous(limits = c(-3, 3)) +
  scale_color_manual(values = custom_colors("treatments")) +
  #scale_y_continuous(limits = c(0, 1.5)) +
  geom_label_repel(data = subset(pd1_v_spt, protein %in% proteins_of_interest), aes(label = protein, color = log2FC_signif),
                   label.size = NA,
                   label.padding = 0.2, 
                   na.rm = TRUE,
                   fill = alpha(c("white"), 0.8), 
                   size = 7, 
                   max.overlaps = 50) 
```

![](geomx_analysis_files/figure-gfm/Protein%20expression%20between%20PD1%20and%20Spt-1.png)<!-- -->
