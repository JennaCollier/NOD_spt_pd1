## utilities for the NOD scRNA-seq analysis

## custom color palettes for various plots
custom_colors <- function(pal) {
  
  if (pal == "tissues") {
    return(c("blood" = "#e41a1c", 
              "pLN" = "#4daf4a", 
              "pancreas" = "#377eb8"))
  }
  
  if (pal == "treatments") {
    return(c("IgG" = "black", 
              "Spt" = rgb(253,128,8, maxColorValue = 255), 
              "PD1" = rgb(128,0,128, maxColorValue = 255)))
  }
  
  if (pal == "pm") {
    return(c("PM" = "#1f78b4",
              "non-PM" = "#33a02c"))
  }
  
  if (pal == "islets") {
    return(c("Insulitis" = "darkred",
             "Peri.insulitis" = "royalblue3",
             "Normal" = "seagreen4"))
  }
  
  if (pal == "colitis") {
    return(c("healthy" = "black",
             "ICI" = "#6B1937",
             "ICI+colitis" = "#D95D8B",
             "colitis" = "#502F2C"))
  }
  
  if (pal == "tetramer") {
    return(c("unknown" = "grey90",
             "NRPv7" = "Firebrick"))
  }
  
  if (pal == "celltypes") {
    return(c("CD8" = "deeppink2",
             "Tcon" = "#011f4b",
             "Treg" = "skyblue"))
  }
  
  if (pal == "clusters") {
    return(c("CD8_cm" = "#00CC99FF",
             "CD8_mem" = "#749B58FF",
             "CD8_pexh" = "#CE3D32FF",
             "CD8_effmem" = "#466983FF",
             "CD8_eff" = "#CC9900FF",
             "CD8_slec" = "#5DB1DDFF",
             "CD8_texh" = "#1A0099FF",
             "Tcon_cm" = "#6BD76BFF",
             "Tcon_mem" = "#D595A7FF",
             "Tcon_effmem" = "#924822FF",
             "Tcon_Tfh" = "#5050FFFF",
             "Tcon_Th21" = "#C75127FF",
             "Tcon_eff" = "#339900FF",
             "Tcon_prog" = "#FFC20AFF",
             "Treg_resting" = "#4775FFFF",
             "Treg_eff" = "#D60047FF",
             "Recently_activated" = "#7A65A5FF",
             "Acinar_contaminated" = "#CDDEB7FF",
             "Interferon_sensing" = "#612A79FF",
             "Proliferating" = "#AE1F63FF"))
  } else {
    
  stop("Invalid palette. Options: clusters, pm, treatments, tissues")
  }
}

## conversion of human gene symbols to mouse gene symbols
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
  mouse = useMart(biomart="ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

## conversion of mouse gene symbols to human gene symbols
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
  mouse = useMart(biomart="ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  mousex <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

## Amino acid functions for TCR analyses
aminos <- function(){
  return(c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))
}

amino_color_scheme <- function(){
  return(make_col_scheme(chars=aminos(), groups=c('neutral', 'hydrophobic', 'acidic', 'acidic', 'very hydrophobic', 'neutral', 'basic', 'hydrophobic', 'basic', 'very hydrophobic', 'hydrophobic', 'neutral', 'neutral', 'neutral', 'basic', 'neutral', 'neutral', 'neutral', 'very hydrophobic', 'very hydrophobic'),
                         
                         cols=c('black', 'darkgoldenrod', 'navyblue', 'navyblue', 'darkgoldenrod2', 'black', 'coral2', 'darkgoldenrod', 'coral2', 'darkgoldenrod2', 'darkgoldenrod', 'black', 'black', 'black', 'coral2', 'black', 'black', 'black', 'darkgoldenrod2', 'darkgoldenrod2')))
  
}