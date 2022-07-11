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
    return(c("pancreas only" = "#33a02c",
              "pancreas & periphery" = "#1f78b4",
              "periphery only" = "#BB5566"))
  }
  
  if (pal == "celltypes") {
    return(c("CD8" = "#fe4a49",
             "Tcon" = "#011f4b",
             "Treg" = "#6497b1"))
  }
  
  if (pal == "clusters") {
    return(c("CD8_cm" = "#837B8DFF",
             "CD8_cm_Gzmm" = "#CE3D32FF",
             "CD8_mem" = "#749B58FF",
             "CD8_pexh" = "#F0E685FF",
             "CD8_effmem" = "#466983FF",
             "CD8_eff_Cxcr6" = "#BA6338FF",
             "CD8_slec" = "#5DB1DDFF",
             "CD8_texh" = "#802268FF",
             "Tcon_cm" = "#6BD76BFF",
             "Tcon_mem" = "#D595A7FF",
             "Tcon_effmem" = "#924822FF",
             "Tcon_fh" = "#5050FFFF",
             "Tcon_Nrn1" = "#C75127FF",
             "Tcon_eff_Cxcr6" = "#D58F5CFF",
             "Treg_Tigit_lo" = "#7A65A5FF",
             "Treg_Tigit_hi" = "#E4AF69FF",
             "Recently_activated" = "#3B1B53FF",
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
  human = useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
  mouse = useMart(biomart="ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = cell_cycle, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}