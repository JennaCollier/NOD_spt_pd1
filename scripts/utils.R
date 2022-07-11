## utilities for the NOD scRNA-seq analysis

#custom color palettes for various plots
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
  if (pal == "clusters") {
    return(c("Tcon_mem" = "#5A655EFF",
             "Tcon_eff_Cd200" = "#E7C76FFF",
             "Tcon_cm" = "#AE1F63FF",
             "CD8_cm_Gzmm" = "#612A79FF",
             "Treg_Tigit_lo" = "#CDDEB7FF",
             "Tcon_effmem" = "#3B1B53FF",
             "Treg_Tigit_hi" = "#E4AF69FF",
             "Recently_activated" = "#7A65A5FF",
             "Tcon_fh" = "#D58F5CFF",
             "Tcon_eff_Ifngr1" = "#C75127FF",
             "Integrin_associated" = "#837B8DFF",
             "Tcon_eff_Lgals1" = "#924822FF",
             "Acinar_contaminated" = "#D595A7FF",
             "Tcon_eff_Cxcr6" = "#6BD76BFF",
             "Interferon_sensing" = "#802268FF",
             "CD8_mem" = "#5DB1DDFF",
             "CD8_slec" = "#BA6338FF",
             "CD8_cm" = "#466983FF",
             "CD8_eff" = "#F0E685FF",
             "CD8_pexh" = "#749B58FF",
             "CD8_texh" = "#CE3D32FF",
             "Proliferating" = "#5050FFFF"))
  } else {
  stop("Invalid palette. Options: clusters, pm, treatments, tissues")
  }
}