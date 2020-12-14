# Produces violin plots in order for user to determine which clusters are which BC subtype
# 
# Inputs:
#   dt - Name of variable containing imputed data set with cluster assignments
# 
# Outputs:
#   Creates violin plot for each gene in PAM50 and stores them in corresponding subtype folder  
# 
# Requirement to run function is the orca command-line utility. Will need to install orca and set
# appropriate path before using this function.
# Installation information: https://github.com/plotly/orca#installation
#
# Example function call: identifyingSubtype(dot)  <- dot is the variable assigned to clusterDT output

identifyingSubtype <- function(dt)
{
  # BC subtypes and their associated genes according to PAM50
  basal <- c("FOXC1","MIA","NDC80","CEP55","ANLN","MELK","ESR1","FOXA1")
  her2 <- c("ERBB2","GRB7","BLVRA","ACTR3B","MYC","SFRP1")
  normal <- c("KRT14","KRT17","KRT5","MLPH","CCNB1","TYMS","UBE2T","RRM2","MMP11")
  luminalB <- c("CXXC5","MKI67","BCL2","EGFR","PHGDH","CDH3")
  luminalA <- c("NAT1","SLC39A6","MAPT","UBE2C","PTTG1","CENPF","NUF2","BIRC5")
  pam50 <- c(basal,her2,luminalB,luminalA)
  
  
  # Removes all genes except PAM50 defined genes
  pam50_dt <- select(dt,c("cluster",all_of(pam50)))
  
  
  # Groups patients by cluster number
  cluster1 <- filter(pam50_dt,cluster == 1)
  cluster2 <- filter(pam50_dt,cluster == 2)
  cluster3 <- filter(pam50_dt,cluster == 3)
  cluster4 <- filter(pam50_dt,cluster == 4)
  
  
  # Produces violin plots for each gene in PAM50. Plots of gene for each subtype will be stored
  # in the same subtype folder. All subtype folders will be stored in "PAM50 Plots" folder
  
  # Creates folders to store violin plots
  dir.create("PAM50 Plots")
  dir.create("./PAM50 Plots/Basal")
  dir.create("./PAM50 Plots/HER2")
  dir.create("./PAM50 Plots/Luminal B")
  dir.create("./PAM50 Plots/Luminal A")
  
  
  # Creates violin plots for Basal genes ----
  for(i in 1:length(basal))
  {
    c1 <- select(cluster1, c("cluster",basal[i]))
    c2 <- select(cluster2, c("cluster",basal[i]))
    c3 <- select(cluster3, c("cluster",basal[i]))
    c4 <- select(cluster4, c("cluster",basal[i]))
    
    all_clusters <- rbind(c1,c2,c3,c4)
    colnames(all_clusters)[2] <- "gene"
    
    fig <- plot_ly(all_clusters,
                   x= ~cluster,
                   y= ~gene,
                   split= ~cluster,
                   type="violin",
                   box=list(visible=T))
    fig <- fig %>% layout(title=basal[i],
                          xaxis=list(title="Cluster"),
                          yaxis=list(title="Relative Protein Abundance",zeroline=F))
    
    # Creates png image of violin plot
    orca(fig,file=basal[i],format="png")
    
    # Moves produced violin plot to Basal folder
    file.copy(paste(basal[i],".png",sep=""),
              "./PAM50 Plots/Basal")
    file.remove(paste(basal[i],".png",sep=""))
  }
  
  
  
  # Creates violin plots for HER2 genes ----
  for(i in 1:length(her2))
  {
    c1 <- select(cluster1, c("cluster",her2[i]))
    c2 <- select(cluster2, c("cluster",her2[i]))
    c3 <- select(cluster3, c("cluster",her2[i]))
    c4 <- select(cluster4, c("cluster",her2[i]))
    
    all_clusters <- rbind(c1,c2,c3,c4)
    colnames(all_clusters)[2] <- "gene"
    
    fig <- plot_ly(all_clusters,
                   x= ~cluster,
                   y= ~gene,
                   split= ~cluster,
                   type="violin",
                   box=list(visible=T))
    fig <- fig %>% layout(title=her2[i],
                          xaxis=list(title="Cluster"),
                          yaxis=list(title="Relative Protein Abundance",zeroline=F))
    
    # Creates png image of violin plot
    orca(fig,file=her2[i],format="png")
    
    # Moves produced violin plot to HER2 folder
    file.copy(paste(her2[i],".png",sep=""),
              "./PAM50 Plots/HER2")
    file.remove(paste(her2[i],".png",sep=""))
  }
  
  
  
  # Creates violin plots for LuminalB genes ----
  for(i in 1:length(luminalB))
  {
    c1 <- select(cluster1, c("cluster",luminalB[i]))
    c2 <- select(cluster2, c("cluster",luminalB[i]))
    c3 <- select(cluster3, c("cluster",luminalB[i]))
    c4 <- select(cluster4, c("cluster",luminalB[i]))
    
    all_clusters <- rbind(c1,c2,c3,c4)
    colnames(all_clusters)[2] <- "gene"
    
    fig <- plot_ly(all_clusters,
                   x= ~cluster,
                   y= ~gene,
                   split= ~cluster,
                   type="violin",
                   box=list(visible=T))
    fig <- fig %>% layout(title=luminalB[i],
                          xaxis=list(title="Cluster"),
                          yaxis=list(title="Relative Protein Abundance",zeroline=F))
    
    # Creates png image of violin plot
    orca(fig,file=luminalB[i],format="png")
    
    # Moves produced violin plot to LuminalB folder
    file.copy(paste(luminalB[i],".png",sep=""),
              "./PAM50 Plots/Luminal B")
    file.remove(paste(luminalB[i],".png",sep=""))
  }
  
  
  
  # Creates violin plots for LuminalA genes ----
  for(i in 1:length(luminalA))
  {
    c1 <- select(cluster1, c("cluster",luminalA[i]))
    c2 <- select(cluster2, c("cluster",luminalA[i]))
    c3 <- select(cluster3, c("cluster",luminalA[i]))
    c4 <- select(cluster4, c("cluster",luminalA[i]))
    
    all_clusters <- rbind(c1,c2,c3,c4)
    colnames(all_clusters)[2] <- "gene"
    
    fig <- plot_ly(all_clusters,
                   x= ~cluster,
                   y= ~gene,
                   split= ~cluster,
                   type="violin",
                   box=list(visible=T))
    fig <- fig %>% layout(title=luminalA[i],
                          xaxis=list(title="Cluster"),
                          yaxis=list(title="Relative Protein Abundance",zeroline=F))
    
    # Creates png image of violin plot
    orca(fig,file=luminalA[i],format="png")
    
    # Moves produced violin plot to LuminalA folder
    file.copy(paste(luminalA[i],".png",sep=""),
              "./PAM50 Plots/Luminal A")
    file.remove(paste(luminalA[i],".png",sep=""))
  }
  
  
  # ----
  
  
 
  
  
}


