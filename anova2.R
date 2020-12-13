# Performs 2 factor ANOVA & produces Estimate Marginal Means plots 
# in order for user to determine which clusters are which BC subtype
# 
# Inputs: 
#   - dt: Name of variable containing imputed data set with cluster number assignments
#
# Outpus:
#   Creates EMM plots for each subtype & stores them in ANOVA folder
#   
# Example function call: anova2(dot)  <- dot is the variable assigned to clusterDT output

anova2 <- function(dt)
{
  # Subtype genes
  basal <- c("FOXC1","MIA","NDC80","CEP55","ANLN","MELK","ESR1","FOXA1")
  her2 <- c("ERBB2","GRB7","BLVRA","ACTR3B","MYC","SFRP1")
  luminalB <- c("CXXC5","MKI67","BCL2","EGFR","PHGDH","CDH3")
  luminalA <- c("NAT1","SLC39A6","MAPT","UBE2C","PTTG1","CENPF","NUF2","BIRC5")
  pam50 <- c(basal,her2,luminalB,luminalA)
  
  dt2 <- select(dt, c("cluster",all_of(pam50)))
  
  # Creates folders to store EMM plot
  dir.create("ANOVA Plots")
  
  # Performs 2 Factor ANOVA for given subtype
  # Input :
  #   - subtype: Name of vector containing genes associated with desired subtype
  #   - total_gene: Number of genes associated with the desired subtype
  #   - subtype_string: Name of subtype as a string
  #   
  # Output:
  #   - Creates Estimated Marginal Means plot & stores it in ANOVA folder
  # 
  # Example function call: subtype_anova(basal,8,"basal")

  subtype_anova <- function(subtype,total_gene,subtype_string)
  {
    subtype_cluster <- select(dt2, all_of(subtype))
    
    # Lists abundances in a single column ----
    abundance <- gather(subtype_cluster)
    abundance <- abundance[-1]
    names(abundance)[1] <- "abundance"
    
    # Lists cluster number for each abundance value ----
    cluster <- dt2["cluster"]
    cluster <- cbind(cluster,replicate((total_gene-1),cluster$cluster)) %>% gather()
    cluster <- cluster[-1]
    names(cluster)[1] <- "cluster"
    cluster <- as.data.frame(sapply(cluster,as.character))
    
    # Lists gene name for each abundance value ----
    gene <- as.data.frame(subtype) %>% transpose()
    gene <- rbind(gene,gene[rep(1,152),]) %>% gather()
    gene <- gene[-1]
    names(gene)[1] <- "gene"
    
    total <- cbind(cluster,gene,abundance)
    
    # Performs 2 Factor ANOVA ----
    anova <- aov(abundance ~ cluster+gene,data=total)
    
    # Estimated Marginal Means ----
    emm <- emmeans(anova,specs="cluster")
    fig <- plot(emm,comparisons=F) + 
      theme_bw() + 
      labs(x="Estimated marginal mean",title=paste("ANOVA for",subtype_string,"subtype",sep=" "))
    
    
    # Creates png image of plot ----
    orca(fig,file=subtype_string,format="png")
    
    # Moves produced violin plot to ANOVA folder
    file.copy(paste(subtype_string,".png",sep=""),
              "./ANOVA Plots")
    file.remove(paste(subtype_string,".png",sep=""))
    
  }
  
  subtype_anova(basal,8,"basal")
  subtype_anova(luminalA,8,"luminal A")
  subtype_anova(luminalB,6,"luminal B")
  subtype_anova(her2,6,"her2")
  
}









