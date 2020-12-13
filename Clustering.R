# Determines which cluster patients are in
#
# Inputs:
#   dt_name - String of the name of xlsx file that has already gone through imputation. 
#             Should not contain gene names. 1st row should not contain patient ID.
#             
#   og_dt - String of the name of original TMT_cleaned file that has not been imputed.
#   
#   pc - Number of principal components to select before clustering
#   
#   cluster_type - String of a type of clustering (options: kmeans, pam, hclust)
# 
# Output: Returns the imputed dataframe with clustering information for each patient.
#         Also creates a csv file of the outputted dataframe called "dt_clusters.csv"
# 
# Example function call: dot <- clustering("knn_OG.xlsx","TMT_cleaned.xlsx",24,"kmeans")

clustering <- function(imputed_dt,og_dt,pc,cluster_type)
{
  # Rows: patients
  # Columns:genes
  dt <- read_excel(imputed_dt,col_names=F)
  
  genes <- read_excel(og_dt,col_names=F)
  genes <- genes[2:nrow(genes),1]     # Selects the column containing gene names
  
  # Assigns the column names in dt with the corresponding gene name
  for(i in 1:ncol(dt))
  {
    colnames(dt)[i] <- genes[i,1]
  }
  
  # Performs PCA on dt
  pca_dt <- prcomp(dt,center=T,scale.=T)
  pc_dt <- as.data.frame(pca_dt$x)
  pc_dt <- pc_dt[,1:pc]    # Selects number of principal components specified in function call
  
  
  # Performs type of clustering specified in function call
  matrix_pc <- as.matrix(pc_dt)
  cluster_result <- PerturbationClustering(data=matrix_pc,
                                           k=4,
                                           clusteringMethod=cluster_type)
  clusters <- as.data.frame(cluster_result$cluster)
  
  
  # Attaches assigned cluster number onto each patient
  # Rows: patients
  # Columns: genes (column name already set to gene name)
  clusters_dt <- cbind(clusters,dt)
  colnames(clusters_dt)[1] <- "cluster"    # Names the column containing cluster info "cluster"
  
  
  # Creates csv file of imputed dataset with cluster assignment for each patient
  write.csv(clusters_dt,"dt_clusters.csv")
  
  return(clusters_dt)
  
}