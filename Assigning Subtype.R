# Assigning cluster numbers with their corresponding subtype
# 
# Inputs:
#   dt_name - String name of xlsx file that has outputted from identifyingSubtype function.
#             File should contain all protein info and cluster number for each patient.
#             
#   basal - Cluster number which has been determined to be basal subtype
#   
#   her2 - Cluster number which has been determined to be HER2 subtype
#   
#   lumA - Cluster number which has been determined to be luminal A subtype
#   
#   lumB - Cluster number which has been determined to be luminal B subtype
#   
# Output:
#   Returns the original input dataframe but instead of containing cluster numbers, each patient
#   has a subtype annotation. Also creates a csv file of the outputted dataframe called "dt_subtypes.csv"
#   
# Example function call: spot <- assigningSubtype("dt_clusters.xlsx",1,2,4,3)
# 
# Note: TMT_cleaned.xlsx file must be in same folder as current working directory

assigningSubtype <- function(dt_name,basal,her2,lumA,lumB)
{
  patientID <- read_excel("TMT_cleaned.xlsx",col_names=F)
  patientID <- patientID[1,c(-1,-2)] %>% transpose()
  colnames(patientID) <- "patientID"
  
  dt <- read_excel(dt_name) %>% cbind(patientID,.)
  
  basal_type <- filter(dt,cluster == basal)
  her2_type <- filter(dt, cluster == her2)
  luminalA_type <- filter(dt,cluster == lumA)
  luminalB_type <- filter(dt,cluster == lumB)
  
  
  # Replaces cluster numbers with corresponding subtype name
  numberToSubtype <- function(type,subtype_name)
  {
    type$subtype <- subtype_name
    subtype_type2 <- cbind(type$patientID,
                           type$subtype,
                           type[,3:(ncol(type)-1)])
    colnames(subtype_type2)[1] <- "patientID"
    colnames(subtype_type2)[2] <- "subtype"
    return(subtype_type2)
  }
  
  basal_type2 <- numberToSubtype(basal_type,"basal")
  her2_type2 <- numberToSubtype(her2_type,"her2")
  luminalA_type2 <- numberToSubtype(luminalA_type,"luminal A")
  luminalB_type2 <- numberToSubtype(luminalB_type,"luminal B")
  
  final_dt <- rbind(basal_type2,her2_type2,luminalA_type2,luminalB_type2)
  
  # Creates csv file of input dataset with subtype assignment for each patient
  write.csv(final_dt,"dt_subtypes.csv")
  
  return(final_dt)
  
}