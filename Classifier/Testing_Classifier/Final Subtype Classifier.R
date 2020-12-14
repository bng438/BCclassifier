#load packages
library(e1071)
library(readxl)

#import cluster subtypes
tmt_subtypes <- read_excel("dt_subtypes-24PC_common18.xlsx")

########################## CREATING KEY and dataset ###########################

#order patients alphabetically
tmt_subtypes_ordered <- tmt_subtypes[order(tmt_subtypes[,1]),]
#isolate subtypes
subtype_key <- subset(tmt_subtypes_ordered, select = c(subtype))
#turn subtypes into numbers
subtype_key_numeric <- transform(subtype_key, subtype=as.numeric(factor(subtype)))
#convert key to factor
subtype_key_factor = factor(subtype_key_numeric[,1])

#get dataset
tmt_data <- tmt_subtypes_ordered[,-1]
tmt_data <- tmt_data[,-1]
#convert into matrix
tmt_matrix <-data.matrix(tmt_data)


######################## Train svm and test on itself #########################

subtype_classifier <- svm(tmt_matrix,subtype_key_factor,
                               type="C-classification",
                               kernel = "radial",
                               gamma=1E-6,cost=100)


##################### Get iTRAQ key and data ##################################

#import iTRAQ key
iTRAQ_key <- read_excel("iTRAQ Subtype & Staging Key.xlsx",col_names = FALSE)
names(iTRAQ_key) <- lapply(iTRAQ_key[1,], as.character)
iTRAQ_key <- iTRAQ_key[-1,]
# order key by patient name
iTRAQ_key <- iTRAQ_key[order(iTRAQ_key[,1]),]
#Isolate subtype column of key
iTRAQ_subtype_key <- subset(iTRAQ_key[,3])
#Rename column
names(iTRAQ_subtype_key)[1] <- "subtype"
#Give iTRAQ key numeric values
iTRAQ_key_numeric <- transform(iTRAQ_subtype_key, subtype=as.numeric(factor(subtype)))
#Turn iTRAQ key into a factor
iTRAQ_key_factor = factor(iTRAQ_key_numeric[,1])

#import iTRAQ data
iTRAQ <- read_excel("iTRAQ_common_genes.xlsx")
#order alphabetically by patient
iTRAQ_ordered <- iTRAQ[order(iTRAQ[,1]),]
# trim patient names to match
iTRAQ_ordered$Gene <- substr(iTRAQ_ordered$Gene, 1, 7)
#remove duplicates from iTRAQ data
iTRAQ_no_dupes <- iTRAQ_ordered[!duplicated(iTRAQ_ordered$Gene), ]
# data handling on iTRAQ data before SVM
iTRAQ_matrix <- data.matrix(iTRAQ_no_dupes)
iTRAQ_matrix <- iTRAQ_matrix[,-1] 
#iTRAQ_matrix <- apply(iTRAQ_matrix, 1, as.numeric)
#iTRAQ_matrix <- t(iTRAQ_matrix)




#test on best and worst parameters
subtype_prediction_iTRAQ_best<- predict(subtype_classifier, iTRAQ_matrix)

i<-0
correct_iTRAQ_subtype_best <- 0

Stage1Count <- 0 #number of subtype 1's in iTRAQ
Stage1Correct <- 0 #number of above guessed correctly

Stage2Count <- 0
Stage2Correct <- 0

Stage3Count <- 0
Stage3Correct <- 0

Stage4Count <- 0
Stage4Correct <- 0

Pred1Count <- 0 #number of times classifier predicted Subtype 1
Pred2Count <- 0
Pred3Count <- 0
Pred4Count <- 0


for(class in iTRAQ_key_factor){
  i<-i+1
  pred <- subtype_prediction_iTRAQ_best[i]
  if (pred == class) {
    correct_iTRAQ_subtype_best <- correct_iTRAQ_subtype_best + 1
  }
  if (class == 1){
    Stage1Count <- Stage1Count + 1
  }
  if (class == 2){
    Stage2Count <- Stage2Count + 1
  }
  if (class == 3){
    Stage3Count <- Stage3Count + 1
  }
  if (class == 4){
    Stage4Count <- Stage4Count + 1
  }
  if (class == 1 & pred == class) {
    Stage1Correct <- Stage1Correct +1
  }
  if (class == 2 & pred == class) {
    Stage2Correct <- Stage2Correct +1
  }
  if (class == 3 & pred == class) {
    Stage3Correct <- Stage3Correct +1
  }
  if (class == 4 & pred == class) {
    Stage4Correct <- Stage4Correct +1
  }
  if (pred == 1){
    Pred1Count <- Pred1Count + 1
  }
  if (pred == 2){
    Pred2Count <- Pred2Count + 1
  }
  if (pred == 3){
    Pred3Count <- Pred3Count + 1
  }
  if (pred == 4){
    Pred4Count <- Pred4Count + 1
  }
}

