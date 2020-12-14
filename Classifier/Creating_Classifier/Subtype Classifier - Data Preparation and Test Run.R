#SUBTYPE CLASSIFIER
#This code includes filtering and data wrangling steps taken to prep data
#The SVM is tested on training data, as well as a preliminary run on ITRAQ

#load packages
library(e1071)
library(readxl)

#import cluster subtypes
tmt_subtypes <- read_excel("dt_subtypes-24PC_common.xlsx")

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

#Create classifier using TMT data and defult parameters
subtype_classifier <- svm(tmt_matrix,subtype_key_factor, kernel = "radial")

#Preict subtypes on training data
subtype_prediction<- predict(subtype_classifier, tmt_matrix)

#Correct predictions from SVM on training data
i<-0
correct_subtype <- 0
for(class in subtype_key_factor){
  i<-i+1
  pred <- subtype_prediction[i]
  if (pred == class) {
    correct_subtype <- correct_subtype + 1
  }
}

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

######  TEST SVM ON ITRAQ  #######

#Use SVM to predict on ITRAQ data using defult parameters
subtype_prediction_iTRAQ<- predict(subtype_classifier, iTRAQ_matrix)

#Get number of correct predictions
i<-0
correct_iTRAQ_subtype <- 0
for(class in iTRAQ_key_factor){
  i<-i+1
  pred <- subtype_prediction_iTRAQ[i]
  if (pred == class) {
    correct_iTRAQ_subtype <- correct_iTRAQ_subtype + 1
  }
}






