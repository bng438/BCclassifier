#SUBTYPE CLASSIFIER
#This code contains data wrangling and filtering steps needed to train the SVM
#A default parameter run of the SVM is included for testing.

#load packages
library(e1071)
library(readxl)

#import TMT/knned/common genes list
kNN_common <- read_excel("knn_OG_common.xlsx",col_names = FALSE)

#import TMT staging key
staging_key <- read_excel("TMT staging key.xlsx")

#remove first row with genes
kNN_common <- kNN_common[-1,] 

#order patient names alphabetically
kNN_common_ordered <- kNN_common[order(kNN_common[,1]),]
staging_key_ordered <- staging_key[order(staging_key[,3]),]

#Remove Not available stages
staging_key_trimmed <-subset(staging_key_ordered, Pathologic_Stage!="Not available")

#Remove kNN rows not included in staging key
kNN_common_trimmed <- subset(kNN_common_ordered, kNN_common_ordered[[1]] %in% c(staging_key_trimmed[[3]]))
'%notin%' <- Negate('%in%')

#people in key not kNN
missing_people <- subset(staging_key_trimmed, staging_key_trimmed[[3]] %notin% c(kNN_common_trimmed[[1]]))

#remove missing people
key_missingremoved <- subset(staging_key_trimmed, staging_key_trimmed[[3]] %in% c(kNN_common_trimmed[[1]]))

#staging keys to integers
numbered_staging_key <- transform(key_missingremoved, id=as.numeric(factor(Pathologic_Stage)))

#order alphabetically by stage
numbered_staging_key <- numbered_staging_key[order(numbered_staging_key[,2]),]
numbered_staging_key[c(1:4), 4] = 1
numbered_staging_key[c(5:73), 4] = 2
numbered_staging_key[c(74:105), 4] = 3

# order by pateint name again
numbered_staging_key <- numbered_staging_key[order(numbered_staging_key[,3]),]

#staging key only numbers
numbers_key <- subset(numbered_staging_key, select = c(id))

#Convert data frame to matrix
kNN_common_trimmed <- kNN_common_trimmed[,-1] 
kNN_common_matrix <-as.matrix(kNN_common_trimmed)
kNN_common_matrix <- apply(kNN_common_matrix, 1, as.numeric)
kNN_common_matrix <- t(kNN_common_matrix)

#Convert key to factor
numbers_key_factor = factor(numbers_key[,1])

#run svm. Predict on training data.
classifier_common <- svm(kNN_common_matrix,numbers_key_factor, kernel = "radial")
predictions_common <- predict(classifier_common, kNN_common_matrix)

#Correct predictions on training data
i<-0
correct_common <- 0
for(class in numbers_key_factor){
  i<-i+1
  pred <- predictions_common[i]
  if (pred == class) {
    correct_common <- correct_common + 1
  }
}

#import iTRAQ
iTRAQ <- read_excel("iTRAQ_common_genes.xlsx")

#order alphabetically by patient
iTRAQ_ordered <- iTRAQ[order(iTRAQ[,1]),]

#import iTRAQ key
iTRAQ_key <- read_excel("iTRAQ Subtype & Staging Key.xlsx",col_names = FALSE)
names(iTRAQ_key) <- lapply(iTRAQ_key[1,], as.character)
iTRAQ_key <- iTRAQ_key[-1,]

# order by pateint name to match iTRAQ ordering
iTRAQ_key <- iTRAQ_key[order(iTRAQ_key[,1]),]

# trim patient names to match
names(iTRAQ_key)[1] <- "PATIENTID"
iTRAQ_key$PATIENTID <- substr(iTRAQ_key$PATIENTID, 6, 12)
iTRAQ_ordered$Gene <- substr(iTRAQ_ordered$Gene, 1, 7)

#Remove iTRAQ rows not included in staging key
iTRAQ_common_trimmed <- subset(iTRAQ_ordered, iTRAQ_ordered[[1]] %in% c(iTRAQ_key[[1]]))
'%notin%' <- Negate('%in%')

#people in key not iTRAQ
missing_iTRAQ <- subset(iTRAQ_key, iTRAQ_key[[1]] %notin% c(iTRAQ_common_trimmed[[1]]))

#remove missing people
iTRAQ_key_missingremoved <- subset(iTRAQ_key, iTRAQ_key[[1]] %in% c(iTRAQ_common_trimmed[[1]]))

# remove duplicates from iTRAQ data
iTRAQ_no_dupes <- iTRAQ_common_trimmed[!duplicated(iTRAQ_common_trimmed$Gene), ]

#order alphabetically by stage
iTRAQ_key <- iTRAQ_key[order(iTRAQ_key[,2]),]

#create factor column of the subtype key
iTRAQ_key <- transform(iTRAQ_key, id=as.numeric(factor("AJCC Stage")))

#relabel each stage such that it is a numeric
iTRAQ_key[c(1:12), 4] = 1
iTRAQ_key[c(13:76), 4] = 2
iTRAQ_key[c(77:103), 4] = 3
iTRAQ_key[c(104:105), 4] = 4

# order by pateint name to match iTRAQ ordering
iTRAQ_key <- iTRAQ_key[order(iTRAQ_key[,1]),]

#remove unnecessary columns
iTRAQ_key = subset(iTRAQ_key, select = c(id))

# convert to factor
iTRAQ_factor = factor(iTRAQ_key[,1])

# Convert iTRAQ to a data matrix, remove dupplicate patient
iTRAQ_matrix <-data.matrix(iTRAQ_no_dupes)
iTRAQ_matrix <- iTRAQ_matrix[,-1] 
iTRAQ_matrix <- apply(iTRAQ_matrix, 1, as.numeric)
iTRAQ_matrix <- t(iTRAQ_matrix)

#run preidction svm on iTRAQ data, trained on KNN data 
predictions_iTRAQ <- predict(classifier_common, iTRAQ_matrix)

#Correct predictions on testing data (iTRAQ)
i<-0
correct_iTRAQ <- 0
for(class in iTRAQ_factor){
  i<-i+1
  pred <- predictions_iTRAQ[i]
  if (pred == class) {
    correct_iTRAQ <- correct_iTRAQ + 1
  }
}

#run svm with default parameters. Predict on training data (iTRAQ)
classifier_common_iTRAQ <- svm(iTRAQ_matrix,iTRAQ_factor, kernel = "radial")

#run preidction svm on iTRAQ data, trained on iTRAQ data
predictions_iTRAQ <- predict(classifier_common_iTRAQ, iTRAQ_matrix)

#Get number of correct predictions on training data
i<-0
correct_iTRAQ <- 0
for(class in iTRAQ_factor){
  i<-i+1
  pred <- predictions_iTRAQ[i]
  if (pred == class) {
    correct_iTRAQ <- correct_iTRAQ + 1
  }
}

