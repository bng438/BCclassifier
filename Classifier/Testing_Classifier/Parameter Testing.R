install.packages('e1071')

library(e1071)
library(readxl)

#import kNN1 for TMT_cleaned (patients are rows / proteins are columns)
kNN1 <- read_excel("knn_OG_label.xlsx",col_names = FALSE)

#import TMT staging key
staging_key <- read_excel("TMT staging key.xlsx")

#order patient names alphabetically
kNN1_ordered <- kNN1[order(kNN1[,1]),]
staging_key_ordered <- staging_key[order(staging_key[,3]),]

#Remove Not available stages
staging_key_trimmed <-subset(staging_key_ordered, Pathologic_Stage!="Not available")

#Remove kNN rows not included in staging key
kNN1_trimmed <- subset(kNN1_ordered, kNN1_ordered[[1]] %in% c(staging_key_trimmed[[3]]))
'%notin%' <- Negate('%in%')

#people in key not kNN
missing_people <- subset(staging_key_trimmed, staging_key_trimmed[[3]] %notin% c(kNN1_trimmed[[1]]))

#remove missing people
key_missingremoved <- subset(staging_key_trimmed, staging_key_trimmed[[3]] %in% c(kNN1_trimmed[[1]]))

#staging keys to integers
numbered_staging_key <- transform(key_missingremoved, id=as.numeric(factor(Pathologic_Stage)))

#staging key only numbers
numbers_key <- subset(numbered_staging_key, select = c(id))

#Convert data frame to matrix
kNN_matrix <-data.matrix(kNN1_trimmed, rownames.force = NA)
numbers_key_factor = factor(numbers_key[,1])

#classifier



iseqc <- function(from=.0000000001, to=100000000000, length.out=11){
  exp(seq(log(from),log(to),length.out = length.out))
}


svm_list <- list()
for(g in iseqc())
{
  for(c in iseqc())
  {
    svm_list[[length(svm_list)+1]] <- svm(kNN_matrix,
                                          numbers_key_factor,
                                          type="C-classification",
                                          kernel="radial",
                                          gamma=g,cost=c,cross=4)
    names(svm_list)[[length(svm_list)]] <- paste("G,C:",g,",",c,sep="")
  }
}
accuracy_list = matrix(nrow=1, ncol=121)

for(i in seq(from=1, to=121, length.out=121))
{
  accuracy_list[1,i] = svm_list[[i]][[29]]
}