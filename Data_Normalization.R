#Reading an input cancer data file
data_file=read.csv("E:/Samiksha/M.Sc_Bioinformatics/SEM_3/Cancer_Genomics/Cancer_genomics_practical/cancer.csv",sep=",",header=T,row.names = 1)

#Create a count per matrix
cpmatrix=data_file
for(i in 1:ncol(data_file)){
  cpmatrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000
}

#Calculate a log of cpm
logcpm=log2(cpmatrix+1)

mat=readRDS('log2cpm_matrix.rds')

summary(logcpm)
saveRDS(logcpm, file='logcpm.rds')
#Calculate a z score
library(matrixStats)
z_score = (logcpm - rowMeans(logcpm))/rowSds(as.matrix(logcpm))[row(logcpm)]


