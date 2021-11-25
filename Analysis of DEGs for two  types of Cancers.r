library(randomForest)
#first we assign files' paths to variables to call and read it.
Lunnorpath <- "~\\lusc-rsem-fpkm-tcga_paired.txt"
Luncanpath <- "~\\lusc-rsem-fpkm-tcga-t_paired.txt"
Luncnvpath <- "~\\lusc_CNV_core.txt"
Rennorpath <- "~\\kirc-rsem-fpkm-tcga_paired.txt"
Rencanpath <- "~\\kirc-rsem-fpkm-tcga-t_paired.txt"
Rencnvpath <- "~\\kirc_CNV_core.txt"
#second, we read files and assign it into dataframes with regard to the file
#seperation type, presence of header, presence of row names.
#then, we trim the unwanted column of Gene_Id to avoid data contamination ;) using subset()
# and we transpose the CNV tabel in the process to match the format of other tables
# also, we remove any NA values if present by using the randomForest library
# using the na.roughfix function
Lunnordata <- na.roughfix(subset(read.table(Lunnorpath, sep = "\t", header = T, row.names = 1), select = -c(Entrez_Gene_Id)))
Luncandata <- na.roughfix(subset(read.table(Luncanpath, sep = "\t", header = T, row.names = 1), select = -c(Entrez_Gene_Id)))
Luncnvdata <- t(as.matrix(read.table(Luncnvpath, sep = "\t", header = T, row.names = 1)))
Rennordata <- na.roughfix(subset(read.table(Rennorpath, sep = "\t", header = T, row.names = 1), select = -c(Entrez_Gene_Id)))
Rencandata <- na.roughfix(subset(read.table(Rencanpath, sep = "\t", header = T, row.names = 1), select = -c(Entrez_Gene_Id)))
Rencnvdata <- t(as.matrix(read.table(Rencnvpath, sep = "\t", header = T, row.names = 1)))


##### LUNG Cancer Samples as PAIRED (DEPENDENT) samples ######

#Performing the stat. hypothesis testing considering paired samples
pout2 <- matrix(NA, nrow=1, ncol=2) #building an intial matrix to add on it
for(i in 1:nrow(Lunnordata)){
  #get the value of t.test and go to the variable called [['p.value']] and put it in the variable pvalue 
  tOutput = t.test(Lunnordata[i,], Luncandata[i,])
  #if P value not nan and its value is less than or equal 0.05 then store this value
  # in pout2
  if(!is.nan(tOutput[['p.value']]) && tOutput[['p.value']] <=0.05){
  
    
    pout2 <- rbind(pout2, c(rownames(Lunnordata[i,]),tOutput[['p.value']]))
    
  }
}

#remove the unwanted NA row we got when we first intiated the matrix variable
pout2 = pout2[-1,]
#perform the correction using False Discovery Rate (FDR) and
#add the p adjusted value as a column to the selected p values
pout2 <- cbind(pout2, p.adjust(pout2[ ,2], method = 'fdr'))
#select only the accepted p adjusted values (<=0.05)
LPadj = pout2[which(pout2[ ,3] <= 0.05),]
#sort the matrix according the the p adj value (least first)
LunDegsSort <- LPadj[order(LPadj[,3]),]
#select the top 5 DEGs with the least p adj. values
LunTopDeg <- LunDegsSort[1:5,]
#remove the original p values (keeping only the Gene_name and the p adj. value)
LunTopDeg <- LunTopDeg[,-2]


#####LUNG Cancer Sampels as INDEPENDANT sampels######


#Intiate a matrix with 2 col to collect the results
pout3 <- matrix(NA, nrow=1, ncol=2)
for(i in 1:nrow(Lunnordata)){
  #get the value of t.test and go to the value of the variable 
  #called 'p.value' and put it in the variable pvalue 
  tOutputind = t.test(x = Luncandata[i,], y = Lunnordata[i,], paired = FALSE, alternative = 'two.sided', var.equal = FALSE)
  #if P value not nan and its value is less than or equal 0.05 then store this value
  # in pout3
  if(!is.nan(tOutputind[['p.value']]) && tOutputind[['p.value']] <=0.05){
    
    
    pout3 <- rbind(pout3, c(rownames(Lunnordata[i,]),tOutputind[['p.value']]))
    
  }
}

#remove the unwanted NA row we got when we first intiated the matrix variable
pout3 = pout3[-1,]
#preform the correction using fdr and add the p adjusted value
#as a column to the selected p values
pout3 <- cbind(pout3, p.adjust(pout3[ ,2], method = 'fdr'))
#select only the accepted p adjusted values (<=0.05)
LPadjind = pout3[which(pout3[ ,3] <= 0.05),]
#sort the matrix according the the p adj value (least first)
LunDegsSortind <- LPadjind[order(LPadjind[,3]),]
#select the top 5 DEGs with the least p adj. values
LunTopDegind <- LunDegsSortind[1:5,]
#omit the original p values (keeping only the Gene_name and the p adj. value)
LunTopDegind <- LunTopDegind[,-2]


####  REGRESSION with CNA for PAIRED LUNG cancer samples ####

#match the patients' codes format 
colnames(Luncandata) <- gsub('.', '-', colnames(Luncandata), fixed=TRUE)
#collect the intesected patients' codes
Lccn <- intersect(colnames(Luncnvdata), colnames(Luncandata))
#fetch the whole CNAs of the intersected patients
Lcnv <- Luncnvdata[ ,Lccn]
#fetch the whole GE (in case of cancer) of the intersected patients
Lcan <- Luncandata[ ,Lccn]
#short listing only the Top 5 highly expressed Genes of the latter list 
LcanTop <- Lcan[LunTopDeg[,1],]

#converting CNAs of the intersected pts to atomic vectors to perform lm
  y <- matrix (NA, nrow=1, ncol=13)

  for (o in 1:nrow(LcanTop)){
  
  y <- rbind( y , unlist(LcanTop[o,], recursive = TRUE, use.names = FALSE))
  }
  y <- y[-1,]
#converting the Top 5 highly exp. genes of the intersected pts to atomic vectors to perform lm
  x <- matrix (NA, nrow=1, ncol=13)
  
  for (p in 1:nrow(Lcnv)){
    
    x <- rbind( x , unlist(Lcnv[o,], recursive = TRUE, use.names = FALSE))
  }
  x <- x[-1,]

#perform Linear regression bet. every highly exp. gene and the whole CNAs 
#to check the causality
LtotReg <- list()
LtotReg1 <- list()
 for (i in 1:nrow(LcanTop)){
  
   for (j in 1:nrow(Lcnv))
   {
     Lunmodel = lm(y[i,] ~ x[j,])
    
     LtotReg1[j] <- summary(Lunmodel)
     #LtotReg1 <- rbind(LtotReg1, Lunmodel)
   }
  
  LtotReg <- rbind(LtotReg, LtotReg1)
}


#####Renal Cancer Samples as PAIRED samples######

#Performing the stat. hypothesis testing considering paired samples
pout4 <- matrix(NA, nrow=1, ncol=2)
for(i in 1:nrow(Rennordata)){
  #get the value of t.test and go to the variable called [['p.value']] and put it in the variable pvalue 
  tOutput = t.test(Rennordata[i,], Rencandata[i,])
  #if P value not nan and its value is less than or equal 0.05 then store this value
  # in pout4
  if(!is.nan(tOutput[['p.value']]) && tOutput[['p.value']] <=0.05){
    
    
    pout4 <- rbind(pout4, c(rownames(Rennordata[i,]),tOutput[['p.value']]))
    
  }
}

#remove the unwanted NA row we got when we first initiated the matrix variable
pout4 = pout4[-1,]
#preform the correction using fdr and add the p adjusted value
#as a column to the selected p values
pout4 <- cbind(pout4, p.adjust(pout4[ ,2], method = 'fdr'))
#select only the accepted p adjusted values (<=0.05)
RPadj = pout4[which(pout4[ ,3] <= 0.05),]
#sort the matrix according the the p adj value (least first)
RenDegsSort <- RPadj[order(RPadj[,3]),]
#select the top 5 DEGs with the least p adj. values
RenTopDeg <- RenDegsSort[1:5,]
#remove the original p values (keeping only the Gene_name and the p adj. value)
RenTopDeg <- RenTopDeg[,-2]


#####RENAL Cancer Sampels as INDEPENDANT sampels######


#Intiate a matrix with 2 col to collect the results
pout5 <- matrix(NA, nrow=1, ncol=2)
for(i in 1:nrow(Rennordata)){
  #get the value of t.test and go to the value of the variable 
  #called 'p.value' and put it in the variable pvalue 
  tOutputind = t.test(x = Rencandata[i,], y = Rennordata[i,], paired = FALSE, alternative = 'two.sided', var.equal = FALSE)
  #if P value not nan and its value is less than or equal 0.05 then store this value
  # in pout5
  if(!is.nan(tOutputind[['p.value']]) && tOutputind[['p.value']] <=0.05){
    
    
    pout5 <- rbind(pout5, c(rownames(Rennordata[i,]),tOutputind[['p.value']]))
    
  }
}

#remove the unwanted NA row we got when we first intiated the matrix variable
pout5 = pout5[-1,]
#preform the correction using fdr and add the p adjusted value
#as a column to the selected p values
pout5 <- cbind(pout5, p.adjust(pout5[ ,2], method = 'fdr'))
#select only the accepted p adjusted values (<=0.05)
RPadjind = pout5[which(pout5[ ,3] <= 0.05),]
#sort the matrix according the the p adj value (least first)
RenDegsSortind <- RPadjind[order(RPadjind[,3]),]
#select the top 5 DEGs with the least p adj. values
RenTopDegind <- RenDegsSortind[1:5,]
#omit the original p values (keeping only the Gene_name and the p adj. value)
RenTopDegind <- RenTopDegind[,-2]


####  REGRESSION with CNA for PAIRED RENAL cancer samples ####

#match the patients' codes format 
colnames(Rencandata) <- gsub('.', '-', colnames(Rencandata), fixed=TRUE)
#collect the intesected patients' codes
Rccn <- intersect(colnames(Rencnvdata), colnames(Rencandata))
#fetch the whole CNAs of the intersected patients
Rcnv <- Rencnvdata[ ,Rccn]
#fetch the whole GE (in case of cancer) of the intersected patients
Rcan <- Rencandata[ ,Rccn]
#short listing only the Top 5 highly expressed Genes of the latter list 
RcanTop <- Rcan[RenTopDeg[,1],]

#converting CNAs of the intersected pts to atomic vectors to perform lm
b <- matrix (NA, nrow=1, ncol=ncol(RcanTop))

for (o in 1:nrow(RcanTop)){
  
  b <- rbind( b , unlist(RcanTop[o,], recursive = TRUE, use.names = FALSE))
}
b <- b[-1,]
#converting the Top 5 highly exp. genes of the intersected pts to atomic vectors to perform lm
a <- matrix (NA, nrow=1, ncol=ncol(Rcnv))

for (p in 1:nrow(Rcnv)){
  
  a <- rbind( a , unlist(Rcnv[o,], recursive = TRUE, use.names = FALSE))
}
a <- a[-1,]

#perform Linear regression bet. every highly exp. gene and the whole CNAs 
#to check the causality
RtotReg <- list()
RtotReg1 <- list()
for (i in 1:nrow(RcanTop)){
  
  for (j in 1:nrow(Rcnv))
  {
    Renmodel = lm(b[i,] ~ a[j,])
    
    RtotReg1[j] <- Renmodel
    #LtotReg1 <- rbind(LtotReg1, Lunmodel)
  }
  
  RtotReg <- rbind(RtotReg, RtotReg1)
}

###########################################

#writing output Files:

LunDegsSort0 <- LunDegsSort
RenDegsSort0 <- RenDegsSort
LunTopDeg0 <- LunTopDeg
RenTopDeg0 <- RenTopDeg
colnames(LunDegsSort0) <- c('Gene_name', 'p_value', 'p_adj_value')
colnames(RenDegsSort0) <- c('Gene_name', 'p_value', 'p_adj_value')
colnames(LunTopDeg0) <- c('Gene_name', 'p_adj_value')
colnames(RenTopDeg0) <- c('Gene_name', 'p_adj_value')

write.table(LunDegsSort0, '~\\Lung_DEGs_sorted.txt', sep = "\t", row.names = F, col.names = T, quote =F)
write.table(RenDegsSort0, '~\\Renal_DEGs_sorted.txt', sep = "\t", row.names = F, col.names = T, quote =F)
write.table(LunTopDeg0, '~\\Lung_Top_DEGs.txt', sep = "\t", row.names = F, col.names = T, quote =F)
write.table(RenTopDeg0, '~\\Renal_Top_DEGs.txt', sep = "\t", row.names = F, col.names = T, quote =F)


#Plotting the Box Plot for the Top 5 DEGs of the two Cancers
boxplot(y[1,], y[2,], y[3,], y[4,], y[5,], at=c(1,2,3,4,5), main = 'Top DEGs in Lung Cancer', names = c(rownames(LcanTop[])), xlab = 'Genes', ylab = 'Relative Expression',ylim = c(min(min(y[1,]),min(y[2,]),min(y[3,]),min(y[4,]),min(y[5,]))-1000, max(max(y[1,]),max(y[2,]),max(y[3,]),max(y[4,]),max(y[5,]))))

boxplot(b[1,], b[2,], b[3,], b[4,], b[5,], at=c(1,2,3,4,5), main = 'Top DEGs in Renal Cancer', names = c(rownames(RcanTop[])), xlab = 'Genes', ylab = 'Relative Expression',ylim = c(min(min(b[1,]),min(b[2,]),min(b[3,]),min(b[4,]),min(b[5,]))-100, max(max(b[1,]),max(b[2,]),max(b[3,]),max(b[4,]),max(b[5,]))))




#Geting CV for DEGs of the two cancers with diff. no. of cases
LunDEGs <- Luncandata[LPadj[,1],]
RenDEGs <- Rencandata[RPadj[,1],]
LunCV = sd(as.matrix(LunDEGs))/mean(as.matrix(LunDEGs))
RenCV = sd(as.matrix(RenDEGs))/mean(as.matrix(RenDEGs))


#Bar Plot for Lung
par(mfrow = c(1,2))
barplot(c(nrow(pout2), nrow(Luncandata)), main = 'Lung Cancer Genes', xlab = 'Type of Genes', ylab = 'No. of Genes', ylim = c(0, nrow(Luncandata)), names=c('Total DEGs', 'Total tested genes'))
abline(h=0)
barplot(c(nrow(pout2), nrow(LPadj)), xlab = 'Type of Genes', ylab = 'No. of Genes', ylim = c(0, nrow(Luncandata)), names=c('Total DEGs', 'Corrected DEGs'))
abline(h=0)

#Bar Polt for Renal
par(mfrow = c(1,2))
barplot(c(nrow(pout4), nrow(Rencandata)), main = '        Renal Cancer Genes', xlab = 'Type of Genes', ylab = 'No. of Genes', ylim = c(0, nrow(Rencandata)), names=c('Total DEGs', 'Total tested genes'))
abline(h=0)
barplot(c(nrow(pout4), nrow(RPadj)), xlab = 'Type of Genes', ylab = 'No. of Genes', ylim = c(0, nrow(Rencandata)), names=c('Total DEGs', 'Corrected DEGs'))
abline(h=0)

