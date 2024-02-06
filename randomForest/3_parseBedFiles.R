### parse test data 

library(randomForest)

load("/hpc/pmc_vanboxtel/projects/pksClassifier/1_Input/forest1.Rdata")
load("/hpc/pmc_vanboxtel/projects/pksClassifier/1_Input/forest.Rdata")


files_snv_test_all <- list.files(path="/hpc/pmc_vanboxtel/projects/pksClassifier/1_Input/Data/HMF/snvs/",pattern="merged.bed",full.names=T, recursive = T)


intervals = list()
count = 1
for(i in 1:4853){

  intervals[[i]] = seq(count,(count+3), 1)
  count = count + 4

}

for(f in 1:length(intervals)){

  files_snv_test = files_snv_test_all[intervals[[f]]]
  test.data <- list()
  val <- c()
  samples <- c()
  for ( PKS_file in files_snv_test) {
    PKS.data <- read.table(PKS_file)
    sample <- gsub(".+/(.+)_snvs.+merged.bed","\\1",PKS_file)
    feature <- gsub(".+sorted.(.+).merged.bed","\\1",PKS_file)
    if ( feature == "TSB" ) { 
      colnames(PKS.data) <- c("CHROM","START","END","NAME","SCORE","STRAND","GENE",feature)
    } else {
      colnames(PKS.data) <- c("CHROM","START","END","NAME","SCORE","STRAND",feature)  
    }
    PKS.data <- PKS.data[PKS.data$CHROM != 'M',]
  
    if ( ! sample %in% names(test.data) ) {
      test.data[[sample]] <- PKS.data
      test.data[[sample]]$SAMPLE <- sample
      test.data[[sample]]$VAL <- "PKS"
    } else{
      test.data[[sample]] <- merge(test.data[[sample]], PKS.data, by=c("CHROM","START","END","NAME","SCORE","STRAND"))
    }
  }


  test.data <- do.call("rbind",test.data)
  test.data <- test.data[!duplicated(test.data[,"START"]),]
  test.data$TSB2 <- "Unknown"
  test.data[test.data$STRAND == "+" & test.data$TSB == "+",]$TSB2 <- "Untranscribed"
  test.data[test.data$STRAND == "+" & test.data$TSB == "-",]$TSB2 <- "Transcribed"
  test.data[test.data$STRAND == "-" & test.data$TSB == "+",]$TSB2 <- "Transcribed"
  test.data[test.data$STRAND == "-" & test.data$TSB == "-",]$TSB2 <- "Untranscribed"

  df.context <- data.frame(do.call('rbind', strsplit(as.character(test.data$NAME),'',fixed=TRUE)))
  colnames(df.context) <- c(paste("POSm",c(10:1),sep=""), "POS0",  paste("POSp",c(1:10),sep=""))
  test.data <- cbind(test.data, df.context)

  test.data.sub <- test.data[,c("CHROM", "START","GENE", "GENEBODY","VAL","REPLISEQ","SIMPLEREPEAT","TSB2",paste("POSm",c(10:1),sep=""),"POS0",paste("POSp",c(1:10),sep=""))]
  if(any(test.data.sub == "N")){
    test.data.sub = test.data.sub[-unique(which(test.data.sub == "N", arr.ind = T)[,1]),]
  }
  test.data.sub$GENEBODY <- as.numeric(as.vector(test.data.sub$GENEBODY))
  test.data.sub$REPLISEQ <- as.numeric(as.vector(test.data.sub$REPLISEQ))
  test.data.sub$SIMPLEREPEAT <- as.numeric(as.vector(test.data.sub$SIMPLEREPEAT))

  test.data.sub$VAL <- as.factor(test.data.sub$VAL)
  test.data.sub$TSB2 <- as.factor(test.data.sub$TSB2)
  test.data.sub$POSm10 <- as.factor(test.data.sub$POSm10)
  test.data.sub$POSm10 = droplevels(test.data.sub$POSm10)
  test.data.sub$POSm9 <- as.factor(test.data.sub$POSm9)
  test.data.sub$POSm9 = droplevels(test.data.sub$POSm9)

  test.data.sub$POSm8 <- as.factor(test.data.sub$POSm8)
  test.data.sub$POSm8 = droplevels(test.data.sub$POSm8)
  
  test.data.sub$POSm7 <- as.factor(test.data.sub$POSm7)
  test.data.sub$POSm7 = droplevels(test.data.sub$POSm7)

  test.data.sub$POSm6 <- as.factor(test.data.sub$POSm6)
  test.data.sub$POSm6 = droplevels(test.data.sub$POSm6)

  test.data.sub$POSm5 <- as.factor(test.data.sub$POSm5)
  test.data.sub$POSm5 = droplevels(test.data.sub$POSm5)
  
  test.data.sub$POSm4 <- as.factor(test.data.sub$POSm4)
  test.data.sub$POSm4 = droplevels(test.data.sub$POSm4)

  test.data.sub$POSm3 <- as.factor(test.data.sub$POSm3)
  test.data.sub$POSm3 = droplevels(test.data.sub$POSm3)

  test.data.sub$POSm2 <- as.factor(test.data.sub$POSm2)
  test.data.sub$POSm2 = droplevels(test.data.sub$POSm2)

  test.data.sub$POSm1 <- as.factor(test.data.sub$POSm1)
  test.data.sub$POSm1 = droplevels(test.data.sub$POSm1)

  test.data.sub$POS0 <- as.factor(test.data.sub$POS0)
  
  test.data.sub$POSp10 <- as.factor(test.data.sub$POSp10)
  test.data.sub$POSp10 = droplevels(test.data.sub$POSp10)

  test.data.sub$POSp9 <- as.factor(test.data.sub$POSp9)
  test.data.sub$POSp9 = droplevels(test.data.sub$POSp9)

  test.data.sub$POSp8 <- as.factor(test.data.sub$POSp8)
  test.data.sub$POSp8 = droplevels(test.data.sub$POSp8)

  test.data.sub$POSp7 <- as.factor(test.data.sub$POSp7)
  test.data.sub$POSp7 = droplevels(test.data.sub$POSp7)

  test.data.sub$POSp6 <- as.factor(test.data.sub$POSp6)
  test.data.sub$POSp6 = droplevels(test.data.sub$POSp6)

  test.data.sub$POSp5 <- as.factor(test.data.sub$POSp5)
  test.data.sub$POSp5 = droplevels(test.data.sub$POSp5)

  test.data.sub$POSp4 <- as.factor(test.data.sub$POSp4)
  test.data.sub$POSp4 = droplevels(test.data.sub$POSp4)

  test.data.sub$POSp3 <- as.factor(test.data.sub$POSp3)
  test.data.sub$POSp3 = droplevels(test.data.sub$POSp3)

  test.data.sub$POSp2 <- as.factor(test.data.sub$POSp2)
  test.data.sub$POSp2 = droplevels(test.data.sub$POSp2)

  test.data.sub$POSp1 <- as.factor(test.data.sub$POSp1)
  test.data.sub$POSp1 = droplevels(test.data.sub$POSp1)

  
  
  predictNoM1 = predict(output.forestNoM, test.data.sub, type = "prob")[,2]
  predictNoM2 = predict(forest2NoM, test.data.sub, type = "prob")[,2]

  finalProbNoM = predictNoM1*predictNoM2
  result = cbind(test.data.sub[,c("CHROM","START", "GENE")], finalProbNoM)
  write.table(result, file = paste("/hpc/pmc_vanboxtel/projects/pksClassifier/3_Output/HMF/snvs/", sample, "_pred.txt", sep = ""), quote = FALSE, row.names= F)
}
