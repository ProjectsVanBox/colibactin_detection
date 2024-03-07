#### 4A, comparison of probability in driver genes 
driverGenes = c("APC", "TP53", "KRAS", "BRAF", "PIK3CA", "SMAD4", "FBXW7", "TCF7L2","FAT4", "ATM")

#### data frame with prediction by random Forest
load("colonCancerScore.Rdata")

### get sampleNames of patients with colibactin-induced damage or not
posSamples = colnames(colonCancerScore)[which(colonCancerScore[1,] > 0.1)]
negSamples = colnames(colonCancerScore)[which(colonCancerScore[1,] <= 0.1)]


### make objects for comparison
diffPerGene = matrix(NA, 2,length(driverGenes))
set.seed(220718)

pValuesDriver = rep(NA, length(driverGenes))
pValuesDriver_neg = rep(NA, length(driverGenes))
sampledPos = list()
sampledNeg = list()
probPos = list()
probNeg = list()
differencesPos = list()
differencesNeg = list()

#### get non-driver genes to sample
genesToSample = intersect(names(which(table(posPatients$GENE) > 0 & table(posPatients$GENE) < 4)),names(which(table(negPatients$GENE) >0  & table(negPatients$GENE) < 20)) )

randomGenes = sample(genesToSample, 1000, replace = F)
meanPos = rep(NA, length(randomGenes))
meanNeg = rep(NA, length(randomGenes))

for(s in 1:length(randomGenes)){
  
  meanPos[s] = mean(posPatients$finalProbNoM[which(posPatients$GENE == randomGenes[s])], na.rm = T)
  meanNeg[s] = mean(negPatients$finalProbNoM[which(negPatients$GENE == randomGenes[s])], na.rm = T)
  
}

diffRandom = meanPos - meanNeg


#### compare probability for driver genes in positive and negative sample
for(g in 1:length(driverGenes)){
    diffPos = rep(NA,length(which(posPatients$GENE == driverGenes[g])))
    diffNeg = rep(NA,length(which(negPatients$GENE == driverGenes[g])))
    
    
    
    meanPosUse = mean(meanPos, na.rm = T)
    meanNegUse = mean(meanNeg, na.rm = T)
    
    probPos[[g]] = posPatients$finalProbability[which(posPatients$GENE == driverGenes[g])]
    probNeg[[g]] = negPatients$finalProbability[which(negPatients$GENE == driverGenes[g])]
    
    for(i in 1:length(which(posPatients$GENE == driverGenes[g]))){
      diffPos[i] = posPatients$finalProbability[which(posPatients$GENE == driverGenes[g])[i]] - meanPosUse
    }
    
    for(i in 1:length(which(negPatients$GENE == driverGenes[g]))){
      diffNeg[i] = negPatients$finalProbability[which(negPatients$GENE == driverGenes[g])[i]] - meanNegUse
    }
    
  
    differencesPos[[g]] = diffPos
    differencesNeg[[g]] = diffNeg
    
    
    diffPerGene[1,g] = mean(differencesPos[[g]])
    diffPerGene[2,g] = mean(differencesNeg[[g]])
    
    
}


#### calculatate significance of difference

pValues_pos = rep(NA, length(driverGenes))
pValues_neg = rep(NA, length(driverGenes))


for(i in 1:length(driverGenes)){
    pValues_t[i] = t.test(posPatients$finalProb[which(posPatients$GENE == driverGenes[i])], meanPos)$p.value
    pValues_neg[i] = t.test(negPatients$finalProb[which(negPatients$GENE == driverGenes[i])], meanNeg)$p.value
    
}

names(pValues_t) = driverGenes

diffPerGeneT = matrix(NA, 2, length(driverGenes))

for(i in 1:length(driverGenes)){
  diffPerGeneT[1,i] = mean(diffPerGene[1,i], na.rm = T)
  diffPerGeneT[2,i] = mean(diffPerGene[2,i], na.rm = T)
}

diffPerGeneT_sd = matrix(NA, 2, length(driverGenes))

for(i in 1:length(driverGenes)){
  diffPerGeneT_sd[1,i] = sd(diffPerGene[1,i], na.rm = T)
  diffPerGeneT_sd[2,i] = sd(diffPerGene[2,i], na.rm = T)
}


####  make a boxplot of the probability
library(ggplot2)

genes = NULL
label = NULL
prob = NULL

diffProb = NULL

for(i in 1:length(driverGenes)){
  diffProb[i] = mean(probPos[[i]]) - mean(probNeg[[i]])
}

driverGenes_ordered = driverGenes[sort(diffProb, decreasing = T, index = T)$ix]
for(i in 1:length(driverGenes_ordered)){
  genes = c(genes, rep(driverGenes_ordered[i], length(which(posPatients$GENE == driverGenes_ordered[i]))), rep(driverGenes_ordered[i], length(which(negPatients$GENE == driverGenes_ordered[i]))))
  label = c(label, rep("Positive", length(which(posPatients$GENE == driverGenes_ordered[i]))), rep("Negative",length(which(negPatients$GENE == driverGenes_ordered[i]))))
  prob = c(prob, posPatients$finalProbNoM[which(posPatients$GENE == driverGenes_ordered[i])], negPatients$finalProbNoM[which(negPatients$GENE == driverGenes_ordered[i])])
  
}

genes = c(genes[1:1403], rep("Random", 2000), genes[1404:2720])
label = c(label[1:1403], rep("Positive", 1000), rep("Negative", 1000), label[1404:2720])
prob = c(prob[1:1403], meanPos, meanNeg, prob[1404:2720])

level_order = factor(c(driverGenes_ordered[1:6], "Random", driverGenes_ordered[7:10]))
 dataToPlot = data.frame(genes, label, prob)

 pdf("Boxplot_genes_prob.pdf", width = 12)
 ggplot(dataToPlot, aes(x=factor(genes, level = level_order), y=prob, fill=label)) + ylab("Posterior probability") + xlab(NA)+
   geom_boxplot() +  scale_fill_manual(values = c( "azure3","aquamarine3")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 16, angle = 90))
 dev.off()
 
 
#### B, plot probability distribution
  
pdf("Density_APC.pdf", height = 2.5)
ggplot(dataToPlot[which(dataToPlot$genes == "APC"),], aes(x=prob, fill=label)) +
  geom_density(alpha = 0.5)  + scale_fill_manual(values = c( "azure3", "aquamarine3")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0,1)
dev.off()

pdf("Density_TP53.pdf", height = 2.5)
ggplot(dataToPlot[which(dataToPlot$genes == "TP53"),], aes(x=prob, fill=label)) +
  geom_density(alpha = 0.5)  + scale_fill_manual(values = c( "azure3", "aquamarine3")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0,1)
dev.off()


#### C, D, E plot boxplot for age at biopsy/diagnosis

pksStatus = rep("Negative", ncol(colonCancerScore))
pksStatus[which(colonCancerScore[1,] > 0.1)] = "Positive"

age_PKSstatus = data.frame(age, pksStatus )

pdf("Age_status.pdf", width = 4)
ggplot(age_PKSstatus, aes(x = age, y = pksStatus, fill = pksStatus)) + geom_violin()  + geom_boxplot(width = 0.1, fill = c("azure2", "aquamarine2"))+ coord_flip() +
  scale_fill_manual(values = c( "azure3","aquamarine3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 16), text = element_text(size = 16))
dev.off()
 
#### F

predict1 = predict(output.forestNoM, test.data, type = "prob")[,2]
predict2 = predict(forest2NoM, test.data, type = "prob")[,2]

finalProb = predict1*predict2


pksClass = rep("noPKS", nrow(test.data))

pksClass[which(finalProb > 0.5)] = "PKS"


test.data$pks = pksClass
test.data$prob = finalProb

ratios = NULL


for(i in 1:length(samples)){
  ratios[i] = table(test.data$pks[which(newNames== samples[i])])[2]/length(which(newNames == samples[i]))
}

names(ratios) = samples

pksStatus = rep("Negative", length(samples))
pksStatus[which(ratios > 0.1)] = "Positive"


#### clinical Info from LeeSix et al., Nature (2019)
rownames(clinInfo) = clinInfo$patient_label_in_files
clinInfo = clinInfo[names(ratios),]
clinInfo$ratio = ratios

clinInfoUse = clinInfo[which(clinInfo$cohort == "bowel_cancer_screening_programme_cohort"),]


pdf("ageVSratio_screeningCohort.pdf")
plot(clinInfoUse$ratio, clinInfoUse$age, bty = "n", lwd = 2, xlim = c(0,0.4),  xlab = "Fraction colibactin-induced", ylab = "Age", cex.lab = 1.5, cex.axis = 1.5)
dev.off()



 
