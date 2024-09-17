###############################################################################
#          Spectra preprocessing and Partial least squares model                                      
#           Rafaela P. Graciano, 2024                                           
###############################################################################

rm(list=ls())


##### Load packages and data #####   
##################################


### Load packages

require(prospectr); require(tidyverse); require(AGHmatrix); require(BGLR); require(dismo); require(pls)

### Load data ##
#NIRS
S19.NIRS= read.csv("NIRS_S19.csv")
S20.NIRS= read.csv("NIRS_S20.csv")

#Phenotypic
C21.Traits <- read.csv("BLUE_C21.csv")
C22.Traits <- read.csv("BLUE_C22.csv")
C19.Traits= read.csv("BLUE_C19.csv")

##### 1. Preprocessing NIRS #########
#####################################
#Average across years
S19_20.NIRS= rbind(S20.NIRS,S19.NIRS)
S19_20.NIRS= aggregate(x = S19_20.NIRS[,2:771], by = list(S19_20.NIRS$Genotype),mean, na.rm = TRUE)

## Tranforming data for next steps
N.19= as.matrix(S19.NIRS[,c(2:771)])
row.names(N.19)= S19.NIRS$Genotype

N.20=as.matrix(S20.NIRS[,2:771])
rownames(N.20)= S20.NIRS$Genotype

N.19_20= as.matrix(S19_20.NIRS[,2:771])
rownames(N.19_20)= S19_20.NIRS$Group.1

#### Plotting raw spectra
wavelengths<-seq(910,1679,by=1)
# matplot(wavelengths,t(N.19_20),lty=1,pch=","
#         ,xlab="Wavelength (nm)",ylab="Absorbance",
#         main="Raw Spectra"
# )

#### Standard normal variate (center and scale)
N19.SNV=t(scale(t(N.19)))
N20.SNV=t(scale(t(N.20)))
N19_20.SNV= t(scale(t(N.19_20)))
### 1st derivative    
N19.SNV.1D <- savitzkyGolay(N19.SNV, m=1, p=2, w=37) 
N20.SNV.1D <- savitzkyGolay(N20.SNV, m=1, p=2, w=37) 
N19_20.SNV.1D <- savitzkyGolay(N19_20.SNV, m=1, p=2, w=37)  

##### 2.formating data for pls  #########
#########################################################
####choose NIRS data source and phenotypic year data 
NIRS20_P21<-inner_join(C21.Traits, data.frame(Genotype= rownames(N20.SNV.1D),N20.SNV.1D)) #traits from 2021 and NIRS 2020
data_N20_P21<-data.frame(NIRS20_P21[,c(1:21)], NIR = I(NIRS20_P21[22:755])) #

##If predictive ability of pls will be compared with PS remember to add a step to match individuals used in both analysis for a fare comparison.

######### 3. Prediction models  ###########
###########################################
nFold = 5
s=20
Traits= colnames(data_N20_P21[,-c(1,22)])##choosing traits
accuracy_intervals = matrix(NA,s,length(Traits))

# Transforming in factor
data_N20_P21= transform(data_N20_P21, Genotype= as.factor(Genotype))
for (t in 1:length(Traits)) {
  
  #### Number of components#####
  comp <- plsr(eval(sym(Traits[t])) ~ as.matrix(data_N20_P21$NIR), ncomp = 25, data = data_N20_P21, validation ="CV") 
  ncomp.permut <- selectNcomp(comp,alpha = 0.05, method = "onesigma", plot = F) ## choosing the number of components with cv and permutation method
  if(ncomp.permut==0){
    ncomp.permut=1 ## Some traits resulted in a "optimal" number of 0 and we used 1 in this cases
  }
  
  geno = data_N20_P21$Genotype
  predictions = matrix(NA,length(geno),s)
  
  for (k in 1:s) {
    seed_value = as.integer(paste("2001",k,sep=""))     # Fix the seed to each iterations
    set.seed(seed_value)                      # Setting the seed
    kGroup = kfold(length(geno), k=5)
    # model
    for(fold in 1:nFold){
      cat("fold number:" , fold, "\n")
      val_geno = geno[which(kGroup==fold)] ## testing individuals
      
      tst = data_N20_P21[which(data_N20_P21$Genotype%in%val_geno),] #testing 
      trn = data_N20_P21[which(!data_N20_P21$Genotype%in%val_geno),] #training 
      
      train.pls <- plsr(eval(sym(Traits[t])) ~ as.matrix(NIR), ncomp = ncomp.permut, data = trn, validation ="CV")
      #summary(train.pls)
      pred = predict(train.pls, tst, ncomp=ncomp.permut)
      predictions[which(data_N20_P21$Genotype%in%val_geno),k] = as.matrix(pred)
      
    }
    # Geno model
    accuracy_intervals[k,t]= cor(predictions[,k],data_N20_P21[,t+1], use="complete")
    
  }
}  
colnames(accuracy_intervals)= Traits

}

accuracy_mean= data.frame( Traits= Traits, colMeans(accuracy_intervals))
write.csv(accuracy_means, "test_accuracy_means.csv")
save( write.csv(accuracy_intervals, "accuracy_intervals_pls.csv")
  
