###############################################################################
#          Spectra preprocessing and phenomic prediction 
#           Rafaela P. Graciano, 2024                                            
###############################################################################

rm(list=ls())


##### Load packages and data #####   
##################################

### Load packages

require(prospectr); require(tidyverse); require(AGHmatrix); require(BGLR); require(dismo)

### Load data ##
#NIRS
S19.NIRS= read.csv("NIRS_S19.csv")
S20.NIRS= read.csv("NIRS_S20.csv")
DH.NIRS= read.csv("NIRS_DH.csv")

#Phenotypic
C22.Traits <- read.csv("BLUE_C22.csv") ## choose the dataset to calibrate the model

##### 1. Preprocessing NIRS #########
#####################################

#Average across years
S19_20.NIRS= rbind(S20.NIRS,S19.NIRS)
S19_20.NIRS= aggregate(x = S19_20.NIRS[,2:771], by = list(S19_20.NIRS$Genotype),mean, na.rm = TRUE)
S19_20.NIRS=S19_20.NIRS%>%rename(Genotype=Group.1)

##combining data sets 
NIRS.comb=rbind(S19_20.NIRS,DH.NIRS[,c(1,8:777)])

## Tranforming data for next steps
N= as.matrix(NIRS.comb[,c(2:771)])
row.names(N)= NIRS.comb$Genotype

#### Plotting raw spectra
wavelengths<-seq(910,1679,by=1)
# matplot(wavelengths,t(N),lty=1,pch=","
#         ,xlab="Wavelength (nm)",ylab="Absorbance",
#         main="Raw Spectra"
# )


#### Standard normal variate (center and scale)
N.SNV=t(scale(t(N)))
### 1st derivative    
N.SNV.1D <- savitzkyGolay(N.SNV, m=1, p=2, w=37) 

##### 2. computing  matrix  #########
#########################################################
###Phenomic
spct <- N.SNV.1D       # Change for other pre-treatments and kernel source
spct2 <- scale(spct, center=T, scale=T) # scale absorbance at each wavelength 
P_matrix<- tcrossprod(spct2)/ncol(spct2)
#superheat(P_matrix , legend = FALSE, title = "SNV") 


######### 3. Prediction models  ###########
###########################################
data=C22.Traits  ## running for blues from 2022, change to run other years
### Matching individuals  ####
Match.names=as.factor(Reduce(intersect, list(row.names(P_matrix[1:675,]), unique(data$Genotype)))) 

## filtering and ordering
data.1=data[which(data$Genotype%in%Match.names),]
data.1= data.1[order(data.1$Genotype),]
n.inb= length(data.1$Genotype) ##number of lines in the diversity panel


n= which(!rownames(P_matrix[1:675,])%in%Match.names)
P_matrix.1=P_matrix[-n,-n]
###
phenos <- colnames(data.1)[-1]  ##choosing traits

# Transforming in factor
data.1= transform(data.1, Genotype= as.factor(Genotype))

### parameters
n_id = length(data.1$Genotype)
n_model = 1 ## number of models
n_dh=length(P_matrix.1[1,])-n.inb #number of dhs
n_total= length(P_matrix.1[1,])
y_pred_PS = matrix(NA,n_dh,1)
y_pred= list()

for(p in 1:length(phenos)){
  y.trait= c(data.1[,p+1],rep(NA,n_dh))
  
      ##PS
      mod1 = BGLR(y = y.trait,
                  ETA=list(G=list(K=P_matrix.1,model='RKHS')),
                  nIter = 10000,
                  burnIn = 1000,
                  thin = 10)
      y_pred_PS = mod1$yHat[(n.inb+1):n_total]# getting prediction for dhs 
      
      y_pred[[p]]=data.frame(Genotype= colnames(P_matrix.1)[(n.inb+1):n_total],
                             y_pred=y_pred_PS )
      

  
}

names(y_pred)= phenos

write.csv(y_pred, "BLUPS_dh.csv")
