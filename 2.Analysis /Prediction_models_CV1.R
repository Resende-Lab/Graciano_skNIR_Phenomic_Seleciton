###############################################################################
#          Spectra preprocessing and prediction models                                       
#           Rafaela P. Graciano, 2024                                           
###############################################################################

rm(list=ls())


##### Load packages and data #####   
##################################


### Load packages

require(prospectr); require(tidyverse); require(AGHmatrix); require(BGLR); require(dismo)

### Load data ##
#Genomic
load("SNP_data.Rdata")

#NIRS
S19.NIRS= read.csv("NIRS_S19.csv")
S20.NIRS= read.csv("NIRS_S20.csv")

#Phenotypic
C21.Traits <- read.csv("BLUE_BLUP_C21.csv")
C22.Traits <- read.csv("BLUE_BLUP_C22.csv")
C19.Traits= read.csv("BLUE_BLUP_C19.csv")
colnames(C19.Traits) =c( "Genotype",sub("^.*_(.*)_.*$", "\\1", colnames(C19.Traits[,2:13])))
colnames(C21.Traits) =c( "Genotype",sub("^.*_(.*)_.*$", "\\1", colnames(C21.Traits[,2:21])))
colnames(C22.Traits) =c( "Genotype",sub("^.*_(.*)_.*$", "\\1", colnames(C22.Traits[,2:24])))

##### 1. Preprocessing NIRS #########
#####################################

# Average spectra per genotype 
S19.NIRS= S19.NIRS %>% filter(outliers!=1)
S19.NIRS.avg=aggregate(x = S19.NIRS[,4:773], by = list(S19.NIRS$genotype),mean, na.rm = TRUE)

S20.NIRS= S20.NIRS %>% filter(outliers!=1)
S20.NIRS.avg=aggregate(x = S20.NIRS[,4:773], by = list(S20.NIRS$Genotype),mean, na.rm = TRUE)

#Average across years
S19_20.NIRS= rbind(S20.NIRS.avg,S19.NIRS.avg)
S19_20.NIRS= aggregate(x = S19_20.NIRS[,2:771], by = list(S19_20.NIRS$Group.1),mean, na.rm = TRUE)


## Tranforming data for next steps
N.19= as.matrix(S19.NIRS.avg[,c(2:771)])
row.names(N.19)= S19.NIRS.avg$Group.1

N.20=as.matrix(S20.NIRS.avg[,2:771])
rownames(N.20)= S20.NIRS.avg$Group.1

N.19_20= as.matrix(S19_20.NIRS[,2:771])
rownames(N.19_20)= S19_20.NIRS$Group.1

#### Plotting raw spectra
wavelengths<-seq(910,1679,by=1)
# matplot(wavelengths,t(N.19_20),lty=1,pch=","
#         ,xlab="Wavelength (nm)",ylab="Absorbance",
#         main="Raw Spectra"
# )


###### Matrix Normalization #######
#### Multiple Scatter Plot (MSC)
N19.msc=msc(N.19)
N20.msc=msc(N.20)
N19_20.msc=msc(N.19_20)

#### Standard normal variate (center and scale)
N19.SNV=t(scale(t(N.19)))
N20.SNV=t(scale(t(N.20)))
N19_20.SNV= t(scale(t(N.19_20)))

##detrand
N19.dt= detrend(X =N.19, wav = wavelengths)
N20.dt= detrend(X =N.20, wav = wavelengths)
N19_20.dt= detrend(X =N.19_20, wav = wavelengths)

### 1st derivative    
N19.SNV.1D <- savitzkyGolay(N19.SNV, m=1, p=2, w=37) 
N20.SNV.1D <- savitzkyGolay(N20.SNV, m=1, p=2, w=37) 
N19_20.SNV.1D <- savitzkyGolay(N19_20.SNV, m=1, p=2, w=37)  

N19.MSC.1D <- savitzkyGolay(N19.msc, m=1, p=2, w=37) 
N20.MSC.1D <- savitzkyGolay(N20.msc, m=1, p=2, w=37) 
N19_20.MSC.1D <- savitzkyGolay(N19_20.msc, m=1, p=2, w=37)  

#### spectra treatments into a list ####

spectra19 <- list("raw" = N.19, "norm_SNV" = N19.SNV, "norm_MSC" = N19.msc,
                  "norm_der1" = N19.SNV.1D,"msc_der1" = N19.MSC.1D,
                  "dt"= N19.dt )
spectra20 <- list("raw" = N.20, "norm_SNV" = N20.SNV, "norm_MSC" = N20.msc,
                  "norm_der1" = N20.SNV.1D, "msc_der1" = N20.MSC.1D,
                  "dt"= N20.dt)
spectra19_20 <- list("raw" = N.19_20, "norm_SNV" = N19_20.SNV, "norm_MSC" = N19_20.msc,
                     "norm_der1" = N19_20.SNV.1D,  "msc_der1" = N19_20.MSC.1D,
                     "dt"= N19_20.dt)

rm(N.19,N19.SNV,N19.msc,N19.SNV.1D,N.20,N20.SNV,N20.msc,
   N20.SNV.1D,N.19_20,N19_20.SNV,N19_20.msc,N19_20.SNV.1D,N19_20.MSC.1D,N19.MSC.1D)



##### 2. computing phenomic and genomic matrix  #########
#########################################################
###Phenomic
spct <- spectra19_20$norm_SNV         # Change for other pre-treatments and kernel source
spct2 <- scale(spct, center=T, scale=T) # scale absorbance at each wavelength 
P_matrix<- tcrossprod(spct2)/ncol(spct2)
#superheat(P_matrix , row.dendrogram = T,col.dendrogram = T, legend = FALSE, title = "SNV") 

#####Genomic
G_matrix <- Gmatrix( snp_mat1, missingValue = "-1",
                 thresh.missing = .3,
                 maf = 0.01)



######### 3. Prediction models  ###########
###########################################
data=C21.Traits  ## running for blues from 2021, change to run other years
### Matching individuals  ####
Match.names=as.factor(Reduce(intersect, list(row.names(G_matrix), row.names(P_matrix), unique(data$Genotype)))) 

## filtering and ordering
data.1=data[which(data$Genotype%in%Match.names),]
data.1= data.1[order(data.1$Genotype),]

G_matrix.1= G_matrix[which(rownames(G_matrix)%in%Match.names), which(rownames(G_matrix)%in%Match.names)]
G_matrix.1=G_matrix.1[order(colnames(G_matrix.1)),order(colnames(G_matrix.1))]

P_matrix.1=P_matrix[which(rownames(P_matrix)%in%Match.names), which(rownames(P_matrix)%in%Match.names)]
P_matrix.1=P_matrix.1[order(colnames(P_matrix.1)),order(colnames(P_matrix.1))]

###
phenos <- colnames(data.1)[2:21]  ##choosing traits

# Transforming in factor
data.1= transform(data.1, Genotype= as.factor(Genotype))

### parameters
n_id = length(data.1$Genotype)
n_model = 3 ## number of models, change if testing more pre-treatments 
s=20  #reps for the cross validation
y_pred_PS = matrix(NA,n_id,s)
y_pred_GS = matrix(NA,n_id,s)
y_pred_PG= matrix(NA,n_id,s)
accuracy_intervals= list()
MSE_intervals= list()

for(p in 1:length(phenos)){
  y.trait= data.1[,p+1]

  for(k in 1:s){
    ## Sets: ####
    seed_value = as.integer(paste("2001",k,sep=""))     # Fix the seed to each iterations
    set.seed(seed_value)                      # Setting the seed
    cv= kfold(n_id, k=5)                # Obtaining cross-validation subsets

    ## Models: ####
    for(i in 1:5){  ##for all sets

      ### defining test and training population
      trn = y.trait
      trn[cv==i]=NA

      ##PS
      mod1 = BGLR(y = trn,
                  ETA=list(G=list(K=P_matrix.1,model='RKHS')),
                  nIter = 10000,
                  burnIn = 1000,
                  thin = 10)
      y_pred_set<-mod1$yHat
      y_pred_PS[which(cv==i),k] = y_pred_set[which(cv==i)] # getting prediction for testing population
  

      ## Genomic
      mod2 = BGLR(y = trn,
                  ETA=list(G=list(K=G_matrix.1,model='RKHS')),
                  nIter = 10000,
                  burnIn = 1000,
                  thin = 10)
      y_pred_set<-mod2$yHat
      y_pred_GS[which(cv==i),k] = y_pred_set[which(cv==i)]

      ## Genomic + phenomic
            mod3 = BGLR(y = trn,
                        ETA=list(G=list(K=G_matrix.1,model='RKHS'),
                                 P=list(K=P_matrix.1,model='RKHS')),
                        nIter = 10000,
                        burnIn = 1000,
                        thin = 10)
            y_pred_set<-mod3$yHat
            y_pred_PG[which(cv==i),k] = y_pred_set[which(cv==i)]
    }
  }

  for(k in 1:s){
    if(k==1){
      accuracy_intervals[[p]] = matrix(NA,s,n_model)
      MSE_intervals[[p]]=matrix(NA,s,n_model)
      colnames(accuracy_intervals[[p]]) <- c("PS","GS","G+P")
      colnames(MSE_intervals[[p]]) <- c("PS","GS","G+P")}

    ##accurracy
    accuracy_intervals[[p]][k,1] = cor(y_pred_PS[,k],y.trait,use="complete")
    accuracy_intervals[[p]][k,2] = cor(y_pred_GS[,k],y.trait,use="complete")
    accuracy_intervals[[p]][k,3] = cor(y_pred_PG[,k],y.trait,use="complete")
    ### Means squared error
    MSE_intervals[[p]][k,1] = mean((y.trait - y_pred_PS[,k])^2,na.rm=T)
    MSE_intervals[[p]][k,2] = mean((y.trait - y_pred_GS[,k])^2,na.rm=T)
    MSE_intervals[[p]][k,3] = mean((y.trait - y_pred_PG[,k])^2,na.rm=T)

 }
 }

##
accuracy_means= matrix(NA,p,n_model)

for(p in 1:length(phenos)){
  accuracy_means[p,]= colMeans(accuracy_intervals[[p]])
}
colnames(accuracy_means)= c("PS","GS","G+P")
accuracy_means= data.frame( Traits= phenos, accuracy_means)

##
MSE_means= matrix(NA,p,n_model)

for(p in 1:length(phenos)){
  MSE_means[p,]= colMeans(MSE_intervals[[p]])
}
colnames(MSE_means)= c("PS","GS","G+P")
MSE_means= data.frame( Traits= phenos, MSE_means)

##saving results
save(accuracy_intervals, file="accuracy_intervals.Rdata")
write.csv(accuracy_means, file="accuracy_means.csv")
save(MSE_intervals, file="MSE_intervals.Rdata")
write.csv(MSE_means, file="MSE_means.csv")
