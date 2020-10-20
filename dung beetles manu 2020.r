## Manu dung beetles
## Laura  Allen - 10/10/2020

## Key questions:
## 1) How does the alpha diversity profile change across disturbance gradient
## 2) How does species composition change (RDA) 
## 3) Does ratio of abundance of large: small beetles change?
## 4) does ratio of abundance of rollers:tunnellers change?
## 5) does the dispersal of small+medium seeds change along the gradient (large seeds = control?).

## other necessary points:
## veg pca to show disturbance
## morans i for spatial autocorrelation?

# libraries and sources ----
rm(list=ls()) ##Clear memory
setwd("C:/Data/PhD/Analysis/Dungbeetles")

source("C://Data/PhD/Analysis/Butterflies/general/diversity.r")
require(ggplot2)
require(reshape2)
require(BiodiversityR)
require(vegan)
require(plotrix)
library(iNEXT)
library(tidyverse)
library(stringr)


## Read in all datasets and combine ----

# #
weather <- read.csv("C://Data/PhD/Processed_data/Dungbeetles/DB_div_weather.csv")
site_data <- read.csv("C://Data/PhD/Processed_data/Site_data/sites_elev_dist_coords.csv")
soil <- read.csv("C://Data/PhD/Processed_data/Soil/Soil_data.csv") ## not much point using soil data as not enough samples

head(weather)
head(site_data)
head(soil)
colnames(weather)[1] <- "Site"
# edit site names so they match
levels(weather$Site)[levels(weather$Site)=="T6-850"] <- "MIN-B"
levels(weather$Site)[levels(weather$Site)=="T7-1350"] <- "MIN-A"
levels(weather$Site)[levels(weather$Site)=="CH-400"] <- "MIN-C"
levels(weather$Site)[levels(weather$Site)=="T1-650"] <- "CCR-B"
levels(weather$Site)[levels(weather$Site)=="T2-800"] <- "CCR-A"
levels(weather$Site)[levels(weather$Site)=="T5-250"] <- "CCR-C"
levels(weather$Site)[levels(weather$Site)=="T2-2150"] <- "MXD-A"
levels(weather$Site)[levels(weather$Site)=="T10-100"] <- "MXD-B"
levels(weather$Site)[levels(weather$Site)=="T3-1800"] <- "MXD-C"


DB_merge1 <- merge(site_data,weather,by="Site")
DB_merge2 <- merge(DB_merge1,soil,by="Site")
head(DB_merge2)
## as we prepare other datasets later, they can be added to dataframe ----
# Edited alpha file in Excel to separate Site_rank into rank and site columns
# Also edited richness output to replace spaces with underscores and remove % and () from colnames

# colnames(veg)[1] <- "Site"
# colnames(rich)[1] <- "Site"
# colnames(alpha)[1] <- "Site"

# newalphasite <- str_split_fixed(alpha$Site, ".",3) # split Site column, to get site values without Rank. infront
# alpha$Site <- newalphasite[,3]
#
# alpha$Site <- as.factor(alpha$Site)
# levels(alpha$Site)
# levels(alpha$Site)[levels(alpha$Site)=="MinD-B"] <- "MIN-B"
# levels(alpha$Site)[levels(alpha$Site)=="MinD-A"] <- "MIN-A"
# levels(alpha$Site)[levels(alpha$Site)=="MinD-C"] <- "MIN-C"


#
# #
# DB_merge1 <- merge(alpha,site_data,by="Site")
# DB_merge2 <- merge(DB_merge1,veg,by="Site")
# DB_merge3 <- merge(DB_merge2,weather,by="Site")
# DB_merge4 <- merge(DB_merge3,abund,by="Site")
# DB_merge5 <- merge(DB_merge4,soil,by="Site")
# DB_merge6 <- merge(DB_merge5,Estqs,by="Site")
# DB_allvars <- merge(DB_merge6,rich,by="Site")
# 
# head(DB_allvars)
##write.csv(DB_allvars,"C://Data/PhD/Processed_data/Dungbeetles/DB_div_allvars.csv")

# 
# table1 <-data.frame(as.character("text"),as.character("text"),0,0,0)
# names(table1) <- c("Variable 1","Variable 2", "Rho", "Rho 95% LCI", "Rho 95% UCI")
# str(table1)
# table1 <- rbind(table1,c(as.character("test1"),as.character("test2"),0.23,0.2,0.4))
# calculating dung beetle diversity


## Code for plots -----
## diversity profile
cf1 <- c('darkred','firebrick2','darksalmon','darkslategray3','dodgerblue2','navyblue') ## colour palette for gradient
cf <- c('darkred','darkred','darkred','firebrick2','firebrick2','firebrick2','darksalmon','darksalmon','darksalmon','darkslategray3','darkslategray3','darkslategray3','dodgerblue2','dodgerblue2','dodgerblue2','navyblue','navyblue','navyblue')
# plot function
gradplot <- function(measure,y,colgrad=cf,legcol=cf1,legpos="topright",title="",ys=c(min(measure),max(measure))){
  plot(rep(0,length(qs1))~qs1,type="l",col="white",ylim=ys,xlim=c(0,2.8),bty='L',xlab=c("q"),
       ylab=y,cex.lab=1.2,cex.axis=1,xaxt="n",main=title)
  axis(1, at=c(log(1),log(2),log(4),log(8),log(16)), labels=c("0","0.5","1","2",expression(infinity)), las=0, cex.axis=1)
  axis.break(axis=1,breakpos=0.35,pos=NA,bgcol="white",breakcol="black",
             style="zigzag",brw=0.05)
  axis.break(axis=1,breakpos=2.5,pos=NA,bgcol="white",breakcol="black",
             style="zigzag",brw=0.05)
  for(rw in 1:18){
    lines(measure[rw,c(1:18)]~log(qs2[1:18]),col=colgrad[rw],lwd=3)
    #points(measure[rw,]~log(qs2),col=colgrad[rw],pch=16,cex=0.7)
  }
  abline(v=c(log(1),log(2),log(4),log(8),log(16)),lty="dotted",col="black")
  legend(legpos,legend= c("1.Banana (most disturbed)","2.Agroforestry","3.Disturbed secondary","4.Cleared regenerating","5.Mixed history","6.Minimally disturbed"),
         cex=1,pch=21,pt.bg=(col=legcol),pt.cex=1.7,bty="n") #Added key later (after submission) to explain subgroups.
}

##
## import data ----
## species by site matrix
db <- read.csv("C://Data/PhD/Processed_data/Dungbeetles/db_spXsite.csv",row.names=1) 
db.mat <- t(as.matrix(db))
head(db.mat)

#checking for singletons and doubletons
length(which(rowSums(db.mat)=="1")) # overall species only found once or twice
length(which(rowSums(db.mat)=="2"))
length(rowSums((db.mat)))
length(which(db.mat=="1")) #in each site
length(which(db.mat=="2"))
length(which(db.mat>"2"))

## 1) How does the alpha diversity profile change across disturbance gradient ----
## alpha diversity
p1 <- db.mat/sum(db.mat) #turned counts into proportions
qs1 <- c(0,seq(from = 0.5, to = 2,by=0.1),Inf)# for calculating
qs2 <- c(1,seq(from = 2, to = 8,by=0.4),16) # for plotting
# qs1 <- c(0,0.5,1,2,Inf)# for calculating
# qs2 <- c(1,2,4,8,16) # for plotting
alpha<- subcommunity.alpha.bar(populations=p1,qs=qs1) # calculate alpha diversity profile (raw values at total sample size)
str(alpha)
#write.csv(alpha,"C://Data/PhD/Processed_data/Dungbeetles/Diversity/alpha.csv")

# plot diversity profile
#tiff(file="C://Data/PhD/Outputs/Dungbeetles/Biotropica figures/Fig1_alphaprofile.tiff",width=190,height=110,units="mm",res=1000, pointsize=9)  
par(mar=c(4.5,4,1,1))
gradplot(measure=alpha,y=c("Effective number of species"),colgrad=cf,legcol=cf1)
dev.off()

## add alpha data main table ----

alphadf <- as.data.frame(alpha)
newalphasite <- str_split_fixed(rownames(alphadf), ".",3) # split Site column, to get site values without Rank. infront
alphadf$Site <- newalphasite[,3]
alphadf$Site <- as.factor(alphadf$Site)
levels(alphadf$Site)
levels(alphadf$Site)[levels(alphadf$Site)=="MinD-B"] <- "MIN-B"
levels(alphadf$Site)[levels(alphadf$Site)=="MinD-A"] <- "MIN-A"
levels(alphadf$Site)[levels(alphadf$Site)=="MinD-C"] <- "MIN-C"

alphadf <- alphadf[,c("q0","q0.5","q1","q2","qInf","Site")] ## select only columns with qs of main interest 
#DB_merge3 <- merge(DB_merge2,alphadf,by="Site") # join with main df

## to get estimate alpha diversity at equal sample sizes and get estimate uncertainty used iNEXT with bootstrapping
m <- c(50, 100, 200,300,400,500,600) # series of sample sizes for extrapolation 
iN123 <- iNEXT(db.mat, q=c(0,0.5,1,2,Inf), datatype="abundance", size=NULL, endpoint=300, knots=80, se=TRUE, conf=0.95, nboot=50)

## extract rows for equal sample size for each size and combine into a table ----
#extract the dataframes for each site
estq123_BAA <- as.data.frame(iN123$iNextEst$`1.BA-A`)
estq123_BAB <- as.data.frame(iN123$iNextEst$`1.BA-B`)
estq123_BAC <- as.data.frame(iN123$iNextEst$`1.BA-C`)
estq123_AFA <- as.data.frame(iN123$iNextEst$`2.AF-A`)
estq123_AFB <- as.data.frame(iN123$iNextEst$`2.AF-B`)
estq123_AFC <- as.data.frame(iN123$iNextEst$`2.AF-C`)
estq123_SFA <- as.data.frame(iN123$iNextEst$`3.SF-A`)
estq123_SFB <- as.data.frame(iN123$iNextEst$`3.SF-B`)
estq123_SFC <- as.data.frame(iN123$iNextEst$`3.SF-C`)
estq123_CCRA <- as.data.frame(iN123$iNextEst$`4.CCR-A`)
estq123_CCRB <- as.data.frame(iN123$iNextEst$`4.CCR-B`)
estq123_CCRC <- as.data.frame(iN123$iNextEst$`4.CCR-C`)
estq123_MXDA <- as.data.frame(iN123$iNextEst$`5.MXD-A`)
estq123_MXDB <- as.data.frame(iN123$iNextEst$`5.MXD-B`)
estq123_MXDC <- as.data.frame(iN123$iNextEst$`5.MXD-C`)
estq123_MINA <- as.data.frame(iN123$iNextEst$`6.MinD-A`)
estq123_MINB <- as.data.frame(iN123$iNextEst$`6.MinD-B`)
estq123_MINC <- as.data.frame(iN123$iNextEst$`6.MinD-C`)

r <- which(estq123_BAA$m==300) #select the end rows for equal sample size estimates
BAA <- estq123_BAA[r,]
r <- which(estq123_BAB$m==300) #select the end rows for equal sample size estimates
BAB <- estq123_BAB[r,]
r <- which(estq123_BAC$m==300) #select the end rows for equal sample size estimates
BAC <- estq123_BAC[r,]
r <- which(estq123_AFA$m==300) #select the end rows for equal sample size estimates
AFA <- estq123_AFA[r,]
r <- which(estq123_AFB$m==300) #select the end rows for equal sample size estimates
AFB<- estq123_AFB[r,]
r <- which(estq123_AFC$m==300) #select the end rows for equal sample size estimates
AFC<- estq123_AFC[r,]
r <- which(estq123_SFA$m==300) #select the end rows for equal sample size estimates
SFA <- estq123_SFA[r,]
r <- which(estq123_SFB$m==300) #select the end rows for equal sample size estimates
SFB<- estq123_SFB[r,]
r <- which(estq123_SFC$m==300) #select the end rows for equal sample size estimates
SFC<- estq123_SFC[r,]
r <- which(estq123_CCRA$m==300) #select the end rows for equal sample size estimates
CCRA <- estq123_CCRA[r,]
r <- which(estq123_CCRB$m==300) #select the end rows for equal sample size estimates
CCRB <- estq123_CCRB[r,]
r <- which(estq123_CCRC$m==300) #select the end rows for equal sample size estimates
CCRC<- estq123_CCRC[r,]
r <- which(estq123_MXDA$m==300) #select the end rows for equal sample size estimates
MXDA <- estq123_MXDA[r,]
r <- which(estq123_MXDB$m==300) #select the end rows for equal sample size estimates
MXDB <- estq123_MXDB[r,]
r <- which(estq123_MXDC$m==300) #select the end rows for equal sample size estimates
MXDC<- estq123_MXDC[r,]
r <- which(estq123_MINA$m==300) #select the end rows for equal sample size estimates
MINA <- estq123_MINA[r,]
r <- which(estq123_MINB$m==300) #select the end rows for equal sample size estimates
MINB <- estq123_MINB[r,]
r <- which(estq123_MINC$m==300) #select the end rows for equal sample size estimates
MINC<- estq123_MINC[r,]

Estqs <- rbind(BAA,BAB,BAC,AFA,AFB,AFC,SFA,SFB,SFC,CCRA,CCRB,CCRC,MXDA,MXDB,MXDC,MINA,MINB,MINC)
Site <- c(rep("BA-A",5),rep("BA-B",5),rep("BA-C",5),rep("AF-A",5),rep("AF-B",5),rep("AF-C",5),rep("SF-A",5),rep("SF-B",5),rep("SF-C",5),rep("CCR-A",5),rep("CCR-B",5),rep("CCR-C",5),rep("MXD-A",5),rep("MXD-B",5),rep("MXD-C",5),rep("MIN-A",5),rep("MIN-B",5),rep("MIN-C",5))
Ests <- cbind(Estqs,Site)
Site <- c("BA-A","BA-B","BA-C","AF-A","AF-B","AF-C","SF-A","SF-B","SF-C","CCR-A","CCR-B","CCR-C","MXD-A","MXD-B","MXD-C","MIN-A","MIN-B","MIN-C")

Est0 <- Ests[Ests$order==0,4]
Est1 <- Ests[Ests$order==1,4]
Est2 <- Ests[Ests$order==2,4]
EstIn <- Ests[Ests$order==Inf,4]
Est_eqSS <- as.data.frame(cbind(Site=as.character(Site),Est0_ss=as.numeric(Est0),Est1_ss=as.numeric(Est1),Est2_ss=as.numeric(Est2),EstInf_ss=as.numeric(EstIn)))
# 
#write.csv(Est_eqSS,"C://Data/PhD/Processed_data/Dungbeetles/db_inext_Ests_EqSS.csv")
#DB_merge4 <- merge(DB_merge3,Est_eqSS,by="Site") # join with main df
#write.csv(DB_merge4,"C://Data/PhD/Processed_data/Dungbeetles/db_mergedata4.csv")
DB_merge4 <- read.csv("C://Data/PhD/Processed_data/Dungbeetles/db_mergedata4.csv")
head(DB_merge4)

## Permutation for estimated diversity ----
## rescramble data and see what proportion of runs the correlation will be found

# function for correlation test
co <- function(ia,ib){
  res <- cor.test(ia, ib, method = "spearman",exact=F)
  cbind(res$p.value)}

# get data into right format
Est_eqSS$Rank <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)) ## add a column for disturbance rank
Est_eqSS$Est0_ss <- as.numeric(as.character(Est_eqSS$Est0_ss))
Est_eqSS$Est1_ss <- as.numeric(as.character(Est_eqSS$Est1_ss))
Est_eqSS$Est2_ss <- as.numeric(as.character(Est_eqSS$Est2_ss))

#  re-ordered ranking, re-run correlation test 1000 times.
perms <- c(NULL)
for(s in 1:10000){
  rs <- sample(rep(1:6,3)) ## artificial reordered ranking
  cr0 <- co(rs,Est_eqSS$Est0_ss)
  cr1 <- co(rs,Est_eqSS$Est1_ss)
  cr2 <- co(rs,Est_eqSS$Est2_ss)
  logcr <- log(cr0)+log(cr1)+log(cr2)# combined p vals
  run <- cbind(s,cr0,cr1,cr2,logcr)
  perms <- rbind(perms,run)
}
head(perms)
perms <- as.data.frame(perms)
names(perms) <- c("Run","p_q0","p_q1","p_q2","p_logcr")

# compare permuted results with real correlation
cr0 <- co(Est_eqSS$Rank,Est_eqSS$Est0_ss)
cr1 <- co(Est_eqSS$Rank,Est_eqSS$Est1_ss)
cr2 <- co(Est_eqSS$Rank,Est_eqSS$Est2_ss)
logcr <- log(cr0)+log(cr1)+log(cr2)
real<- as.data.frame(cbind(0,cr0,cr1,cr2,logcr)) #Run 0 is the real values
names(real) <- c("Run","p_q0","p_q1","p_q2","p_logcr")
comb <- rbind(perms,real)
head(comb)
tail(comb)

# check if real (Run 0) is in bottom 5% of pvals for each
pq0 <- rank(comb$p_q0)[10001]/length(comb[,1]) #divide the rank of the real p val by the total number to ge the p value of the probability of getting that result by chance
pq1 <- rank(comb$p_q1)[10001]/length(comb[,1])
pq2 <- rank(comb$p_q2)[10001]/length(comb[,1])
pqlog <- rank(comb$p_logcr)[10001]/length(comb[,1])

pvals_perm <- cbind(pq0,pq1,pq2,pqlog)
pvals_perm
# pq0        pq1       pq2      pqlog
# 0.01029897 0.07649235 0.2748725 0.03559644

###///////////
## Vegetation data ----

# Clean dataset 
options(stringsAsFactors = FALSE) # import not as factors
veg=read.csv("C://Data/PhD/Rawdata/Vegetation/Vegmap_gradient_means.csv",header=TRUE)

# rename sites so in order
veg$Site[veg$Site=="CH-400"] <- "xf.MIN-C"
veg$Site[veg$Site=="T1-650"] <- "xd.CCR-B"
veg$Site[veg$Site=="T10-100"] <- "xe.MXD-B"
veg$Site[veg$Site=="T2-2150"] <- "xe.MXD-A"
veg$Site[veg$Site=="T2-800"] <- "xd.CCR-A"
veg$Site[veg$Site=="T3-1800"] <- "xe.MXD-C"
veg$Site[veg$Site=="T5-250"] <- "xd.CCR-C"
veg$Site[veg$Site=="T6-850"] <- "xf.MIN-B"
veg$Site[veg$Site=="T7-1350"] <- "xf.MIN-A"
veg$Site[veg$Site=="AF-A"] <- "xb.AF-A"
veg$Site[veg$Site=="AF-B"] <- "xb.AF-B"
veg$Site[veg$Site=="AF-C"] <- "xb.AF-C"
veg$Site[veg$Site=="BA-A"] <- "xa.BA-A"
veg$Site[veg$Site=="BA-B"] <- "xa.BA-B"
veg$Site[veg$Site=="BA-C"] <- "xa.BA-C"
veg$Site[veg$Site=="SF-A"] <- "xc.SF-A"
veg$Site[veg$Site=="SF-B"] <- "xc.SF-B"
veg$Site[veg$Site=="SF-C"] <- "xc.SF-C"

# calculate means per site from quadrat data
mean.veg <- aggregate(veg,list(veg$Site),mean,na.rm=T)
sort(mean.veg$Group.1)
mean.veg <- mean.veg[order(mean.veg$Group.1),]
str(mean.veg)
Sitenames <- mean.veg$Group.1
veg1 <- mean.veg[,-c(1,2,3,10,17)] #remove unwanted columns
for(r in 2:12){
veg1[,r] <- as.numeric(veg1[,r])
}
veg.mat <- as.matrix(veg1)
rownames(veg.mat) <- str_sub(Sitenames, 4)

## Vegetation PCA
# calculate rda of veg.mat
rda.out <- rda(veg.mat,scale=T) # default scale option =F 
# the rda (redundancy analysis) function calculates the PCA and can also be used for canonical RDA
# scale = F means the data should be centred by columns but not standardised
summary(rda.out) # shows the two axes of the PCA for each site
# scaling = 1 preserves Euclidean distance among objects
# scaling = 2 preserves correlation between descriptors

# PCA w scaling 2
veg_pca2 <- summary(rda.out, scaling=1) # Eigenvalues, species scores, site scores 
biplot(rda.out, scaling=2) # Function biplot.rda() of vegan is used 
scrs <- scores(rda.out,display = c("sites", "species"), scaling =1) #extract sepecies and site scores from rda summary
sort(scrs$species[,1])
veg_pc1 <- veg_pca2$sites[,1] #PC1 axis
veg_pc2 <- veg_pca2$sites[,2] #PC2 axis
veg_pca_out <- as.data.frame(cbind("Site"=rownames(veg.mat),veg_pc1,veg_pc2))
# write.csv(veg_pca_out,"C://Data/PhD/Processed_data/Vegetation/veg_pca_outOct18.csv")
DB_merge5 <- merge(DB_merge4,veg_pca_out,by="Site") # join with main df

# custom PCA plot ----
str(veg.mat)
rda.out <- rda(veg.mat,scale=T)
scrs <- scores(rda.out,display = c("sites", "species"), scaling =2) #extract sepecies and site scores from rda summary
summary(rda.out)
max(scrs$sites)
#tiff(file="C://Data/PhD/Outputs/Dungbeetles/Biotropica figures/FigS6_PCA_veg.tiff",width=140,height=140,units="mm",res=1000, pointsize=10)  
par(mar=c(1,1,1,1))
par(mfrow=c(1,1))
plot.new()
plot.window(xlim = c(-2,2), ylim = c(-2.5,2), asp = 1)
axis(side=1, at=c(-2,-1.5,-1,-0.5,0.5,1,1.5,2), lty="solid", pos=c(0,0),padj=-1,cex.axis=0.7)
axis(side=2, at=c(-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2), lty="solid", pos=c(0,0),padj=1,cex.axis=0.7)
colvec <- c('#d73027','#d73027','#d73027','#fc8d59','#fc8d59','#fc8d59','#fee090','#fee090','#fee090','#e0f3f8','#e0f3f8','#e0f3f8','#91bfdb','#91bfdb','#91bfdb','#4575b4','#4575b4','#4575b4')
symvec <- c(rep("1",3),rep("2",3),rep("3",3),rep("4",3),rep("5",3),rep("6",3))
points(scrs$sites, col=colvec,pch = 1,cex=4, lwd=4.5)
points(scrs$sites, col=colvec, bg= "grey",pch = symvec,cex=1.5,lwd=1)

# plot arrows for each environental variable
arrows(x0=0, y0=0, scrs$species[1,1], scrs$species[1,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[2,1], scrs$species[2,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[3,1], scrs$species[3,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[4,1], scrs$species[4,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[5,1], scrs$species[5,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[6,1], scrs$species[6,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[7,1], scrs$species[7,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[8,1], scrs$species[8,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[9,1], scrs$species[9,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[10,1], scrs$species[10,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[11,1], scrs$species[11,2], code=2, lwd=1, col="black",length=0.1)
arrows(x0=0, y0=0, scrs$species[12,1], scrs$species[12,2], code=2, lwd=1, col="black",length=0.1)

text(x=scrs$species[1,1]+0.4, y= scrs$species[1,2]-0.05, labels = "Canopy height", cex=1,col="black")
text(x=scrs$species[2,1], y= scrs$species[2,2]+0.1, labels = "Canopy cover", cex=1,col="black")
text(x=scrs$species[3,1]+0.3, y= scrs$species[3,2], labels = "Mid-canopy height", cex=1,col="black")
text(x=scrs$species[4,1], y= scrs$species[4,2]+0.1, labels = "Mid-canopy cover", cex=1,col="black")
text(x=scrs$species[5,1], y= scrs$species[5,2], labels = "Leaf-litter depth", cex=1,col="black")
text(x=scrs$species[6,1]-0.2, y= scrs$species[6,2]+0.2, labels = "DBH largest trees", cex=1,col="black")
text(x=scrs$species[7,1], y= scrs$species[7,2], labels = "n trees DBH>5", cex=1,col="black")
text(x=scrs$species[8,1]+0.1, y= scrs$species[8,2]+0.1, labels = "Shrub layer", cex=1,col="black")
text(x=scrs$species[9,1]+0.2, y= scrs$species[9,2]+0.1, labels = "Herb layer", cex=1,col="black")
text(x=scrs$species[10,1], y= scrs$species[10,2], labels = "Herb abundance", cex=1,col="black")
text(x=scrs$species[11,1], y= scrs$species[11,2], labels = "Bare ground", cex=1,col="black")
text(x=scrs$species[12,1], y= scrs$species[12,2], labels = "Woody debris", cex=1,col="black")
text(x=1.8, y=0.1,labels = c("PCA1"),cex=1)
text(x=0.1, y=-2,labels = c("PCA2"),srt=90,cex=1)
legend(x=0.5,y=-1.3,bty="n",pch=21,pt.cex=1,cex=1,pt.bg=c('#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4'),legend=c("1.Banana (most disturbed)","2.Agroforestry","3.Disturbed secondary","4.Cleared regenerating","5.Mixed history","6.Minimally disturbed"))
dev.off()


###########################
##//////////
## Makes sense up to this line. ----
##//////////

#iNEXT resampling plots -----
inr0 <- read.csv("C://Data/PhD/Outputs/Dungbeetles/DBq0ests_rsampled.csv") ## no idea where these files were created
inr1 <- read.csv("C://Data/PhD/Outputs/Dungbeetles/DBq1ests_rsampled.csv")
inr2 <- read.csv("C://Data/PhD/Outputs/Dungbeetles/DBq2ests_rsampled.csv")
inrI <- read.csv("C://Data/PhD/Outputs/Dungbeetles/DBqInfests_rsampled.csv")

tiff(file="C://Data/PhD/Outputs/Dungbeetles/Biotropica figures/Fig2_alphaests_inextresampling.tiff",width=190,height=70,units="mm",res=1000, pointsize=12)  
par(mfrow=c(1,3))
par(mar=c(5.1,5.1,4.1,2.1))
plot(inr0$meandiffs~inr0$rankdiff,1,pch=21,cex=1,cex.axis=1,cex.lab=1.2,col='black',bty="n",ylab="Diversity estimates:
     Less disturbed > More disturbed",xlab="Difference in disturbance rank",ylim=c(0,1),xlim=c(1,5))
abline(lm(inr0$meandiffs~inr0$rankdiff))
mtext(side = 3, line = 0,adj=0,"a",font=2,cex=1.2)
plot(inr1$meandiffs~inr1$rankdiff,1,pch=21,cex=1,cex.axis=1,cex.lab=1.2,col='black',bty="n",ylab="Diversity estimates:
     Less disturbed > More disturbed",xlab="Difference in disturbance rank",ylim=c(0,1),xlim=c(1,5))
abline(lm(inr1$meandiffs~inr1$rankdiff))
mtext(side = 3, line = 0,adj=0,"b",font=2,cex=1.2)
plot(inr2$meandiffs~inr2$rankdiff,1,pch=21,cex=1,cex.axis=1,cex.lab=1.2,col='black',bty="n",ylab="Diversity estimates:
     Less disturbed > More disturbed",xlab="Difference in disturbance rank",ylim=c(0,1),xlim=c(1,5))
abline(lm(inr2$meandiffs~inr2$rankdiff))
mtext(side = 3, line = 0,adj=0,"c",font=2,cex=1.2)
dev.off()





## Models ----

