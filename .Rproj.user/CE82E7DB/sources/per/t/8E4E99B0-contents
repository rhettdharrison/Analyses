setwd("~/Documents/Projects/BMZ trees on farms Zambia/Analyses")
library(vegan)
library(BIOMASS)
##Changes from Tree level data tk 3.csv to Tree level data tk 4.csv
# trees$ALIVE replace "D" with "d"
# trees$PLANTED replace "a" with "n"
# delete duplicate plot 3163 data
# trees$STEM_NO replace "a" with 1
# tree$PLANT_ID remove extra ids where trees$STEM_NO > 1
#index <- c(rep(0,dim(trees)[1]))
#for(i in 1:dim(trees)[1]){
#  if(trees$STEM_NO[i]==1) {
#    next
#  } else {
#    if(trees$PLANT_NO[i] == trees$PLANT_NO[i-1]) {
#      next
#    } else {
#      index[i] <- 1
#    }
#    
#  }
#}
#which(index==1)
# [1]   99  200  206  238 2295 2299 2708
#99 PLANT_ID -> 15, HEIGHT -> NA
#200 PLANT_ID -> 4, HEIGHT -> NA
#206 STEM_NO -> 1
#238 STEM_NO -> 1
#2295 STEM_NO -> 1, 2296 STEM_NO -> 2
#2299 STEM_NO -> 1
#2708 PLANT_ID -> 4, HEIGHT -> NA

trees <- read.csv("Data/Tree level data tk 4.csv",sep=",",header=T)
head(trees)
str(trees)
trees$HEIGHT[trees$STEM_NO>1] <- NA

plots <- read.csv("Data/Plot_Level_Data Tk2.csv",sep=",",header=T)
head(plots)
str(plots)
plots <- plots[plots$PLOT_ID %in% unique(trees$PLOT_ID),]
dim(plots)
plots <- plots[order(plots$PLOT_ID),]


trees$BA <- pi*(trees$GBH/200)^2

plot_BA <- tapply(trees$BA,trees$PLOT_ID,sum)

#calculate plot BA per ha
plot_BA <- (plot_BA/(pi*plots$PLOT_SIZE^2))*10000
plot(log(plot_BA))

png("Basal Area.png")
boxplot(plot_BA~droplevels(plots$RECLASSIFICATION),ylab="Basal Area")
dev.off()

#need to calculate biomass
head(trees)
trees$GENUS <- unlist(strsplit(as.character(trees$SPECIES)," "))[(1:length(trees$SPECIES))*2 - 1]
trees$SPP <- unlist(strsplit(as.character(trees$SPECIES)," "))[(seq(2,length(trees$SPECIES)*2,2))]

WSG <- getWoodDensity(trees$GENUS,trees$SPP)
head(WSG)

trees$PLANTS <- paste(trees$PLOT_ID,trees$PLANT_NO,sep="_")

DBH_sq <- tapply(trees$BA,trees$PLANTS,sum)/pi

AGB <- 0.0673*(WSG$meanWD[!(is.na(trees$HEIGHT))]*trees$HEIGHT[!(is.na(trees$HEIGHT))]*DBH_sq)^0.976

plot_AGB <- tapply(AGB,trees$PLOT_ID[!(is.na(trees$HEIGHT))],sum)
png("AGB.png")
boxplot(plot_AGB~droplevels(plots$RECLASSIFICATION),ylab="Above Ground Biomass")
dev.off()

## Classification accuracy
LUC_accuracy <- c(rep(NA,dim(plots)[1]))
for(i in 1:dim(plots)[1]){
    if(!(is.na(match(plots$RECLASSIFICATION[i],plots$CATEGORY[i])))){
        LUC_accuracy[i] <-1
    } else {
        LUC_accuracy[i] <- 0
    }
}

LUC_accuracy_pr <- tapply(LUC_accuracy,plots$RECLASSIFICATION,sum)/table(plots$RECLASSIFICATION)
LUC_accuracy_pr <- c(LUC_accuracy_pr,sum(LUC_accuracy)/dim(plots)[1])
names(LUC_accuracy_pr) <- c("Crop","Forest","Grass","Wet","Total")

png("assign.png")
barplot(LUC_accuracy_pr, ylab="Proportion of plots correctly assigned")
dev.off()

#Create community matrix
trees_2 <- trees[trees$STEM_NO==1,]
trees_2 <- trees_2[!(is.na(match(trees_2$ALIVE,"a"))),]

abundance <- tapply(trees_2$STEM_NO,list(trees_2$PLOT_ID,trees_2$SPECIES),sum)
abundance[is.na(abundance)] <- 0
dim(abundance)


boxplot(diversity(abundance,index="simpson")~plots$RECLASSIFICATION,ylab="Simpson Diversity")

png("Outputs/diversity.png")
boxplot(diversity(abundance,index="shannon")~plots$RECLASSIFICATION, ylab="Shannon Diversity")
dev.off()

png("Outputs/abundance.png")
barplot(sort(apply(abundance,2,sum),decreasing=T),names.arg = "",xlab="SPECIES",ylab="ABUNDANCE")
dev.off()

# remove species occuring in one plot
abundance <- abundance[,apply(sign(abundance),2,sum)>1]
dim(abundance)

abundance <- abundance[apply(sign(abundance),1,sum)>0,]
dim(abundance)

meta_plots <- metaMDS(abundance, k=3, try = 40) #uses abundance data SOLUTION REACHED
meta_plots_j <- metaMDS(abundance, dist="jaccard", k=4, try = 40) #presence absence NO CONVERGENCE
##metaMDS is not converging. Many plots are identical (which seems strange but I do not know the anaswer yet)


trees_2 <- trees[trees$STEM_NO==1,]
dim(trees_2)
trees_2 <- cbind(trees_2,AGB)

species_AGB <- tapply(trees_2$AGB,list(trees_2$PLOT_ID,trees_2$SPECIES),sum,na.rm=T)
species_AGB[is.na(species_AGB)] <- 0

# remove species occuring in one plot
species_AGB <- species_AGB[,apply(sign(species_AGB),2,sum)>1]
dim(species_AGB)

species_AGB <- species_AGB[apply(sign(species_AGB),1,sum)>0,]
dim(species_AGB)

meta_plots_AGB <- metaMDS(species_AGB,distance = "mahalanobis",try =40, k = 4) # NO CONVERGENCE

##ploting 

png("Outputs/NMDS.png")
par(pty='s')
par(mar=c(8,8,0.5,0.5))
ordiplot(meta_plots,display="sites", type="none",)
points(meta_plots$points[plots$RECLASSIFICATION[-60]=="Forest",1:2],col="dark green",pch=20, cex=3)
points(meta_plots$points[plots$RECLASSIFICATION[-60]=="Cropland",1:2],col="yellow",pch=20, cex=3)
points(meta_plots$points[plots$RECLASSIFICATION[-60]=="Wetland",1:2],col="blue",pch=20, cex=3)
points(meta_plots$points[plots$RECLASSIFICATION[-60]=="Grassland",1:2],col="green",pch=20, cex=3)
dev.off()

