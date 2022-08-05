###################### NMP Filter by Signal Intensity  Script ##################

library(dplyr)
library(data.table)
library(ggplot2)
library(plotly)
library(tidyverse)
library(factoextra)
library(viridisLite)
library(viridis)
library(reshape2)
library(gplots)
library(gatepoints)

###############  Formatting Data - adapted from ASCRIBE ################ 

DMSO <- list() # creates a list
listcsv <- dir(pattern = "*DMSO.csv") # creates the list of all the csv files in the directory

for (k in 1:length(listcsv)){
  
  DMSO[[k]] <- read.csv(listcsv[k])
  
  ########### Fixing Voxel Length - Co-ordinates in pixels not microns
  DMSO[[k]]$X.microns<-DMSO[[k]]$X*0.3787879
  DMSO[[k]]$Y.microns<-DMSO[[k]]$Y*0.3787879
  DMSO[[k]]$Z.microns<-DMSO[[k]]$Z*0.9999284
  
  DMSO[[k]]$Reduced.X <- DMSO[[k]]$X.microns - min(DMSO[[k]]$X.microns)
  DMSO[[k]]$Normalised.X <- DMSO[[k]]$Reduced.X / max(DMSO[[k]]$Reduced.X)
  
  DMSO[[k]]$Reduced.Y <- DMSO[[k]]$Y.microns - min(DMSO[[k]]$Y.microns)
  DMSO[[k]]$Normalised.Y <- DMSO[[k]]$Reduced.Y / max(DMSO[[k]]$Reduced.Y)
  
  DMSO[[k]]$Reduced.Z <- DMSO[[k]]$Z.microns - min(DMSO[[k]]$Z.microns)
  DMSO[[k]]$Normalised.Z <- DMSO[[k]]$Reduced.Z / max(DMSO[[k]]$Reduced.Z)
  
  DMSO[[k]]$Embryo <- k
  
}

######## Takes sum of all embryos 
sumDMSO <- rbindlist(DMSO)

####### Calculates Mean across different co-ordinates 
MeanYDMSO <- mean(max(DMSO[[k]]$Y.microns)) - mean(min(DMSO[[k]]$Y.microns))
MeanXDMSO <- mean(max(DMSO[[k]]$X.microns)) - mean(min(DMSO[[k]]$X.microns))
MeanZDMSO <- mean(max(DMSO[[k]]$Z.microns)) - mean(min(DMSO[[k]]$Z.microns))

sumDMSO$Scaled.X <- sumDMSO$Normalised.X * MeanXDMSO

sumDMSO$Scaled.Y <- sumDMSO$Normalised.Y * MeanYDMSO
sumDMSO$Scaled.Y <- (max(sumDMSO$Scaled.Y) - sumDMSO$Scaled.Y) 

sumDMSO$Scaled.Z <- sumDMSO$Normalised.Z * MeanZDMSO

sumDMSO$Condition <- "DMSO"

sumDMSO


#############################################################################################

DAPT <- list() # creates a list
listcsv <- dir(pattern = "*DAPT.csv") # creates the list of all the csv files in the directory
for (k in 1:length(listcsv)){
  
  DAPT[[k]] <- read.csv(listcsv[k])
  
  
  #######Convert to microns 
  DAPT[[k]]$X.microns<-DAPT[[k]]$X*0.3787879
  DAPT[[k]]$Y.microns<-DAPT[[k]]$Y*0.3787879
  DAPT[[k]]$Z.microns<-DAPT[[k]]$Z*0.9999284
  
  DAPT[[k]]$Reduced.X <- DAPT[[k]]$X.microns - min(DAPT[[k]]$X.microns)
  DAPT[[k]]$Normalised.X <- DAPT[[k]]$Reduced.X / max(DAPT[[k]]$Reduced.X)
  
  DAPT[[k]]$Reduced.Y <- DAPT[[k]]$Y.microns - min(DAPT[[k]]$Y.microns)
  DAPT[[k]]$Normalised.Y <- DAPT[[k]]$Reduced.Y / max(DAPT[[k]]$Reduced.Y)
  
  DAPT[[k]]$Reduced.Z <- DAPT[[k]]$Z.microns - min(DAPT[[k]]$Z.microns)
  DAPT[[k]]$Normalised.Z <- DAPT[[k]]$Reduced.Z / max(DAPT[[k]]$Reduced.Z)
  
  DAPT[[k]]$Embryo <- k
  
}

sumDAPT <- rbindlist(DAPT)

MeanYDAPT <- mean(max(DAPT[[k]]$Y.microns)) - mean(min(DAPT[[k]]$Y.microns))
MeanXDAPT <- mean(max(DAPT[[k]]$X.microns)) - mean(min(DAPT[[k]]$X.microns))
MeanZDAPT <- mean(max(DAPT[[k]]$Z.microns)) - mean(min(DAPT[[k]]$Z.microns))

sumDAPT$Scaled.X <- sumDAPT$Normalised.X * MeanXDAPT

sumDAPT$Scaled.Y <- sumDAPT$Normalised.Y * MeanYDAPT

sumDAPT$Scaled.Z <- sumDAPT$Normalised.Z * MeanZDAPT

sumDAPT$Condition <- "DAPT"

sumDAPT


##############################################################################################


sumintensities <- rbind(sumDMSO, sumDAPT)

sumintensities$ID <- seq.int(nrow(sumintensities))


##############################################################################################

##### Normalisation of Signal Intensities - DMSO 
sumDMSO$Reduced4 <- sumDMSO$signal4 - 0
sumDMSO$Reduced3 <- sumDMSO$signal3 - 0
sumDMSO$Reduced2 <- sumDMSO$signal2 - 0
sumDMSO$Reduced1 <- sumDMSO$signal1 - 0

sumDMSO$Reduced4[sumDMSO$Reduced4<0] <- 0
sumDMSO$Reduced3[sumDMSO$Reduced3<0] <- 0
sumDMSO$Reduced2[sumDMSO$Reduced2<0] <- 0
sumDMSO$Reduced1[sumDMSO$Reduced1<0] <- 0

sumDMSO$Normalised4 <- sumDMSO$Reduced4 / max(sumDMSO$Reduced4)
sumDMSO$Normalised3 <- sumDMSO$Reduced3 / max(sumDMSO$Reduced3)
sumDMSO$Normalised2 <- sumDMSO$Reduced2 / max(sumDMSO$Reduced2)
sumDMSO$Normalised1 <- sumDMSO$Reduced1 / max(sumDMSO$Reduced1)

sumDMSO$Normalised4 <- sumDMSO$Normalised4 / max(sumDMSO$Normalised4)
sumDMSO$Normalised3 <- sumDMSO$Normalised3 / max(sumDMSO$Normalised3)
sumDMSO$Normalised2 <- sumDMSO$Normalised2 / max(sumDMSO$Normalised2)
sumDMSO$Normalised1 <- sumDMSO$Normalised1 / max(sumDMSO$Normalised1)

##### Normalisation of Signal Intensities - DAPT 
sumDAPT$Reduced4 <- sumDAPT$signal4 - 0
sumDAPT$Reduced3 <- sumDAPT$signal3 - 0
sumDAPT$Reduced2 <- sumDAPT$signal2 - 0
sumDAPT$Reduced1 <- sumDAPT$signal1 - 0

sumDAPT$Reduced4[sumDAPT$Reduced4<0] <- 0
sumDAPT$Reduced3[sumDAPT$Reduced3<0] <- 0
sumDAPT$Reduced2[sumDAPT$Reduced2<0] <- 0
sumDAPT$Reduced1[sumDAPT$Reduced1<0] <- 0

sumDAPT$Normalised4 <- sumDAPT$Reduced4 / max(sumDAPT$Reduced4)
sumDAPT$Normalised3 <- sumDAPT$Reduced3 / max(sumDAPT$Reduced3)
sumDAPT$Normalised2 <- sumDAPT$Reduced2 / max(sumDAPT$Reduced2)
sumDAPT$Normalised1 <- sumDAPT$Reduced1 / max(sumDAPT$Reduced1)

sumDAPT$Normalised4 <- sumDAPT$Normalised4 / max(sumDAPT$Normalised4)
sumDAPT$Normalised3 <- sumDAPT$Normalised3 / max(sumDAPT$Normalised3)
sumDAPT$Normalised2 <- sumDAPT$Normalised2 / max(sumDAPT$Normalised2)
sumDAPT$Normalised1 <- sumDAPT$Normalised1 / max(sumDAPT$Normalised1)

############ NMP Overlay - Signal Intensity Range Histograms ###################

### Soxb1b/c 
ggplot(sumintensities, aes(x = Normalised1)) + geom_histogram()+
  labs(x="Normalised Signal Intensity Soxb1b/c",
       y= "Count")

### Bra1/2 
ggplot(sumintensities, aes(x = Normalised2)) + geom_histogram()+
  labs(x="Normalised Signal Intensity Bra1/2",
       y= "Count")
### EDU
ggplot(sumintensities, aes(x = Normalised4)) + geom_histogram()+
  labs(x="Normalised Signal Intensity EdU",
       y= "Count")


library(dplyr)

###### Signal Names for Different Embryo Groups 

####### 8_10 signal 1 = Sox / All others = Sox
#       8_10 signal 2 = Bra / All others = Tbx 
#       8_10 signal 3 = Tbx / All others = Bra 
#       8_10 signal 4 = Edu / All others = Sox


############ NMP Only - DMSO  ################

sumDMSO$binary <- 0
sumDMSO[sumDMSO$Normalised1 > 0.1 & sumDMSO$Normalised3 > 0.09]$binary <- 1
names(sumDMSO)[names(sumDMSO)=="X.1"] <- 'Label'
Embryo1 <- subset(sumDMSO, Embryo == 1)
ordered_Embryo1 <- Embryo1[order(Embryo1$Label),]
write.csv(ordered_Embryo1, 'Embryo1_NMP_DMSO.csv')

######### NMP Only - DAPT  ############### 

sumDAPT$binary <- 0
sumDAPT[sumDAPT$Normalised1 > 0.1 & sumDAPT$Normalised3 > 0.09]$binary <- 1
names(sumDAPT)[names(sumDAPT)=="X.1"] <- 'Label'
Embryo1 <- subset(sumDAPT, Embryo == 1)
ordered_Embryo1 <- Embryo1[order(Embryo1$Label),]
write.csv(ordered_Embryo1, 'Embryo1_NMP_DAPT.csv')

######## Proliferating NMP DMSO #########

sumDMSO$binary <- 0
sumDMSO[sumDMSO$Normalised1 > 0.1 & sumDMSO$Normalised3 > 0.09 & sumDMSO$Normalised4 >0.2]$binary <- 1
names(sumDMSO)[names(sumDMSO)=="X.1"] <- 'Label'
Embryo1 <- subset(sumDMSO, Embryo == 1)
ordered_Embryo1 <- Embryo1[order(Embryo1$Label),]
write.csv(ordered_Embryo1, 'Embryo1_NMP_DMSO_Proliferating.csv')

######## Proliferating NMP DAPT ####################

sumDAPT$binary <- 0
sumDAPT[sumDAPT$Normalised1 > 0.1 & sumDAPT$Normalised3 > 0.09 & sumDAPT$Normalised4> 0.2]$binary <- 1
names(sumDAPT)[names(sumDAPT)=="X.1"] <- 'Label'
Embryo1 <- subset(sumDAPT, Embryo == 1)
ordered_Embryo1 <- Embryo1[order(Embryo1$Label),]
write.csv(ordered_Embryo1, 'Embryo1_NMP_DAPT_Proliferating.csv')


############ Signal Intensity Graphs for NMPS ###############

intensity<- read.csv("data/sox.csv")
head(intensity)

sox<-subset(intensity,gene=="sox")
bra<-subset(intensity,gene=="bra")
tbx<-subset(intensity,gene=="tbx")

####### Sox Graph 
graph_sox <- ggplot(sox,aes(x=tail, y=count, fill=tail))
graph_sox + geom_boxplot()+
  scale_fill_manual(values =c("#00CC33", "#99FF99"))+ 
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("A", "B"),
                   labels=c("Trunk NMPs", "Tailbud NMPs"))+
  labs(x ="Location",
       y= "Normalised Soxb1b/c Intensity")

###### Bra Graph 
graph_bra <- ggplot(bra,aes(x=tail, y=count, fill=tail))
graph_bra + geom_boxplot()+
  scale_fill_manual(values =c("#CC0066", "#FF66CC"))+ 
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("A", "B"),
                   labels=c("Trunk NMPs", "Tailbud NMPs"))+
  labs(x ="Location",
       y= "Normalised Bra1/2 Intensity")

##### Tbx Graph 
graph_tbx <- ggplot(tbx,aes(x=tail, y=count, fill=tail))
graph_tbx + geom_boxplot()+
  scale_fill_manual(values =c("#FFCC33", "#FFFF66"))+ 
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("A", "B"),
                   labels=c("Trunk NMPs", "Tailbud NMPs"))+
  labs(x ="Location",
       y= "Normalised Tbx6/16 Intensity")

###### Stats 
lm_sox <- lm(count ~ tail, data = sox)
lm_bra <- lm(count ~ tail, data = bra)
lm_tbx <- lm(count ~ tail, data = tbx) 

##### Testing for Normality 
resid_count1<-residuals(lm_sox)
resid_count2<-residuals(lm_bra)
resid_count3<-residuals(lm_tbx)

shapiro.test(resid_count1)
shapiro.test(resid_count2)
shapiro.test(resid_count3)

# Data not normally distributed 

library(carData)
library(car)

### Testing for equal variance 
leveneTest(count ~ group, data = sox)
leveneTest(count ~ group, data = sox)
leveneTest(count ~ group, data = sox)

##### Significance 
kruskal.test(count ~ tail, data = sox)
kruskal.test(count ~ tail, data = bra)
kruskal.test(count ~ tail, data = tbx)

############## NMP Graphs #############

nmps <- read.csv("data/final_workbook.csv")
head(nmps)
library(ggplot2)

count$group<-as.factor(count$group)

####### Manual Count by group 
graph1 <- ggplot(nmps,aes(x=group, y=count, fill=nmp))
graph1 + geom_boxplot()+
  scale_fill_brewer(palette = "Greys")+
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("A", "B", "C", "D"),
                   labels=c("8-10", "10-12", "12-14", "14-16"))+
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8))+
  labs(x="Edu Labelled Time Period (hpf)",
       y= "Proliferating NMP Count")

########## Automatic NMP Count Graphs   
###### NMP Count 

graph2 <- ggplot(count,aes(x=group, y=total, fill=group ))
graph2 + geom_boxplot()+
  scale_fill_brewer(palette = "Greys")+
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("A", "B", "C", "D"),
                   labels=c("8-10", "10-12", "12-14", "14-16"))+
  labs(x="Edu Labelled Time Period (hpf)",
       y= "NMP Count")  

graph3 <- ggplot(count,aes(x=group, y=prol, fill=group ))
graph3 + geom_boxplot()+
  scale_fill_brewer(palette = "Greys")+
  theme(legend.position="none")+
  scale_x_discrete(breaks=c("A", "B", "C", "D"),
                   labels=c("8-10", "10-12", "12-14", "14-16"))+
  labs(x="Edu Labelled Time Period (hpf)",
       y= "Proliferating NMP Count")  

lm_count1 <- lm(count~ group, data = nmps)
lm_count2 <- lm(total~group, data = nmps)
lm_count3 <- lm(nmp~group, data= nmps)

shapiro.test(nmps$count)
shapiro.test(nmps$total)
shapiro.test(nmps$nmp)

# Data not normally distributed 

library(carData)
library(car)

## Testing for equal variance 
leveneTest(count ~ group, data = nmps)
leveneTest(total ~ group, data = nmps)
leveneTest(nmp ~ group, data = nmps)

## Significance 
kruskal.test(count ~ group, data = nmps)
kruskal.test(total ~ group, data = nmps)
kruskal.test(nmp ~ group, data = nmps)

####### NMP Count with trunk and tail  

trunk <- read.csv("data/final_trunk.csv")
head(trunk)

graph2 <- ggplot(trunk,aes(x=group, y=count, fill=tail))
graph2 + geom_boxplot()+
  scale_fill_brewer(palette = "Greys", labels=c("Trunk", "Tail"))+
  scale_x_discrete(breaks=c("A", "B", "C", "D"),
                   labels=c("8-10", "10-12", "12-14", "14-16"))+
  labs(fill="Position",
       x="Edu Labelled Time Period (hpf)",
       y= "Proliferating NMP Count") 
lm_count1 <- lm(count~ group, data = trunk)

resid_count1<-residuals(lm_count1)
shapiro.test(resid_count1)

#### Not Normally distributed 

library(carData)
library(car)

leveneTest(count ~ group, data = trunk)
##### Has equal variance 

kruskal.test(count ~ group, data= trunk)

####### Comparison of count #########

avg<- read.csv("data/count_avg.csv")
head(avg)

count<-as.numeric(avg$count)
library(ggplot2)

ggplot(avg, aes(x=group, y=count, fill=type))+
  geom_bar(position ="dodge",stat="identity")+
  scale_fill_brewer(palette = "Paired", labels=c("Manual Count", "Automated"))+
  scale_x_discrete(guide=guide_axis(angle=45),breaks=c("A", "B","C","D"),
                   labels=c("8-10", "10-12","12-14","14-16"))+
  labs(fill="Count Type",
       x="EdU labelled (hpf)",
       y="Average Proliferating NMP Count")

total<- read.csv("data/final_count.csv")
head(total)

A<-subset(total,group=="A")
B<-subset(total,group=="B")
C<-subset(total,group=="C")
D<-subset(total,group=="D")

#### Data normally distributed 

t.test(A$count, A$total,
       alternative = "two.sided", paired = TRUE)

t.test(B$count, B$total,
       alternative = "two.sided", paired = TRUE)

t.test(C$count, C$total,
       alternative = "two.sided", paired = TRUE)

t.test(D$count, D$total,
       alternative = "two.sided", paired = TRUE)