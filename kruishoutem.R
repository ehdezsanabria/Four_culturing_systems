##### THIS IS THE CODE FOR THE FINAL ANALYSES OF OLIVER#####

library(stringr)

# Provide your data path
data <- read.csv2("kruis18.csv")

# If you want Genus instead of family,
#change accordingly in data$ and ask Ruben for further assistance 
data.sorted <- data[with(data, order(str_trim(data$Genus))),]

# Provide number of samples
number_samples = 152

#### Provide starting column of samples ####
start_col = 2

#### row number of last sample !!!!
end_col = number_samples+start_col-1

#### str_trim function removes all spaces to avoid classifying same OTU 
###as diff Genera ####
Genus <-unique(str_trim(data.sorted$Genus))
positions <- matrix(nrow=length(Genus), ncol=3)

#### Find the start and end position of all Genera in the dataframe ####
for(i in 1:length(Genus)){
  positions[i,1] = Genus[i]
  temp = which(sapply(data.sorted$Genus, 
                      function(x) any(str_trim(x) == str_trim(paste(Genus[i])))))
  positions[i,2] = temp[1]
  positions[i,3] = temp[(length(temp))] 
}
Samples.0 = colnames(data[,start_col:(end_col)])


#### Add all the rows (from start to end) of the same Genus together ####
#### end_col + as.integer() function to force integer values of positions-> this I added!!!!
data.pool<-data.frame()
for(i in 1: length(Genus)){
  temp2 = colSums(data.sorted[as.integer(positions[i,2]):as.integer(positions[i,3]),
                              (start_col:(end_col))])
  data.pool = rbind(data.pool,temp2)
}
data.pool = cbind(data.pool, Genus)
Samples = c(Samples.0, "Genus")
colnames(data.pool) = Samples

#Write into an Excel file
write.csv2(file="pooled kruis18.csv",data.pool)

####BIODIVERSITY INDICES####

library(vegan)
library(xtable)

#load data for all samples
mydata <- read.csv2('indiceskruis.csv',
                    stringsAsFactors = FALSE)

#set up data as data matrix
counts <- as.matrix(mydata[-1])

#alpha diversity indices
shannon <- diversity(counts,index="shannon", MARGIN=1, base=exp(1))
simpson <- diversity(counts,index="simpson", MARGIN=1, base=exp(1))
inversesimp <- diversity(counts,index="inv", MARGIN=1, base=exp(1))

#species richness and eveness indices
totalspecies<-specnumber(counts)
Pielou <- shannon/log(totalspecies)
alpha <- fisher.alpha(counts)

#Export indices to Excel

library(xlsx)
library(rJava)
write.xlsx(Pielou, file="indicesKruis.xlsx",
           sheetName="Pielou", append=FALSE)
write.xlsx(shannon, file="indicesKruis.xlsx", sheetName="Shannon", 
           append=TRUE)
write.xlsx(simpson, file="indicesKruis.xlsx", sheetName="Simpson", 
           append=TRUE)
write.xlsx(alpha, file="indicesKruis.xlsx", sheetName="alpha", 
           append=TRUE)
write.xlsx(inversesimp, file="indicesKruis.xlsx", sheetName="InverseSimpson", 
           append=TRUE)
write.xlsx(totalspecies, file="indicesKruis.xlsx", sheetName="TotalSpecies", 
           append=TRUE)

####testing dissimilarities in the distance matrix of the samples####

library(vegan)
library(phyloseq)
library(doParallel)
library(ggplot2)
library(foreach)
library(picante)


# read in data of raw counts without taxonomy
mydata <- read.csv2('permanovakruis.csv',
                    stringsAsFactors=FALSE)
counts <- as.matrix(mydata[-1])


#read data of counts with taxonomical identification, 
#this is better to use
forbarplot <- read.csv2("barplotkruis.csv",
                        stringsAsFactors = FALSE)

keeps <- c("Genus",
           grep("\\d",names(forbarplot), value = TRUE))
famdata <- forbarplot[keeps]


mydata$type <- factor(paste(substr(mydata$Sample,1,3)))

distmatrix <- vegdist(counts, method="bray")
anova(betadisper(distmatrix, mydata$type))


# Prepare family data
fammatrix <- do.call(cbind,
                     by(famdata[-1],famdata$Genus,
                        colSums))

famdist <- vegdist(fammatrix, method = "bray")

anova(betadisper(famdist, mydata$type))


# Check where's the problem

plot(meta11<-metaMDS(distmatrix),type="n",
     main="Bray-Curtis distances among bacterial genera")

points(meta11,select=which(mydata$type=="GF1"),col="mistyrose",pch=16)
points(meta11,select=which(mydata$type=="GF2"),col="lightpink",pch=16)
points(meta11,select=which(mydata$type=="GF3"),col="coral",pch=16)
points(meta11,select=which(mydata$type=="GF4"),col="magenta",pch=16)
points(meta11,select=which(mydata$type=="GF5"),col="purple",pch=16)
points(meta11,select=which(mydata$type=="GF6"),col="slateblue",pch=16)
points(meta11,select=which(mydata$type=="GF7"),col="blue",pch=16)
points(meta11,select=which(mydata$type=="GF8"),col="navy",pch=16)

points(meta11,select=which(mydata$type=="GO1"),col="mistyrose",pch=17)
points(meta11,select=which(mydata$type=="GO2"),col="lightpink",pch=17)
points(meta11,select=which(mydata$type=="GO3"),col="coral",pch=17)
points(meta11,select=which(mydata$type=="GO4"),col="magenta",pch=17)
points(meta11,select=which(mydata$type=="GO5"),col="purple",pch=17)
points(meta11,select=which(mydata$type=="GO6"),col="slateblue",pch=17)
points(meta11,select=which(mydata$type=="GO7"),col="blue",pch=17)
points(meta11,select=which(mydata$type=="GO8"),col="navy",pch=17)

points(meta11,select=which(mydata$type=="SA1"),col="mistyrose",pch=18)
points(meta11,select=which(mydata$type=="SA2"),col="lightpink",pch=18)
points(meta11,select=which(mydata$type=="SA3"),col="coral",pch=18)
points(meta11,select=which(mydata$type=="SA4"),col="magenta",pch=18)
points(meta11,select=which(mydata$type=="SA5"),col="purple",pch=18)
points(meta11,select=which(mydata$type=="SA6"),col="slateblue",pch=18)
points(meta11,select=which(mydata$type=="SA7"),col="blue",pch=18)
points(meta11,select=which(mydata$type=="SA8"),col="navy",pch=18)

points(meta11,select=which(mydata$type=="SP1"),col="mistyrose",pch=15)
points(meta11,select=which(mydata$type=="SP2"),col="lightpink",pch=15)
points(meta11,select=which(mydata$type=="SP3"),col="coral",pch=15)
points(meta11,select=which(mydata$type=="SP4"),col="magenta",pch=15)
points(meta11,select=which(mydata$type=="SP5"),col="purple",pch=15)
points(meta11,select=which(mydata$type=="SP6"),col="slateblue",pch=15)
points(meta11,select=which(mydata$type=="SP7"),col="blue",pch=15)
points(meta11,select=which(mydata$type=="SP8"),col="navy",pch=15)

adonis(distmatrix ~ mydata$type)

# Check where's the problem, THIS IS THE CODE FOR THE PERMANOVA WITH LABELS
plot(meta11<-metaMDS(famdist),type="n",
     main="Chao distances by family")


points(meta11,select=which(mydata$type=="GF1"),col="mistyrose",pch=16)
points(meta11,select=which(mydata$type=="GF2"),col="lightpink",pch=16)
points(meta11,select=which(mydata$type=="GF3"),col="coral",pch=16)
points(meta11,select=which(mydata$type=="GF4"),col="magenta",pch=16)
points(meta11,select=which(mydata$type=="GF5"),col="purple",pch=16)
points(meta11,select=which(mydata$type=="GF6"),col="slateblue",pch=16)
points(meta11,select=which(mydata$type=="GF7"),col="blue",pch=16)
points(meta11,select=which(mydata$type=="GF8"),col="navy",pch=16)

points(meta11,select=which(mydata$type=="GO1"),col="mistyrose",pch=17)
points(meta11,select=which(mydata$type=="GO2"),col="lightpink",pch=17)
points(meta11,select=which(mydata$type=="GO3"),col="coral",pch=17)
points(meta11,select=which(mydata$type=="GO4"),col="magenta",pch=17)
points(meta11,select=which(mydata$type=="GO5"),col="purple",pch=17)
points(meta11,select=which(mydata$type=="GO6"),col="slateblue",pch=17)
points(meta11,select=which(mydata$type=="GO7"),col="blue",pch=17)
points(meta11,select=which(mydata$type=="GO8"),col="navy",pch=17)

points(meta11,select=which(mydata$type=="SA1"),col="mistyrose",pch=18)
points(meta11,select=which(mydata$type=="SA2"),col="lightpink",pch=18)
points(meta11,select=which(mydata$type=="SA3"),col="coral",pch=18)
points(meta11,select=which(mydata$type=="SA4"),col="magenta",pch=18)
points(meta11,select=which(mydata$type=="SA5"),col="purple",pch=18)
points(meta11,select=which(mydata$type=="SA6"),col="slateblue",pch=18)
points(meta11,select=which(mydata$type=="SA7"),col="blue",pch=18)
points(meta11,select=which(mydata$type=="SA8"),col="navy",pch=18)

points(meta11,select=which(mydata$type=="SP1"),col="mistyrose",pch=15)
points(meta11,select=which(mydata$type=="SP2"),col="lightpink",pch=15)
points(meta11,select=which(mydata$type=="SP3"),col="coral",pch=15)
points(meta11,select=which(mydata$type=="SP4"),col="magenta",pch=15)
points(meta11,select=which(mydata$type=="SP5"),col="purple",pch=15)
points(meta11,select=which(mydata$type=="SP6"),col="slateblue",pch=15)
points(meta11,select=which(mydata$type=="SP7"),col="blue",pch=15)
points(meta11,select=which(mydata$type=="SP8"),col="navy",pch=15)


# And for the family

adonis(famdist ~ mydata$type)

####the following are for the PERMANOVA and betadispersion####
library(vegan)
### Read data
kruis <- read.csv2("betadispersion.csv",header=TRUE,stringsAsFactors = FALSE,dec=".")

### Check for homogeneity of variance
betaDisp<-betadisper(dist(kruis[,-(1:3)]),gr=kruis$TRT)
TukeyHSD(betaDisp)

betaDisp<-betadisper(dist(kruis[,-(1:3)]),gr=kruis$TPT)
TukeyHSD(betaDisp)

### Test if each factor is important for the relative abundances
adonis2(kruis[,-(1:3)]~kruis$TRT,method="bray",perm=999)

adonis(kruis[,-(1:3)]~kruis$TPT,method="bray")

### Test all interactions assuming that each factor 
### are also accounted for
adonis(kruis[,-(1:3)]~kruis$TRT*kruis$TPT,method="bray")

library(FactoMineR)
library(missMDA)

# Read in data for bulk soil
kruis<- read.csv2('mfabag.csv',header=TRUE,row.names = 1, 
                   sep=";",dec='.') 
#names <-colnames(kruis[1:198])
# get properties of the data

Genus <- colnames(kruis[-(1:2)]) # First  columns are descriptors
ngenus <- length(Genus)       # the number of families
#This -1 is just in case that the number of columns doesn't match. The number 3 indicates
#the starting column
supvar <- colnames(kruis)[1:2] #This -1 is just in case 
nvar <- ncol(kruis)

nsupvar <- length(supvar)

# Convert to factors
kruis$TRT <- factor(kruis$TRT)
kruis$TPT <- factor(kruis$TPT)


# Do the analysis weighing physical and chemical characteristics
res <- MFA(kruis,
           group = c(1,1,1,1,2,7,1,187),ncp=5,
           type = c("n","n","s","s","s","s","s","s"),
           name.group = c("Treatment","Time","pH", "EC",
                          "Nitrogen","Minerals","plant length","Bacteria"),
           num.group.sup = 1:2)

#Do the analysis taking separately all the factors and the physical variables
res <- MFA(kruis,
           group = c(1, nvar),
           type = c(rep('n',nsupvar),rep('s', ngenus)),
           name.group = names(kruis),
           num.group.sup = 1:2)

# Select the variables that contribute the most to total variance 
#This is the code where I selected the most important for each dimension
#Dim 1 was the lupine and Dim 3 was the tomato
plot(res, axes=c(1,2), choix="var", habillage="group", cex=0.8, 
     select= "contrib 20")
plot(res, axes=c(1,2), choix="var", habillage="group", 
     cex=0.8, select= c("pH","Pyrinomonas","Bythopirellula",
                        "Acidobact_Gp6","Vulgatibacter","Acidobact_Gp16",
                        "Actinoplanes","Hypomicrobium","Gemmatimonas",
                        "Unclas_Bacteria","Dongia","plantlength","Ca","P",
                        "SO4","Na","Mg","Cl","NO3.N","NH4.N","EC",
                        "K","Anammoximicrobium","Rhizobium"))

# Ellipses can be plotted around the samples that are significantly similar 95% confidence
plotellipses(res,keepvar=c(1:3),axes=c(1,2),keepnames=TRUE,pch = 15:19)
plot(res, axes=c(1,2),choix="group", cex=0.8,ylim=c(0,1.0),xlim=c(0,0.5))

# Create an object containing the P value and the correlation value of each variable
# to the corresponding axis
stats<-dimdesc(res,axes=1:5,proba=0.8)
summary(res, nbelements=Inf,ncp=5)
write.infile(stats, file="statskruisbag.csv", sep=" ", nb.dec=7)


library(FactoMineR)
library(missMDA)

# Read in data for bulk soil
kruis<- read.csv2('mfasoil.csv',header=TRUE,row.names = 1, 
                  sep=";",dec='.') 
#names <-colnames(kruis[1:198])
# get properties of the data

Genus <- colnames(kruis[-(1:2)]) # First  columns are descriptors
ngenus <- length(Genus)       # the number of families
#This -1 is just in case that the number of columns doesn't match. The number 3 indicates
#the starting column
supvar <- colnames(kruis)[1:2] #This -1 is just in case 
nvar <- ncol(kruis)

nsupvar <- length(supvar)

# Convert to factors
kruis$TRT <- factor(kruis$TRT)
kruis$TPT <- factor(kruis$TPT)


# Do the analysis weighing physical and chemical characteristics
res <- MFA(kruis,
           group = c(1,1,1,1,2,7,1,196),ncp=5,
           type = c("n","n","s","s","s","s","s","s"),
           name.group = c("Treatment","Time","pH", "EC",
                          "Nitrogen","Minerals","plant length","Bacteria"),
           num.group.sup = 1:2)

#Do the analysis taking separately all the factors and the physical variables
res <- MFA(kruis,
           group = c(1, nvar),
           type = c(rep('n',nsupvar),rep('s', ngenus)),
           name.group = names(kruis),
           num.group.sup = 1:2)

# Select the variables that contribute the most to total variance 
#This is the code where I selected the most important for each dimension
#Dim 1 was the lupine and Dim 3 was the tomato
plot(res, axes=c(1,2), choix="var", habillage="group", cex=0.8, 
     select= "contrib 20")
plot(res, axes=c(1,2), choix="var", habillage="group", 
     cex=0.8, select= c("pH","Pyrinomonas","Bythopirellula",
                        "Acidobact_Gp6","Vulgatibacter","Acidobact_Gp16",
                        "Actinoplanes","Hypomicrobium","Gemmatimonas",
                        "Unclas_Bacteria","Dongia","plantlength","Ca","P",
                        "SO4","Na","Mg","Cl","NO3.N","NH4.N","EC",
                        "K","Anammoximicrobium","Rhizobium"))

# Ellipses can be plotted around the samples that are significantly similar 95% confidence
plotellipses(res,keepvar=c(1:2),axes=c(1,2),keepnames=TRUE,pch = 15:19)
plot(res, axes=c(1,2),choix="group", cex=0.8,ylim=c(0,1.0),xlim=c(0,0.5))

# Create an object containing the P value and the correlation value of each variable
# to the corresponding axis
stats<-dimdesc(res,axes=1:5,proba=0.8)
summary(res, nbelements=Inf,ncp=5)
write.infile(stats, file="statskruissoil.csv", sep=" ", nb.dec=7)



library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)




library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)

#### Start with Spiec-Easi for the total network both treatments####
mydata <- read.csv2('GO8network.csv',
                    stringsAsFactors = FALSE)
name.sorted <- mydata[-1]

#Used SparCC, the graph is awuful but better than using the other code with the
#glasso

sp.est <- sparcc(name.sorted)
sp.net <- abs(sp.est$Cor)>0.9
diag(sp.net) <-0

###Select the following all together to run###
ig.sp <-adj2igraph(sp.est$Cor*sp.net,rmEmptyNodes=TRUE,diag=TRUE,
                   edge.attr=list(weights=factor(weight.bin)),
                   vertex.attr=list(Label=colnames(name.sorted),
                                    names=colnames(name.sorted)))

weight.bin <- c(); weight.bin[E(ig.sp)$weight>0] <- "Positive"; 
weight.bin[E(ig.sp)$weight<0] <- "Negative"
ig.sp <-adj2igraph(sp.est$Cor*sp.net,rmEmptyNodes=TRUE,diag=TRUE,
                   edge.attr=list(weights=weight.bin),
                   vertex.attr=list(Label=colnames(name.sorted),
                                    names=colnames(name.sorted)))
weight.bin <- c(); weight.bin[E(ig.sp)$weight>0] <- "Positive"; 
weight.bin[E(ig.sp)$weight<0] <- "Negative"

## Change weights (correlation coef) to absolute values so that negative edges are displayed 
## with equal thickness as positive correlations

E(ig.sp)$weight <- abs(E(ig.sp)$weight)

###Export for gephi or cytoscape
write_graph(ig.sp,file=" network_Got8.sp3.graphml",
            format="graphml")

##PHYSICAL AND CHEMICAL CHARACTERISTICS MFA
library (lattice)
library (gplots)
library (heatmap.plus)
library(plotrix)
library (gdata)
library (stats)
library (ggplot2)
library (RColorBrewer)
library(FactoMineR)
library(missMDA)


####figure MFA with the time and treatment as supplementary categories

#Load your file
test<-read.csv2('pFLA data2.csv',header=TRUE,dec='.')

#This code is for running the MFA, change at the end which are the 
#categorical variables
res<- MFA(test,group=c(1,1,1,3,5,1,2,2,7), type=c("n","n","s","s","s",
                                                     "s","s","s","s"),
          ncp=5,name.group=c("Trt","time","Total","Bacterial", "Fungi","Protozoa",
                             "pH","Nitrogen","Minerals"),
          num.group.sup=c(1:2))

#Plot each data point and label according to treatment
plot(res, choix="ind",habillage="Trt", cex=0.5,autoLab = "auto")

#Select the variables that contribute the most (six in this case)
plot(res, choix="var", habillage="group", cex=0.8, select= "contrib 6")

#Summary of the statistics
summary(res)

#Values of the correlations with each axis
stats<-dimdesc(res, axes=1:5,proba=0.8)
write.infile(stats, file="kruishoutem_stats.csv", sep=" ", nb.dec=7)

#Plot 95% confidence intervals around categories (treatments or time points)
plotellipses(res,keepvar=c(1,2),keepnames=TRUE)

#####using HAEMERLINCK####
#Load your file
test<-read.csv2('pFLA data.csv',header=TRUE,dec='.')

#This code is for running the MFA, change at the end which are the 
#categorical variables
res<- MFA(test,group=c(1,1,1,3,5,1,2,2,7), type=c("n","n","s","s","s",
                                                  "s","s","s","s"),
          ncp=5,name.group=c("Trt","time","Total","Bacterial", "Fungi","Protozoa",
                             "pH","Nitrogen","Minerals"),
          num.group.sup=c(1:2))

#Plot each data point and label according to treatment
plot(res, choix="ind",habillage="Trt", cex=0.8,autoLab = "auto",axes=c(1,4),
     lab.ind=F, lab.var=F)

#Select the variables that contribute the most (six in this case)
plot(res, choix="var", habillage="group", cex=0.8, select= "contrib 6")

#Summary of the statistics
summary(res)

#Values of the correlations with each axis
stats<-dimdesc(res, axes=1:5,proba=0.8)
write.infile(stats, file="kruishoutemHAMERLINCK_stats.csv", sep=" ", nb.dec=7)

#Plot 95% confidence intervals around categories (treatments or time points)
plotellipses(res,keepvar=c(1,2),axes=c(1,4),keepnames=TRUE,pch = 15:19)


####RCCA###
library(mixOmics)

X <- read.csv2('bacterial.csv',dec='.') #file with bacterial variables
Y <-read.csv2('chemical.csv',dec='.') #file with chemical variables

#create a grid to plot the network, these values don't need to be changed
grid1 <- seq(0, 0.2, length = 51)
grid2 <- seq(0.0001, 0.2, length = 51)
cv.score<-tune.rcc(X,Y, grid1=grid1,grid2=grid2, plt= FALSE,validation = "loo")

#This command is to open an additional window to show the plots
windows()

#this is the code that worked with more links and changing the threshold
oliver.res <- rcc(X,Y,ncomp=5,lambda1 = 0.064,lambda2 = 0.008)
network(oliver.res, comp = 1:5, threshold = 0.5)

#thisis the threshold I sent for publication in the crazy roots paper
network(oliver.res, comp = 1:3,threshold = 0.65,
        color.node = c("mistyrose", "lightcyan"),
        shape.node = c("rectangle", "rectangle"), 
        color.edge = color.jet(100),
        lty.edge = "solid", lwd.edge = 2, 
        show.edge.labels = FALSE)

#Close the window of the plots
dev.off()
