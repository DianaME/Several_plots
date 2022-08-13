#################################
###Preparing the GWAS analysis####
##################################
#gen is genetic data
#dat phenotype information

############creating PCA plots #################
#########PRINCIPAL COMPONENT ANALYSIS
library(NAM)
Adjusted_genotypes= snpQC(gen= gen, MAF=0.05, impute=TRUE, remove = TRUE)   
gen<- Adjusted_genotypes
gen<- as.matrix(gen)

K<- GRM(gen)

Spectral = eigen(K,symmetric = TRUE) ##espectral decomposition of a matrix
PCs = Spectral$vectors
write.csv(PCs, file="PCs_popstructure.csv")

Var_explained<- as.data.frame(Spectral$values)
Var_explained<- Var_explained %>% mutate(var.exp = `Spectral$values`/ sum(`Spectral$values`))
colnames(Var_explained)[1]<- "eigenvalues"
sum(Var_explained$var.exp) ##the sum should be equal to 1 
PC<- seq(1:329)
Var_explained<- cbind (PC, Var_explained)
write.csv(Var_explained, file="Varex_PCA_Accessions.csv")

##barplot of the variance explained of the first 10 PCs 
data<- Var_explained[1:10,-(2)]
data$var.exp<- data$var.exp*100
png("Variance_explained_byPC_.png", width = 700, height = 700)
barplot(data$var.exp, main="Explained Variance by PCs", names.arg = 1:10, xlab = "PC", ylab = "Percent variance", col="salmon")
dev.off()


# extract pc scores for first two component and add to dat dataframe
Genetic$pc1 <- Spectral$vectors[, 1] # indexing the first column

Genetic$pc2 <- Spectral$vectors[, 2]  # indexing the second column
list<- rownames(Genetic)
Genetic<- cbind(list,Genetic)
Genetic<- Genetic[,c(31691:31693,1:31690)]

colnames(Genetic)[4]<-"Accession"
library(RColorBrewer)
#cbbPalette <- c( "#8DD3C7", "#BEBADA","#FB8072", "#80B1D3",
Genetic$CLUSTER<- dat$Continent ##vector of the grouping factor

#Genetic$CLUSTER<- dat$Maturity_Group

library(RColorBrewer)
cbbPalette <- c( "green4",
                 "#6A3D9A", # purple
                 "gold1",
                 "skyblue2",  # lt pink
                 "#FDBF6F", # lt orange
                 "maroon", "orchid1", "deeppink1",  "steelblue4")
png("PCA_continent.png", width = 700, height = 500)

ggplot(data = Genetic, aes(x = pc1, y = pc2, color = CLUSTER)) + 
  
  # scale_fill_discrete()+
  scale_colour_manual(values=cbbPalette)+
  
  geom_hline(yintercept = 0, lty = 2) +
  
  geom_vline(xintercept = 0, lty = 2) +
  
  guides(color = guide_legend(title = "Continent"), shape = guide_legend(title = "Continent"))+
  
  #scale_shape_manual(values = c(15, 16)) +
  
  geom_point(alpha = 0.8, size = 4) +
  
  #stat_ellipse(geom="polygon", aes(fill = CLUSTER), 
  
  #            alpha = 0.2, 
  
  #           show.legend = FALSE, 
  
  #          level = 0.95) +
  
  xlab("PC 1 (20.09%)") + 
  
  ylab("PC 2 (7.4%)") +
  
  
  theme_minimal() +
  
  theme(panel.grid = element_blank(), 
        
        panel.border = element_rect(fill= "transparent"), text = element_text(size=16, face="bold"))



dev.off()


##genomic relationship matrix for 1096 subset of USDA germplasm#######################################################################
library(ComplexHeatmap)

library(scales)

color<-c("#FDE725FF", "#EDE51BFF" ,"#DDE318FF", "#CBE11EFF", "#BADE28FF", "#A9DB33FF", "#97D83FFF", "#87D549FF",
"#76D153FF", "#67CC5CFF", "#59C864FF", "#4CC26CFF", "#40BC72FF", "#35B779FF", "#2DB17EFF" ,"#25AC82FF",
 "#21A585FF", "#1FA088FF" ,"#1F998AFF" ,"#20938CFF" ,"white", "#24878EFF" ,"#26818EFF" ,"#297B8EFF",
"#2B748EFF", "#2E6F8EFF" ,"#31688EFF", "#34618DFF", "#375B8DFF", "#3A538BFF", "#3D4D8AFF", "#404588FF",
 "#433E85FF", "#453581FF", "#472E7CFF" ,"#482576FF", "#481D6FFF" ,"#481467FF", "#460B5DFF" ,"#440154FF")

f1 = colorRamp2(seq(min(K), max(K), length = 11), c("#FDE725FF",  "#C2DF23FF", "#51C56AFF", "#2BB07FFF", "#1E9B8AFF", "white", "#25858EFF", "#2D708EFF", "#38598CFF", "#433E85FF" ,"#482173FF" ))



png(filename="heatmap1.png",height=6,width=9,res=200,units="in")
par(mar=c(2,3,3,2))
Heatmap(K,
        name = "Genomic Relationship", #title of legend,
        row_names_gp = gpar(fontsize = 7),
        column_title = "GRM heatmap", 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        col=rev(viridis_pal()(12)))
dev.off()   

####cluster ###############################################################
show_col( )
library(NAM)
Genetic<-gen
GeneticDistance = Gdist(Genetic,method=6)

Tree = hclust(GeneticDistance,method = 'ward.D2')
png(filename="cluster.png",height=6,width=9,res=200,units="in")
plot(Tree,labels = FALSE)

rect.hclust(Tree, k=9, border="red")

dev.off() 

library(factoextra)
library(fpc)
library(NbClust)
group<- cutree(Tree, k=9) #33.86
group<- cutree(Tree, k=8) #34.77
group<- cutree(Tree, k=7) #35.93
group<- cutree(Tree, k=6) #37.18
group<- cutree(Tree, k=10) #33.08
group<- cutree(Tree, k=5) #38.79
group<- cutree(Tree, k=4) #41.055
group<- cutree(Tree, k=3) #43.97
group<- cutree(Tree, k=2) #48.45
group<- cutree(Tree, k=11) #32.45
group<- cutree(Tree, k=12) #31.85
group<- cutree(Tree, k=13) #31.31
group<- cutree(Tree, k=14) #30.77
cluster.stats(GeneticDistance, group)

x<- seq(2,14, by=1)
ss<- c(48.45,43.97,41.05,38.79,37.18,35.93,34.77,33.86,33.08,32.45,31.85,31.31,30.77)

png("Ball_Hall.png", width = 700, height = 500)
plot(x = x, y = ss, xlab="Cluster number", ylab="Ball Hall Criterion (within-group SS)")


dev.off()



##historgram of kinship values#################
#making the k matrix a vector

K_table <- as.matrix(K)
K_table[upper.tri(K_table, diag = TRUE)] <- NA
K_table<- c(K_table)
K_table<- na.omit(K_table)
K_table<- as.data.frame(K_table)

a<- seq(1:length(K_table$K_table))
K_table<- cbind(a, K_table)
colnames(K_table)<- c("number", "b")
#K_table<- as.data.frame(K)
#K_table$b<- abs(K_tabl$b)

library(dplyr)

# sample data
# condition

K_table<-filter(K_table, b<0.7)

png("Relative_K_coeff.png", width = 800, height = 600)
ggplot(K_table, aes(x= b)) +
  geom_density(aes(y = after_stat(count)), color="plum3", fill="plum3") +  xlab("Relative Kinship Coefficients") + scale_x_continuous(breaks = c(seq(-0.4,0.6, by=0.1)))+
  theme_minimal() +
  
  theme(panel.grid = element_blank(), 
        
        panel.border = element_rect(fill= "transparent"), text = element_text(size=16))
dev.off()

##summary of data ################################
country<- dat%>% group_by(Maturity_Group)%>% summarise(n())
country$perc<- country$`n()`/sum(country$`n()`)

write.csv(country, file="maturity_groupperc.csv")

#########
country<- dat%>% group_by(Country)%>% summarise(n())
country$perc<- country$`n()`/sum(country$`n()`)

write.csv(country, file="country_groupperc.csv")

#########
country<- dat%>% group_by(Year_acquisition)%>% summarise(n())
country$perc<- country$`n()`/sum(country$`n()`)

write.csv(country, file="yearofacqui_groupperc.csv")




# converting SNP in -1,0,1 system
#gen<-gen-1
##doing quality control in the genotypes
library(NAM)


#Subsetting pheno
##Pheno = Phenotype file 

#Converting first column in index #first column is usually acession or line name
Pheno <- data.frame(Pheno[,-1], row.names = Pheno[,1])
write.csv(Pheno, file="Pheno.csv")

#Extracting the number of SNP per chromosome
SNP<-colnames(gen)
Chr<-numeric(20)
SNP_name<-c("^Gm_01", "^Gm_02", "^Gm_03","^Gm_04", "^Gm_05", "^Gm_06", "^Gm_07", "^Gm_08", "^Gm_09",
            "^Gm_10", "^Gm_11", "^Gm_12","^Gm_13", "^Gm_14", "^Gm_15", "^Gm_16", "^Gm_17", "^Gm_18",
            "^Gm_19", "^Gm_20")
for (j in 1:20)
{
  Chr[j]<-length(grep(SNP_name[j], SNP))
}
sum(Chr)
write.csv(Chr, file="SNP_Chr.csv")

# number of SNP not matching. Scaffold presented ####what does this do?
Chrom = substr(SNP, 4, 5)
SNP<-as.data.frame(SNP)
B<-subset(Chrom,Chrom>20)

#Eliminating scaffold from SNPs list
#gen <- gen[,!grepl("^scaf", colnames(gen))]

#A <- SNP[grep("Gm_01|Gm_02|Gm_03|Gm_04|Gm_05|Gm_06|Gm_07|Gm_08|Gm_09|Gm_10|Gm_11|Gm_12|Gm_13|Gm_14|Gm_15|Gm_16|Gm_17|Gm_18|Gm_19|Gm_20",SNP), "Type"] <- "Chrom"

#Converting Geno in a matrix from data frame
gen= data.matrix(gen)
save(gen, file="Geno_ok.Rdata")

#Cleaning the memory
rm(list=ls())

#importing genome and genome
# importing genome
library(readxl)
Geno <- load("./Geno_ok.Rdata")
Pheno <- read.csv("Pheno.csv")
Pheno <- data.frame(Pheno[,-1], row.names = Pheno[,1])
Chr <- read.csv("SNP_Chr.csv")
#group<-read.csv("cluster.csv")
group<-read.csv("cluster.csv") #
group<- group[,-1]
group <- data.frame(group[,-1], row.names = group[,1])
colnames(group)<-"group"

##for other years


#Eigenvalue decomposition
library(NAM)
spectral = eigen(GRM(gen,T),T)
###Unique segments based on Eigenvalues
#No of unique segments
unseg<- sum(spectral$values>1)
#threshold
thres<- 0.01/unseg
-(log10(thres))
##Benjamini-Hochberg FDR: ??=0.05m×(1???FDR)
alpha=0.05/(31689*.95)
alpha= 1.660879e-06



################GWAs TRAIT 1
GWASR1= gwas2(y=Pheno[,1], gen=gen, chr= Chr, fam=NULL, cov=PCs) #cov = Chromosome18$Gm_18_54531027_T_C)
saveRDS(GWASR1, "GWASR1_1920.rds")
stat<- GWASR1[["PolyTest"]]
write.csv(stat, file = "Stat_gwas_R1_1920.csv")
#stat<- read.csv("./Stat_gwas_R1_17.CSV")
SNP_name<- paste("SNP", seq(1:nrow(stat)),"")
a<- data.frame(rownames(stat))
colnames(a)<-"name"
info<- a%>% separate(name,c("x","chr","Position", "y","z"))

Data<- as.data.frame(cbind(SNP_name, info$chr, info$Position, stat$pval))
colnames(Data)<- c("SNP", "chromosome","Position","logpval")
Data$Position<- as.numeric(as.character(Data$Position))
Data$logpval<-as.numeric(as.character(Data$logpval))
Data$p<- 10^(-Data$logpval) ##to obtain pval
Data<- Data[,-4]



library(CMplot)

CMplot(Data, plot.type="m", LOG10 = TRUE, ylim=NULL, threshold=c(thres, thres),threshold.lty=c(1,2),
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),signal.cex=c(1,1),
       signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       width=14,height=6)



###getting the significant SNPs
stat<-read.csv("Stat_gwas_R1.csv")
thr<- -log10(thres)
df<- stat[order(stat$pval, decreasing = TRUE),]
df<- df[df$pval >thr,c(1,12,14)]
df<- df%>% separate(X,c("x","chr","Position", "y","z"))  
df<-df[,c(2,3,6,7)]
Trait<-rep("R1", nrow(df))
df<-cbind(Trait,df)

###GWAs for  TRAIT 2
GWASR8= gwas2(y=Pheno$R8, gen=gen, chr= Chr, fam=NULL, cov=PCs) #cov = Chromosome18$Gm_18_54531027_T_C)
saveRDS(GWASR8, "GWASR8_1920.rds")
stat<-GWASR8[["PolyTest"]]
write.csv(stat, file = "Stat_gwas_R8_1920.csv")


qqman::qq(Data$p)

SNP_name<- paste("SNP", seq(1:nrow(stat)),"")
a<- data.frame(rownames(stat))
colnames(a)<-"name"
info<- a%>% separate(name,c("x","chr","Position", "y","z"))

Data<- as.data.frame(cbind(SNP_name, info$chr, info$Position, stat$pval))
colnames(Data)<- c("SNP", "chromosome","Position","logpval")
Data$Position<- as.numeric(as.character(Data$Position))
Data$logpval<-as.numeric(as.character(Data$logpval))
Data$p<- 10^(-Data$logpval) ##to obtain pval
Data<- Data[,-4]


CMplot(Data, plot.type="m", LOG10 = TRUE, ylim=NULL, threshold=c(thres, thres),threshold.lty=c(1,2),
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),signal.cex=c(1,1),
       signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       width=14,height=6)


stat<- read.csv("Stat_gwas_R8.csv")
thr<- -log10(thres)
df1<- stat[order(stat$pval, decreasing = TRUE),]
df1<- df1[df1$pval >thr,c(1,12,14)]
df1<- df1%>% separate(X,c("x","chr","Position", "y","z"))  
df1<-df1[,c(2,3,6,7)]
Trait<-rep("R8", nrow(df1))
df1<-cbind(Trait,df1)

####TRAIT 3#######################################################################

GWASRL= gwas2(y=Pheno$RL, gen=gen, chr= Chr, fam=NULL, cov=PC) #cov = Chromosome18$Gm_18_54531027_T_C)
saveRDS(GWASRL, "GWASRL_1920.rds")
stat<-GWASRL[["PolyTest"]]
write.csv(stat, file = "Stat_GWASRL_1920.csv")


SNP_name<- paste("SNP", seq(1:nrow(stat)),"")
a<- data.frame(rownames(stat))
colnames(a)<-"name"
info<- a%>% separate(name,c("x","chr","Position", "y","z"))

Data<- as.data.frame(cbind(SNP_name, info$chr, info$Position, stat$pval))
colnames(Data)<- c("SNP", "chromosome","Position","logpval")
Data$Position<- as.numeric(as.character(Data$Position))
Data$logpval<-as.numeric(as.character(Data$logpval))
Data$p<- 10^(-Data$logpval) ##to obtain pval
Data<- Data[,-4]


CMplot(Data, plot.type="m", LOG10 = TRUE, ylim=NULL, threshold=c(thres, thres),threshold.lty=c(1,2),
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),signal.cex=c(1,1),
       signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       width=14,height=6)

stat<- read.csv("Stat_GWASRP.csv")
thr<- -log10(thres)
df2<- stat[order(stat$pval, decreasing = TRUE),]
df2<- df2[df2$pval >thr,c(1,12,14)]
df2<- df2%>% separate(X,c("x","chr","Position", "y","z"))  
df2<-df2[,c(2,3,6,7)]
Trait<-rep("RL", nrow(df2))
df2<-cbind(Trait,df2)


##########significant SNPs for all three traits
TopSNP<- rbind(df,df1,df2)
write.csv(TopSNP, file="sig_SNPs.csv")


#####################################################################table to summary GWAs or all years###################
library(plyr)
library(readr)
library(readxl)
library(lubridate)
library(dplyr)
setwd("C:/folderwith_GWAs_data")
mydir = "C:/folderwith_GWAs_data"
myfiles = list.files(path=mydir, pattern="sig_SNPs*", full.names=TRUE)

dat = ldply(myfiles, read.csv)
b<-paste(dat$chr,dat$Position, sep="_")
dat<-cbind(dat,b)
dat<- dat[,c(1,2,3,4,7,5,6)]
colnames(dat)[5]<-"SNP"

write.csv(dat, file="alltopSNPs.csv")

dat<-read.csv("alltopSNPs.csv")
library(ggplot2)
dat$SNP <- factor(dat$SNP)
dat$Trait <- factor(dat$Trait)
dat$year <- factor(dat$year)
unique(dat$year)
levels(dat$year)<-c("2017", "2018",  "2019", "2020", "2017-18","2018-19", "2019-20",  "2017-19","2018-20", "2017-20")

colnames(dat)[c(1,6,8)]<- c("year", "Chr & Position", "SNP Effect")

png("sigSNPs.png", width = 1000, height = 1000)
ggplot(dat, aes(`year`, `Chr & Position`)) + geom_tile(aes(fill =`SNP Effect`)) + 
  scale_fill_viridis_c(na.value="transparent") + facet_grid(Trait~., scales = "free", space = "free") +
  theme_linedraw(base_size = 18)+ xlab("Years")+ theme(axis.text.x = element_text(angle = 90))

dev.off()


#############circular GWAs plot three traits ##############3
setwd("C:/all_traits ")
stat<-read.csv("Stat_gwas_Trait1.csv")
stat <- data.frame(stat[,-1], row.names = stat[,1])

SNP_name<- paste("SNP", seq(1:nrow(stat)),"")
a<- data.frame(rownames(stat))
colnames(a)<-"name"
info<- a%>% separate(name,c("x","chr","Position", "y","z"))

Data<- as.data.frame(cbind(SNP_name, info$chr, info$Position, stat$pval))
colnames(Data)<- c("SNP", "chromosome","Position","logpval")
Data$Position<- as.numeric(as.character(Data$Position))
Data$logpval<-as.numeric(as.character(Data$logpval))
Data$p<- 10^(-Data$logpval) 
#Data$p[Data$pl<0] = 0##to obtain pval
Data<- Data[,-5]
colnames(Data)[4]<-"trait1"


stat1<-read.csv("Stat_gwas_Trait2.csv")
stat1 <- data.frame(stat1[,-1], row.names = stat1[,1])

SNP_name<- paste("SNP", seq(1:nrow(stat1)),"")
a<- data.frame(rownames(stat1))
colnames(a)<-"name"
info<- a%>% separate(name,c("x","chr","Position", "y","z"))

Data1<- as.data.frame(cbind(SNP_name, info$chr, info$Position, stat1$pval))
colnames(Data1)<- c("SNP", "chromosome","Position","logpval")
Data1$Position<- as.numeric(as.character(Data1$Position))
Data1$logpval<-as.numeric(as.character(Data1$logpval))
Data1$p<- 10^(-Data1$logpval) ##to obtain pval
Data1<- Data1[,-5]
colnames(Data1)[4]<-"trait2"


stat2<-read.csv("Stat_GWAS_trait3.csv")
stat2 <- data.frame(stat2[,-1], row.names = stat2[,1])

SNP_name<- paste("SNP", seq(1:nrow(stat2)),"")
a<- data.frame(rownames(stat2))
colnames(a)<-"name"
info<- a%>% separate(name,c("x","chr","Position", "y","z"))
colnames(info)[2]<-"chromosome"

Data2<- as.data.frame(cbind(SNP_name, info$chr, info$Position, stat2$pval))
colnames(Data2)<- c("SNP", "chromosome","Position","logpval")
Data2$Position<- as.numeric(as.character(Data2$Position))
Data2$logpval<-as.numeric(as.character(Data2$logpval))
Data2$p<- 10^(-Data2$logpval) ##to obtain pval
Data2<- Data2[,-4]
colnames(Data2) [4]<-"trait3"

#final_Dat<- merge(info,dat, by=c("chromosome","Position"))

final_Dat<- merge(Data,Data1, by=c("chromosome","Position"))
final_Dat<-merge(final_Dat,Data2,by=c("chromosome","Position") )
#final_Dat<-final_Dat[,c(7,1,2,4,6,8)]
#final_Dat$SNP<-stat$X

write.csv(final_Dat, file="GWAsThreeTraits.csv")

final_Dat<-read.csv("GWAsThreeTraits.csv")
final_Dat<- final_Dat[,-1]
SNPs <- list(
  final_Dat$SNP[final_Dat$trait1<thres],
  final_Dat$SNP[final_Dat$trait2<thres],
  final_Dat$SNP[final_Dat$trait3<thres])

#SNPs[[3]]<- SNPs[[3]][-3]
####SNP density
#var(Final_DB$R8,na.rm = TRUE)
library(CMplot)
###circular manhattan plot
CMplot(final_Dat, type="p",plot.type="c", chr.labels=paste("Chr",c(1:20),sep = ""),r=0.4,cir.legend = TRUE, LOG10 = TRUE, ylim=NULL, threshold=c(thres),threshold.lty=c(1,2),
       threshold.lwd=c(1,1), cir.chr.h=1.5,amplify=TRUE,threshold.col=c("black","grey"),signal.line=1,signal.col=c("red","green"),chr.den.col="black",signal.cex=c(1,1),
       signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       width=14,height=6)



############multitrack plot
library(broman)
c<-brocolors("crayons")[c("Yellow Orange","Royal Purple","Forest Green","Blush","Olive Green","Pacific Blue","Atomic Tangerine","Orchid", "Chestnut","Cornflower")]
CMplot(final_Dat, plot.type="m",multracks=TRUE,col=c, threshold=c(thres,thres),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green","blue"),
       signal.cex=1, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1)#,
       #highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4)


### qqplot
CMplot(final_Dat,plot.type="q",col=c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),threshold=thres,
       ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",conf.int=TRUE,box=FALSE,multracks=
         TRUE,cex.axis=2,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,ylim=c(0,8),width=5,height=5)



################Local GWAs plots#################
Data$SNP<-final_Dat$SNP
library(devtools) 
library(ggstar)
####close up chromosome 4
chr4<-subset(final_Dat, chromosome=="4", select =col1:Lastcol)
chr4$Trait1<-as.numeric(as.character(chr4$Trait1))
chr4$Trait2<-as.numeric(as.character(chr4$Trait2))
chr4<- gather(chr4,Trait,logpval,5:6)#gathe the trait columns
chr4<-chr4[,-4]

#####add new positions
SNPID<-read.table("SNPID.txt")##inforamtion with the SNP positions
colnames(SNPID)<-c("SNP","Long_ID","SNP_name","chr","position2")

chr4<-merge(chr4,SNPID,by="SNP")
chr4<-chr4[,c(1,2,3,4,5,9)]
chr4$position2<- chr4$position2/1000000
ylim<- max(chr4$logpval)
write.csv(chr4,file="plot1.csv") ###add the positions of the genes that you want to plot
chr4<-read.csv("plot1.csv")
chr4<-chr4[,-1]
th<- -log10(thres)


#genes<-read.csv("GEN_IDsPlot.csv")
#genes$start<- genes$start/1000000
#genes$end<- genes$end/1000000

#genes$point<- (genes$start+genes$end)/2
#chr4gen<-subset(genes,Chr== 4)
library(ggrepel)

chr4$Trait<-as.factor(as.character(chr4$Trait))
colnames(chr4)[4]<-"Category"

c<-brocolors("crayons")[c("Yellow Orange","Royal Purple","Forest Green","Blush","Olive Green","Pacific Blue","Atomic Tangerine","Orchid", "Chestnut","Cornflower")]
c<-brocolors("crayons")[c("Atomic Tangerine","Cornflower","Pacific Blue","Yellow Orange")]

png("chr4zoom.png", width = 1000, height = 600)
ggplot(chr4, aes(x = position2, y = logpval, color = Category)) +
  geom_point(alpha = 0.9, size=4) +
  geom_text_repel(data=subset(chr4, logpval>th),aes(label=SNP),size=4)+
  geom_hline(yintercept = th, color = "black", linetype = "dashed") +
     scale_color_manual(values = c("R8"="#ffae42","RL"= "#1ca9c9","Glyma Gene"="#7851a9","E genes"="#de5d83" ) )+
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim+0.5)) +
    labs(x = "Chromosome 4 Position Mb", 
       y = "-log10(p-value)") + 
      theme_minimal() +
  theme(text = element_text(size = 18))
dev.off()



###chromosome 10 #ohter chromosome

chr10<-subset(final_Dat, chromosome=="10", select =SNP:R1)

SNPID<-read.table("SNPID.txt")
colnames(SNPID)<-c("SNP","Long_ID","SNP_name","chr","position2")

chr10<-merge(chr10,SNPID,by="SNP")
chr10<-chr10[,c(1,2,3,4,8)]
chr10$position2<- chr10$position2/1000000
colnames(chr10)[4]<-"logpval"
ylim<- max(chr10$logpval)
write.csv(chr10,file="plot2.csv")
chr10<-read.csv("plot2.csv")
chr10<-chr10[,-1]
th<- -log10(thres)


library(ggrepel)

chr10$Trait<-as.factor(as.character(chr10$Trait))
colnames(chr10)[4]<-"Category"


c<-brocolors("crayons")[c("Atomic Tangerine","Cornflower","Pacific Blue","Yellow Orange")]

png("chr10zoom.png", width = 900, height = 600)
ggplot(chr10, aes(x = position2, y = logpval, color = Category)) +
  geom_point(alpha = 0.9, size=4) +
  geom_text_repel(data=subset(chr10, logpval>th),aes(label=SNP),size=4)+
  geom_hline(yintercept = th, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("R1"=  "#ffae42","Glyma Gene"="#7851a9","E genes"="#de5d83" ) )+
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim+2)) +
  labs(x = "Chromosome 10 Position Mb", 
       y = "-log10(p-value)") + 
  theme_minimal() +
  theme(text = element_text(size = 18))
dev.off()



#################histograms#################################################################################
library(ggplot2)
library(factoextra)
library(ggpubr)                                                                                                                            
##load blues all years

library(plyr)
Final_DB<- rbind.fill(year1,year2,year3,year4)
var(Final_DB$RL, na.rm = TRUE)

BLUP_ALL<- Final_DB
a<- BLUP_ALL[,c(1,2,3,4)]##getting columns for each trait individually
a$Trait<- rep("Trait1",nrow(a))
colnames(a)[4]<-"value"

b<- BLUP_ALL[,c(1,2,3,5)]
b$Trait<- rep("Trait2",nrow(b))
colnames(b)[4]<-"value"

c<- BLUP_ALL[,c(1,2,3,6)]
c$Trait<- rep("Trait3",nrow(c))
colnames(c)[4]<-"value"

BLUP_ALL<-rbind(a,b,c)
BLUP_ALL$value<-as.numeric(as.character(BLUP_ALL$value))
BLUP_ALL<- na.omit(BLUP_ALL)
min(BLUP_ALL$value)
max(BLUP_ALL$value)
BLUP_ALL$Year<- as.factor(BLUP_ALL$Year)
BLUP_ALL$Trait<- as.factor(BLUP_ALL$Trait)
BLUP_ALL$Location<- as.factor(BLUP_ALL$Location)

write.csv(BLUP_ALL, file="filehistall_final.csv")
BLUP_ALL<- read.csv('filehistalls_final.csv')#check plots


png("histograms1allwithlocation.png", width = 1000, height = 900)

ggplot(BLUP_ALL, aes(`value`, fill= `Location`)) + geom_histogram(color="Black", alpha=0.5, position = "identity", binwidth = 3) + 
  scale_fill_viridis_d()+
  labs(x="Days After Planting (DAP)")+
  #scale_x_continuous(breaks = seq(31,168,by=20))+
  facet_grid(Year~Trait, scales='free') +
  theme_bw(base_size = 18)+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 14), strip.text.y = element_text(size = 14))
  
dev.off()


#########barplot of the most significant SNPs in four and 10 ##########
#########boxplor of the most significant SNPs in four and 10 ##########
SNP_list<-c("Gm_03_38448001_C_T","Gm_10_40908884_G_A","Gm_04_16183920_T_C","Gm_04_37126858_A_G","Gm_04_16031274_G_A","Gm_04_36874657_C_T", "Gm_04_37010886_T_C", "Gm_04_37078558_G_A")
gen<-as.data.frame(gen)
colnames(gen)[1:5]

col.num <- which(colnames(gen) %in% SNP_list)

Genetic<- gen[,c(col.num)]
Genetic$Accession<- rownames(Genetic)
sum(Genetic$Gm_10_44734696_C_T)

Accessions<-read.csv("Alleles_candidateG.csv")
Final<- merge(Accessions,Genetic, by="Accession")


write.csv(Final,"AllelesandSNPs.csv")
final<- read.csv("AllelesandSNPS.csv")
final<- read.csv("AllelesandSNPs.csv")
final<-final[,c(2:14)]
final<- merge(Final_DB,final, by="Accession")
Final$Gm_04_37078558_G_A<- final$Gm_04_37078558_G_A
l$Gm_04_37126858_A_G<-replace(Final$Gm_04_37126858_A_G, Final$Gm_04_37126858_A_G==0, "Major")
Final<-Final[,c(4:18)]
colnames(Final)[c(10:15)]<-c("04_16031274","04_16183920","04_36874657","04_37010886","04_37078558","04_37126858")


Long<- gather(Final,SNP,Allele,4:15)
#Long<- gather(Long,Trait,Value,4:6)
Long<- Long[,c(4:8)]
Long<-reshape2::melt(Long, id.vars=c("SNP","Allele"))

write.csv(Long, file="candidat_Long.csv")
Long<-read.csv("candidateAlle_Long.csv")
Long<-Long[,-1]
Long$SNP<-replace(Long$SNP, Long$SNP== "Gm_04_16031274_G_A", "04_16031274_G_A")

library(dplyr)
library(tidyverse)
library(rstatix)
library(ggpubr)
Long$Allele<-as.factor(as.character(Long$Allele))
#Long$Allele<-replace(Long$Allele, Long$Allele== "Major", "Ref_Allele")
compare_means(value~Allele, data=Long, method = "t.test")
Long$Trait<-as.factor(as.character(Long$Trait))
mycomparisons<- list( c("Ref_Allele", "Alt_Allele") )
levels(Long$SNP)<- c("E8.SNP1", "E8.SNP2","E8.SNP3","Glyma.04g125700","Glyma.04g159300","Glyma.04g163100",
                      "04_16031274" ,"04_16183920", "04_36874657","04_37010886", "04_37078558", "04_37126858")
levels(Long$Allele)<-c("Ref_Allele","Alt_Allele")
Long<-na.omit(Long)
Long$variable<-as.factor(as.character(Long$variable))

###R1#######################################################################
LongR1<-subset(Long,variable=="R1",SNP:value)
#LongR1<-subset(LongR1,SNP=="Gm_10_44972284_T_C",SNP:value)

stat.test <- LongR1 %>%
  group_by(SNP) %>%
  t_test(value ~ Allele) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% add_xy_position(fun="mean_sd",x = "Allele")
stat.test

# Create a bar plot with error bars (mean +/- sd)
bp <- ggbarplot(LongR1, x = "Allele", y = "value", add = "mean_sd", 
                fill = "Allele", facet = c("variable","SNP")
)+
  scale_fill_brewer(palette = "Dark2") +
  labs(y="DAP")+
  theme(axis.text.x=element_blank())+
  stat_pvalue_manual(stat.test, hide.ns = TRUE, tip.length = 0, step.increase = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))


###R8######################################################################
LongR8<-subset(Long,variable=="R8",SNP:value)
LongR8<-filter(LongR8,LongR8$SNP != "Gm_10_44972284_T_C")

stat.test <- LongR8 %>%
  group_by(SNP) %>%
  t_test(value ~ Allele) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% add_xy_position(fun="mean_sd",x = "Allele")
stat.test

# Create a bar plot with error bars (mean +/- sd)
bp1<- ggbarplot(LongR8, x = "Allele", y = "value", add = "mean_sd", 
                fill = "Allele", facet = c("variable","SNP")
)+
  scale_fill_brewer(palette = "Dark2") +
  labs(y="DAP")+
  theme(axis.text.x=element_blank())+
  stat_pvalue_manual(stat.test, hide.ns = TRUE, tip.length = 0, step.increase = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

###RL
LongRL<-subset(Long,variable=="RL",SNP:value)
LongRL<-filter(LongRL,LongRL$SNP != "Gm_10_44972284_T_C")

stat.test <- LongRL %>%
  group_by(SNP) %>%
  t_test(value ~ Allele) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% add_xy_position(fun="mean_sd",x = "Allele")
stat.test

# Create a bar plot with error bars (mean +/- sd)
bp3<- ggbarplot(LongRL, x = "Allele", y = "value", add = "mean_sd", 
                fill = "Allele", facet = c("variable","SNP")
)+
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_blank())+
  labs(y="DAP")+
  stat_pvalue_manual(stat.test, hide.ns = TRUE, tip.length = 0, step.increase = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

png("TtestEgnes_SNPs.png", width = 1200, height = 900)
ggarrange(bp,bp1,bp3,labels = c( "A","B", "C"), nrow = 3)

dev.off()


