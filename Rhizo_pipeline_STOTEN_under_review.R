#Data preparation
#Read table
library(data.table)
library(vegan)
#KEGG table input
kegg_raw_new<-fread(input="J://shangzhuang/kegg/kegg_raw.csv", sep=',', header=TRUE) 
#check colsums and select the value of minumun sum for rarefy
colSums(kegg_raw_new[,7:30])
#Rarefy
kegg_table<-kegg_raw_new[,7:30]
kegg_table_rarefy<-as.data.frame(t(rrarefy(t(kegg_otu),110000)))
rarefy_kegg<-cbind(kegg_raw_new[,1:6],kegg_otu_rarefy)
#check colsums after rarefy
colSums(rarefy_kegg[,7:30])
#MN2P2-1    MN2P2-2    MN2P2-3    MN2P2-4     N0P0-1     N0P0-2     N0P0-3     N0P0-4     N2P2-1     N2P2-2     N2P2-3 
#110000     110000     110000     110000     110000     110000     110000     110000     110000     110000     110000 
#N2P2-4 RS-MN2P2-1 RS-MN2P2-2 RS-MN2P2-3 RS-MN2P2-4   RSN0P0-1   RSN0P0-2   RSN0P0-3   RSN0P0-4  RS-N2P2-1  RS-N2P2-2 
#110000     110000     110000     110000     110000     110000     110000     110000     110000     110000     110000 
#RS-N2P2-3  RS-N2P2-4 
#110000     110000 
#write.csv(rarefy_kegg,file="J://shangzhuang/kegg/kegg_rarefy.csv")


#Alpha diversity, PCoA based on CAZyme genes in bluk soil and rhizosphere
cazyme<-read.csv("J://shangzhuang/cazyme/GH_table_rarefied.csv",row.names = 2)
design<-read.csv("J://shangzhuang/design.csv")
#Fig.1a alpha diversity
library(vegan)
#alpha diversity
index<-diversity(t(cazyme[,2:25]), index= "shannon", MARGIN = 1)
index<-as.data.frame(index)
index<-cbind(design,index)
library(ggplot2)
p=ggplot(index,aes(x=Group, y=index,fill=Compartment))+
  geom_violin()+geom_boxplot(width=0.15,fill="white")+
  labs(x="", y="Gene diversity")+theme_bw()+
  scale_fill_manual(values=c("#9B8261","#5C8E78"))+
  #  scale_fill_brewer()+
  # scale_y_continuous(limit = c(2.15, 2.35))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=9))

#Beta diversity of CAZyme genes
library(ape)
data <- vegdist(t(cazyme[,2:25]), method = "bray")
#PCoA analysis
pcoa<- pcoa(data, correction = "none", rn = NULL)
#Generate PCoA table with the first and second axis
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
aa<- data.frame(rownames(pcoa$vectors),PC1,PC2)
colnames(aa) <-c("sample","PC1","PC2")
points<-cbind(design,aa)
#PCoA tables visualization
p=ggplot(data=points,aes(x=PC1,y=PC2))+geom_point(aes(colour=Compartment,shape=Treatment),size=3.5)+
  scale_colour_manual(values=c("#9B8261","#5C8E78"))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  xlab(paste("PCoA ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PCoA ( ",pc2,"%"," )",sep=""))+ 
  stat_ellipse(aes(fill=Compartment),data=points,geom ="polygon",alpha=0.12)+
  scale_fill_manual(values =c("#9B8261","#5C8E78"))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour='black',size=9)) 

#Beta diversity of motility genes
#PCoA of motility genes
library(vegan)
library(ape)
#beta diversity
motility<-read.csv("J://shangzhuang/kegg/motility_gene.csv",row.names = 6)
design<-read.csv("J://shangzhuang/design.csv")
data <- vegdist(t(motility[,7:30]), method = "bray")
#PCoA analysis
pcoa<- pcoa(data, correction = "none", rn = NULL)
#Generate PCoA table with the first and second axis
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
aa<- data.frame(rownames(pcoa$vectors),PC1,PC2)
colnames(aa) <-c("sample","PC1","PC2")
points<-cbind(design,aa)
#PCoA tables visualization
p=ggplot(data=points,aes(x=PC1,y=PC2))+geom_point(aes(colour=Compartment,shape=Treatment),size=3.5)+
  scale_colour_manual(values=c("#5C8E78","#9B8261"))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+ 
  stat_ellipse(aes(fill=Compartment),data=points,geom ="polygon",alpha=0.12)+
  scale_fill_manual(values =c("#5C8E78","#9B8261"))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour='black',size=9))

#Taxonomic composition of CAZyme and motility taxa
library(car)
cazyme<-read.csv("J://shangzhuang/cazyme/cazyme_taxonomy.csv")
design<-read.csv("J://shangzhuang/design.csv")
se_phy<-aggregate(cazyme[,8:31],by=list(cazyme$phylum),FUN=sum)
rownames(se_phy)<-se_phy[,1]
mean<-se_phy[,-1]
mean<-mean[order(rowSums(mean),decreasing=T),]
se<-rbind(colSums(mean[14:39,]),mean[13:1,])
rownames(se)[1]<-"Others"
aa<-cbind(t(mean),design)
se_phy<-aggregate(aa[,1:39],by=list(aa$Compartment),FUN=mean)
rownames(se_phy)<-se_phy[,1]
mean<-se_phy[,-1]
#write.csv(t(mean),file="J://shangzhuang/motility_taxa/CAZyme_phyla.csv")
#motility
motility<-read.csv("J://shangzhuang/motility_taxa/motility_taxonomy.csv")
design<-read.csv("J://shangzhuang/design.csv")
se_phy<-aggregate(motility[,9:32],by=list(motility$Level_2),FUN=sum)
rownames(se_phy)<-se_phy[,1]
mean<-se_phy[,-1]
mean<-mean[order(rowSums(mean),decreasing=T),]
se<-rbind(colSums(mean[14:26,]),mean[13:1,])
rownames(se)[1]<-"Others"
aa<-cbind(t(mean),design)
se_phy<-aggregate(aa[,1:26],by=list(aa$Compartment),FUN=mean)
rownames(se_phy)<-se_phy[,1]
mean<-se_phy[,-1]
#write.csv(t(mean),file="J://shangzhuang/motility_taxa/motility_phyla.csv")

#Figure generation for CAZyme taxa and motility taxa
library(ggplot2)
library(reshape2)
library(ggalluvial)
library(vegan)
aa<-read.csv("J://shangzhuang/motility_taxa/barplot_tax.csv",row.names = 1)
tax<-cbind(rownames(aa),aa[,1:2])
dat <- reshape2::melt(tax, id = 'rownames(aa)')
Group<-as.data.frame(c(rep("CAZyme taxonomy",18)))
colnames(Group)<-"Group"
dat<-cbind(dat,Group)
#color
color <-c( "#FBD7C7","#B07873","#D8BFD8","#F5DEB3","#C67171","#5E98A1","#8FB28D","#7381B0","#7D9EC0")
#figure
p=ggplot(dat, aes(x =dat$variable, y = value, 
                fill =dat$`rownames(aa)`,stratum =dat$`rownames(aa)`, alluvium =dat$`rownames(aa)`)) +
  geom_stratum(width = 0.6,alpha=1) +geom_flow(alpha = 0.3) + 
  scale_fill_manual(values=color)+
  theme(axis.text=element_text(colour='black',size=9))+
  labs(x = '', y = 'Relative abundance',fill=" ")+theme_bw()+
  theme(axis.text=element_text(colour='black',size=11))
#motility taxa
aa<-read.csv("J://shangzhuang/motility_taxa/barplot_tax.csv",row.names = 1)
tax<-cbind(rownames(aa),aa[,4:5])
dat <- reshape2::melt(tax, id = 'rownames(aa)')
Group<-as.data.frame(c(rep("Motility taxonomy",18)))
colnames(Group)<-"Group"
dat<-cbind(dat,Group)
#color
color <-c( "#FBD7C7","#B07873","#D8BFD8","#F5DEB3","#C67171","#5E98A1","#8FB28D","#7381B0","#7D9EC0")
p=ggplot(dat, aes(x =dat$variable, y = value, 
                fill =dat$`rownames(aa)`,stratum =dat$`rownames(aa)`, alluvium =dat$`rownames(aa)`)) +
  geom_stratum(width = 0.6,alpha=1) +geom_flow(alpha = 0.3) + 
  scale_fill_manual(values=color)+
  theme(axis.text=element_text(colour='black',size=9))+
  labs(x = '', y = 'Relative abundance',fill=" ")+theme_bw()+
  theme(axis.text=element_text(colour='black',size=11))




