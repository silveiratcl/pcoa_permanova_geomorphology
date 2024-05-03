library(readxl)
library(ape)
library(vegan)
library(FD)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(ggforce)
library(RColorBrewer)

maindf <- read_excel("~/Documents/IMBRSea/Thesis Brazil/Data Anlaysis/R/3dvsdafor_version2.xlsx")
maintranf<-maindf[,-c(1,2)] #remove the columns with sites and observers
site<-maindf[,-(2:7)] #extracting the sites column
observers<-maindf[,-c(1,3:7)] #extracting the observers column
gowerdf<-gowdis(maintranf) #calculating the gower dissimilarity since the data is categorical/ordinal
gowerpcoa<-cmdscale(gowerdf,k=(nrow(maindf)-1),eig=TRUE,add=FALSE) #creating the dissimilarity matrix using cmdscale funtion with 2 dimensions
head(gowerpcoa$points)[,1:2]

gowerpcoa$eig[1:5]
(explic <- gowerpcoa$eig[1:2]/sum(gowerpcoa$eig)*100) #calculating the variation explained by both the axis

pcoadf<-as.data.frame(gowerpcoa$points[,1:2]) #creating a dataframe to input the gower dissimilarity into columns
colnames(pcoadf)<-c("PC1","PC2") #column names changed to PC1 and PC2
head(pcoadf)

merged_df<-cbind(site,observers,pcoadf[,-c(3:4)]) #merging the columns from the pcoa analysis and the observers and sites column
merged_df
attach(merged_df)

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

# Define the number of sites
num_sites <- 10
# Choose a color palette with 10 distinct colors
site_colors1 <- brewer.pal(n = num_sites, name = "Set1")
site_colors_dark2 <- brewer.pal(n = num_sites, name = "Dark2")
custom_site_colors <- c("red","skyblue", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E","yellow")
custom_site_colors_paired <- brewer.pal(n = num_sites, name = "Paired")

plotgg2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("PCo1(48.1%)") +
  ylab("PCo2(42.7%)") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_mark_ellipse(data=merged_df, 
                    aes(x=PC1, y=PC2, colour=Site, 
                        fill=after_scale(alpha(colour, 0.2))), 
                    expand=0, size=0.2, show.legend=FALSE) +
  
  geom_point(data=merged_df, 
             aes(x=PC1, y=PC2, colour=Site, shape=observers), 
             size=5) +
  BioR.theme +
  scale_colour_manual(values = custom_site_colors) +
  coord_fixed(ratio=1) +
  ggtitle("PCoA")+
  theme_minimal()
plotgg2


###############################PERMANOVA#######################################

comm <- maindf[,3:7]
beta.comm <- betadisper(gowerdf, group=maindf$observers)
permutest(beta.comm)
permanova1 <- adonis2(comm ~ observers, data=maindf, method="gower", permutations=9999) 
View(permanova1)

#restricted permanova blocking each of the site to test the difference between observers
comm <- maindf[,3:7]
beta.comm <- betadisper(gowerdf, group=maindf$observers)
permutest(beta.comm)
perm.ctrl <- how(blocks=maindf$Site, nperm=9999)
permanova2 <- adonis2(comm ~ observers, data=maindf, method="gower", permutations= perm.ctrl) 
View(permanova2)

####permanova version2#######
attach(pv)
comm <- pv[,3:7]
beta.comm <- betadisper(gowerdf, group=pv$observers)
permutest(beta.comm)
permanova1 <- adonis2(comm ~ observers, data=pv, method="gower", permutations=9999) 
View(permanova1)

perm.ctrl <- how(blocks=pv$Site, nperm=9999)
permanova2 <- adonis2(comm ~ observers, data=pv, method="gower", permutations= perm.ctrl)
View(permanova2)
