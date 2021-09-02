#############################################################################################################
### analysis for the paper, Consistent patterns of microbial diversity, prediction accuracy and associated
#   bacterial organisms between pediatric ulcerative colitis and healthy children using 16S rRNA and 
#   metagenomic shotgun sequencing data
### remove uc with therapies
### including the rarefaction, alpha, beta diversity analysis, and differential analysis
### 2021/08/31
#############################################################################################################



#=================================================================================================#
# packages needed in the analysis
library(vegan)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ape)
library(edgeR)
library(ranacapa)



#=================================================================================================#
# metadata
metadata <- read.csv("metadata.csv",row.names=1)

# filter samples, remove UC with therapies 
metadata <- metadata[metadata$group=="Healthy"|metadata$group=="UC",]



#=================================================================================================#
###### function for reading the amplicon count table at different levels ######
# filename="E:/USC/IBD/20210614_analysis/data_from_server/table_l2.qza";samples=samples_keep$SampleID
read.feature.table <- function(filename,samples){
  df <- as.data.frame(read_qza(filename)$data)
  colnames(df) <- gsub("^s","",colnames(df))
  # filter samples
  df <- df[,which(colnames(df)%in%samples)]
  return(df)
}
###### feature table ######
amplicon_phylum_table  <- read.feature.table("qiime2/table-l2.qza",samples=metadata$SampleID)
amplicon_class_table   <- read.feature.table("qiime2/table-l3.qza",samples=metadata$SampleID)
amplicon_order_table   <- read.feature.table("qiime2/table-l4.qza",samples=metadata$SampleID)
amplicon_family_table  <- read.feature.table("qiime2/table-l5.qza",samples=metadata$SampleID)
amplicon_genus_table   <- read.feature.table("qiime2/table-l6.qza",samples=metadata$SampleID)
amplicon_species_table <- read.feature.table("qiime2/table-l7.qza",samples=metadata$SampleID)




#=================================================================================================#
###### function for reading the shotgun count table at different levels ######
# filename="E:/USC/IBD/20210513_revision/shotgun_taxonomy/shotgun_abundance_L7_m.txt";samples=samples_keep$shotgun.SampleID
read.shotgun.table <- function(filename,samples){
  # reads number for each sample
  count_table <- read.csv("E:/USC/IBD/20210607_sequence_depth/count_table_shotgun.csv")
  count_table$X <- sub("-",".",count_table$X)
  rownames(count_table) <- count_table$X
  count_table <- count_table[,-1]
  count_for_samples <- rowSums(count_table)
  # relative abundance
  df <- read.table(filename,sep="\t",quote = "")
  for(i in 1:ncol(df)){df[,i] <- as.integer(df[,i]*count_for_samples[colnames(df)[i]])}
  # filter samples
  df <- df[,which(colnames(df)%in%samples)]
  # return
  return(df)
}
###### feature table ######
shotgun_phylum_table  <- read.shotgun.table("shotgun_tables/shotgun_abundance_L2.txt",samples=metadata$shotgun.SampleID)
shotgun_class_table   <- read.shotgun.table("shotgun_tables/shotgun_abundance_L3.txt",samples=metadata$shotgun.SampleID)
shotgun_order_table   <- read.shotgun.table("shotgun_tables/shotgun_abundance_L4.txt",samples=metadata$shotgun.SampleID)
shotgun_family_table  <- read.shotgun.table("shotgun_tables/shotgun_abundance_L5.txt",samples=metadata$shotgun.SampleID)
shotgun_genus_table   <- read.shotgun.table("shotgun_tables/shotgun_abundance_L6.txt",samples=metadata$shotgun.SampleID)
shotgun_species_table <- read.shotgun.table("shotgun_tables/shotgun_abundance_L7.txt",samples=metadata$shotgun.SampleID)



#=================================================================================================#
###### shotgun subspecies count table ######
shotgun_subspecies_table <- read.csv("count_table_shotgun.csv",row.names=1)
rownames(shotgun_subspecies_table) <- sub("-",".",rownames(shotgun_subspecies_table))
shotgun_subspecies_table <- shotgun_subspecies_table[which(rownames(shotgun_subspecies_table)%in%metadata$shotgun.SampleID),]
shotgun_subspecies_table <- as.data.frame(t(shotgun_subspecies_table))



#=================================================================================================#
###### rarefaction curve for shotgun and 16s ######
### shotgun subspecies ###
rownames(metadata) <- metadata$shotgun.SampleID
SampleType <- metadata[colnames(shotgun_subspecies_table),c("shotgun.SampleID","SampleType")]
shotgun_table <- shotgun_subspecies_table
shotgun_table$sum.taxonomy <- rownames(shotgun_subspecies_table)
shotgun_subspecies_physeq_ob <- convert_anacapa_to_phyloseq(shotgun_table,SampleType)
shotgun_subspecies_rarefaction <- ggrare(shotgun_subspecies_physeq_ob, step = 100, se = F,color="SampleType")
shotgun_subspecies_rarefaction <- shotgun_subspecies_rarefaction+
  labs(x="Number of sequences", y="Number of sepcies", title="Rarafaction curve for shotgun data")+
  scale_x_continuous(breaks = c(0,1000000,2000000),labels=c(0,1000000,2000000))+
  scale_color_discrete(name="")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(legend.justification=c(0.98,0.02), legend.position=c(0.98,0.02),legend.background = element_rect(fill = NA))+
  theme(plot.caption = element_text(hjust = 0, face="italic"))

### amplicon genus ###
rownames(metadata) <- metadata$SampleID
SampleType <- metadata[colnames(amplicon_genus_table),c("SampleID","SampleType")]
amplicon_table <- amplicon_genus_table
# remove sample "PT.58" because of the large sequencing depth
SampleType <- SampleType[-grep("PT.58",SampleType$SampleID),]
amplicon_table <- amplicon_table[,-grep("PT.58",colnames(amplicon_table))]
# make names for the 16S samples
SampleType$SampleID <- paste0("x",SampleType$SampleID)
colnames(amplicon_table) <- paste0("x",colnames(amplicon_table))
# convert to physeq object
amplicon_table$sum.taxonomy <- rownames(amplicon_genus_table)
amplicon_genus_physeq_ob <- convert_anacapa_to_phyloseq(amplicon_table,SampleType)
amplicon_genus_rarefaction <- ggrare(amplicon_genus_physeq_ob, step=100, se=F, color="SampleType")
amplicon_genus_rarefaction <- amplicon_genus_rarefaction+
  labs(x="Number of sequences", y="Number of genus", title="Rarafaction curve for 16S data")+
  scale_color_discrete(name="")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(legend.justification=c(0.98,0.02), legend.position=c(0.98,0.02),legend.background = element_rect(fill = NA))+
  theme(plot.caption = element_text(hjust = 0, face="italic"))

### arrange the plots ###
ggarrange(amplicon_genus_rarefaction,shotgun_subspecies_rarefaction,nrow=1,ncol=2,labels=c("A","B"))



#=================================================================================================#
##### box plot for Shannon diversities at different rarefaction level ######
plot.shannon.box <- function(filename,metadata,plot_title){
  
  # feature table
  feature_table <- read.csv(filename)
  rownames(feature_table) <- feature_table$X
  feature_table <- feature_table[,-1]
  colnames(feature_table) <- sub("^X","",colnames(feature_table))
  
  # compute shannon index
  alpha_diversities <- diversity(t(feature_table),index="shannon")
  df <- data.frame(SampleID=names(alpha_diversities),diversities=alpha_diversities)
  df <- merge(df,metadata,by="SampleID")
  
  # boxplot
  boxplot <- ggplot(df, aes(x=SampleType,y=diversities)) +
    geom_boxplot(aes(colour=SampleType),width=0.5,size=0.5,outlier.fill="white",outlier.color="white")+ 
    geom_jitter(aes(colour=SampleType,fill=SampleType),width =0.2,shape = 21,size=2)+
    scale_y_continuous(name = "Shannon Index",limits = c(0,max(df$diversities)+0.5))+
    scale_x_discrete(name = "Phenotypes")+ 
    labs(title=plot_title)+
    stat_compare_means(comparisons=list(c("Healthy","UC")),correct=FALSE,label="p.format",method = "wilcox.test")+
    theme_bw()+
    theme(panel.grid.minor = element_blank())+
    theme(legend.position="none")
  
  # return
  return(boxplot)
}

# rarefaction level (6): 1k, 5k, 10k, 30k, 50k, 100k
amplicon_genus_shannon_1k <- plot.shannon.box(filename="amplicon_tables/genus_table_rarefaction_1k.csv",metadata=metadata,plot_title="1k reads per sample")
amplicon_genus_shannon_5k <- plot.shannon.box(filename="amplicon_tables/genus_table_rarefaction_5k.csv",metadata=metadata,plot_title="5k reads per sample")
amplicon_genus_shannon_10k <- plot.shannon.box(filename="amplicon_tables/genus_table_rarefaction_10k.csv",metadata=metadata,plot_title="10k reads per sample")
amplicon_genus_shannon_30k <- plot.shannon.box(filename="amplicon_tables/genus_table_rarefaction_30k.csv",metadata=metadata,plot_title="30k reads per sample")
amplicon_genus_shannon_50k <- plot.shannon.box(filename="amplicon_tables/genus_table_rarefaction_50k.csv",metadata=metadata,plot_title="50k reads per sample")
amplicon_genus_shannon_100k <- plot.shannon.box(filename="amplicon_tables/genus_table_rarefaction_100k.csv",metadata=metadata,plot_title="100k reads per sample")

# arrange the plot
ggarrange(amplicon_genus_shannon_1k,amplicon_genus_shannon_5k,amplicon_genus_shannon_10k,
          amplicon_genus_shannon_30k,amplicon_genus_shannon_50k,amplicon_genus_shannon_100k,
          ncol=3,nrow=2,labels=c("A","B","C","D","E","F"))



#=================================================================================================#
##### PCoA plots at different rarefaction level ######
plot.pcoa <- function(filename,metadata,plot_title){
  
  # feature table
  feature_table <- read.csv(filename,row.names=1)
  colnames(feature_table) <- sub("^X","",colnames(feature_table))
  
  # filter the samples
  feature_table <- feature_table[,colnames(feature_table)%in%metadata$SampleID]
  
  # normalize the feature table
  feature_table <- as.data.frame(apply(feature_table,2,function(x) x/sum(x)))
  
  # compute the bray curtis distance
  beta_diversity <- as.matrix(vegdist(t(feature_table),method = "bray"))
  
  # permanova
  metadata <- metadata[rownames(beta_diversity),]
  permanova <- adonis(beta_diversity~SampleType+Gender+Age, data=metadata, permutations=1000)
  r2 <- permanova$aov.tab["SampleType","R2"]
  p.value <- permanova$aov.tab["SampleType","Pr(>F)"]
  
  # annotate the r2 and p value in the figure
  r2 <- sprintf("italic(R^2) == %.3f",r2)
  p.value <- sprintf("italic(p) == %.3f",p.value) 
  permanova_labels <- data.frame(r2=r2,p.value=p.value,stringsAsFactors = FALSE)
  
  # pcoa plot
  PCOA <- pcoa(as.dist(beta_diversity))
  # data frame for pcoa plot
  pcoa_df <- as.data.frame(PCOA$vectors[,1:2])
  pcoa_df$SampleID <- rownames(pcoa_df)
  pcoa_df <- merge(pcoa_df,metadata,by="SampleID")
  # axis
  pro1 <- as.numeric(sprintf("%.3f",PCOA$values[,"Relative_eig"][1]))*100  
  pro2 <- as.numeric(sprintf("%.3f",PCOA$values[,"Relative_eig"][2]))*100
  # plot
  pcoa_plot <- ggplot(pcoa_df,aes(x=-Axis.1,y=-Axis.2,col=group)) +
    geom_point(size=3) +
    xlim(-0.7,0.6) + ylim(-0.6,0.4) +
    labs(x=paste0("PCOA1(",pro1,"%)"), y=paste0("PCOA2(",pro2,"%)"),title=plot_title) +
    geom_vline(aes(xintercept=0),linetype="dotted") +
    geom_hline(aes(yintercept=0),linetype="dotted") +
    stat_ellipse(aes(group=SampleType,col=SampleType), level = 0.8, show.legend=FALSE)+
    annotate("text",label="Phenotypes:",x=-0.65,y=-0.4,size=4,hjust=0)+
    geom_text(data=permanova_labels,mapping=aes(x=-0.65,y=-0.47,label=r2),parse=TRUE,inherit.aes=FALSE,size=4,hjust=0)+
    geom_text(data=permanova_labels,mapping=aes(x=-0.65,y=-0.54,label=p.value),parse=TRUE,inherit.aes=FALSE,size=4,hjust=0)+
    theme_bw() +
    theme(legend.title=element_blank(), legend.text=element_text(size=12))+
    theme(title = element_text(size = 14))+
    theme(axis.title = element_text(size = 16),axis.text = element_text(size = 12,colour="black"))+
    theme(legend.justification=c(0.02,0.98), legend.position=c(0.02,0.98),legend.background = element_rect(fill = NA))
  
  # return
  return(pcoa_plot)
}

rownames(metadata) <- metadata$SampleID

# rarefaction level (6): 1k, 5k, 10k, 30k, 50k, 100k
amplicon_genus_pcoa_1k <- plot.pcoa(filename="amplicon_tables/genus_table_rarefaction_1k.csv",metadata=metadata,plot_title="1k reads per sample")
amplicon_genus_pcoa_5k <- plot.pcoa(filename="amplicon_tables/genus_table_rarefaction_5k.csv",metadata=metadata,plot_title="5k reads per sample")
amplicon_genus_pcoa_10k <- plot.pcoa(filename="amplicon_tables/genus_table_rarefaction_10k.csv",metadata=metadata,plot_title="10k reads per sample")
amplicon_genus_pcoa_30k <- plot.pcoa(filename="amplicon_tables/genus_table_rarefaction_30k.csv",metadata=metadata,plot_title="30k reads per sample")
amplicon_genus_pcoa_50k <- plot.pcoa(filename="amplicon_tables/genus_table_rarefaction_50k.csv",metadata=metadata,plot_title="50k reads per sample")
amplicon_genus_pcoa_100k <- plot.pcoa(filename="amplicon_tables/genus_table_rarefaction_100k.csv",metadata=metadata,plot_title="100k reads per sample")

# arrange the plot
ggarrange(amplicon_genus_pcoa_1k,amplicon_genus_pcoa_5k,amplicon_genus_pcoa_10k,
          amplicon_genus_pcoa_30k,amplicon_genus_pcoa_50k,amplicon_genus_pcoa_100k,
          ncol=3,nrow=2,labels=c("A","B","C","D","E","F"))



#=================================================================================================#
###### function for ploting the pcoa plot and permanova ######
plot.pcoa <- function(feature_table,metadata,plot_title,samples){
  
  # normalize the feature table
  feature_table <- as.data.frame(apply(feature_table,2,function(x) x/sum(x)))
  
  # compute the bray curtis distance
  beta_diversity <- as.matrix(vegdist(t(feature_table),method = "bray"))
  
  # permanova
  rownames(metadata) <- metadata[,samples]
  metadata <- metadata[rownames(beta_diversity),]
  permanova <- adonis(beta_diversity~SampleType+Gender+Age, data=metadata, permutations=1000)
  r2 <- permanova$aov.tab["SampleType","R2"]
  p.value <- permanova$aov.tab["SampleType","Pr(>F)"]
  
  # annotate the r2 and p value in the figure
  r2 <- sprintf("italic(R^2) == %.3f",r2)
  p.value <- sprintf("italic(p) == %.3f",p.value) 
  permanova_labels <- data.frame(r2=r2,p.value=p.value,stringsAsFactors = FALSE)
  
  # pcoa plot
  PCOA <- pcoa(as.dist(beta_diversity))
  # data frame for pcoa plot
  pcoa_df <- as.data.frame(PCOA$vectors[,1:2])
  pcoa_df[,samples] <- rownames(pcoa_df)
  pcoa_df <- merge(pcoa_df,metadata,by=samples)
  # axis
  pro1 <- as.numeric(sprintf("%.3f",PCOA$values[,"Relative_eig"][1]))*100  
  pro2 <- as.numeric(sprintf("%.3f",PCOA$values[,"Relative_eig"][2]))*100
  # plot
  pcoa_plot <- ggplot(pcoa_df,aes(x=-Axis.1,y=Axis.2,col=group)) +
    geom_point(size=3) +
    xlim(-0.7,0.6) + ylim(-0.65,0.4) +
    labs(x=paste0("PCOA1(",pro1,"%)"), y=paste0("PCOA2(",pro2,"%)"),title=plot_title) +
    geom_vline(aes(xintercept=0),linetype="dotted") +
    geom_hline(aes(yintercept=0),linetype="dotted") +
    stat_ellipse(aes(group=SampleType,col=SampleType), level = 0.8, show.legend=FALSE)+
    annotate("text",label="Phenotypes:",x=-0.65,y=-0.4,size=4,hjust=0)+
    geom_text(data=permanova_labels,mapping=aes(x=-0.65,y=-0.47,label=r2),parse=TRUE,inherit.aes=FALSE,size=4,hjust=0)+
    geom_text(data=permanova_labels,mapping=aes(x=-0.65,y=-0.54,label=p.value),parse=TRUE,inherit.aes=FALSE,size=4,hjust=0)+
    theme_bw() +
    theme(legend.title=element_blank())+
    theme(title = element_text(size = 14))+
    theme(axis.title = element_text(size = 16),axis.text = element_text(size = 12,colour="black"))+
    theme(legend.justification=c(0.02,0.98), legend.position=c(0.02,0.98),legend.background = element_rect(fill = NA))
  
  # return
  return(pcoa_plot)
}

# shotgun
shotgun_subspecies_pcoa <- plot.pcoa(feature_table=shotgun_subspecies_table, metadata=metadata,
                                     plot_title="PCoA plot using shotgun data",
                                     samples="shotgun.SampleID")
# amplicon
amplicon_genus_pcoa <- plot.pcoa(feature_table=amplicon_genus_table, metadata=metadata,
                                 plot_title="PCoA plot using 16S data",
                                 samples="SampleID")

# arrange the plots
ggarrange(amplicon_genus_pcoa,shotgun_subspecies_pcoa,nrow=1,ncol=2,labels=c("A","B"))



#=================================================================================================#
### correlation for shannon index between shotgun and 16s ###
SampleType <- metadata[,c("SampleID","shotgun.SampleID","SampleType")]

# shotgun shannon
shotgun_subspecies_shannon <- diversity(t(shotgun_subspecies_table),index='shannon')
shotgun_subspecies_shannon_df <- data.frame(shotgun.SampleID=names(shotgun_subspecies_shannon),shotgun.shannon=shotgun_subspecies_shannon)
shannon_corr_df <- merge(shotgun_subspecies_shannon_df,SampleType,by="shotgun.SampleID")

# amplicon genus shannon
amplicon_genus_shannon <- diversity(t(amplicon_genus_table),index='shannon')
amplicon_genus_shannon_df <- data.frame(SampleID=names(amplicon_genus_shannon),amplicon.shannon=amplicon_genus_shannon)
shannon_corr_df <- merge(shannon_corr_df,amplicon_genus_shannon_df,by="SampleID")

# linear regression analysis
healthy_lm <- lm(amplicon.shannon~shotgun.shannon,data=shannon_corr_df[shannon_corr_df$SampleType=="Healthy",])
uc_lm <- lm(amplicon.shannon~shotgun.shannon,data=shannon_corr_df[shannon_corr_df$SampleType=="UC",])
overall_lm <- lm(amplicon.shannon~shotgun.shannon,data=shannon_corr_df)

# spearman correlation test
spearman_test <- cor.test(shannon_corr_df$shotgun.shannon,shannon_corr_df$amplicon.shannon,method="spearman")
corr_name <- "Spearman's correlation test:"
corr_coef <- sprintf("italic(rho) == %.3f",spearman_test$estimate)
corr_pvalue <- paste0("p = ", format(spearman_test$p.value,scientific=TRUE,digit=3))
corr_label <- data.frame(title=corr_name,rho=corr_coef,p=corr_pvalue)

# scatter plot
ggplot(shannon_corr_df,aes(y=amplicon.shannon,x=shotgun.shannon,col=SampleType))+
  geom_point(size=3)+
  scale_y_continuous(name="16S genus Shannon Indices")+
  scale_x_continuous(name="shotgun species Shannon Indicies")+
  geom_abline(intercept=coef(healthy_lm)[1],slope = coef(healthy_lm)[2],color="#F8766D",lwd=1)+
  geom_abline(intercept=coef(uc_lm)[1],slope = coef(uc_lm)[2],color="#00BFC4",lwd=1)+ 
  geom_abline(intercept=coef(overall_lm)[1],slope = coef(overall_lm)[2],color="black",lwd=1)+
  geom_text(data=corr_label,mapping=aes(x=1,y=3.2,label=title),inherit.aes=FALSE,size=4,color="black",hjust=0)+ 
  geom_text(data=corr_label,mapping=aes(x=1,y=3.1,label=rho),parse=TRUE,inherit.aes=FALSE,size=4,color="black",hjust=0)+ 
  geom_text(data=corr_label,mapping=aes(x=1,y=3.0,label=p),inherit.aes=FALSE,size=4,color="black",hjust=0)+
  theme_bw()+
  theme(legend.title=element_blank()) + 
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 14,colour="black"))+
  theme(panel.grid.minor = element_blank())+
  theme(legend.justification=c(0.98,0.02),legend.position=c(0.98,0.02),legend.background=element_rect(fill = NA))



#=================================================================================================#
### correlation for shannon index between shotgun and 16s ###
SampleType <- metadata[,c("SampleID","shotgun.SampleID","SampleType")]

# shotgun subspecies pcoa1
shotgun_subspecies_ab <- as.data.frame(apply(shotgun_subspecies_table,2,function(x) x/sum(x)))
shotgun_subspecies_beta_div <- as.matrix(vegdist(t(shotgun_subspecies_ab),method = "bray"))
shotgun_subspecies_pcoa <- pcoa(as.dist(shotgun_subspecies_beta_div))
shotgun_subspecies_pcoa_df <- data.frame(shotgun.Axis1=-shotgun_subspecies_pcoa$vectors[,1])
shotgun_subspecies_pcoa_df$shotgun.SampleID <- rownames(shotgun_subspecies_pcoa_df)
pcoa_corr_df <- merge(shotgun_subspecies_pcoa_df,SampleType,by="shotgun.SampleID")

# amplicon genus pcoa1
amplicon_genus_ab <- as.data.frame(apply(amplicon_genus_table,2,function(x) x/sum(x)))
amplicon_genus_beta_div <- as.matrix(vegdist(t(amplicon_genus_table),method = "bray"))
amplicon_genus_pcoa <- pcoa(as.dist(amplicon_genus_beta_div))
amplicon_genus_pcoa_df <- data.frame(amplicon.Axis1=-amplicon_genus_pcoa$vectors[,1])
amplicon_genus_pcoa_df$SampleID <- rownames(amplicon_genus_pcoa_df)
pcoa_corr_df <- merge(amplicon_genus_pcoa_df,pcoa_corr_df,by="SampleID")

# linear regression analysis
healthy_lm <- lm(amplicon.Axis1~shotgun.Axis1,data=pcoa_corr_df[pcoa_corr_df$SampleType=="Healthy",])
uc_lm <- lm(amplicon.Axis1~shotgun.Axis1,data=pcoa_corr_df[pcoa_corr_df$SampleType=="UC",])
overall_lm <- lm(amplicon.Axis1~shotgun.Axis1,data=pcoa_corr_df)

# spearman correlation test
cor.test(pcoa_corr_df[pcoa_corr_df$SampleType=="Healthy","shotgun.Axis1"],
         pcoa_corr_df[pcoa_corr_df$SampleType=="Healthy","amplicon.Axis1"],method="spearman")
cor.test(pcoa_corr_df[pcoa_corr_df$SampleType=="UC","shotgun.Axis1"],
         pcoa_corr_df[pcoa_corr_df$SampleType=="UC","amplicon.Axis1"],method="spearman")
spearman_test <- cor.test(pcoa_corr_df$shotgun.Axis1,pcoa_corr_df$amplicon.Axis1,method="spearman")
corr_name <- "Spearman's correlation test:"
corr_coef <- sprintf("italic(rho) == %.3f",spearman_test$estimate)
corr_pvalue <- "p < 2.2e-16"
corr_label <- data.frame(title=corr_name,rho=corr_coef,p=corr_pvalue)

# scatter plot
ggplot(pcoa_corr_df,aes(y=amplicon.Axis1,x=shotgun.Axis1,col=SampleType))+
  geom_point(size=3)+
  scale_y_continuous(name="16S genus PCoA1 (30.2%)")+
  scale_x_continuous(name="shotgun species PCoA1 (24.6)")+
  geom_abline(intercept=coef(healthy_lm)[1],slope = coef(healthy_lm)[2],color="#F8766D",lwd=1)+
  geom_abline(intercept=coef(uc_lm)[1],slope = coef(uc_lm)[2],color="#00BFC4",lwd=1)+ 
  geom_abline(intercept=coef(overall_lm)[1],slope = coef(overall_lm)[2],color="black",lwd=1)+
  geom_text(data=corr_label,mapping=aes(x=-0.4,y=0.5,label=title),inherit.aes=FALSE,size=4,color="black",hjust=0)+ 
  geom_text(data=corr_label,mapping=aes(x=-0.4,y=0.45,label=rho),parse=TRUE,inherit.aes=FALSE,size=4,color="black",hjust=0)+ 
  geom_text(data=corr_label,mapping=aes(x=-0.4,y=0.4,label=p),inherit.aes=FALSE,size=4,color="black",hjust=0)+
  theme_bw()+
  theme(legend.title=element_blank()) + 
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 14,colour="black"))+
  theme(panel.grid.minor = element_blank())+
  theme(legend.justification=c(0.98,0.02),legend.position=c(0.98,0.02),legend.background=element_rect(fill = NA))



#=================================================================================================#
###### function for differential analysis ######
#count_table=amplicon_species_table; metadata=metadata
differential.analysis <- function(count_table, metadata){
  
  # filter the features with >90% zeros (reduce the effect of zero inflation)
  count_table <- count_table[!(rowSums(count_table==0) > ncol(count_table)*0.9),]
  # remove features with very low variance (below half the median of all feature wise variances)
  feature_variances <- apply(count_table,1,var)
  count_table <- count_table[!(feature_variances < median(feature_variances)*0.5),]
  
  # group information 
  group <- metadata[colnames(count_table),"SampleType"]
  
  # TMM normalization from edgeR package
  dge <- DGEList(counts=count_table,group=group)
  tmm_dge <- edgeR::calcNormFactors(dge, method = "TMM")
  
  # tagwise dispersion estimation
  design <- model.matrix(~group)
  tmm_dge <- estimateDisp(tmm_dge,design)
  
  # differential analysis
  et <- exactTest(tmm_dge)
  tTag <- topTags(et, n=nrow(count_table))
  tTag <- as.data.frame(tTag)
  
  # return
  return(tTag)
}


### amplicon results
rownames(metadata) <- metadata$SampleID
amplicon_phylum_diff <- differential.analysis(count_table=amplicon_phylum_table, metadata=metadata)
amplicon_class_diff <- differential.analysis(count_table=amplicon_class_table, metadata=metadata)
amplicon_order_diff <- differential.analysis(count_table=amplicon_order_table, metadata=metadata)
amplicon_family_diff <- differential.analysis(count_table=amplicon_family_table, metadata=metadata)
amplicon_genus_diff <- differential.analysis(count_table=amplicon_genus_table, metadata=metadata)
amplicon_species_diff <- differential.analysis(count_table=amplicon_species_table, metadata=metadata)
# merge the results
amplicon_diff_results <- rbind(amplicon_phylum_diff,amplicon_class_diff)
amplicon_diff_results <- rbind(amplicon_diff_results,amplicon_order_diff)
amplicon_diff_results <- rbind(amplicon_diff_results,amplicon_family_diff)
amplicon_diff_results <- rbind(amplicon_diff_results,amplicon_genus_diff)
amplicon_diff_results <- rbind(amplicon_diff_results,amplicon_species_diff)
amplicon_diff_results <- amplicon_diff_results[order(rownames(amplicon_diff_results),decreasing=F),]
# save the results
write.csv(amplicon_species_diff[amplicon_species_diff$FDR<0.25,],"amplicon_tables/rm_therapy_amplicon_species_differential.csv")


### shotgun results
rownames(metadata) <- metadata$shotgun.SampleID
shotgun_phylum_diff <- differential.analysis(count_table=shotgun_phylum_table, metadata=metadata)
shotgun_class_diff <- differential.analysis(count_table=shotgun_class_table, metadata=metadata)
shotgun_order_diff <- differential.analysis(count_table=shotgun_order_table, metadata=metadata)
shotgun_family_diff <- differential.analysis(count_table=shotgun_family_table, metadata=metadata)
shotgun_genus_diff <- differential.analysis(count_table=shotgun_genus_table, metadata=metadata)
shotgun_species_diff <- differential.analysis(count_table=shotgun_species_table, metadata=metadata)
# merge the results
shotgun_diff_results <- rbind(shotgun_phylum_diff,shotgun_class_diff)
shotgun_diff_results <- rbind(shotgun_diff_results,shotgun_order_diff)
shotgun_diff_results <- rbind(shotgun_diff_results,shotgun_family_diff)
shotgun_diff_results <- rbind(shotgun_diff_results,shotgun_genus_diff)
shotgun_diff_results <- rbind(shotgun_diff_results,shotgun_species_diff)
shotgun_diff_results <- shotgun_diff_results[order(rownames(shotgun_diff_results),decreasing=F),]
# save the results
write.csv(shotgun_species_diff[shotgun_species_diff$FDR<0.25,],"shotgun_tables/rm_therapy_shotgun_species_differential.csv")



################### compare shotgun and 16s with vila et. al. ####################
# read the comparison file
comparison_df <- read.csv("amplicon_tables/rm_therapy_vila_comparison.csv")
colnames(comparison_df)[1] <- "family"
comparison_df$number <- as.integer(comparison_df$number)
# plot
require(ggplot2)
ggplot(data = comparison_df, aes(x=family, y=number)) + 
  geom_bar(stat="identity",aes(fill=effect))+
  geom_text(aes(label=abs(number)), hjust = 0.5,size=3)+
  # geom_col(width = 0.5) + 
  coord_flip()+
  scale_y_continuous(name = "number of species in each family")+
  scale_x_discrete(name = "disease associated taxa at family level") +
  theme_bw()+
  # theme(panel.grid =element_blank()) + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  facet_wrap(~methods, ncol=3)





























