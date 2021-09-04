#############################################################################################################
### analysis for the paper, Consistent patterns of microbial diversity, prediction accuracy and associated
#   bacterial organisms between pediatric ulcerative colitis and healthy children using 16S rRNA and 
#   metagenomic shotgun sequencing data
### including the rarefaction, alpha, beta diversity analysis, and differential analysis
### 2021/08/31
#############################################################################################################


#===========================================================================================================#
###### metadata and count table used in the following analysis ######
# metadata
metadata <- read.csv("metadata.csv",row.names=1)

# genus level abundance of 16S rRNA
require(qiime2R)
amplicon_genus_table <- as.data.frame(read_qza("./qiime2/table-l6.qza")$data)
colnames(amplicon_genus_table) <- gsub("^s","",colnames(amplicon_genus_table))
amplicon_genus_table <- amplicon_genus_table[,which(colnames(amplicon_genus_table)%in%metadata$SampleID)]

# species level abundance of shotgun data
shotgun_subspecies_table <- as.data.frame(t(read.csv("count_table_shotgun.csv",row.names=1)))
colnames(shotgun_subspecies_table) <- sub("-",".",colnames(shotgun_subspecies_table))
shotgun_subspecies_table <- shotgun_subspecies_table[,which(colnames(shotgun_subspecies_table)%in%metadata$shotgun.SampleID)]



#===========================================================================================================#
###### rarefaction curve for shotgun and 16s ######
require(ranacapa)
### shotgun subspecies ###
rownames(metadata) <- metadata$shotgun.SampleID
SampleType <- metadata[colnames(shotgun_subspecies_table),c("shotgun.SampleID","SampleType")]
shotgun_table <- shotgun_subspecies_table
# convert to physeq object
shotgun_table$sum.taxonomy <- rownames(shotgun_subspecies_table)
shotgun_subspecies_physeq_ob <- convert_anacapa_to_phyloseq(shotgun_table,SampleType)
# rarefaction curve
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
# rarefaction curve
amplicon_genus_rarefaction <- ggrare(amplicon_genus_physeq_ob, step=100, se=F, color="SampleType")
amplicon_genus_rarefaction <- amplicon_genus_rarefaction+
  labs(x="Number of sequences", y="Number of genera", title="Rarafaction curve for 16S data")+
  scale_color_discrete(name="")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(legend.justification=c(0.98,0.02), legend.position=c(0.98,0.02),legend.background = element_rect(fill = NA))+
  theme(plot.caption = element_text(hjust = 0, face="italic"))

### arrange the plots ###
ggarrange(amplicon_genus_rarefaction,shotgun_subspecies_rarefaction,nrow=1,ncol=2,labels=c("A","B"))



#===========================================================================================================#
###### rarefy the samples at different rarefaction levels ######
# randomly rarefied our samples at different rarefaction levels 100 times. 
# If the number of reads in a sample is less than the rarefaction level, we directly use all the reads.
# This was done using the 'feature-table rarefy' command in qiime2.

###### function for getting the rarefaction genus count table ######
compute.rarefaction.table <- function(original_table,rarefaction_level){
  
  # features in the original feature table
  features_all <- rownames(original_table)
  # samples in the original feature table
  samples_all <- colnames(original_table)
  
  # one of the rarefaction table
  filename <- paste0("./qiime2/table_l6_rarefaction/table_l6_rarefy_",rarefaction_level,"_1.qza")
  rarefied_table_example <- read_qza(filename)$data
  # samples in the special rarefaction level
  samples <- colnames(rarefied_table_example)
  
  # table for saving the rarefied count table (all the features, samples in that rarefaction level)
  rarefied_count_table <- as.data.frame(matrix(0,nrow=length(features_all),ncol=length(samples),
                                               dimnames=list(features_all[order(features_all)],samples[order(samples)])))
  
  # perform rarefaction 100 times
  for(i in 1:100){
    filename <- paste0("./qiime2/table_l6_rarefaction/table_l6_rarefy_",rarefaction_level,"_",i,".qza")
    df <- as.data.frame(read_qza(filename)$data)
    # uniform the format of the rarefied table, make them have the same features
    features_not_included <- features_all[-which(features_all%in%rownames(df))]
    features_added <- matrix(0,nrow=length(features_not_included),ncol=length(samples),dimnames=list(features_not_included,samples))
    df <- rbind(df,features_added)
    df <- df[order(rownames(df)),order(colnames(df))]
    # add the rarefaction table together
    rarefied_count_table <- rarefied_count_table+df
  }
  rarefied_count_table <- rarefied_count_table/100
  colnames(rarefied_count_table) <- sub("^s","",colnames(rarefied_count_table))
  
  # add the missing samples from the original table
  samples_not_included <- samples_all[-which(samples_all%in%colnames(rarefied_count_table))]
  samples_added <- original_table[order(rownames(original_table)),samples_not_included]
  rarefied_table_final <- cbind(rarefied_count_table,samples_added)
  # remove extra samples, keep the same samples with original table
  rarefied_table_final <- rarefied_table_final[,samples_all]
  
  # return
  return(rarefied_table_final)
}

# rarefaction table for 1k per reads
genus_table_rarefaction_1k <- compute.rarefaction.table(original_table=amplicon_genus_table,rarefaction_level="1k")
write.csv(genus_table_rarefaction_1k, "amplicon_tables/genus_table_rarefaction_1k.csv")

# rarefaction table for 5k per reads
genus_table_rarefaction_5k <- compute.rarefaction.table(original_table=amplicon_genus_table,rarefaction_level="5k")
write.csv(genus_table_rarefaction_5k, "amplicon_tables/genus_table_rarefaction_5k.csv")

# rarefaction table for 10k per reads
genus_table_rarefaction_10k <- compute.rarefaction.table(original_table=amplicon_genus_table,rarefaction_level="10k")
write.csv(genus_table_rarefaction_10k, "amplicon_tables/genus_table_rarefaction_10k.csv")

# rarefaction table for 30k per reads
genus_table_rarefaction_30k <- compute.rarefaction.table(original_table=amplicon_genus_table,rarefaction_level="30k")
write.csv(genus_table_rarefaction_30k, "amplicon_tables/genus_table_rarefaction_30k.csv")

# rarefaction table for 50k per reads
genus_table_rarefaction_50k <- compute.rarefaction.table(original_table=amplicon_genus_table,rarefaction_level="50k")
write.csv(genus_table_rarefaction_50k, "amplicon_tables/genus_table_rarefaction_50k.csv")

# rarefaction table for 100k per reads
genus_table_rarefaction_100k <- compute.rarefaction.table(original_table=amplicon_genus_table,rarefaction_level="100k")
write.csv(genus_table_rarefaction_100k, "amplicon_tables/genus_table_rarefaction_100k.csv")



#===========================================================================================================#
##### PCoA plots at different rarefaction level ######
require(vegan)
require(ape)
require(ggplot2)
require(ggpubr)
plot.pcoa <- function(filename,metadata,plot_title){
  
  # feature table
  feature_table <- read.csv(filename,row.names=1)
  colnames(feature_table) <- sub("^X","",colnames(feature_table))
  
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
  pcoa_plot <- ggplot(pcoa_df,aes(x=-Axis.1,y=-Axis.2,col=group,shape=group)) +
    geom_point(size=3) +
    xlim(-0.6,0.7) + ylim(-0.6,0.6) +
    labs(x=paste0("PCOA1(",pro1,"%)"), y=paste0("PCOA2(",pro2,"%)"),title=plot_title) +
    scale_color_manual(breaks=c("Healthy","UC","UC with A","UC with B","UC with BI","UC with BIS","UC with S"),
                       values=c("#F8766D","#00BFC4","#00BFC4","#00BFC4","#00BFC4","#00BFC4","#00BFC4"))+
    scale_shape_manual(breaks=c("Healthy","UC","UC with A","UC with B","UC with BI","UC with BIS","UC with S"),
                       values=c(16,16,1,2,3,4,7))+
    geom_vline(aes(xintercept=0),linetype="dotted") +
    geom_hline(aes(yintercept=0),linetype="dotted") +
    stat_ellipse(aes(group=SampleType,col=SampleType), level = 0.8, show.legend=FALSE)+
    annotate("text",label="Phenotypes:",x=-0.6,y=-0.4,size=4,hjust=0)+
    geom_text(data=permanova_labels,mapping=aes(x=-0.6,y=-0.47,label=r2),parse=TRUE,inherit.aes=FALSE,size=4,hjust=0)+
    geom_text(data=permanova_labels,mapping=aes(x=-0.6,y=-0.54,label=p.value),parse=TRUE,inherit.aes=FALSE,size=4,hjust=0)+
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
amplicon_genus_pcoa_1k <- plot.pcoa("amplicon_tables/genus_table_rarefaction_1k.csv",metadata,"1k reads per sample")
amplicon_genus_pcoa_5k <- plot.pcoa("amplicon_tables/genus_table_rarefaction_5k.csv",metadata,"5k reads per sample")
amplicon_genus_pcoa_10k <- plot.pcoa("amplicon_tables/genus_table_rarefaction_10k.csv",metadata,"10k reads per sample")
amplicon_genus_pcoa_30k <- plot.pcoa("amplicon_tables/genus_table_rarefaction_30k.csv",metadata,"30k reads per sample")
amplicon_genus_pcoa_50k <- plot.pcoa("amplicon_tables/genus_table_rarefaction_50k.csv",metadata,"50k reads per sample")
amplicon_genus_pcoa_100k <- plot.pcoa("amplicon_tables/genus_table_rarefaction_100k.csv",metadata,"100k reads per sample")

# arrange the plot
ggarrange(amplicon_genus_pcoa_1k,amplicon_genus_pcoa_5k,amplicon_genus_pcoa_10k,
          amplicon_genus_pcoa_30k,amplicon_genus_pcoa_50k,amplicon_genus_pcoa_100k,
          ncol=3,nrow=2,labels=c("A","B","C","D","E","F"))



#===========================================================================================================#
##### box plot for Shannon diversities at different rarefaction level ######
require(vegan)
require(ggplot2)
require(ggpubr)
plot.shannon.box <- function(filename,metadata,plot_title){
  
  # feature table
  feature_table <- read.csv(filename,row.names=1)
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
amplicon_genus_shannon_1k <- plot.shannon.box("amplicon_tables/genus_table_rarefaction_1k.csv",metadata,"1k reads per sample")
amplicon_genus_shannon_5k <- plot.shannon.box("amplicon_tables/genus_table_rarefaction_5k.csv",metadata,"5k reads per sample")
amplicon_genus_shannon_10k <- plot.shannon.box("amplicon_tables/genus_table_rarefaction_10k.csv",metadata,"10k reads per sample")
amplicon_genus_shannon_30k <- plot.shannon.box("amplicon_tables/genus_table_rarefaction_30k.csv",metadata,"30k reads per sample")
amplicon_genus_shannon_50k <- plot.shannon.box("amplicon_tables/genus_table_rarefaction_50k.csv",metadata,"50k reads per sample")
amplicon_genus_shannon_100k <- plot.shannon.box("amplicon_tables/genus_table_rarefaction_100k.csv",metadata,"100k reads per sample")

# arrange the plot
ggarrange(amplicon_genus_shannon_1k,amplicon_genus_shannon_5k,amplicon_genus_shannon_10k,
          amplicon_genus_shannon_30k,amplicon_genus_shannon_50k,amplicon_genus_shannon_100k,
          ncol=3,nrow=2,labels=c("A","B","C","D","E","F"))



#===========================================================================================================#
# sequencing depth had little effect on alpha and beta divesities according to the above results
# in the following analysis, we directly use all the data

###### ploting the pcoa plot and permanova ######
require(vegan)
require(ape)
require(ggplot2)
require(ggpubr)

plot.pcoa <- function(feature_table,metadata,plot_title){
  
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
  pcoa_df$shotgun.SampleID <- rownames(pcoa_df)
  pcoa_df <- merge(pcoa_df,metadata,by="shotgun.SampleID")
  # axis
  pro1 <- as.numeric(sprintf("%.3f",PCOA$values[,"Relative_eig"][1]))*100  
  pro2 <- as.numeric(sprintf("%.3f",PCOA$values[,"Relative_eig"][2]))*100
  # plot
  pcoa_plot <- ggplot(pcoa_df,aes(x=Axis.1,y=Axis.2,col=group,shape=group)) +
    geom_point(size=3) +
    xlim(-0.6,0.6) + ylim(-0.6,0.6) +
    labs(x=paste0("PCOA1(",pro1,"%)"), y=paste0("PCOA2(",pro2,"%)"),title=plot_title) +
    scale_color_manual(breaks=c("Healthy","UC","UC with A","UC with B","UC with BI","UC with BIS","UC with S"),
                       values=c("#F8766D","#00BFC4","#00BFC4","#00BFC4","#00BFC4","#00BFC4","#00BFC4"))+
    scale_shape_manual(breaks=c("Healthy","UC","UC with A","UC with B","UC with BI","UC with BIS","UC with S"),
                       values=c(16,16,1,2,3,4,7))+
    geom_vline(aes(xintercept=0),linetype="dotted") +
    geom_hline(aes(yintercept=0),linetype="dotted") +
    stat_ellipse(aes(group=SampleType,col=SampleType), level = 0.8, show.legend=FALSE)+
    annotate("text",label="Phenotypes:",x=-0.6,y=-0.4,size=4,hjust=0)+
    geom_text(data=permanova_labels,mapping=aes(x=-0.6,y=-0.47,label=r2),parse=TRUE,inherit.aes=FALSE,size=4,hjust=0)+
    geom_text(data=permanova_labels,mapping=aes(x=-0.6,y=-0.54,label=p.value),parse=TRUE,inherit.aes=FALSE,size=4,hjust=0)+
    theme_bw() +
    theme(legend.title=element_blank())+
    theme(title = element_text(size = 14))+
    theme(axis.title = element_text(size = 16),axis.text = element_text(size = 12,colour="black"))+
    theme(legend.justification=c(0.02,0.98), legend.position=c(0.02,0.98),legend.background = element_rect(fill = NA))
  
  # return
  return(pcoa_plot)
}

# shotgun
rownames(metadata) <- metadata$shotgun.SampleID
shotgun_subspecies_pcoa <- plot.pcoa(shotgun_subspecies_table,metadata,"PCoA plot using shotgun data")

# amplicon
rownames(metadata) <- metadata$SampleID
amplicon_genus_pcoa <- plot.pcoa(amplicon_genus_table,metadata,"PCoA plot using 16S data")

# arrange the plots
ggarrange(amplicon_genus_pcoa,shotgun_subspecies_pcoa,nrow=1,ncol=2,labels=c("A","B"))



#===========================================================================================================#
### correlation for PCA1 between shotgun and 16s ###
require(vegan)
require(ape)
require(ggplot2)
require(ggpubr)

# shotgun subspecies pcoa1
shotgun_subspecies_ab <- as.data.frame(apply(shotgun_subspecies_table,2,function(x) x/sum(x)))
shotgun_subspecies_beta_div <- as.matrix(vegdist(t(shotgun_subspecies_ab),method = "bray"))
shotgun_subspecies_pcoa <- pcoa(as.dist(shotgun_subspecies_beta_div))
shotgun_subspecies_pcoa_df <- data.frame(shotgun.Axis1=shotgun_subspecies_pcoa$vectors[,1])
shotgun_subspecies_pcoa_df$shotgun.SampleID <- rownames(shotgun_subspecies_pcoa_df)
pcoa_corr_df <- merge(shotgun_subspecies_pcoa_df,metadata,by="shotgun.SampleID")

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



#===========================================================================================================#
### correlation for shannon index between shotgun and 16s ###
require(vegan)
require(ggplot2)
require(ggpubr)

# shotgun shannon
shotgun_subspecies_shannon <- diversity(t(shotgun_subspecies_table),index='shannon')
shotgun_subspecies_shannon_df <- data.frame(shotgun.SampleID=names(shotgun_subspecies_shannon),shotgun.shannon=shotgun_subspecies_shannon)
shannon_corr_df <- merge(shotgun_subspecies_shannon_df,metadata,by="shotgun.SampleID")

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



#===========================================================================================================#
###### differential analysis of 16S rRNA ######
require(qiime2R)
require(edgeR)

### function for reading the amplicon count table at different levels
read.feature.table <- function(filename,samples){
  df <- as.data.frame(read_qza(filename)$data)
  colnames(df) <- gsub("^s","",colnames(df))
  # filter samples
  df <- df[,which(colnames(df)%in%samples)]
  return(df)
}
### feature tables
amplicon_phylum_table  <- read.feature.table("qiime2/table-l2.qza",samples=metadata$SampleID)
amplicon_class_table   <- read.feature.table("qiime2/table-l3.qza",samples=metadata$SampleID)
amplicon_order_table   <- read.feature.table("qiime2/table-l4.qza",samples=metadata$SampleID)
amplicon_family_table  <- read.feature.table("qiime2/table-l5.qza",samples=metadata$SampleID)
amplicon_genus_table   <- read.feature.table("qiime2/table-l6.qza",samples=metadata$SampleID)
amplicon_species_table <- read.feature.table("qiime2/table-l7.qza",samples=metadata$SampleID)

### function for differential analysis
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
write.csv(amplicon_species_diff[amplicon_species_diff$FDR<0.25,],"amplicon_tables/amplicon_species_differential.csv")



#===========================================================================================================#
###### preparation for the differential analysis of shotgun data ######
require(taxonomizr)

# relative abundance table for shotgun data
# difference between shotgun_abundance.csv and count_table_shotgun.csv: 
# shotgun_abundance.csv is relative abundance for samples and the feature name is taxon id rather than species name.
# to obtain the full linkage and the abundance tables at different level, we need to use shotgun_abundance.csv
shotgun_ab <- read.csv("shotgun_tables/shotgun_abundance.csv",row.names=1)
# save as the txt file for converting to biom format
write.table(shotgun_ab,"shotgun_tables/shotgun_abundance.txt",sep="\t",quote=F)

### linkage for the shotgun species according to the taxon id
shotgun_linkage <- data.frame(taxon_id=NA,superkingdom=NA,phylum=NA,class=NA,order=NA,family=NA,genus=NA,species=NA)
for(i in 1:nrow(shotgun_ab)){
  shotgun_linkage[i,"taxon_id"] <- sub("^X","",rownames(shotgun_ab)[i])
  shotgun_linkage[i,2:8] <- getTaxonomy(shotgun_linkage[i,"taxon_id"],"accessionTaxa.sql")
}

### correct for the NA values mannully
# taxon id 330: Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Pseudomonadaceae; Pseudomonas; Pseudomonas oleovorans 
# taxon id 219334 (merged into 1423732): Bacteria; Firmicutes; Bacilli; Lactobacillales; Lactobacillaceae; Lacticaseibacillus; Lacticaseibacillus casei 
# taxon id 471821 (merged into 1408204): Bacteria; Elusimicrobia; Endomicrobia; Endomicrobiales; Endomicrobiaceae; Endomicrobium; Candidatus Endomicrobium trichonymphae
# taxon id 1476577 (merged into 2093824): Bacteria; Candidatus Saccharibacteria; Candidatus Saccharimonia; Candidatus Nanosynbacterales; Candidatus Nanosynbacteraceae; Candidatus Nanosynbacter; Candidatus Nanosynbacter lyticus
# taxon id 1759437 (merged into 2741499): Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Yersiniaceae; Serratia; Serratia surfactantfaciens
# taxon id 1849491 (merged into 1817405): Bacteria; Firmicutes; Bacilli; Bacillales; Staphylococcaceae; Abyssicoccus; Abyssicoccus albus

# taxon id 330, couldn't find the linkage, add it mannually
shotgun_linkage[grep("^330$",shotgun_linkage$taxon_id),2:8] <- c("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas", "Pseudomonas oleovorans")
shotgun_linkage[grep("^219334$",shotgun_linkage$taxon_id),2:8] <- c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lacticaseibacillus", "Lacticaseibacillus casei")
shotgun_linkage[grep("^471821$",shotgun_linkage$taxon_id),2:8] <- c("Bacteria", "Elusimicrobia", "Endomicrobia", "Endomicrobiales", "Endomicrobiaceae", "Endomicrobium", "Candidatus Endomicrobium trichonymphae")
shotgun_linkage[grep("^1476577$",shotgun_linkage$taxon_id),2:8] <- c("Bacteria", "Candidatus Saccharibacteria", "Candidatus Saccharimonia", "Candidatus Nanosynbacterales", "Candidatus Nanosynbacteraceae", "Candidatus Nanosynbacter", "Candidatus Nanosynbacter lyticus")
shotgun_linkage[grep("^1759437$",shotgun_linkage$taxon_id),2:8] <- c("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Yersiniaceae", "Serratia", "Serratia surfactantfaciens")
shotgun_linkage[grep("^1849491$",shotgun_linkage$taxon_id),2:8] <- c("Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Abyssicoccus", "Abyssicoccus albus")

# save the results
write.csv(shotgun_linkage,"shotgun_tables/shotgun_linkage.csv")

### taxonomy file needed by qiime
shotgun_taxonomy <- data.frame(Feature_ID=shotgun_linkage$taxon_id,
                               Taxon=paste0("k__",shotgun_linkage$superkingdom,
                                               "; p__",shotgun_linkage$phylum,
                                               "; c__",shotgun_linkage$class,
                                               "; o__",shotgun_linkage$order,
                                               "; f__",shotgun_linkage$family,
                                               "; g__",shotgun_linkage$genus,
                                               "; s__",shotgun_linkage$species))
write.table(shotgun_taxonomy,"shotgun_tables/shotgun_taxonomy.tsv",quote=F,row.names=F,sep="\t")

# import the taxonomy file and abundance file to qiime and summarize the abundance at different taxonomy levels
#biom convert -i shotgun_original_abundance.txt -o shotgun_original_abundance.biom --table-type="OTU table" --to-hdf5
#biom add-metadata -i shotgun_original_abundance.biom --observation-metadata-fp  shotgun_taxonomy.tsv -o shotgun_abundance.biom --sc-separated taxonomy --observation-header OTUID,taxonomy 
#biom convert -i shotgun_abundance.biom -o shotgun_abundance.txt --table-type="OTU table" --to-tsv
#summarize_taxa.py -i shotgun_abundance.biom -L 7 -o ./



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
###### differential analysis for shotgun data ######
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
write.csv(shotgun_species_diff[shotgun_species_diff$FDR<0.25,],"shotgun_tables/shotgun_species_differential.csv")



#=================================================================================================#
###### compare shotgun and 16s with vila et. al. ######
# read the comparison file
comparison_df <- read.csv("amplicon_tables/vila_comparison.csv")
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










