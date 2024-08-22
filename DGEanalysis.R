#This code is for Differential gene expression between mir483 overexpressed/knockdown livers vs controls.

library(edgeR)
data<-read.delim("raw_count.tsv",row.names=1,sep="\t")
meta<-read.delim("metadata.txt",row.names=2,sep="\t")
genes<-read.delim("mouse_gene_names.ens100.txt",sep="\t",row.names=1)
genes<-genes[order(rownames(genes)),]
y<-DGEList(data)
y$genes<-genes[rownames(y),]
y$samples$group<-as.factor(paste(substring(rownames(meta),1,6),meta$Genotype,sep="_"))
plotMDS(y,col=as.numeric(c(y$samples$group)))

#overexpression of Mir483 in 13 weeks
x<-y[,1:12]
keep1<-rowSums(cpm(x)>=1)>=6
table(keep1)
# keep1
# FALSE  TRUE 
# 42941 12460 
x<-x[keep1,,keep.lib.sizes=F]
x<-calcNormFactors(x)
plotMDS(cpm(x),col=as.numeric(c(x$samples$group)))
#looks well seperated on 1st dim (86%)
x_design<-model.matrix(~0+x$samples$group)
colnames(x_design)<-c("mutant","wt")
x<-estimateDisp(x,x_design)
plotBCV(x)
#fit
x_fit<-glmQLFit(x,x_design)
#contrast
x_qlftest<-glmQLFTest(x_fit,contrast=c(1,-1))
topTags(x_qlftest)
table(decideTestsDGE(x_qlftest,lfc = 0))
write.table(topTags(x_qlftest,n="inf"),"13w_miR483overexpression_liver_comparison_mutvswt.txt",sep="\t",quote=F,col.names=NA)
library(pheatmap)
x_cpm<-cpm(x,prior.count = 5)

#draw heatmap to check
pheatmap(x_cpm[decideTestsDGE(x_qlftest)!=0,],scale = "row")
#export cpm to table
x_cpm<-as.data.frame(cbind(x_cpm,x$genes))
write.table(x_cpm,"normalised_cpm_13w_samples_only.txt",sep="\t",col.names=NA,quote=F)


#KD in E19
z<-y[,13:24]
keep2<-rowSums(cpm(z)>=1)>=6
table(keep2)
# keep2
# FALSE  TRUE 
# 41753 13648 
z<-z[keep2,,keep.lib.sizes=F]
z<-calcNormFactors(z)
plotMDS(cpm(z),col=as.numeric(c(z$samples$group)))
#looks well seperated on 1st dim (86%)
z_design<-model.matrix(~0+z$samples$group)
colnames(z_design)<-c("mutant","wt")
z<-estimateDisp(z,z_design)
plotBCV(z)
#fit
z_fit<-glmQLFit(z,z_design)
#contrast
z_qlftest<-glmQLFTest(z_fit,contrast=c(1,-1))
topTags(z_qlftest)
table(decideTestsDGE(z_qlftest,lfc = 0))
# 0 
# 13648 

write.table(topTags(z_qlftest,n="inf"),"e19_miR483KD_liver_comparison_mutvswt.txt",sep="\t",quote=F,col.names=NA)

library(pheatmap)
z_cpm<-cpm(z,prior.count = 5)

#draw heatmap to check
pheatmap(z_cpm[decideTestsDGE(z_qlftest,adjust.method = "none",lfc=0.5)!=0,],scale = "row")
#58 and 59 are outliers

#export cpm to table
z_cpm<-as.data.frame(cbind(z_cpm,z$genes))
write.table(z_cpm,"normalised_cpm_e19_samples_only.txt",sep="\t",col.names=NA,quote=F)

#removing 58 and 59
z<-z[,-c(6,7)]
colnames(z)
z_design<-model.matrix(~0+z$samples$group)
colnames(z_design)<-c("mutant","wt")
z_design
z<-estimateDisp(z,z_design)
plotBCV(z)
#fit
z_fit<-glmQLFit(z,z_design)
#contrast
z_qlftest<-glmQLFTest(z_fit,contrast=c(1,-1))
topTags(z_qlftest)
table(decideTestsDGE(z_qlftest,lfc = 0))
# -1     0     1 
# 395 12902   351 
write.table(topTags(z_qlftest,n="inf"),"e19_miR483KD_liver_comparison_removed_418_423_mutvswt.txt",sep="\t",quote=F,col.names=NA)
z2_cpm<-cpm(z,prior.count = 5)
#looks a lot better
pheatmap(z2_cpm[decideTestsDGE(z_qlftest)!=0,],scale = "row")


#to export miR483 normalised counts
keep<-rowSums(cpm(y)>=1)>=5
table(keep)
y<-y[keep,,keep.lib.sizes=F]
y<-calcNormFactors(y)
y_cpm<-cpm(y,prior.count=5)
y_cpm<-as.data.frame(cbind(y_cpm,y$genes))
write.table(y_cpm,"normalised_cpm_all_samples.txt",sep="\t",col.names=NA,quote=F)
