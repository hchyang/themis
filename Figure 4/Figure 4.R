### Figure 4A. tSNE plot for individual biomarkers
library(Rtsne)
group.df=read.table('group.txt',header = T)
fem.df=read.table('EndMotif_frequency_matrix.txt',header = T)
group.df$PID=fem.df$PID[match(group.df$SampleID,fem.df$Sample)]
set.seed(1)
library(ggplot2)
library(RColorBrewer)
#FEM
group_h.df=subset(group.df,Cancer_type=='HEALTHY')
fem.df=read.table('EndMotif_frequency_matrix.txt',header = T)
fem.m=as.matrix(fem.df[,-c(1:2)])
rownames(fem.m)=fem.df$Sample
fem_h.m=fem.m[which(rownames(fem.m)%in%group_h.df$SampleID),]
fem_h.m=t(fem_h.m)
library(preprocessCore)
fem_t.m=t(fem.m)
for (i in 1:ncol(fem_t.m)){
  temp.m=cbind(fem_h.m,fem_t.m[,i])
  temp.m=normalize.quantiles(temp.m)
  fem.m[i,]=temp.m[,ncol(temp.m)]
}
rownames(fem.m)=fem.df$Sample

tsne<-Rtsne(fem.m,dims = 2,pca = T,max_iter = 2000,theta = 0.05,perplexity = 10,verbose = T,pca_center = F,normalize =F)
tsne_out=as.data.frame(tsne$Y)
tsne_out$type=group.df$Cancer_type[match(rownames(fem.m),group.df$SampleID)]
tsne_out$type=factor(tsne_out$type,levels=c('HEALTHY','BRCA','NSCLC','COREAD','STAD','ESCA','LIHC','PACA'))
png('FEM_tSNE.png',width = 2000,height = 2000,res = 600)
p<-ggplot(tsne_out,aes(V1,V2,color=type,shape=type))+geom_point()+xlab('tSNE 1')+ylab('tSNE 2')+ggtitle('FEM')+theme(plot.title = element_text(hjust = 0.5,face='bold'),panel.background = element_rect(fill = "white",color = 'black'),legend.position = 'none')+scale_shape_manual(values = c(17,rep(20,7)))+scale_color_manual(values = c('darkgray',brewer.pal(8,'Dark2')))
plot(p)
dev.off()
#CAFF
caff.df=read.table('CAFF_arm_level_zscore_matrix.txt',header = T)
caff.m=as.matrix(caff.df[,-c(1:2)])
rownames(caff.m)=caff.df$Sample
for (i in 1:nrow(caff.m)){
  caff.m[i,]=(caff.m[i,]-mean(caff.m[i,]))/sd(caff.m[i,])
}
tsne<-Rtsne(caff.m,dims = 2,pca = F,max_iter = 2000,theta = 0.05,perplexity = 10,verbose = T,pca_center = F,normalize = F)
tsne_out=as.data.frame(tsne$Y)
tsne_out$type=group.df$Cancer_type[match(rownames(caff.m),group.df$SampleID)]
tsne_out$type=factor(tsne_out$type,levels=c('HEALTHY','BRCA','NSCLC','COREAD','STAD','ESCA','LIHC','PACA'))
png('CAFF_tSNE.png',width = 2000,height = 2000,res = 600)
p<-ggplot(tsne_out,aes(V1,V2,color=type,shape=type))+geom_point()+xlab('tSNE 1')+ylab('tSNE 2')+ggtitle('CAFF')+theme(plot.title = element_text(hjust = 0.5,face='bold'),panel.background = element_rect(fill = "white",color = 'black'),legend.position = 'none')+scale_shape_manual(values = c(17,rep(20,7)))+scale_color_manual(values = c('darkgray',brewer.pal(8,'Dark2')))
plot(p)
dev.off()

#MFR
mfr.df=read.table('fragment.stat_Methy.Ratio.txt',header = T,check.names = F)
mfr.m=t(as.matrix(mfr.df[,-c(1)]))
for (i in 1:nrow(mfr.m)){
  mfr.m[i,]=(mfr.m[i,]-mean(mfr.m[i,]))/sd(mfr.m[i,])
}
tsne<-Rtsne(mfr.m,dims = 2,pca = F,max_iter = 2000,theta = 0.05,perplexity = 10,verbose = T,pca_center = F,normalize = F)
tsne_out=as.data.frame(tsne$Y)

tsne_out$type=group.df$Cancer_type[match(rownames(mfr.m),group.df$PID)]

tsne_out$type=factor(tsne_out$type,levels=c('HEALTHY','BRCA','NSCLC','COREAD','STAD','ESCA','LIHC','PACA'))
png('MFR_tSNE.png',width = 2000,height = 2000,res = 600)
p<-ggplot(tsne_out,aes(V1,V2,color=type,shape=type))+geom_point()+xlab('tSNE 1')+ylab('tSNE 2')+ggtitle('MFR')+theme(plot.title = element_text(hjust = 0.5,face='bold'),panel.background = element_rect(fill = "white",color = 'black'),legend.position = 'none')+scale_shape_manual(values = c(17,rep(20,7)))+scale_color_manual(values = c('darkgray',brewer.pal(8,'Dark2')))
plot(p)
dev.off()

#FSI
fsi.df=read.table('05.first_ratio.corrected.matrix_tsne_zscore.txt',header = T,check.names = F)
fsi.m=t(as.matrix(fsi.df[,-c(1)]))
tsne<-Rtsne(fsi.m,dims = 2,pca = F,max_iter = 2000,theta = 0.05,perplexity = 10,verbose = T,pca_center = F,normalize = F)
tsne_out=as.data.frame(tsne$Y)
tsne_out$type=group.df$Cancer_type[match(rownames(fsi.m),group.df$SampleID)]
tsne_out$type=factor(tsne_out$type,levels=c('HEALTHY','BRCA','NSCLC','COREAD','STAD','ESCA','LIHC','PACA'))
png('FSI_tSNE.png',width = 2000,height = 2000,res = 600)
p<-ggplot(tsne_out,aes(V1,V2,color=type,shape=type))+geom_point()+xlab('tSNE 1')+ylab('tSNE 2')+ggtitle('FSI')+theme(plot.title = element_text(hjust = 0.5,face='bold'),panel.background = element_rect(fill = "white",color = 'black'),legend.position = 'none')+scale_shape_manual(values = c(17,rep(20,7)))+scale_color_manual(values = c('darkgray',brewer.pal(8,'Dark2')))
plot(p)
dev.off()

### Figure 4B. ROC plots of individual models.
data.df=read.csv('prediction_ scores.csv',header = T,check.names = F,na.strings = 'AAA')
data.df$Diagnosis=factor(data.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))

validation.df=subset(data.df,Cohort=='Test')
validation.df$Diagnosis=factor(validation.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))

discovery.df=subset(data.df,Cohort=='Training')
discovery.df$Diagnosis=factor(discovery.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))

cancer_type.v=c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD')

plot_roc_overall_combine_omics_with_ensemble<-function(comp,cohort){
  if(cohort=='Training'){data.df=discovery.df}
  if(cohort=='Test'){data.df=validation.df}
  op_names=c('MFR','FSI','CAFF','FEM','THEMIS')
  data.df=cbind(data.df,cancer='OVERALL')
  data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
  data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
  tiff(paste0('Overall_',cohort,'_omics_with_ensemble_ROC.tif'),width = 3400,height = 3400,res = 600,compression = 'lzw',pointsize = 13)
  library(pROC)
  library(RColorBrewer)
  cols=c(brewer.pal(4,'Dark2'),'black')
  plot.roc(data.df$cancer,data.df[,op_names[1]],levels=c('HEALTHY','OVERALL'),ci=T,direction = comp,print.auc=T,col=cols[1],type='S',xlim=c(1,0),ylim=c(0,1),add=F,print.auc.y=0.3-0.05*1,print.auc.x=0.6,lwd=3,main=cohort)
  text(0.61,0.3-0.05*1,labels = paste0(op_names[1],' '),adj = c(1,1),col=cols[1])
  for (i in 2:5){
    plot.roc(data.df$cancer,data.df[,op_names[i]],levels=c('HEALTHY','OVERALL'),ci=T,direction = comp,print.auc=T,col=cols[i],type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.3-0.05*i,print.auc.x=0.6,lwd=3)
    text(0.61,0.3-0.05*i,labels = paste0(op_names[i],' '),adj = c(1,1),col=cols[i])
    
  }
  dev.off()
}
plot_roc_overall_combine_omics_with_ensemble('<','Training')
plot_roc_overall_combine_omics_with_ensemble('<','Test')

### Figure 4C. Venn diagram depicting overlapping true and false positive.
data.df=read.csv('prediction_ scores.csv',header = T,check.names = F,na.strings = 'AAA')
colnames(data.df)[c(6:10)]=c('CAFF','FSI','MFR','FEM','THEMIS')
validation.df=subset(data.df,Cohort=='Test')
validation.df$Diagnosis=factor(validation.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','STAD','LIHC','NSCLC','PACA'))

discovery.df=subset(data.df,Cohort=='Training')
discovery.df$Diagnosis=factor(discovery.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','STAD','LIHC','NSCLC','PACA'))
data.df=cbind(data.df,cancer='OVERALL')
data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
cutoff.v=rep(NA,5)
for (i in 1:5){
  h_train.df=subset(discovery.df,Diagnosis=='HEALTHY')
  idx=round(nrow(h_train.df)*0.01)
  cutoff.v[i]=sort(h_train.df[,5+i],decreasing = T)[idx]
}
names(cutoff.v)=colnames(discovery.df[6:10])
# False positives.
pos.l=list()
for (i in 1:5){
  data_cancer.df=subset(data.df,Diagnosis!='HEALTHY')
  pos.l[[i]]=data_cancer.df[which(data_cancer.df[,5+i]>=cutoff.v[i]),]$PID
}
names(pos.l)=colnames(discovery.df)[6:10]
library(RColorBrewer)
library(venn)
library(ggplot2)
library(ggpolypath)
tiff('Overlap_0.99_spec_cancer_training_test.tif',res=600,compression = 'lzw',height = 3500,width = 3500,pointsize = 12)
venn(list(`MFR`=pos.l[[3]],`FSI`=pos.l[[2]],`CAFF`=pos.l[[1]],`FEM`=pos.l[[4]],`THEMIS`=pos.l[[5]]),zcolor = brewer.pal(5,'Dark2')[c(1,2,3,4,5)],sncs = 1.4,ilcs=1.1,box = F)
dev.off()

# True positives.
pos.l=list()
for (i in 1:5){
  data_cancer.df=subset(data.df,Diagnosis=='HEALTHY')
  pos.l[[i]]=data_cancer.df[which(data_cancer.df[,5+i]>=cutoff.v[i]),]$PID
}
names(pos.l)=colnames(discovery.df)[6:10]

library(RColorBrewer)
library(venn)
library(ggplot2)
library(ggpolypath)
tiff('Overlap_0.99_spec_healthy_training_test.tif',res=600,compression = 'lzw',height = 3500,width = 3500,pointsize = 12)
venn(list(`MFR`=pos.l[[3]],`FSI`=pos.l[[2]],`CAFF`=pos.l[[1]],`FEM`=pos.l[[4]],`THEMIS`=pos.l[[5]]),zcolor = brewer.pal(5,'Dark2')[c(1,2,3,4,5)],sncs = 1.4,ilcs=1.1,box = F)
dev.off()