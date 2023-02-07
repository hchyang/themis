### Figure 6B. Confusion matrices for CSO prediction
library(pheatmap)
library(RColorBrewer)
# Training cohort
mat<-read.delim('merge.MD.first_short.first_long.baseline.txt.result_train.txt',header=T,stringsAsFactors = F)
Group_top1<-dcast(as.data.frame(table(mat$Group,mat$TOP1)),Var1 ~ Var2)
class(Group_top1)
Group_top1
rownames(Group_top1)<-Group_top1$Var1
Group_top1<-Group_top1[,-1]
Group_top1<-Group_top1[Group_sort,Group_sort]
value<-Group_top1
Group_top1$sum<-apply(Group_top1,1,function(x){sum(x)})
Group_top1
tmp<-apply(Group_top1,2,function(x){x/Group_top1$sum})
tmp<-tmp[,-8]
bk<-seq(0,1,by=0.25)
pheatmap(tmp,display_numbers = value ,fontsize_number = 12,cluster_rows = F,cluster_cols = F,angle_col = '45',cellheight = 30,cellwidth = 30,fontsize=10,filename = 'training_TOO_6_groups.pdf',breaks=seq(0,1,length.out = 100),legend_breaks = seq(0,1.25,0.25),legend_labels = seq(0,1.25,0.25)*100,main='Training')
# Test cohort
mat<-read.delim('merge.MD.first_short.first_long.baseline.txt.result_test.txt',header=T,stringsAsFactors = F)
Group_top1<-dcast(as.data.frame(table(mat$Group,mat$TOP1)),Var1 ~ Var2)
class(Group_top1)
Group_top1
rownames(Group_top1)<-Group_top1$Var1
Group_top1<-Group_top1[,-1]
Group_top1<-Group_top1[Group_sort,Group_sort]
value<-Group_top1
Group_top1$sum<-apply(Group_top1,1,function(x){sum(x)})
Group_top1
tmp<-apply(Group_top1,2,function(x){x/Group_top1$sum})
tmp<-tmp[,-8]
bk<-seq(0,1,by=0.25)
pheatmap(tmp,display_numbers = value ,fontsize_number = 12,cluster_rows = F,cluster_cols = F,angle_col = '45',cellheight = 30,cellwidth = 30,fontsize=10,filename = 'testing_TOO_6_groups.pdf',breaks=seq(0,1,length.out = 100),legend_breaks = seq(0,1.25,0.25),legend_labels = seq(0,1.25,0.25)*100,main='Test')