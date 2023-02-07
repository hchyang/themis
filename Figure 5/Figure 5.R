### Figure 5A. ROC plots of THEMIS model performance by cancer types.
data.df=read.csv('../Figure 4/prediction_ scores.csv',header = T,check.names = F,na.strings = 'AAA')
data.df$Diagnosis=factor(data.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))

validation.df=subset(data.df,Cohort=='Test')
validation.df$Diagnosis=factor(validation.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))

discovery.df=subset(data.df,Cohort=='Training')
discovery.df$Diagnosis=factor(discovery.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))

cancer_type.v=c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD')
plot_roc_by_cancer<-function(op_num,comp,cohort){
  if(cohort=='Training'){data.df=discovery.df}
  if(cohort=='Test'){data.df=validation.df}
  op_names=c('FSI','CAFF','MFR','FEM','THEMIS')
  data.df$Diagnosis=factor(data.df$Diagnosis,levels = c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD'))
  cancer_type.v=c('HEALTHY','BRCA','COREAD','ESCA','LIHC','NSCLC','PACA','STAD')
  tiff(paste0(op_names[op_num],'_',cohort,'_by_cancer_ROC.tif'),width = 3950,height = 3950,res = 600,pointsize = 14,compression = 'lzw')
  library(pROC)
  library(RColorBrewer)
  plot.roc(data.df$Diagnosis,data.df[,op_names[op_num]],levels=c('HEALTHY',cancer_type.v[2]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(cancer_type.v),'Dark2')[1],type='S',xlim=c(1,0),ylim=c(0,1),add=F,print.auc.y=0.4-0.05*1,print.auc.x=0.5,main=cohort)
  text(0.51,0.4-0.05*1,labels = paste0(cancer_type.v[2],' '),adj = c(1,1),col=brewer.pal(length(cancer_type.v),'Dark2')[1])
  for (i in 3:length(cancer_type.v)){
    plot.roc(data.df$Diagnosis,data.df[,op_names[op_num]],levels=c('HEALTHY',cancer_type.v[i]),ci=T,direction = comp,print.auc=T,col=brewer.pal(length(cancer_type.v),'Dark2')[i-1],type='S',xlim=c(1,0),ylim=c(0,1),add=T,print.auc.y=0.4-0.05*(i-1),print.auc.x=0.5)
    text(0.51,0.4-0.05*(i-1),labels = paste0(cancer_type.v[i],' '),adj = c(1,1),col=brewer.pal(length(cancer_type.v),'Dark2')[i-1])
  }
  data.df=cbind(data.df,cancer='OVERALL')
  data.df$cancer[which(data.df$Diagnosis=='HEALTHY')]='HEALTHY'
  data.df$cancer[which(data.df$Diagnosis!='HEALTHY')]='OVERALL'
  dev.off()
}
plot_roc_by_cancer(5,'<','Training')
plot_roc_by_cancer(5,'<','Test')

### Figure 5B. THEMIS sensitivity by cancer types and stages
library(tidyverse)
library(stats)
library(openxlsx)
ci_plot <- function(file_in){
  file_in <- file_in
  file_out <- paste0(file_in, '.hist.pdf')
  df <- read.table(file_in, header = T)
  df$Cohort<-factor(df$Cohort,levels = c('Training','Test'))
  df=df[order(df$Cohort,decreasing = F),]
  df <- df[df$Group != 'OV',]
  df <- df[df$Group != 'HEALTHY',]
  #    df <- df[df$Group != 'Cancer',]
  df <- df[df$Stage != 'unknown',]
  df <- df[!is.na(df$Stage),]
  df$Group <- as.character(df$Group)
  #    df$Group <- ifelse(df$Group == 'Cancer', 'Pan-Cancer', df$Group)
  df$Group <- ifelse(df$Group == 'Cancer', 'OVERALL', df$Group)
  #df$Group <- ifelse(df$Group == 'Cancer' , 'Pan-Cancer', df$Group)
  df$CI_low <- as.numeric(gsub("-.*", "", df$CI))
  df$CI_high <- as.numeric(gsub(".*-", "", df$CI))
  df$Sens <- df$Sens
  df$CI_low <- df$CI_low
  df$CI_high <- df$CI_high
  df <- df %>% group_by(Stage, Group) %>% mutate(group_number = paste0(total, collapse=' | '))
  df$Stage_new <- paste0(df$Stage, "\n", '(', df$group_number ,')')
  df <- transform(df, Group = factor(Group, levels = c('BRCA', 'COREAD', 'ESCA', 'LIHC', 'NSCLC', 'PACA','STAD','OVERALL')))
  #    write.table(df, 'df', sep = "\t", quote =FALSE, row.names =FALSE)
  #    df <- df %>% mutate(across(Group, factor, levels = c('BRCA', 'COREAD', 'ESCA', 'STAD', 'LIHC', 'NSCLC', 'PACA', 'Ovreall')))
  Stage_new <- df[df$Stage_new == 'Training',]
  #df <- df[order(df$Group),]
  pdf(file_out, width = 8, height = 9,pointsize = 11)
  p <- ggplot(data = df, aes(x = Stage_new, y = Sens, ymin = CI_low, ymax = CI_high, colour = Cohort, fill = Cohort)) +
    #         geom_histogram(stat = "identity", binwidth = 0.1,  position = position_dodge(0.8)) + 
    geom_histogram(stat = "identity", width = 0.6,  position = position_dodge(0.7)) + 
    geom_errorbar(position = position_dodge(width = 0.7), width = 0.2, color = 'black',size=0.4) + 
    scale_color_manual(values = c("#27B7B7", "#FF66B2")) + 
    scale_fill_manual(values = c("#27B7B7", "#FF66B2")) +
    #         scale_color_manual(values = c("#008080", "#800000")) +
    #         scale_fill_manual(values = c("#008080", "#800000")) +
    #         scale_x_discrete(labels = Stage_new) + 
    #         theme_bw() + 
    #         theme(axis.text = element_text(face = 'bold')) + 
    geom_point(position = position_dodge(width = 0.7), size = 0.5, color = 'black') + 
    theme_bw() +
    facet_wrap(vars(Group), ncol = 3, scales = "free_x") +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),legend.position = 'none') +
    scale_y_continuous(limits  = c(0, 100), labels = function(x) paste0(x, "%"), breaks = seq(0, 100, by = 25)) + 
    ylab('Sensitivity')  + 
    xlab('Stage')
  print(p)
  dev.off()
}
ci_plot('Sen_by_Cancer_and_Stage_spf0.99')
