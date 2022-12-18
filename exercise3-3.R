#code for computer exercise 3

install.packages("ggplot2")
install.packages("reshape2")

library(ggplot2)
library(reshape2)

x = read.table('counts_matrix.txt')
metadata = read.table('metadata.txt',sep='\t',header=TRUE) 

#Exercise 1
table(metadata$Sex) #get how many patients of each gender are part of the study
table(metadata$diagnosis) #get how many patients have which diagnosis

rownames(x)
colnames(x)
metadata$patient.id

colnames(x) == metadata$patient.id


count.nonzero=function(counts){
  i=0
  for (j in 1:80)
    if (counts[j]>0){i=i+1}
  return(i)
}

nonzero = vector(length=58037)
for(j in 1:58037){
  nonzero[j]=count.nonzero(as.vector(x[j, ]))
}

#this takes very long -> maybe there' a better way?
x.filtered = x
for(j in 58037:1){
  if(nonzero[j] < 20){
    x.filtered = x.filtered[-c(j),]
  }
}

#Exercise 2
cpm = x.filtered+1 #add pseudocount
sum = apply(cpm, 2, sum) #calculate the total mapped reads in each sample
cpm = cpm/(sum/1000000) 
logCPM = log(cpm)

#Exercise 3
#cleaning the data for using ggplot for plotting boxplot
logCPM2 = logCPM
col_names = names(logCPM)
logCPM2 <- cbind(Gene = rownames(logCPM), logCPM)
melt_data = melt(logCPM2, id = c("Gene"))
colnames(melt_data) = c("Gene", "Sample", "Expression")

#plotting boxplot for normalized data using ggplot
p1 = (ggplot(melt_data, mapping = aes(x = factor(Sample), y = Expression))
      + geom_boxplot(color = "blue4", fill="skyblue")
      + labs(y = "Gene expression", x = "Sample",
             title = " The distribution of gene expression in each sample in normalized data")
      + theme_bw()
      + theme(axis.text.x = element_text(angle = 90),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5)
      ))
p1

#logtransform the non-normalized data:
logx = log(x.filtered+0.01) #add a small number to avoid taking log(0)

col_names = names(logx)
logx <- cbind(Gene = rownames(logx), logx)
melt_non_norm = melt(logx, id = c("Gene"))
colnames(melt_non_norm) = c("Gene", "Sample", "Expression")
#plotting boxplot for non_normalized data using ggplot
p2 = (ggplot(melt_non_norm, mapping = aes(x = factor(Sample), y = Expression))
      + geom_boxplot(color = "blue4", fill="skyblue")
      + labs(y = "Gene expression", x = "Sample",
             title = " The distribution of gene expression in each sample in non-normalized data")
      + theme_bw()
      + theme(axis.text.x = element_text(angle = 90),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5)
      ))
p2

#plot SRR1782694 against SRR1782687 using scatter plot
p3 = (ggplot(logCPM, 
             mapping = aes(x = logCPM[,'SRR1782694'], y = logCPM[,'SRR1782687'])
)
+ geom_point(alpha=0.5, size=2)
+ geom_abline(slope=1, color = "red", linewidth = 1.25)
+ labs(x = "Sample SRR1782694", y = "Sample SRR1782687",
       title = "Comparing gene expressions between two different samples")
+ theme_bw()
+ theme( panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5))
)

p3
#Exercise 4
tspan6 = merge(t(logCPM[1,]), metadata,  all = FALSE, by.x=0, by.y=1)

p4 = (ggplot(tspan6, mapping = aes(x = factor(diagnosis), y = tspan6[,2]))
      + geom_boxplot(color = "blue4", fill="skyblue")
      + labs(y = "Gene expression", x = "Samples",
             title = " The gene ENSG00000000003 expression")
      + theme_bw()
      + theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))
)

p4
#linear model:
diagnosis=as.factor(tspan6[,7]) 
diagnosis = relevel(diagnosis, "Not IBD") #set healthy as "standard" value
fit1=lm(tspan6[,2]~diagnosis)
summary(fit1)

gender=as.factor(tspan6[,4])
fit2=lm(tspan6[,2]~diagnosis+gender+tspan6[,5])
summary(fit2)

#Exercise 5
diagnosis=as.factor(metadata[,"diagnosis"])
diagnosis = relevel(diagnosis, "Not IBD")
gender=as.factor(metadata[,"Sex"])

p1 = vector(length=38558)
param1 = vector(length=38558)

p2_diag = vector(length=38558)
p2_gend = vector(length=38558)
p2_age = vector(length=38558)
param2_diag = vector(length=38558)
param2_gend = vector(length=38558)
param2_age = vector(length=38558)


for(j in 1:38558){
  fit1_all = lm(t(logCPM[j,])~diagnosis)
  fit2_all = lm(t(logCPM[j,])~diagnosis+gender+metadata[,4])
  
  p1[j] = summary(fit1_all)$coefficients[2,4]
  param1[j] = summary(fit1_all)$coefficients[2,1]
  
  p2_diag[j] = summary(fit2_all)$coefficients[2,4]
  p2_gend[j] = summary(fit2_all)$coefficients[3,4]
  p2_age[j] = summary(fit2_all)$coefficients[4,4]
  
  param2_diag[j] = summary(fit2_all)$coefficients[2,1]
  param2_gend[j] = summary(fit2_all)$coefficients[3,1]
  param2_age[j] = summary(fit2_all)$coefficients[4,1]
}

p1_adj = p.adjust(p1, "fdr")
p2_diag_adj = p.adjust(p2_diag, "fdr")
p2_gend_adj = p.adjust(p2_gend, "fdr")
p2_age_adj = p.adjust(p2_age, "fdr")

fit1_coeff = cbind(logCPM[,0], p1, p1_adj, param1)
fit2_coeff = cbind(logCPM[,0], p2_diag, p2_gend, p2_age, p2_diag_adj, p2_gend_adj, p2_age_adj, param2_diag, param2_gend, param2_age)

#have a look at the smallest p values
fit1_coeff[order(fit1_coeff[,"p1_adj"], decreasing = F)[1:100],]
fit2_coeff[order(fit2_coeff[,"p2_diag_adj"], decreasing = F)[1:100],]

#count significant genes for both models and the up-regulated sig. genes
i_p1=0
i_p2=0
i_gen=0
i_age=0
i_up=0
for (j in 1:38558){
  if (p1_adj[j] < 0.05){i_p1=i_p1+1}
  if (p2_diag_adj[j] < 0.05){
    i_p2=i_p2+1
    if (param2_diag[j]>0){i_up=i_up+1}
  }
  if (p2_age_adj[j] < 0.05){i_age=i_age+1}
  if (p2_gend_adj[j] < 0.05){i_gen=i_gen+1}
}

#Exercise 6
library(gplots)
library(RColorBrewer)
xsig = as.matrix(logCPM[order(p2_diag_adj, decreasing = F)[1:100],])

mypalette = brewer.pal(11,"RdYlBu") 
morecols = colorRampPalette(mypalette) 
mycols=rev(morecols(255)) 
column.cols=c("purple","orange")[as.factor(metadata$diagnosis)] 


pdf("top100sigGenesHeatmap.pdf",height=10,width=8) 
heatmap.2(xsig,trace='none',col=mycols,main='The 100 most significant 
genes',ColSideColors=column.cols) 
dev.off() 

#Exercise 7
pca = prcomp(t(logCPM))
summary(pca)
pca$x[,1]

p5 = (ggplot(as.data.frame(pca$x), mapping = aes(x=as.vector(pca$x[,1]), y = as.vector(pca$x[,2]), color = diagnosis))
      + geom_point(alpha=1, size=3)
      + labs(x = "PCA1", y = "PCA2",
             title = "Comparing PCA1 and PCA2")
      + theme_bw()
      + scale_color_manual(values = c("CD" = "blue", "Not IBD" = "chocolate1"))
      + theme (panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size = 1),
               plot.title = element_text(hjust = 0.5))
)
p5

p6 = (ggplot(as.data.frame(pca$x), mapping = aes(x=as.vector(pca$x[,1]), y = as.vector(pca$x[,3]), color = diagnosis))
      + geom_point(alpha=1, size=3)
      + labs(x = "PCA1", y = "PCA3",
             title = "Comparing PCA1 and PCA3")
      + theme_bw()
      + scale_color_manual(values = c("CD" = "blue", "Not IBD" = "chocolate1"))
      + theme (panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size = 1),
               plot.title = element_text(hjust = 0.5))
)

p6

p7 = (ggplot(as.data.frame(pca$x), mapping = aes(x=as.vector(pca$x[,2]), y = as.vector(pca$x[,3]), color = diagnosis))
      + geom_point(alpha=1, size=3)
      + labs(x = "PCA2", y = "PCA3",
             title = "Comparing PCA2 and PCA3")
      + theme_bw()
      + scale_color_manual(values = c("CD" = "blue", "Not IBD" = "chocolate1"))
      + theme (panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size = 1),
               plot.title = element_text(hjust = 0.5))
)

p7

#color the sample based on their gender
p8 = (ggplot(as.data.frame(pca$x), mapping = aes(x=as.vector(pca$x[,1]), y = as.vector(pca$x[,2]), color = gender))
      + geom_point(alpha=1, size=3)
      + labs(x = "PCA1", y = "PCA2",
             title = "Comparing PCA1 and PCA2")
      + theme_bw()
      + scale_color_manual(values = c("Female" = "purple3", "Male" = "goldenrod1"))
      + theme (panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size = 1),
               plot.title = element_text(hjust = 0.5))
)
p8

p9 = (ggplot(as.data.frame(pca$x), mapping = aes(x=as.vector(pca$x[,1]), y = as.vector(pca$x[,3]), color = gender))
      + geom_point(alpha=1, size=3)
      + labs(x = "PCA1", y = "PCA2",
             title = "Comparing PCA1 and PCA3")
      + theme_bw()
      + scale_color_manual(values = c("Female" = "purple3", "Male" = "goldenrod1"))
      + theme (panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size = 1),
               plot.title = element_text(hjust = 0.5))
)
p9

p10 = (ggplot(as.data.frame(pca$x), mapping = aes(x=as.vector(pca$x[,2]), y = as.vector(pca$x[,3]), color = gender))
       + geom_point(alpha=1, size=3)
       + labs(x = "PCA2", y = "PCA3",
              title = "Comparing PCA2 and PCA3")
       + theme_bw()
       + scale_color_manual(values = c("Female" = "purple3", "Male" = "goldenrod1"))
       + theme (panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size = 1),
                plot.title = element_text(hjust = 0.5))
)
p10

