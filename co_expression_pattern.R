#Rscript
#Author: Yunlong Ma
#Co-expression patterns by using the Pearson correlation analysis
#Usage: it is designed to uncover the co-expression patterns among genes associated with Childhood-onset asthma (CoA)

#Set up the working directory depended on user own need.
setwd("C:\\Users\\Administrator\\Desktop\\corplot")
setwd("F:\\Desktop\\corplot")

#Part I install package

#Installing corrplot package of R platform 
if(!require("corrplot"))install.packages("corrplot")
if(!require("corrr"))install.packages("corrr")
if(!require("dplyr"))install.packages("dplyr")
if(!require("ggplot2"))install.packages("ggplot2")
if(!require("reshape2"))install.packages("reshape2")



#Load the corrplot package
library(corrplot)
library(corrr)
library(dplyr)
library(ggplot2)
library(reshape2)



set.seed(141079)


#Take a example for using the corrplot package
#y <- cor(mtcars)
#corrplot(y,method="circle")
#corrplot(y,method="circle",tl.col = "black")
#corrplot.mixed(y,tl.col = "black")
#corrplot.mixed(y, lower.col = "black", number.cex = .7,tl.col = "black")
#corrplot(y, method = "pie",tl.col = "black")


#Part II read relevant RNA expression data of co-expression analysis

#Read and re-organized CoA-related expression data
pearson<- read.delim("asthma.txt",header = T)
pearson_new2<-pearson[,c(-1,-2)] 
mat_pearson<-as.matrix(pearson_new2)
row.names(mat_pearson)<- pearson[,2]
#write.table(mat_pearson,file = "CoA_gene_expression.txt",append=FALSE, row.names=TRUE,col.names = FALSE)
data_for_CoA <-t(mat_pearson)
correlation_CoA <- cor(data_for_CoA )


#Read and re-organized Control-related expression data
pearson_new<- read.delim("control.txt",header = T)
pearson_new3<-pearson_new[,c(-1,-2)]
mat_pearson1<-as.matrix(pearson_new3)
row.names(mat_pearson1)<- pearson_new[,2]
#write.table(mat_pearson1,file = "control_gene_expression.txt",append=FALSE, row.names=TRUE,col.names = FALSE)
data_for_Control <-t(mat_pearson1)
correlation_Control <- cor(data_for_Control)



#Part III Visualization

#Ploting the identified genes co-expression patterns in the CoA group
corrplot(correlation_CoA, method = "pie",tl.col = "black",tl.cex=0.8)
write.csv(correlation_CoA, file = "correlation_CoA.csv", 
            row.names=TRUE)


#Ploting the identified genes co-expression patterns in the control group
corrplot(correlation_Control, method = "pie",tl.col = "black",tl.cex=0.8,number.cex = 0.8)
write.csv(correlation_Control, file = "correlation_Control.csv", 
       row.names=TRUE)


#Part IV Use corrr package to visualize gene-gene co-expression patterns in a network


#Co-expression patterns for CoA
data_CoA <- data_for_CoA %>% correlate() %>%  rearrange() %>%  shave()
ydata_CoA <- fashion(data_CoA)
data_for_CoA %>% correlate() %>%  rearrange() %>% rplot()

#Visualization for co-expression network
#Set the correlation efficient >=0.3 to show
data_for_CoA %>% correlate() %>%  rearrange() %>% network_plot(min_cor = 0.3) 

#Co-expression patterns for Control
data_control <- data_for_Control %>% correlate() %>%  rearrange() %>%  shave()
ydata_control <- fashion(data_control)
data_for_Control %>% correlate() %>%  rearrange() %>% rplot()
#Visualization for co-expression network
data_for_Control %>% correlate() %>%  rearrange() %>% network_plot(min_cor = 0.3)


#Part V Paired Student's t-test

# Create correlation data frame (cor_CoA)
temp_CoA<- data_for_CoA %>% correlate() %>%  rearrange() %>% stretch()
cor_data_CoA_file <- na.omit(temp_CoA)
write.csv(cor_data_CoA_file, file = "cor_data_CoA_file.csv", 
          row.names=TRUE)
#cor_data_CoA <- na.omit(temp_CoA$r)

# Create correlation data frame (cor_Control)
temp_control <- data_for_Control %>% correlate() %>%  rearrange() %>% stretch()
cor_data_Control_file <- na.omit(temp_control)
write.csv(cor_data_Control_file, file = "cor_data_Control_file.csv", 
          row.names=TRUE)
#cor_data_Control <- na.omit(temp_control$r)


#define a function for match the correlation data betwen CoA and control
Cor_reshape <- function(x,y){
   len = length(x[,1])
   for (i in seq(1,len)){
     for(j in seq(1,len)){
       if (x[i,1]==y[j,1] & x[i,2]==y[j,2]){
           x[i,4] <- y[j,3]
       }
     }
   }
  colnames(x) <- c("Gene1","Gene2","CoA","Control")
  final_data <- na.omit(x)
  return(final_data)
}

#use the Cor_reshape for data organization
data <- Cor_reshape(temp_CoA,temp_control)

#Paired t-test
tt <- t.test(data$CoA,data$Control, paired=TRUE) 
Pvalue <- tt[3]

#Plotting the density figure
md <- melt(data,id=c("Gene1","Gene2"))
p<-ggplot(md, aes(x = value))
p + geom_density(color = "black", fill = "gray")
p + geom_density(aes(color =variable)) 
p + geom_density(aes(fill = variable), alpha=0.6)



#End

