#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Rscript to remove outliers from a table
# Load the package
library(ggstatsplot)

#Most used methods : Z-score method and the Interquartile Range (IQR)
#Use IQR method because it does not depend on the mean and standard deviation of a dataset 

#Concept
#The interquartile range is the central 50% or the area between the 75th and the 25th percentile of a distribution. 
#A point is an outlier if it is above the 75th or below the 25th percentile by a factor of 1.5 times the IQR.

#Load data : By default 1st column is considered as "Name" & 2nd column as "property" to be quality checked
df <- read.table(args[1], header=as.logical(args[2]))

#quantile function to find the 25th and the 75th percentile of the dataset
Q <- quantile(as.matrix(df[2]), probs=c(.25, .75), na.rm = FALSE)

#IQR function gives the difference of the 75th and 25th percentiles.
iqr <- IQR(as.matrix(df[2]))

#find the cut-off ranges
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range﻿

#subset function to extract the part of your dataset between the upper and lower ranges leaving out the outliers. 
filtered<-subset(df, as.matrix(df[2]) > (Q[1] - 1.5*iqr) & as.matrix(df[2]) < (Q[2]+1.5*iqr))

#If there are not outliers just headings would be printed
print(filtered) 

#link
#https://www.reneshbedre.com/blog/find-outliers.html
