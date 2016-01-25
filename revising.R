# Andrew Fairless, February 2011
# modified May 2015 for posting onto Github
# This script analyzes whether litter size affects sociability by running a
# robust OP regression on the residuals of a trimmed means ANOVA with sex and 
# age as grouping factors in the C57BL/6J MICE, as described in Fairless et al 2012
# Fairless et al 2012, doi: 10.1016/j.bbr.2011.12.001, PMID:  22178318, PMCID:  PMC3474345

# The fictional data in "altereddata.txt" were modified from the original 
# empirical data used in Fairless et al 2012.
# I am using fictional data instead of the original data because I do not have 
# permission of my co-authors to release the data into the public domain.  
# NOTE:  Because these data are fictional, several important characteristics of
# these data may be different from those of the original data.

# Each row is a separate mouse.
# The left-most 3 columns are quasi-independent variables (mouse strain, sex, and age).
# The right-most 5 columns are dependent variables describing behaviors of the
# mice during the Social Approach/Choice Test.

# analysis as described in Fairless et al 2012:
# "This result was confirmed by a simple OP regression of litter size and the
# residuals of a 2 × 5 (sex by age) trimmed means ANOVA of social cylinder 
# investigation, which yielded a slope significantly less than zero, D = 0.004, 
# p < 0.001 (Fig. 3c)."

source("Rallfun-v13.txt")

alldata = read.table("altereddata.txt", header = TRUE)      
strainsplit = split(alldata, alldata$strain)
mousestrain = 2
data = strainsplit[[mousestrain]]                                               # C57BL/6J mice only
depvarcol = 18                                                                  # dependent variable column; specifies 'soc.3rd5min'
datatemp = data.matrix(cbind(data[ , 2], data[ , 3], data[ , depvarcol]))       # creates matrix with only 'sex', 'age', and 'soc.3rd5min' variables

model = t2way(2, 5, datatemp, MAT = T)            # robust trimmed means ANOVA
mns = model$means                                 # table of means for each group defined by sex and age
rownames(mns) = c("Female", "Male")
colnames(mns) = c("20", "24", "32", "43", "71")

# adds column with group mean of 'soc.3rd5min' for each mouse
data[ , dim(data)[[2]] + 1] = NA                  
colnames(data)[dim(data)[[2]]] = "groupmeans"
for(iter in 1:dim(datatemp)[[1]]) {
     if(datatemp[iter, 2] == 20) data[iter, dim(data)[[2]]] = mns[datatemp[iter, 1], 1]
     if(datatemp[iter, 2] == 24) data[iter, dim(data)[[2]]] = mns[datatemp[iter, 1], 2]
     if(datatemp[iter, 2] == 32) data[iter, dim(data)[[2]]] = mns[datatemp[iter, 1], 3]
     if(datatemp[iter, 2] == 43) data[iter, dim(data)[[2]]] = mns[datatemp[iter, 1], 4]
     if(datatemp[iter, 2] == 71) data[iter, dim(data)[[2]]] = mns[datatemp[iter, 1], 5]
}

# adds column with residual for each mouse (i.e., difference between each mouse's 
# score for 'soc.3rd5min' and the mouse's group mean as defined by sex and age)
lastcol = dim(data)[[2]]
data[ , dim(data)[[2]] + 1] = data[ , depvarcol] - data[ , lastcol]
colnames(data)[dim(data)[[2]]] = "resids"

strainname = gsub("[[:punct:]]", "", names(strainsplit)[mousestrain])   # removes "/" from strain name
sink(file = paste("littersize,", strainname, ".txt", sep = ""))

paste("This analysis covers the", names(strainsplit)[mousestrain], "mice")
print("robust trimmed means ANOVA")
paste("The dependent variable is:  ", colnames(data)[depvarcol])
paste("Factor 1 or A is:  ", colnames(data)[2])
paste("Factor 2 or B is:  ", colnames(data)[3])
t2way(2, 5, datatemp, MAT = T)

# robust OP regression on litter size and the residuals of the trimmed means ANOVA
opreg(data$littersize, data$resids)
# tests whether slope of regression line is significantly different from zero
regtest(data$littersize, data$resids, regfun = opreg, nboot = 1000, plotit = TRUE, grp = c(1))
print("'$test' denotes test statistic value")
print("'$crit' denotes critical value for test statistic")
print("'$p.value' denotes p value for test statistic")

sink(file = NULL)
