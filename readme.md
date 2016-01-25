Andrew Fairless, February 2011
modified May 2015 for posting onto Github
This script analyzes whether litter size affects sociability by running a
robust OP regression on the residuals of a trimmed means ANOVA with sex and 
age as grouping factors in the C57BL/6J MICE, as described in Fairless et al 2012
Fairless et al 2012, doi: 10.1016/j.bbr.2011.12.001, PMID:  22178318, PMCID:  PMC3474345

The fictional data in "altereddata.txt" were modified from the original 
empirical data used in Fairless et al 2012.
I am using fictional data instead of the original data because I do not have 
permission of my co-authors to release the data into the public domain.  
NOTE:  Because these data are fictional, several important characteristics of
these data may be different from those of the original data.

Each row is a separate mouse.
The left-most 3 columns are quasi-independent variables (mouse strain, sex, and age).
The right-most 5 columns are dependent variables describing behaviors of the
mice during the Social Approach/Choice Test.

analysis as described in Fairless et al 2012:
"This result was confirmed by a simple OP regression of litter size and the
residuals of a 2 × 5 (sex by age) trimmed means ANOVA of social cylinder 
investigation, which yielded a slope significantly less than zero, D = 0.004, 
p < 0.001 (Fig. 3c)."
