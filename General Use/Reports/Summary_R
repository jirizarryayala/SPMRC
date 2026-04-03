

#------1
#LINE1 SPMRC analysis consult starter code 
# for this code, you need some packages for data wrangling and some for statistics 
#Hmisc and rms for statistics 
#data.table and dplyr for data wrangling 
#------
library(Hmisc) # use install.packages("Hmisc")
library(rms)   # use install.packages("rms")
library(dplyr) # it should be lazy loaded in the new version of the R 4/4.4.1
library(tidyr) #datawrangling ; i guess tidyr includes dplyr so either one is fine 
# to learn more about each packages use help(Hmisc); help(rms); help(dplyr)
#to learn more about function of each package use help again : help(describe); help(summary)


#-----2
#just reading the data
#make sure to change the directory to read your file 
#-----
dir <- "Users/Kam/google_drive/tulane/SPMRC/"   # this is where you access your file 
line1 <- read.csv(paste(dir/line.csv", sep=""), header=T, as.is=T) # reading the data



#----data cleaning 3
#using base function for cbind() , rbind() 
# or using tidyr function for bind_row(), and creating ad, and ctrl
# this was my cleanin approach 
# you can use cbind() to bind columns create a case and control column 
# or you can do rbind() to add case and control row for patients data 
# as long as anything allows you to collapse case and control in its own file, is fine 
#---
# a. Add group labels and pivot to long format
ad <- line1 %>%
  pivot_longer(cols = everything(), values_to = "Expression", names_to = "Patient") %>%
  mutate(Condition = factor("AD"))

ctrl <- line1 %>%
  pivot_longer(cols = everything(), values_to = "Expression", names_to = "Patient") %>%
  mutate(Condition = factor("Control"))

#b. combine data 
data <- bind_rows(ad, ctrl) # this will stack up rows together 


#----
#looking at the data pre and post 
# describe() and summary()function from Hmisc
#describe() function shows an overall view of data : min, max, median , mean 
#latex() and sink() are how you can save your data nicely #install.packages("xtable")
sink() with txt can save output as a text file to read in the text editor 
#---

sink(file=paste(dir, "this_is_where_your_latex_file_is_going/line_describe.txt", sep=""))
d1 <- describe(line1[, -1]) 
d2 <- describe(ad[, -1])
d3 <- describe(ctrl[, -1])

print(d1); print(d2); print(d3)
latex(d1, file=paste(dir, "this_is_where_your_latex_file_is_going/line_describe.tex", sep=""))
latex(d2, file=paste(dir, "this_is_where_your_latex_file_is_going/line_describe.tex", sep=""))
latex(d3, file=paste(dir, "this_is_where_your_latex_file_is_going/line_describe.tex", sep=""))
sink()


#-----
#use Summary function in Hmisc for pre vs post 
#summary()
#------

# the stat formula is Dependent ~ Independent (Expression ~ Condition)
summary <- summary(Expression ~ Condition, 
                        data = combined_data, 
                        continuous = 2, 
                        overall = TRUE, 
                        test = TRUE) # test=TRUE is appropriate here as you have two groups

print(summary)
