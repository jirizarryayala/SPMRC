#--------------------------
# BIMI 6100
# Crosslin 20230907
# Introduction to R
#--------------------------
rm(list=objects()); options(stringsAsFactors = FALSE)
library(Hmisc)
library(rms)
options(digits=3)
#--------------------------
dir <- "/Users/davidcrosslin/google_drive/tulane/graduate_program/classes/bimi_6100/"
#/Volumes/GoogleDrive/My Drive/tulane/graduate_program/classes/bimi_6100/"

# read in data
rhc <- read.csv(paste(dir,"2023/r/data/rhc.csv",sep=""),header=T,as.is=T);
# the does the same as above
#rhc <- read.table(paste(dir,"2023/r/data/rhc.csv",sep=""),header=T,as.is=T,sep=",");
#help(read.csv )

# using describe function in Hmisc
d1 <- describe(rhc[,-1])
help(describe)

# output to latex
latex(d1,file=paste(dir,"2023/r/latex/subfiles/rhc.tex",sep=""))

#print to text file
sink(file=paste(dir,"2023/r/latex/subfiles/rhc.txt",sep=""))
print(d1)
sink()

# using summary function in Hmisc---------------------------

#---must make categorical as factors when using summary 
rhc$sex<- as.factor(rhc$sex)
rhc$death<- as.factor(rhc$death)
rhc$race <- as.factor(rhc$race)

s1 <- summary(sex~death+race+age+bili1+crea1, data=rhc,continuous = 2,overall=TRUE, method='reverse', test=TRUE,   prmsd=TRUE)      
print(s1)

# output to latex
latex(s1,file=paste(dir,"2023/r/latex/subfiles/rhc_by_sex.tex",sep=""), prn=TRUE)

#print to text file
sink(file=paste(dir,"2023/r/latex/subfiles/rhc_by_sex.txt",sep=""))
print(s1)
sink()

#------------------------------------------------
# look at correlation
cor(rhc$meanbp1,rhc$wblc1)
plot(rhc$meanbp1,rhc$wblc1)

cor(rhc$meanbp1,rhc$hrt1)
plot(rhc$meanbp1,rhc$hrt1)

cor(rhc$age,rhc$hrt1)
plot(rhc$age,rhc$hrt1)

cor(rhc$age,rhc$wtkilo1)
plot(rhc$age,rhc$meanbp1)

cor(rhc$sadmdte, rhc$dschdte)
cor(rhc$sadmdte, rhc$dschdte,  use ="complete.obs")
plot(rhc$sadmdte, rhc$dschdte)

#--------------------------
# BIMI 6100
# Crosslin 20231019
# Data Management 
#--------------------------

# SQLite database
library(RSQLite)

#setting the driver....can be MySQL, Oracle etc.
drv <- dbDriver("SQLirte")  

#set connection
con <- dbConnect(drv, dbname = paste(dir,"2023/r/data/bimi_6100.db",sep=""))
# dbDisconnect(con) # to disconnect database driver

dbWriteTable(con, "rhc", rhc, row.names = F)
#dbExecute(con,"drop table rhc")

dbListTables(con)
dbListFields(con = con, name = "rhc")

# get a distinct count of the number of distinct records
count <- dbGetQuery(con, "select distinct count(*) as count from rhc"); print(count) #5735

count_death <- dbGetQuery(con, "select distinct  count(*) as count from rhc where death in ('Yes')"); print(count_death) #3722

age_summary <- dbGetQuery(con, "select death, avg(age) as average_age from rhc group by death"); print(age_summary) 

# dbExecute(con,"drop table age")
dbExecute(con, "create table age as select death, avg(age) as average_age from rhc group by death"); 

temp_age_table <- dbGetQuery(con,"select * from age"); print(temp_age_table) 

dbExecute(con, "create table rhc2 as select * from rhc as a inner join age as b on a.death=b.death"); 

rhc2 <- dbGetQuery(con,"select * from rhc2"); head(rhc2) 

# writing data directly from R
write.csv(rhc, file=paste(dir,"2023/r/data/rhc_output.csv",sep=""),row.names=FALSE,quote=FALSE)
save(rhc, file = paste(dir,"2023/r/data/rhc_output.RData",sep="")) 
rhc_temp <- get(load(file = paste(dir,"2023/r/data/rhc_output.RData",sep=""))); head(rhc_temp ); dim(rhc_temp )

# more tips for slicing and extracting data in R: https://cran.r-project.org/doc/contrib/Short-refcard.pdf


#--------------------------
# BIMI 6100
# Crosslin 20231019
# Regression modeling strategies 
#--------------------------
rhc$sex<- as.factor(rhc$sex)
rhc$cardiohx<- as.factor(rhc$cardiohx)
rhc$race <- as.factor(rhc$race)


attach(rhc)
dd<- datadist(rhc)
detach(rhc)
options(datadist="dd")  # for plotting w hmisc, rms
options(digits=3)

m <- ols(wblc1~race+sex+income+cardiohx, data=rhc)
summary(m)
anova(m)
plot(summary(m))

latex(summary(m),file=paste(dir,"2023/r/latex/subfiles/rhc_ols_summary.tex",sep=""))

options(prType='latex')
latex(anova(m),file=paste(dir,"2023/r/latex/subfiles/rhc_ols_anova.tex",sep=""))


png(file=paste(dir, "2023/r/latex/subfiles/rhc_ols_summary.png",sep=""))
plot(summary(m))
dev.off()

png(file=paste(dir, "2023/r/latex/subfiles/rhc_ols_anova.png",sep=""))
plot(anova(m))
dev.off()

m <- ols(wblc1~race+sex+income+cardiohx, data=rhc)

# -------------------logistic regression   
m <- lrm(cardiohx~wblc1+race+sex+income, data=rhc )
summary(m)
anova(m)
plot(summary(m))
plot(nomogram(m))

options(prType='latex')
latex(m,file=paste(dir,"2023/r/latex/subfiles/rhc_lrm.tex",sep=""))

latex(anova(m),file="/google_drive/bime/classes/mebi_537/R/latex_output/m2_anova.tex")
latex(summary(m),file="/google_drive/bime/classes/mebi_537/R/latex_output/m2_summary.tex")

png(file=paste(dir, "2023/r/latex/subfiles/rhc_ols_summary.png",sep=""))
plot(summary(m))
dev.off()

png(file=paste(dir, "2023/r/latex/subfiles/rhc_ols_anova.png",sep=""))
plot(anova(m))
dev.off()
