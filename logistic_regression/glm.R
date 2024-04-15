#logistic regression in R

glm_table<-read.csv("log_regression.csv",header=TRUE,sep="\t",row.names=1)
str(glm_table)

glm_table$disease<-as.factor(glm_table$disease)
glm_table$Archaea_final<-as.factor(glm_table$Archaea_final)
glm_table$sex<-as.factor(glm_table$sex)
glm_table$endoscopic_biofilm<-as.factor(glm_table$endoscopic_biofilm)

mylogit<-glm(Archaea_final~disease+Abx.yes_no+PPI_yes_no+endoscopic_biofilm,data=glm_table, family="binomial")
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))