library(Maaslin2)
library(data.table)
library(tidyverse)

#fread instead of read.table due to problems with feature names that contain ' and " ect
#convert to rownames afterwards

#use Maaslin2 to generate linear models, for analysis of archaea influence on microbiota composition

#cohort 2

a <- fread(file="maaslin_table_new_final.csv",sep="\t",header=T)  %>% as.data.frame()
row.names(a) <- a$sample
table<-subset(a, select = -c(sample))

meta<-read.table("maaslin_meta_new_final.csv",header=T,row.names=1,sep="\t")

fit_data_ixn = Maaslin2(input_data     = table, 
                        input_metadata = meta, 
                        normalization  = "TSS",
                        transform = "LOG",
                        output         = "maaslin2", 
                        fixed_effects  = c("Sample_classification","Disease",	"Archaea_final"),
                        reference      = c("Sample_classification,biofilm negative","Disease,control",	"Archaea_final,no"),
                        min_prevalence = 0.1,
                        min_abundance  = 0.0001)

#cohort 1

a <- fread(file="maaslin_table_OTU_counts.csv",sep="\t",header=T)  %>% as.data.frame()
row.names(a) <- a$sample
table<-subset(a, select = -c(sample))

meta<-read.table("maaslin_meta.csv",header=T,row.names=1,sep="\t")


fit_data_ixn = Maaslin2(input_data     = table, 
                        input_metadata = meta, 
                        normalization  = "TSS",
                        transform = "LOG",
                        output         = "maaslin2_counts_tss_log", 
                        fixed_effects  = c("Sample_classification","Disease2",	"Archaea_final"),
                        reference      = c("Sample_classification,biofilm negative","Disease2,control",	"Archaea_final,no"),
                        min_prevalence = 0.1,
                        min_abundance  = 0.0001)
