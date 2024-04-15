#script to load ASV table from cohort 1 into phyloseq
#filter and use phyloseq to dseq2 function to convert for DSEq2
#run DSeq2 on each disease cohort, compare archaea-pos vs. neg.
#phyloseq_import
library(phyloseq)
library(ape)
library(phangorn)

otus<-read.table("table.csv",header=T,row.names=1)
meta<-import_qiime_sample_data("meta.csv")


tax<-as.matrix(read.table("taxa.csv",sep="\t",row.names=1))

colnames(tax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

OTU = otu_table(otus, taxa_are_rows = TRUE)
TAX = tax_table(tax)
tree<-read.tree("SINA_iqtree.tree")
tree<-midpoint(tree)


physeq = phyloseq(OTU, meta,TAX,tree)
physeq->ps

#filter according to https://f1000research.com/articles/5-1492

ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


filterPhyla = c("Chloroflexi", "Deinococcota","Gemmatimonadota","Myxococcota","Thermotogota")

# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps1
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)

# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))

ps3 = tax_glom(ps2, "Genus", NArm = TRUE)

#DESeq testing

library(reshape2)
library(ggplot2)
library(DESeq2)


#generate tables, use subset for UC, IBS and controls
#examplary for UC

psCex<-subset_samples(ps3, Disease%in%c("UC"))

ps_dds <- phyloseq_to_deseq2(psCex, ~ Archaea.final)
diagdds = DESeq(ps_dds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps3)[rownames(sigtab), ], "matrix"))
head(sigtab)

write.csv(sigtab,file="UC_Archaea_DESeq.csv")
