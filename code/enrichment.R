setwd('/share/lactase/prostate_studies/SEKJ_work/enrichment')
version
prostate <- read.delim('../meta_analysis/meta-analysis_733_812_1.tbl', stringsAsFactors = F)
LDdata <- read.delim("/share/lactase/uk10k_1000G/BPD_analysis/uk10k_1000G_chr2_with_maf_and_rs.ld_gt_0.2.sorted.tsv",
                     stringsAsFactors = F)
LD.snps <- unlist(strsplit(LDdata$RSID_B, split = ",", fixed = T))

prostate$in.LD <- ifelse(prostate$MarkerName %in% LD.snps, 1, 0) #this pull the SNPs in LD with LP to test
prostate$P.value <- as.numeric(prostate$P.value)

#correct for multiple comparisons
tmp <- table(ifelse(prostate$P.value < 3e-8, 1, 0), prostate$in.LD) # this is approximately the Bonferroni correction for 1954033 SNPS (calculated by [alpha]/[number of SNPs] = Bonferroni Comparison
tmp
#y = if SNP is sig (0=no,1=yes); x = if SNP is in LD block (0=no,1=yes)
prop.table(tmp, 2)

fisher.test(tmp, alternative = "greater")

wilcox.test(P.value ~ in.LD, data = prostate, alternative = "greater") #are the P.values the same for those in LD and no in LD with LP
summary(prostate$P.value[which(prostate$in.LD == 0)])
summary(prostate$P.value[which(prostate$in.LD == 1)])

pca_SNPs_inLD <- prostate[which(prostate$in.LD == 1),]
write.table(pca_SNPs_inLD, file = "PCa_SNPs_inLD.txt", quote = FALSE, sep = "\t", row.names = FALSE)
prostate[which(prostate$in.LD == 1),]

prostate$pval.upper <- pnorm(prostate$Effect/prostate$StdErr, lower.tail = F)
prostate$pval.lower <- pnorm(prostate$Effect/prostate$StdErr)

wilcox.test(pval.upper ~ in.LD, data = prostate, alternative = "greater")
wilcox.test(pval.lower ~ in.LD, data = prostate, alternative = "greater")

summary(prostate$pval.upper[which(prostate$in.LD == 0)])
summary(prostate$pval.upper[which(prostate$in.LD == 1)])

summary(prostate$pval.lower[which(prostate$in.LD == 0)])
summary(prostate$pval.lower[which(prostate$in.LD == 1)])
