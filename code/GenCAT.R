setwd('/share/lactase/prostate_studies/SEKJ_work/enrichment')
version

LDdata <- read.delim("/share/lactase/uk10k_1000G/BPD_analysis/uk10k_1000G_chr2_with_maf_and_rs.ld_gt_0.2.sorted.tsv",
                     stringsAsFactors = F)

tmp <- strsplit(LDdata$RSID_B, split = ",", fixed = T)
LD2 <- data.frame(RSID = unlist(tmp), UK10KID = rep(LDdata$SNP_B, unlist(lapply(tmp, length))),
                  stringsAsFactors = F)
LD2 <- subset(LD2, RSID != "N/A")

prostate <- read.delim('../meta_analysis/meta-analysis_733_812_1.tbl')

# take snps in prostate data with LD > 0.2
prostate.chr2.LD.gt.0.2 <- merge(prostate, LD2[,c("RSID", "UK10KID")], by.x = "MarkerName", by.y = "RSID", all.x = F, all.y = F)

library(snpStats)
library(GenCAT)
SNPdata <- read.plink(bed = "/share/lactase/uk10k_1000G/BPD_analysis/plink.bed",
                      bim = "/share/lactase/uk10k_1000G/BPD_analysis/plink.bim",
                      fam = "/share/lactase/uk10k_1000G/BPD_analysis/plink.fam",
                      select.snps = prostate.chr2.LD.gt.0.2$UK10KID)
geno_dat <- SNPdata$genotypes

# Reorder prostate data AGAIN to match order in geno_dat
tmp <- match(colnames(geno_dat), prostate.chr2.LD.gt.0.2$UK10KID)
prostate3 <- prostate.chr2.LD.gt.0.2[tmp,]

# Modify data for use with GenCAT
prostate.modified <- prostate3
prostate.modified$effect_allele <- toupper(prostate.modified$Allele1)
prostate.modified$other_allele <- toupper(prostate.modified$Allele2)
prostate.modified$testStat <- prostate.modified$Effect
prostate.modified$SNP <- prostate.modified$UK10KID
prostate.modified$chr <- substr(prostate.modified$UK10KID, 1, 1)
prostate.modified$class <- "A"

snpInfo <- SNPdata$map
snpInfo$chr <- snpInfo$chromosome
snpInfo$SNP <- snpInfo$snp.name
snpInfo$A1 <- snpInfo$allele.1
snpInfo$A2 <- snpInfo$allele.2


PCa_GenCAT <- GenCAT(SNPdata = prostate.modified, genoData = geno_dat,
       snpInfo = snpInfo)
write.table(PCa_GenCAT$GenCAT, file = "PCa_GenCAT_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(PCa_GenCAT$Used, file = "PCa_GenCAT_SNPs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(PCa_GenCAT$notFound, file = "PCa_GenCAT_notFound.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(PCa_GenCAT$unMatched, file = "PCa_GenCAT_unMatched.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(PCa_GenCAT$TransStats, file = "PCa_GenCAT_TransStats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
PCa_GenCAT
