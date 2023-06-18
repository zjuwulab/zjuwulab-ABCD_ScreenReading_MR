library(OneSampleMR)
library(data.table)
library(SNPRelate)
# parameters
screens = c("tvm","video","game","text","socialN","Vchat");
cognitions = c("picvocab","flanker","pattern","picture","cread");
behaviors = c("anxdep","withdep","somatic","social","thought","attention","rulebreak","aggressive");
p_thrh = 1e-5 
maf_thrp = 0.005

# prepare data
exp_name = "screens/tvm"
out_name = "behaviors/anxdep";
revse    = 0 # if 1 apply reverse analyses
if (revse == 1){
  temp = exp_name
  exp_name = out_name
  out_name = temp
  rm(temp)
}

# prepare data
geneQC <- "/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_imputed/imputed_QC08/merge"                         # genetic data after preprocessing
gwa_exp <- fread(paste0("/data/limingyang/Projects/ABCD/Projects/GWAS/reports/",exp_name,"/GWAS/ABCD_gwa.fastGWA"))   # GWAS summary data of exposure
phe_exp <- fread(paste0("/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_phen/phenotypes/",exp_name,".phen"))  # phenotype data of exposure

gwa_out <- fread(paste0("/data/limingyang/Projects/ABCD/Projects/GWAS/reports/",out_name,"/GWAS/ABCD_gwa.fastGWA"))  # GWAS summary data of outcome
phe_out <- fread(paste0("/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_phen/phenotypes/",out_name,".phen")) # phenotype data of outcome

# confounders
cov_inc  <- fread("/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_phen/phenotypes/income.phen")   # income
cov_edu  <- fread("/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_phen/phenotypes/edu.phen")      # education
cov_age  <- fread("/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_phen/covs/age.txt")             # age
cov_sex  <- fread("/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_phen/covs/sex_num.txt")         # sex
dat_cov  <- cbind(cov_inc$V3,cov_edu$V3,cov_age$V3,cov_sex$V3)                                            
dat_mri  <- fread("/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_phen/covs/MRIinfo.txt",head = F) # MRI serial number, when involve MRI derivative data (brain volume)
dat_mri  <- dat_mri$V3
dat_sit  <- fread("/data/limingyang/Projects/ABCD/Projects/GWAS/data/ABCD_phen/covs/site.txt",head = F)    # site
dat_sit  <- dat_sit$V3

# keep snps which pass the p_thrh
dat_exp <- subset(gwa_exp, P < p_thrh)
dat_wrt <- dat_exp[,"SNP"]
keeptxt <- paste0(resfder,"/keep.txt") 
write.table(dat_wrt, keeptxt, sep = "\t", row.names = F, col.names = F, quote = F)
filenm <- paste0(resfder,"/gwas_sig")
cmd = paste0("plink --bfile ",geneQC," --extract ",keeptxt," --make-bed --out ",filenm)
system(cmd) 

# LD pruned
setwd(resfder)
snpgdsBED2GDS(bed.fn = "gwas_sig.bed", 
              bim.fn = "gwas_sig.bim", 
              fam.fn = "gwas_sig.fam", 
              out.gdsfn = "gwas_sig.gds")                
gds <- snpgdsOpen("gwas_sig.gds")
maf <- dat_exp$AF1
maf_thr <- which(maf >= maf_thrp)
snp_sel <- dat_exp$SNP[maf_thr]
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=1e6, snp.id = snp_sel,
                          ld.threshold=sqrt(0.001) ,start.pos = "first", verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)  
       
# get pruned ID of SNPs
dat_snp <- dat_exp$SNP
snpids  <- matrix(data = 0, nrow = length(pruned), ncol = 1)
for (i in 1:length(pruned)){
  snpids[i] = which(dat_snp == pruned[i])
}                  
  
# extract significant results and remove snps which are significantly associated with outcome                              
dat_out_rf <- data.frame(SNP = gwa_out$SNP, pval = gwa_out$P)
dat_exp_ld <- dat_exp[snpids,]
pruned_all <- pruned
dat_exp_ld <- merge(dat_exp_ld,dat_out_rf,by = "SNP")
dat_eo <- subset(dat_exp_ld, pval < p_thrh)
if (dim(dat_eo)[1] > 0){
  snpnm = dat_eo$SNP
  pruned_all = subset(pruned_all,pruned_all != snpnm)
}

# code the genotype as 0 , 1 , 2 bsaed on the number of effective alleles
dat_gen <- snpgdsGetGeno(gds, snp.id=pruned_all, verbose=TRUE)
# add exposure and outcome
data_exp <- phe_exp$V3
data_out <- phe_out$V3
# one sample MR - ivreg function
fit <- ivreg::ivreg(data_out ~ dat_cov + data_exp + dat_sit | dat_cov + dat_sit + dat_gen) 

snpgdsClose(gds)
write.table(res_tab,paste0(resfder,"/diff_pvalue.txt"),sep = "\t",row.names = F)