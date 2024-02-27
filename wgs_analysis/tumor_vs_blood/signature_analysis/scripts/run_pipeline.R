# command line args
args = commandArgs(trailingOnly = TRUE)
mutation_file <- args[1]
signature_file <- args[2]
outpath <- args[3]
cores <- as.integer(args[4])

library(pbmcapply)

# load mSigTools and mSigAct
setwd("/workspace/projects/sjd_pediatric_tumors/code/run_msigact/mSigAct")
source("/workspace/projects/sjd_pediatric_tumors/code/run_msigact/mSigAct/mSigAct.v0.10.R")
source("/workspace/projects/sjd_pediatric_tumors/code/run_msigact/mSigAct/mSigTools.v0.13.R")

# TESTING CONFIGURATION
# mutation_file <- '/workspace/users/fmuinos/AML/synthetic/synthetic_data/catalogue_SBS35_100.tsv'
# signature_file <- '/workspace/users/fmuinos/AML/synthetic/synthetic_data/signatures_SBS35.tsv'
# outpath <- '/workspace/users/fmuinos/AML/synthetic/runs/'
# cores <- as.integer("30")

name_file <- basename(mutation_file)

outpath_name <- paste(outpath, name_file, sep ='/')
dir.create(outpath_name, showWarnings = FALSE)

sigs <- as.matrix(read.table(signature_file, sep = '\t', header = T, row.names = 1))
muts <- as.matrix(read.table(mutation_file, sep = '\t', header = T, row.names = 1))

names_sigs <- colnames(sigs)
print(colnames(sigs))
target_signature <- names_sigs[length(names_sigs)]

mSigAct <- process.one.group(muts, sigs,
                             target.sig.name = target_signature,
                             path.root = outpath_name,
                             obj.fun = obj.fun.nbinom.maxlh,
                             nbinom.size=10, ## = dispersion parameter
                             mc.cores=cores) ## = number of cores to use

pval<-mSigAct$pval

apval<-p.adjust(pval,"bonferroni")
exposure<-mSigAct$exposure
df<-t(rbind(pval,apval,exposure))
df<-df[order(df[,1],decreasing = F),]

name_outfile <- paste("results", name_file, "mSigAct", target_signature, "tsv", sep ='.')
write.table(df, file = paste(outpath, name_outfile, sep ="/"), sep ='\t')

name_outfile <- paste("results", name_file, "mSigAct", target_signature, "rds", sep ='.')
saveRDS(mSigAct, file = paste(outpath, name_outfile, sep ="\t"))
