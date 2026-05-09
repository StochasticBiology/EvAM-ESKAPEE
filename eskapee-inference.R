library(arrow)
library(dplyr)
library(tidyr)
library(hyperinf)
library(ggrepel)
library(ggpubr)
library(hypermk2)

expt = 0
run.mash = TRUE
run.hmk2 = TRUE
fname = paste0("eskapee-phylo-fits-", expt, "-", run.hmk2, ".Rdata", collapse="")

# see details
# https://huggingface.co/datasets/ayates/amr_portal/blob/main/README.md?utm_source=chatgpt.com
# important fields here: antibiotic_name, BioSample_ID, resistance_phenotype

#system("wget https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2025-11/phenotype.parquet")
df = read_parquet("phenotype.parquet")

if(expt == 0) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus", 
              "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Enterococcus faecium",
              "Enterobacter")
  to.get = 4
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus", 
              "Acinetobacter baumannii", "Enterococcus faecium",
              "Enterobacter")
}

# simpler set  
if(expt == 1) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae")
  to.get = 12
}

# simpler set, smaller  
if(expt == 4) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae")
  to.get = 5
}

if(expt == 5) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae", 
              "Acinetobacter baumannii", "Pseudomonas aeruginosa")
  to.get = 8
}

if(expt == 2) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae", "Acinetobacter baumannii")
  to.get = 9
}

if(expt == 3) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae", 
              "Acinetobacter baumannii", 
              "Enterobacter")
  to.get = 7
}

df$eskapee = df$species
df$eskapee[grep("[Ee]nterobacter", df$species)] = "Enterobacter"

# 1. Clean + binarize + collapse duplicates
df_bin <- df %>%
  filter(resistance_phenotype %in% c("resistant", "susceptible"),
         eskapee %in% ESKAPEE) %>%
  mutate(value = ifelse(resistance_phenotype == "resistant", 1, 0)) %>%
  group_by(BioSample_ID, assembly_ID, eskapee, antibiotic_name) %>%
  summarise(value = max(value), .groups = "drop")

# summarise counts of drug-bug pairs
tmp = df_bin
counts = tmp %>% count(eskapee, antibiotic_name) 
drugs = c()
ggplot(counts, aes(x=eskapee, y=antibiotic_name, fill=n)) + geom_tile()

# count appearances across pathogens for each drug
appears = data.frame()
for(drug in unique(df_bin$antibiotic_name)) {
  appears = rbind(appears, data.frame(drug=drug,
                                      species=length(which(counts$antibiotic_name==drug)),
                                      min = min(counts$n[counts$antibiotic_name==drug]),
                                      sum = sum(counts$n[counts$antibiotic_name==drug]))
  )
}
appears = appears[appears$species == length(ESKAPEE),]
appears = appears[order(-appears$min),]

# "appears" now summarises appearances across pathogens
drugs = appears$drug[1:to.get]

# pull the resistance profiles for these drugs into a wide dataframe 
wide_all <- tmp[tmp$antibiotic_name %in% drugs,] %>%
  pivot_wider(names_from = antibiotic_name,
              values_from = value)  %>%
  drop_na()

# output specifics
cat(nrow(wide_all))
unique(wide_all$eskapee)
table(wide_all$eskapee)

build_tree = function(filename) {
  #lines <- readLines("mash_distances.tab")
  lines <- readLines(filename)
  
  # Number of genomes
  n <- as.integer(lines[1])
  
  # Initialize distance matrix
  dist_mat <- matrix(0, n, n)
  
  # Extract filenames
  files <- character(n)
  files[1] <- strsplit(lines[2], "/")[[1]][[4]]  # first genome
  for(i in 2:n){
    parts <- strsplit(lines[i + 1], "\\s+")[[1]]
    tmp <- parts[1]
    files[i] <- strsplit(tmp, "/")[[1]][[4]]
  }
  rownames(dist_mat) <- colnames(dist_mat) <- files
  
  # Fill lower triangle
  # first row is zero by definition (dist to self)
  for(i in 2:n){
    parts <- strsplit(lines[i + 1], "\\s+")[[1]]
    vals <- as.numeric(parts[-1])  # skip filename
    dist_mat[i, 1:length(vals)] <- vals
  }
  
  # Mirror to upper triangle
  dist_mat[upper.tri(dist_mat)] <- t(dist_mat)[upper.tri(dist_mat)]
  
  # Convert to 'dist' object
  dist_obj <- as.dist(dist_mat)
  
  # Build NJ tree
  tree <- ape::nj(dist_obj)
  
  return(tree)
}

trees = list()
for(bug in ESKAPEE) {
  ids = unique(wide_all$assembly_ID[wide_all$eskapee == bug])
  idfile = paste0(bug, "-", expt, "-gca_ids.txt", collapse="")
  idfile = gsub(" ", "", idfile)
  sketchdir = paste0("sketches_", idfile, "/", collapse="")
  if(run.mash == TRUE) {
    write.table(ids, file=idfile,
                row.names = FALSE, quote = FALSE, col.names = FALSE)
    cmd.str = paste0("bash -lc './mash-batch-cli.sh ", idfile, "'", collapse="")
    system(cmd.str)
    cmd.str = paste0("mash triangle ", sketchdir, "*.msh > ", idfile, "-mash_distances.tab", collapse="")
    system(cmd.str)
  }
  trees[[bug]] = build_tree(paste0(idfile, "-mash_distances.tab", collapse=""))
}

# fit bug-specific models
final_df = wide_all

fit = fit.hmk2 = fit.hmm = fit.hmm.phy = tree.set = data.set = list()
for(bug in ESKAPEE) {
  df = final_df[final_df$eskapee == bug,c(2,4:ncol(final_df))]
  tree = ape::multi2di(trees[[bug]])
  df = df[match(tree$tip.label, df$assembly_ID), ]
  m = as.matrix(df[,2:ncol(df)])
  if(run.hmk2 == TRUE) {
    fit.hmk2[[bug]] = hyperinf(m, tree, method="hypermk2", reversible=TRUE)
  }
  fit.hmm[[bug]] = hyperinf(df, boot.parallel = 10)
  fit.hmm.phy[[bug]] = hyperinf(df, tree, boot.parallel = 10)
  data.set[[bug]] = df
  tree.set[[bug]] = tree
}

if(run.hmk2 == TRUE) {
  all.fits = list(fit.hmk2=fit.hmk2, fit.hmm=fit.hmm, fit.hmm.phy=fit.hmm.phy, data.set=data.set, tree.set=tree.set)
} else {
  all.fits = list(fit.hmm=fit.hmm, fit.hmm.phy=fit.hmm.phy, data.set=data.set, tree.set=tree.set)
}

save(all.fits, file=fname)
