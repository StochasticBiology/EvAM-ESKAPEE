library(arrow)
library(dplyr)
library(tidyr)
library(hyperinf)
library(ggrepel)
library(ggpubr)
library(phytools)

run.bash = TRUE
run.mash = TRUE
get.premade.tree = FALSE

# pip install ncbi-genome-download
# brew install mash
#system("wget https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2025-11/genotype.parquet")
#gdf = read_parquet("genotype.parquet")
pdf = read_parquet("phenotype.parquet")

ESKAPEE = c("Klebsiella pneumoniae")
ids = unique(pdf$assembly_ID[pdf$species %in% ESKAPEE])

if(get.premade.tree == FALSE) {
  write.table(ids, file="gca_ids.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
  if(run.bash == TRUE) {
    system("bash -lc './mash-batch.sh'")
  }
  
  if(run.mash == TRUE) {
    # Mash distances -- change directory if needed
    system("mash triangle sketches-hpylori/*.msh > mash_distances.tab")
  }
  lines <- readLines("mash_distances.tab")
  
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
  tree <- nj(dist_obj)

  # Save tree
  write.tree(tree, paste0(ESKAPEE, "-tree.nwk", collapse=""))
} else {
  tree = read.tree(paste0(ESKAPEE, "-tree.nwk", collapse=""))
}

#plot(tree)

########

df = pdf
df$eskapee = df$species
df$eskapee[grep("[Ee]nterobacter", df$species)] = "Enterobacter"

# 1. Clean + binarize + collapse duplicates
df_bin <- df[pdf$assembly_ID %in% ids,] %>%
  filter(resistance_phenotype %in% c("resistant", "susceptible"),
         eskapee %in% ESKAPEE) %>%
  mutate(value = ifelse(resistance_phenotype == "resistant", 1, 0)) %>%
  group_by(assembly_ID, eskapee, antibiotic_name) %>%
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

to.get = 8
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

# fit bug-specific models
final_df = wide_all

bug = ESKAPEE[1]
this.df = final_df[final_df$eskapee == bug,c(1,3:ncol(final_df))]
this.df$assembly_ID = as.character(this.df$assembly_ID)

plot.1 = plot_hyperinf_data(this.df, tree, bmargin = 80)
#plot.1

fit = hyperinf(this.df, boot.parallel = 10)
fit.est.phy = hyperinf(this.df, auto.cluster = TRUE, boot.parallel = 10)
fit.phy = hyperinf(this.df, tree, boot.parallel = 10)

expt.names = c("CS", "Est phy", "Genome phy")
fits = list(fit, fit.est.phy, fit.phy)

plot.2 = plot_hyperinf_comparative(fits, expt.names=expt.names)

plot.3 = plot_hyperinf_bubbles(fits, expt.names=expt.names)

plot.4 = plot_hyperinf_ordering_matrices(list(fit, fit.est.phy, fit.phy), expt.names=expt.names ) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

sf = 2
png(paste0(ESKAPEE, "-summaries.png", collapse=""), width=800*sf, height=800*sf, res=72*sf)
ggarrange(plot.1, plot.2, plot.3, plot.4)
dev.off()

