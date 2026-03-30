library(arrow)
library(dplyr)
library(tidyr)
library(hyperinf)
library(ggrepel)
library(ggpubr)

# see details
# https://huggingface.co/datasets/ayates/amr_portal/blob/main/README.md?utm_source=chatgpt.com
# important fields here: antibiotic_name, BioSample_ID, resistance_phenotype

system("wget https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2025-11/phenotype.parquet")
df = read_parquet("phenotype.parquet")
expt = 1

if(expt == 0) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus", 
              "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Enterococcus faecium",
              "Enterobacter")
  to.get = 4
}

# simpler set  
if(expt == 1) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae")
  to.get = 5
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
  group_by(BioSample_ID, eskapee, antibiotic_name) %>%
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

# fit bug-specific models
final_df = wide_all
fit = fit.phy = fit.rev = fit.phy.rev = list()

for(bug in ESKAPEE) {
  set.seed(1)
  this.df = final_df[final_df$eskapee == bug,c(1,3:ncol(final_df))]
  this.df = this.df[sample(1:nrow(this.df), 50),]
  fit[[bug]] = hyperinf(this.df, boot.parallel = 10)
  fit.phy[[bug]] = hyperinf(this.df, auto.cluster = TRUE, boot.parallel = 10)
  fit.rev[[bug]] = hyperinf(this.df, reversible=TRUE)
  fit.phy.rev[[bug]] = hyperinf(this.df, auto.cluster = TRUE, reversible=TRUE)
}

drug.labels = colnames(final_df[3:ncol(final_df)])
drug.labels
plot.drugs = ggtexttable(data.frame(Drug=drug.labels))

################

for(this.bug in 1:2) {
bug = ESKAPEE[this.bug]
set.seed(1)
this.df = final_df[final_df$eskapee == bug,c(1,3:ncol(final_df))]
this.df = this.df[sample(1:nrow(this.df), 50),]
fit.list = list(fit[[bug]]$boots[[1]], fit.phy[[bug]]$boot[[1]],
                fit.rev[[bug]], fit.phy.rev[[bug]])
plot.0 = ggarrange(
plot_hyperinf_data(this.df),
plot_hyperinf_data(this.df, auto.cluster = TRUE, font.size = 2, bmargin = 50),
ncol=2)

plot.1 = plot_hyperinf_bubbles(fit.list, 
                      expt.names = c("CS", "Phy", "CS rev", "Phy rev"),
                          feature.names = substr(drug.labels, start=1, stop=3)) 

plot.2 = plot_hyperinf_ordering_matrices(fit.list, 
                                expt.names = c("CS", "Phy", "CS rev", "Phy rev"),
                                
                                feature.names = substr(drug.labels, start=1, stop=3)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
        
plot.3 = plot_hyperinf_comparative(fit.list, expt.names = c("CS", "Phy", "CS rev", "Phy rev"),
                          feature.names = substr(drug.labels, start=1, stop=3),
                          style="full") 

sf = 2
png(paste0("principle-", bug, "-", to.get, ".png", collapse=""), width=800*sf, height=800*sf, res=72*sf)
print(ggarrange(plot.0, plot.1, plot.2, plot.3, nrow=2, ncol=2))
dev.off()
}
