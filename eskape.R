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
expt = 2

if(expt == 0) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus", 
              "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Enterococcus faecium",
              "Enterobacter")
  to.get = 4
}

# simpler set  
if(expt == 1) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae")
  to.get = 9
}

if(expt == 2) {
  ESKAPEE = c("Escherichia coli", "Klebsiella pneumoniae", "Acinetobacter baumannii")
  to.get = 9
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
              values_from = value) %>%
  drop_na()

# output specifics
cat(nrow(wide_all))
unique(wide_all$eskapee)
table(wide_all$eskapee)

# fit bug-specific models
final_df = wide_all
fit = fit.plot = list()
for(bug in ESKAPEE) {
  this.df = final_df[final_df$eskapee == bug,c(1,3:ncol(final_df))]
  fit[[bug]] = hyperinf(this.df, boot.parallel = 10)
  fit.plot[[bug]] = plot_hyperinf(fit[[bug]]) + theme(legend.position="none")
}
drug.labels = colnames(final_df[3:ncol(final_df)])
drug.labels
plot.drugs = ggtexttable(data.frame(Drug=drug.labels))
fit.plot[[length(fit.plot)+1]] = plot.drugs

png(paste0("eskapee-expt-", expt, "-networks.png", collapse=""),
    width=800*sf, height=800*sf, res=72*sf)
print(
  ggarrange(plotlist=fit.plot, labels=ESKAPEE)
)
dev.off()


#plot_hyperinf_data(this.df[sample(1:nrow(this.df), 100),])

# pull sets of interactions implied from fitting hypertraps model
# not fully tested! check that fitted model has right dynamics
int.df = data.frame()
for(bug in ESKAPEE) {
  l2rep = full_to_squared(fit[[bug]])
  # pull interactions from a reconstructed squared fit
  
  nf = sqrt(ncol(l2rep))
  for(i in 1:nf) {
    for(j in 1:nf) {
      if(i != j) {
        x = l2rep[,(i-1)*nf+(j-1)+1]
        if(mean(x) > 0 & length(which(x < 0)) < length(x)*0.05 & sd(x)/mean(x) < 0.3) {
          int.df = rbind(int.df, data.frame(bug=bug, i=i, j=j,
                                            mean=mean(x),
                                            sd=sd(x)))
        }
        if(mean(x) < 0 & length(which(x > 0)) < length(x)*0.05 & sd(x)/mean(x) < 0.3) {
          int.df = rbind(int.df, data.frame(bug=bug, i=i, j=j,
                                            mean=mean(x),
                                            sd=sd(x)))
        }
      }
    }
  }
}
int.df

##### here are some ways we could regularise these interaction models

# start penalised HyperTraPS with the HyperHMM initial params
trial = l2rep[nrow(l2rep),]
this.m = as.matrix(final_df[final_df$eskapee == bug,3:ncol(final_df)])
ht.fit = HyperTraPS(this.m, initialparams=trial, 
                    length=4, kernel=3, penalty=1)
ht.fit$featurenames = drug.labels
plot_hyperinf(ht.fit)
plotHypercube.lik.trace(ht.fit)
plotHypercube.influencegraph(ht.fit, cv.thresh = 0.5)

## can we -- awkwardly -- pass exactly these params to HyperTraPS and regularise? 
trial = l2rep[nrow(l2rep),]
priors = cbind(matrix(trial, ncol=1), matrix(trial+1e-3, ncol=1))
this.m = as.matrix(final_df[final_df$eskapee == bug,3:ncol(final_df)])
try.fit.0 = HyperTraPS(this.m, priors=priors, length=1)
try.fit = HyperTraPS(this.m, priors=priors, length=1, regularise=1)
plotHypercube.regularisation(try.fit)

ggarrange(fit.plot[[length(fit.plot)-1]],
plot_hyperinf(try.fit.0),
plot_hyperinf(try.fit))

###########

# don't trust these interactions -- overfitting likely from HyperHMM

ggplot(int.df[abs(int.df$mean) > 0.5,], aes(x=i, y=j, fill=mean)) + 
  geom_tile() + geom_abline() + facet_wrap(~ bug)

ggplot(int.df, aes(x=i, y=j, color=factor(bug))) + 
  geom_point(alpha=0.5, size=20) + geom_abline() + 
  scale_y_continuous(breaks = 1:length(drug.labels), labels=drug.labels) +
  scale_x_continuous(breaks = 1:length(drug.labels), labels=drug.labels) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

int.df$buglabel = sapply(int.df$bug, function (x) {which(ESKAPEE==x)})

png(paste0("eskapee-expt-", expt, "-interactions.png", collapse=""),
    width=600*sf, height=600*sf, res=72*sf)
print(
  ggplot(int.df, aes(x=i, y=j, label=buglabel)) + 
    geom_text_repel(max.overlaps=30) + geom_abline() + 
    scale_y_continuous(breaks = 1:length(drug.labels), labels=drug.labels) +
    scale_x_continuous(breaks = 1:length(drug.labels), labels=drug.labels) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
)
dev.off()

fit.ht = plot.ht = plot.comp = list()
for(i in 1:length(ESKAPEE)) {
  fit.ht[[i]] = hypertraps_from_params(full_to_squared(fit[[i]]))
  plot.hct[[i]] = plot_hyperinf(fit.ht[[i]])
  plot.comp[[i]] = ggarrange(plot_hyperinf(fit[[i]]), plot.hct[[i]])
}
plot.comp[[length(plot.comp)+1]] = plot.drugs

png(paste0("eskapee-expt-", expt, "-check-plot.png", collapse=""),
    width=2000*sf, height=1000*sf, res=72*sf)
print(ggarrange(plotlist = plot.comp))
dev.off()
###???? problem

# look just at prevalences
drug.labels[c(2,5,8,6,3,7,4,1)]
prev.df = data.frame()
for(i in 1:length(ESKAPEE)) {
  prev.df = rbind(prev.df, 
                  data.frame(bug=ESKAPEE[i],
                             drug=colnames(final_df[,3:ncol(final_df)]),
                             prevalence=colMeans(final_df[final_df$eskapee==ESKAPEE[i],3:ncol(final_df)]))
  )
}
ggplot(prev.df, aes(x=drug, y=prevalence, fill=bug)) + geom_col(position="dodge")
