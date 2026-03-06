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
df$eskapee = df$species

top.bugs = c("Mycobacterium tuberculosis",
             "Staphylococcus aureus",
             "Klebsiella pneumoniae",
             "Escherichia coli",
             "Pseudomonas aeruginosa")

bug.name = top.bugs[3]
expt = 2

# choose what covariate we want to explore
# for each, we place requirements on the minimum number of distinct levels
# and the minimal number of samples for each of these levels
# and we find as many drugs as we can that allow us complete records for that combination
if(expt == 1) {
  covariate.col = "host_sex"
  reqd.covariate.levels = 2
  reqd.samples = 20
} else if(expt == 2) {
  covariate.col = "country"
  reqd.covariate.levels = 5
  reqd.samples = 5
} else if(expt == 3) {
  covariate.col = "host"
  reqd.covariate.levels = 2
  reqd.samples = 20
} else if(expt == 4) {
  covariate.col = "collection_year"
  reqd.covariate.levels = 5
  reqd.samples = 5
} else if(expt == 5) {
  df$latitude_band = 10*round(as.numeric(df$isolation_latitude)/10)
  covariate.col = "latitude_band"
  reqd.covariate.levels = 3
  reqd.samples = 5
} else if(expt == 6) {
  covariate.col = "geographical_region"
  reqd.covariate.levels = 5
  reqd.samples = 5
} else if(expt == 7) {
  df$age_band = 10*round(as.numeric(df$host_age)/10)
  covariate.col = "age_band"
  reqd.covariate.levels = 3
  reqd.samples = 5
}

# subset our bug of interest
bug.df = as.data.frame(df[df$species==bug.name,])

# pull covariate information into a given column
bug.df$covariate = as.character(bug.df[,c(covariate.col)])

# reframe resistance scores
df_bin <- bug.df %>%
  filter(resistance_phenotype %in% c("resistant", "susceptible")) %>%
  mutate(value = ifelse(resistance_phenotype == "resistant", 1, 0)) %>%
  group_by(BioSample_ID, eskapee, covariate, antibiotic_name) %>%
  summarise(value = max(value), .groups = "drop")

# count how many drug entries we have
tmp = df_bin
counts = tmp %>% count(eskapee, antibiotic_name) 
drugs = c()
ggplot(counts, aes(x=eskapee, y=antibiotic_name, fill=n)) + geom_tile()

# count appearances for each drug and produce an ordered list
appears = data.frame()
for(drug in unique(df_bin$antibiotic_name)) {
  appears = rbind(appears, data.frame(drug=drug,
                                      species=length(which(counts$antibiotic_name==drug)),
                                      min = min(counts$n[counts$antibiotic_name==drug]),
                                      sum = sum(counts$n[counts$antibiotic_name==drug]))
  )
}
appears = appears[order(-appears$min),]

# start our filtering to try and get required levels with required samples
to.get = 12
n.covariate.levels = 0
# loop, decreasing number of drugs each time
while(n.covariate.levels < reqd.covariate.levels) {
  to.get = to.get-1
  if(to.get == 0) {
    message("We don't have samples meeting these requirements!")
    stop()
  }
  drugs = appears$drug[1:to.get]
  
  # pull the resistance profiles for these drugs into a wide dataframe 
  wide_all_unc <- tmp[tmp$antibiotic_name %in% drugs,] %>%
    pivot_wider(names_from = antibiotic_name,
                values_from = value) 
  
  wide_all_unc <- wide_all_unc %>%
    dplyr::select(all_of(c("BioSample_ID", "covariate")), all_of(drugs))
  
  wide_all <- wide_all_unc %>%
    drop_na()
  
  # count the samples for each covariate level
  t.counts = table(wide_all[,c("covariate")]) 
  good.covariates = names(t.counts)[which(as.vector(t.counts)>reqd.samples)]
  n.covariate.levels = length(good.covariates)
}

# subset out just those covariate levels we keep
wide_all = wide_all[wide_all$covariate %in% good.covariates,]

# fit the aggregated data
fit = hyperinf(wide_all[2:ncol(wide_all)])
ggarrange(plot_hyperinf(fit), 
          ggtexttable(data.frame(Drug=colnames(wide_all)[3:ncol(wide_all)])), 
          widths=c(3,1))

# produce a list of model fits for each covariate level
res.set = plot.set = list()
for(this.covariate in unique(wide_all$covariate)) {
  this.sub = wide_all[wide_all$covariate==this.covariate,2:(ncol(wide_all))]
  this.sub$covariate = paste("l-", this.sub$covariate)
  res.set[[this.covariate]] = hyperinf(this.sub, boot.parallel = 10)
  #plot.set[[this.covariate]] = plot_hyperinf(res.set[[this.covariate]])
}

# get labels for drugs and covariates
drug.names = colnames(wide_all)[3:(ncol(wide_all))]
drug.abbrevs = substr(drug.names, start=1,stop=3)
covariate.names = unique(wide_all$covariate)

boots = list()
for(i in 1:length(res.set)) {
  boots = c(boots, res.set[[i]]$boots)
}

plot_hyperinf_comparative(res.set, threshold=0.1,
                          expt.names = covariate.names,
                          feature.names = drug.abbrevs) +theme(legend.position="none")

# produce the comparative plot
comp.plot = ggarrange(
  plot_hyperinf_comparative(res.set, threshold=0.1,
                            expt.names = covariate.names,
                            feature.names = drug.abbrevs) +theme(legend.position="none")
  ,
  plot_hyperinf_bubbles(boots, thetastep=3, p.scale=0.5, 
                        expt.names=rep(covariate.names, each=11),
                        sqrt.trans = TRUE,
                        fill.name="Covariate",
                        feature.names = drug.abbrevs),
  widths=c(1,1)
)

# output to file
sf = 2
png(paste0(bug.name, "-", covariate.col, "-boots.png", collapse=""), width=1400*sf, height=500*sf, res=72*sf)
print(comp.plot)
dev.off()

if(FALSE) {
  ### non-TB-related tests
  
  data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
  fit.1 = hyperinf(data)
  fit.2 = hyperinf(data, reversible=TRUE)
  fit.3 = hyperinf(data, method="hypertraps")
  fit.4 = hyperinf(data, method="hyperlau")
  
  plot_hyperinf_comparative(list(fit.1, fit.2, fit.3, fit.4))
  
  
  data = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
  fit.1 = hyperinf(data)
  fit.2 = hyperinf(1-data)
  fit.3 = hyperinf(data, reversible=TRUE)
  
  plot_hyperinf_comparative(list(fit.1, fit.2, fit.3))
  
  data = matrix(c(0,0,1, 0,1,1, 1,1,1, 1,0,0, 0,1,1), ncol=3, nrow=5)
  fit.1 = hyperinf(data, boot.parallel = 10)
  fit.2 = hyperinf(1-data, boot.parallel = 10)
  
  plot_hyperinf_comparative(list(fit.1, fit.2))
  
}