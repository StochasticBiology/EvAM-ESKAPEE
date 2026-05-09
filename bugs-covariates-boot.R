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

bug.name = top.bugs[4]
expt = 1

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
nboots = 10
for(this.covariate in unique(wide_all$covariate)) {
  this.sub = wide_all[wide_all$covariate==this.covariate,2:(ncol(wide_all))]
  this.sub$covariate = paste("l-", this.sub$covariate)
  if(nrow(this.sub) > 1) {
  res.set[[this.covariate]] = hyperinf(this.sub, boot.parallel = nboots)
  }
  #plot.set[[this.covariate]] = plot_hyperinf(res.set[[this.covariate]])
}

mean.fit = hyperinf(wide_all[2:ncol(wide_all)], boot.parallel = nboots)
res.set[["mean"]] = mean.fit

# get labels for drugs and covariates
drug.names = colnames(wide_all)[3:(ncol(wide_all))]
drug.abbrevs = substr(drug.names, start=1,stop=3)
covariate.names = names(res.set)

plot_hyperinf_comparative(res.set, threshold=0.1,
                          expt.names = covariate.names,
                          feature.names = drug.abbrevs) +theme(legend.position="none")

# produce the comparative plot
comp.plot = ggarrange(
  plot_hyperinf_comparative(res.set, threshold=0.1,
                            expt.names = covariate.names,
                            label_size=3,
                            feature.names = drug.abbrevs) +theme(legend.position="none")
  ,
  plot_hyperinf_bubbles(res.set, thetastep=3, p.scale=0.5, 
                        expt.names=rep(covariate.names),
                        sqrt.trans = TRUE,
                        fill.name="Covariate", 
                        feature.names = gsub("-", "-\n", drug.names)) + 
    xlab("Step") + scale_x_continuous(breaks = 1:10),
  widths=c(0.6,1)
)

comp.plot
compare_orderings(res.set[["Georgia"]], res.set[["Republic of Moldova"]], 
                  type="relative", threshold=0.4)

plot_hyperinf_compare_orderings(res.set[["Georgia"]], res.set[["Republic of Moldova"]])
compare_orderings(res.set[[]], res.set[[2]], type="absolute", threshold=0.2)

plot_hyperinf_bubbles(list(res.set[[1]], res.set[[2]]))
# output to file
sf = 2
png(paste0(bug.name, "-", covariate.col, "-big-boots.png", collapse=""), 
    width=700*sf, height=400*sf, res=72*sf)
print(comp.plot)
dev.off()

covariate.names
# at the moment, this function only works for comparing exactly two model fits
# it should show cases where bootstrap estimates for a given P_ij show complete separation
# across the two model fits
plot_hyperinf_bootstrap(res.set[[1]], res.set[[2]])
