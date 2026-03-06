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
             "Streptococcus pneumoniae",
             "Staphylococcus aureus",
             "Neisseria gonorrhoeae",
             "Klebsiella pneumoniae",
             "Salmonella enterica",
             "Escherichia coli",
             "Pseudomonas aeruginosa")

bug.name = top.bugs[1]

bug.df = df[df$species==bug.name,]

df_bin <- bug.df %>%
  filter(resistance_phenotype %in% c("resistant", "susceptible")) %>%
  mutate(value = ifelse(resistance_phenotype == "resistant", 1, 0)) %>%
  group_by(BioSample_ID, eskapee, country, antibiotic_name) %>%
  summarise(value = max(value), .groups = "drop")

tmp = df_bin
counts = tmp %>% count(eskapee, antibiotic_name) 
drugs = c()
ggplot(counts, aes(x=eskapee, y=antibiotic_name, fill=n)) + geom_tile()

appears = data.frame()
for(drug in unique(df_bin$antibiotic_name)) {
  appears = rbind(appears, data.frame(drug=drug,
                                      species=length(which(counts$antibiotic_name==drug)),
                                      min = min(counts$n[counts$antibiotic_name==drug]),
                                      sum = sum(counts$n[counts$antibiotic_name==drug]))
  )
}
appears = appears[order(-appears$min),]

to.get = 12
n.countries = 0

while(n.countries < 5) {
to.get = to.get-1
drugs = appears$drug[1:to.get]

# pull the resistance profiles for these drugs into a wide dataframe 
wide_all_unc <- tmp[tmp$antibiotic_name %in% drugs,] %>%
  pivot_wider(names_from = antibiotic_name,
              values_from = value) 

wide_all_unc <- wide_all_unc %>%
  dplyr::select(all_of(c("BioSample_ID", "country")), all_of(drugs))

wide_all <- wide_all_unc %>%
  drop_na()

t.counts = table(wide_all$country) 
good.countries = names(t.counts)[which(as.vector(t.counts)>5)]
n.countries = length(good.countries)
}

wide_all = wide_all[wide_all$country %in% good.countries,]

fit = hyperinf(wide_all[2:ncol(wide_all)])
ggarrange(plot_hyperinf(fit), 
          ggtexttable(data.frame(Drug=colnames(wide_all)[3:ncol(wide_all)])), 
          widths=c(3,1))

refs = wide_all$BioSample_ID
bug.sub = bug.df[bug.df$BioSample_ID %in% refs,]
countries = unique(bug.sub[,c("BioSample_ID", "country")]) %>% drop_na()
wide_all$country = ""
for(i in 1:nrow(wide_all)) {
  ref = which(countries$BioSample_ID==wide_all$BioSample_ID[i])
  if(length(ref) > 1) { cat(i,"\n")}
  if(length(ref) == 0) {wide_all$country[i] = "na"}
  else {
    wide_all$country[i] = countries$country[ref] }
}
res.set = plot.set = bubble.set= list()
for(country in unique(wide_all$country)) {
  this.sub = wide_all[wide_all$country==country,2:(ncol(wide_all))]
  res.set[[country]] = hyperinf(this.sub)
  plot.set[[country]] = plot_hyperinf(res.set[[country]])
  bubble.set[[country]] = plot_bubbles(res.set[[country]], formatted=TRUE)
}

drug.names = colnames(wide_all)[3:(ncol(wide_all))]
drug.abbrevs = substr(drug.names, start=1,stop=3)
country.names = unique(wide_all$country)

comp.plot = ggarrange(
  plot_hyperinf_comparative(res.set, threshold=0.1,
                            expt.names = country.names,
                            feature.names = drug.abbrevs) +theme(legend.position="none")
  ,
  plot_hyperinf_bubbles(res.set, thetastep=3, p.scale=0.5, 
                        expt.names=country.names,
                        sqrt.trans = TRUE,
                        fill.name="Country",
                        feature.names = drug.abbrevs),
  widths=c(0.5,1)
)

sf = 2
png(paste0(bug.name, ".png", collapse=""), width=1000*sf, height=500*sf, res=72*sf)
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
}