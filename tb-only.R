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

tb.df = df[df$species=="Mycobacterium tuberculosis",]

df_bin <- tb.df %>%
  filter(resistance_phenotype %in% c("resistant", "susceptible")) %>%
  mutate(value = ifelse(resistance_phenotype == "resistant", 1, 0)) %>%
  group_by(BioSample_ID, eskapee, antibiotic_name) %>%
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
to.get = 8
drugs = appears$drug[1:to.get]

# INH (isoniazid); RIF (rifampicin, rifampin in the United States); PZA (pyrazinamide); EMB (ethambutol); STR (streptomycin); AMI (amikacin); CAP (capreomycin); MOX (moxifloxacin); OFL (ofloxacin); and PRO (prothionamide).

#orig.drugs = c("isoniazid", "rifampin", "pyrazinamide", "ethambutol", "streptomycin", "amikacin", "capreomycin", "moxifloxacin", "ofloxacin", "prothionamide")

#drugs = unique(c(drugs, orig.drugs))
#drugs = orig.drugs

# pull the resistance profiles for these drugs into a wide dataframe 
wide_all_unc <- tmp[tmp$antibiotic_name %in% drugs,] %>%
  pivot_wider(names_from = antibiotic_name,
              values_from = value) 

wide_all_unc <- wide_all_unc %>%
  dplyr::select(all_of(c("BioSample_ID", "eskapee")), all_of(drugs))

wide_all <- wide_all_unc %>%
  drop_na()

fit = hyperinf(wide_all[2:ncol(wide_all)])
ggarrange(plot_hyperinf(fit), 
          ggtexttable(data.frame(Drug=colnames(wide_all)[3:ncol(wide_all)])), 
          widths=c(3,1))

refs = wide_all$BioSample_ID
tb.sub = tb.df[tb.df$BioSample_ID %in% refs,]
countries = unique(tb.sub[,c("BioSample_ID", "country")]) %>% drop_na()
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
  this.sub = wide_all[wide_all$country==country,2:(ncol(wide_all)-1)]
  res.set[[country]] = hyperinf(this.sub)
  plot.set[[country]] = plot_hyperinf(res.set[[country]])
 # bubble.set[[country]] = plot_bubbles(res.set[[country]], formatted=TRUE)
}

drug.names = colnames(amr.test)[3:(ncol(amr.test)-1)]
drug.abbrevs = substr(drug.names, start=1,stop=3)
country.names = unique(wide_all$country)

ggarrange(
  plot_hyperinf_comparative(res.set, threshold=0.1,
                            expt.names = country.names,
                            feature.names = drug.abbrevs) +theme(legend.position="none")
  ,
  plot_hyperinf_bubbles(res.set, thetastep=3, p.scale=0.5, 
                        expt.names=unique(wide_all$country),
                        fill.name="Country",
                        feature.names = drug.abbrevs),
  widths=c(0.5,1)
)

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