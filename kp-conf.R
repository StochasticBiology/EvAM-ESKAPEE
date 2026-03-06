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
kp.df = df[df$species=="Klebsiella pneumoniae",]

df_bin <- kp.df %>%
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
to.get = 10
drugs = appears$drug[1:to.get]
drugs = appears$drug

# INH (isoniazid); RIF (rifampicin, rifampin in the United States); PZA (pyrazinamide); EMB (ethambutol); STR (streptomycin); AMI (amikacin); CAP (capreomycin); MOX (moxifloxacin); OFL (ofloxacin); and PRO (prothionamide).

orig.drugs = c("isoniazid", "rifampin", "pyrazinamide", "ethambutol", "streptomycin", "amikacin", "capreomycin", "moxifloxacin", "ofloxacin", "prothionamide")

#drugs = orig.drugs

# pull the resistance profiles for these drugs into a wide dataframe 
wide_all_unc <- tmp[tmp$antibiotic_name %in% drugs,] %>%
  pivot_wider(names_from = antibiotic_name,
              values_from = value) 

wide_all_unc <- wide_all_unc %>%
  dplyr::select(all_of(c("BioSample_ID", "eskapee")), all_of(drugs))

wide_all <- wide_all_unc %>%
  drop_na()

library(dplyr)

amr_map <- list(
  AGly = c("gentamicin","amikacin","tobramycin"),
  Col  = c("colistin"),
  Flq  = c("ciprofloxacin","levofloxacin"),
  MLS  = c("azithromycin"),
  Phe  = c("chloramphenicol"),
  Sul  = c("trimethoprim-sulfamethoxazole"),
  Tmt  = c("trimethoprim","trimethoprim-sulfamethoxazole"),
  Tet  = c("tetracycline"),
  Tgc  = c("tigecycline"),
  Bla_a = c("ampicillin","cefazolin","cefuroxime"),
  Bla_ESBL = c("ceftriaxone","cefotaxime","ceftazidime","cefepime","aztreonam","cefoxitin"),
  Bla_inhR = c("piperacillin-tazobactam","ampicillin-sulbactam","amoxicillin-clavulanic acid"),
  Bla_ESBL_inhR = c("ceftazidime-avibactam"),
  Bla_Carb = c("meropenem","imipenem","ertapenem")
)

amr_map <- list(
  AGly = c("gentamicin","amikacin","tobramycin"),
  
  Flq  = c("ciprofloxacin","levofloxacin"),
  
  Sul  = c("trimethoprim-sulfamethoxazole"),
  Tmt  = c("trimethoprim-sulfamethoxazole"),
  
  Bla_a = c("ampicillin","cefazolin","cefuroxime"),
  
  Bla_ESBL = c(
    "ceftriaxone",
    "cefotaxime",
    "ceftazidime",
    "cefepime",
    "aztreonam",
    "cefoxitin"
  ),
  
  Bla_inhR = c(
    "piperacillin-tazobactam",
    "ampicillin-sulbactam"
  ),
  
  Bla_Carb = c(
    "meropenem",
    "imipenem",
    "ertapenem"
  )
)

amr_map <- list(
  
  AGly = c(
    "gentamicin","amikacin","tobramycin",
    "kanamycin","streptomycin","plazomicin",
    "netilmicin","neomycin"
  ),
  
  Col = c(
    "colistin","polymyxin B"
  ),
  
  Fcyn = c(
    "fosfomycin"
  ),
  
  Flq = c(
    "ciprofloxacin","levofloxacin","norfloxacin",
    "nalidixic acid","moxifloxacin","ofloxacin",
    "delafloxacin"
  ),
  
  Gly = c(
    "vancomycin","teicoplanin"
  ),
  
  MLS = c(
    "azithromycin","erythromycin","clarithromycin",
    "clindamycin"
  ),
  
  Phe = c(
    "chloramphenicol"
  ),
  
  Rif = c(
    "rifampin"
  ),
  
  Sul = c(
    "trimethoprim-sulfamethoxazole",
    "sulfamethoxazole",
    "sulfisoxazole",
    "trimethoprim-sulfobactam"
  ),
  
  Tet = c(
    "tetracycline","doxycycline","minocycline"
  ),
  
  Tgc = c(
    "tigecycline","eravacycline","omadacycline"
  ),
  
  Tmt = c(
    "trimethoprim",
    "trimethoprim-sulfamethoxazole",
    "trimethoprim-sulfobactam"
  ),
  
  Bla_a = c(
    "ampicillin","amoxicillin","penicillin",
    "carbenicillin","piperacillin",
    "cefazolin","cefuroxime","cephalothin",
    "cephalexin","oxacillin"
  ),
  
  Bla_inhR = c(
    "piperacillin-tazobactam",
    "ampicillin-sulbactam",
    "amoxicillin-clavulanic acid",
    "ticarcillin-clavulanic acid",
    "cefpodoxime-clavulanic acid",
    "ceftazidime-clavulanic acid",
    "cefotaxime-clavulanic acid"
  ),
  
  Bla_ESBL = c(
    "ceftriaxone","cefotaxime","ceftazidime",
    "cefepime","aztreonam","cefoxitin",
    "cefpodoxime","cefotetan","ceftiofur",
    "cefatrizine","cefixime","ceftaroline"
  ),
  
  Bla_ESBL_inhR = c(
    "ceftazidime-avibactam",
    "ceftolozane-tazobactam"
  ),
  
  Bla_Carb = c(
    "meropenem","imipenem","ertapenem","doripenem",
    "meropenem-vaborbactam","imipenem-relebactam"
  ),
  
  SHV = c(),
  
  Bla_chr = c(),
  
  Omp = c()
  
)

wide_all_0 = wide_all_unc
wide_all_0[is.na(wide_all_0)] = 0
tmp = matrix(0, nrow=nrow(wide_all_0), ncol=1+length(names(amr_map)))
colnames(tmp) = c("BioSample_ID", names(amr_map))
tmp = as.data.frame(tmp)
for(i in 1:nrow(wide_all_0)) {
  tmp[i,1] = wide_all_0$BioSample_ID[i]
  for(j in 1:length(names(amr_map))) {
    this.score = sum(wide_all_0[i,amr_map[[j]]])
    if(this.score == 0) { tmp[i,j+1] = 0} else {tmp[i,j+1] = 1}
  }
}
amr.new = as.data.frame(tmp)

colnames(amr.new)
amr.test = amr.new[,c("BioSample_ID", "AGly", "Col", "Fcyn", "Bla_ESBL", "Sul", "Rif", "Flq", "Tgc", "Tet", "MLS", "Gly")]
fit = hyperinf(amr.test)
ggarrange(plot_hyperinf(fit), 
          ggtexttable(data.frame(Drug=colnames(amr.test)[2:(ncol(amr.test)-1)])), 
          widths=c(3,1))

set.seed(1)
amr.small = amr.new[sample(1:nrow(amr.new),20), c("BioSample_ID", "AGly", "Bla_ESBL", "Sul","Flq", "Tet", "Gly")]
fit.rev = hyperinf(amr.small, reversible=TRUE)
fit.nonrev = hyperinf(amr.small)

ggarrange(plot_hyperinf(fit.rev, threshold = 0.0005), plot_hyperinf(fit.nonrev), 
          ggtexttable(data.frame(Drug=colnames(amr.small)[2:(ncol(amr.small))])),
          labels=c("A", "B", ""), nrow=1, widths=c(2,2,1))
refs = amr.test$BioSample_ID
kp.sub = kp.df[kp.df$BioSample_ID %in% refs,]
countries = unique(kp.sub[,c("BioSample_ID", "country")]) %>% drop_na()
amr.test$country = ""
for(i in 1:nrow(amr.test)) {
  ref = which(countries$BioSample_ID==amr.test$BioSample_ID[i])
  if(length(ref) > 1) { cat(i,"\n")}
  if(length(ref) == 0) {amr.test$country[i] = "na"}
  else {
  amr.test$country[i] = countries$country[ref][1] }
}
library(hyperhmm)
res.set = plot.set = bubble.set= list()
for(country in unique(amr.test$country)) {
  this.sub = amr.test[amr.test$country==country,1:(ncol(amr.test)-1)]
  res.set[[country]] = hyperinf(this.sub)
  plot.set[[country]] = plot_hyperinf(res.set[[country]])
  bubble.set[[country]] = plot_bubbles(res.set[[country]], formatted=TRUE)
}
plot_hyperinf_comparative(res.set[1:10], threshold=0.05) 

sf = 2
png("all-kp-estimates.png", width=3000*sf, height=3000*sf, res=72*sf)
ggarrange(plotlist=bubble.set, labels=unique(amr.test$country))
dev.off()

ggarrange(plotlist=bubble.set, labels=unique(amr.test$country))
arr = plot_hyperinf(res.set[[1]], "native")
plot_bubbles(res.set[[1]], formatted=TRUE)

res.set[[1]]$stats
