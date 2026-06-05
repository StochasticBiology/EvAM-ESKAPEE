library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexUpset)
library(ggplot2)
library(igraph)
library(ggpubr)
library(hyperinf)
library(arrow)

run.inference = TRUE
sf = 2

# NB: this list of drugs to ignore is derived from TASK 2 - TASK 1 comparison
prune = c("amoxicillin-clavulanic acid", "cefoxitin", "doxycycline",
          "cefotaxime", "chloramphenicol", "tigecycline", "aztreonam", "fosfomycin")

####### TASK 1 -- process new phenotypes
####### TASK 2 -- fits from CABBAGE data
####### TASK 3 -- try and do predictions
####### TASK 4 -- embed new observations in previous model vis

##### helpers

BinToDec <- function(state) {
  this.ref = 0
  for(j in 1:length(state)) {
    this.ref = this.ref + state[j]*(2**(length(state)-j))
  }
  return(this.ref)
}
DecToBinS <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste0(s, collapse=""))
}
DecToLevel <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(sum(s))
}

####### TASK 1 -- process new phenotypes

df = read_excel("tanzania-phenotypes.xlsx")

length(unique(df$`study ID number`))
nrow(df)
df[,13:ncol(df)]

rdf <- df %>%
  mutate(across(13:ncol(.), ~ replace(., . == 1, "resistant"))) %>%
  mutate(across(13:ncol(.), ~ replace(., . == 0, "sensitive")))

new_df <- rdf %>%
  pivot_longer(
    cols = 13:ncol(.),
    names_to = "drug",
    values_to = "drug_outcome"
  ) %>%
  rename(ID = `study ID number`) %>%
  filter(drug_outcome %in% c("sensitive", "resistant"))

new_df$species = NA
new_df$species[grep("coli", new_df$`name of bacteria`)] = "Ec"
new_df$species[grep("neum", new_df$`name of bacteria`)] = "Kp"

small_df = new_df[!is.na(new_df$species),c("ID", "admission date", "species", "drug", "drug_outcome")]

slashes = grep("/", small_df$`admission date`)
dots = grep("[.]", small_df$`admission date`)
small_df$date = small_df$`admission date`
small_df$date[dots] = sapply(small_df$`admission date`[dots], as.Date, format = "%d.%m.%y")-as.numeric(as.Date("1899-12-30"))
small_df$date[slashes] = sapply(small_df$`admission date`[slashes], as.Date, format = "%d/%m/%y")-as.numeric(as.Date("1899-12-30"))
small_df$date = as.numeric(small_df$date)

unique(small_df$drug)

clean_names <- small_df$drug %>%
  str_to_lower() %>%
  str_replace_all("ecoli|klebpneumo|kleboxytoka|kleb|oxytoka", "") %>%
  str_replace_all("[^a-z]", "")   # remove punctuation, dots, spaces

drug_map <- c(
  "cipro" = "ciprofloxacin",
  "cip"   = "ciprofloxacin",
  
  "doxy"  = "doxycycline",
  
  "caz"   = "ceftazidime",
  "ctx"   = "cefotaxime",
  
  "mem"   = "meropenem",
  "mero"  = "meropenem",
  
  "sxt"   = "trimethoprim-sulfamethoxazole",
  
  "atm"   = "aztreonam",
  
  "fep"   = "ceftazidime",   # XXX borderline (cefepime vs ceftazidime!)
  "fox"   = "cefoxitin",
  
  "fosfo" = "fosfomycin",         # keep as-is (no canonical version present)
  "fosco" = "fosfomycin", # XXX assuming typo??? (cefepime vs ceftazidime!)
  
  "genta" = "gentamicin",          # already canonical
  
  "tige" = "tigecycline",
  
  "bactrim" = "trimethoprim-sulfamethoxazole",
  
  "tzp" = "piperacillin-tazobactam",
  
  "augmentin" = "amoxicillin-clavulanic acid" 
)

small_df$drug = clean_names

for(i in 1:nrow(small_df)) {
  if(small_df$drug[i] %in% names(drug_map)) {
    small_df$drug[i] = drug_map[small_df$drug[i]]
  }
  if(small_df$drug[i] %in% c("c", "", "atmiae", "esbl", "feppmeumo", "ptt")) {
    small_df$drug[i] = NA
  }
}
table(small_df$drug)


clean_df = small_df[!is.na(small_df$drug),]
unique(clean_df$drug)

sens_mat <- clean_df %>%
  filter(drug != "") %>%
  filter(drug_outcome == "resistant") %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = drug,
    values_from = value,
    values_fn = max,   # collapse duplicates safely
    values_fill = 0
  )

kp.upset = upset(
  sens_mat[sens_mat$species=="Kp",],
  intersect = setdiff(colnames(sens_mat), c("species", "date")),  # all drug columns
  n_intersections = 20
)
ec.upset = upset(
  sens_mat[sens_mat$species=="Ec",],
  intersect = setdiff(colnames(sens_mat), c("species", "date")),  # all drug columns
  n_intersections = 20
)

ggarrange(kp.upset, ec.upset, nrow=2)


kp.upset.1 = upset(
  sens_mat[sens_mat$species=="Kp" & sens_mat$date < 41500,],
  intersect = setdiff(colnames(sens_mat), c("species", "date")),  # all drug columns
  n_intersections = 20
)
kp.upset.2 = upset(
  sens_mat[sens_mat$species=="Kp" & sens_mat$date >= 41500,],
  intersect = setdiff(colnames(sens_mat), c("species", "date")),  # all drug columns
  n_intersections = 20
)

ec.upset.1 = upset(
  sens_mat[sens_mat$species=="Ec" & sens_mat$date < 41500,],
  intersect = setdiff(colnames(sens_mat), c("species", "date")),  # all drug columns
  n_intersections = 20
)
ec.upset.2 = upset(
  sens_mat[sens_mat$species=="Ec" & sens_mat$date >= 41500,],
  intersect = setdiff(colnames(sens_mat), c("species", "date")),  # all drug columns
  n_intersections = 20
)

old.new.upset = ggarrange(kp.upset.1, ec.upset.1, 
                          kp.upset.2, ec.upset.2, ncol=2, nrow=2,
          labels = c("Kp old", "Ec old", "Kp new", "Ec new"))

png("old-new-upset.png", width=1000*sf, height=600*sf, res=72*sf)
print(old.new.upset)
dev.off()

fits = this.mat = list()
for(species in unique(sens_mat$species)) {
  this.mat[[species]] = as.matrix(sens_mat[sens_mat$species == species,6:ncol(sens_mat)])
  if(run.inference == TRUE) {
    fits[[species]] = hyperinf(this.mat[[species]], method="hyperhmm", boot.parallel=3)
  }
}

if(run.inference == TRUE) {
  save(fits, file="tanzania-fits.Rdata")
} else {
  load("tanzania-fits.Rdata")
}

plot.comp = plot_hyperinf_comparative(fits, 
                                      expt.names=unique(sens_mat$species),
                                      style="full",
                                      feature.names = substr(fits[[1]]$feature.names, 1, 3))
#plot.ec = plot_hyperinf(fits[[1]], feature.names = substr(fits[[1]]$feature.names, 1, 3))
#plot.kp = plot_hyperinf(fits[[2]], feature.names = substr(fits[[2]]$feature.names, 1, 3))

#ggarrange(plot.ec, plot.kp)

png("tanzania-hyperinf.png", width=1000*sf, height=800*sf, res=72*sf)
print(plot.comp)
dev.off()

####### TASK 2 -- fits from CABBAGE data

drug.names = fits[[1]]$feature.names

qdf = read_parquet("phenotype.parquet")

qdf = qdf[qdf$species %in% c("Klebsiella pneumoniae", "Escherichia coli") &
            #as.numeric(qdf$collection_year) <= 2015 &
            qdf$antibiotic_name %in% drug.names, 
          c("BioSample_ID", "species", "collection_year", "antibiotic_name", "resistance_phenotype")]

qsens_mat <- qdf %>%
  # filter(resistance_phenotype == "resistant") %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = antibiotic_name,
    values_from = value,
    values_fn = max,   # collapse duplicates safely
    values_fill = 0
  )

q.upset = upset(
  qsens_mat,
  intersect = colnames(qsens_mat)[5:ncol(qsens_mat)],  # all drug columns
  n_intersections = 20
)
q.upset

qdf = qdf[!(qdf$antibiotic_name %in% prune),]
table(qdf$antibiotic_name)

# get IDs with complete records for each drug
good_ids <- qdf %>%
  count(BioSample_ID) %>%
  filter(n == n_distinct(qdf$antibiotic_name)) %>%
  pull(BioSample_ID)

good_ids

qdf_filtered <- qdf %>%
  filter(BioSample_ID %in% good_ids)

qsens_mat <- qdf_filtered %>%
  filter(resistance_phenotype == "resistant") %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = antibiotic_name,
    values_from = value,
    values_fn = max,   # collapse duplicates safely
    values_fill = 0
  )

qfits = qthis.mat = list()
for(species in unique(qsens_mat$species)) {
  qthis.mat[[species]] = as.matrix(qsens_mat[qsens_mat$species == species,5:ncol(qsens_mat)])
  if(run.inference == TRUE) {
    qfits[[species]] = hyperinf(qthis.mat[[species]], method="hyperhmm", boot.parallel=3)
  }
}

# only Kp has complete records pre-2010
qthis.mat[["past"]] = as.matrix(qsens_mat[qsens_mat$species == "Klebsiella pneumoniae" &
                                            !is.na(qsens_mat$collection_year) &
                                            qsens_mat$collection_year < 2010,
                                          5:ncol(qsens_mat)])
if(run.inference == TRUE) {
  qfits[["past"]] = hyperinf(qthis.mat[["past"]], method="hyperhmm", boot.parallel=3)
}

if(run.inference == TRUE) {
  save(qfits, file = "qtanzania-fits.Rdata")
} else {
  load("qtanzania-fits.Rdata")
}

qplot.comp = plot_hyperinf_comparative(qfits[1:2], expt.names=c("Ec", "Kp"),
                                       feature.names = substr(qfits[[1]]$feature.names, 1, 3),
                                       style = "full")

qplot.comp

png("both-fits.png", width=800*sf, height=600*sf, res=72*sf)
ggarrange(plot.comp, qplot.comp)
dev.off()

####### TASK 3 -- try and do predictions

for(speciesref in 1:3) {
  if(speciesref == 1) { 
    species = "Ec"
  } else {
    species = "Kp"
  }
  
  cabbage.fit = qfits[[speciesref]]
  inf.trans = cabbage.fit$transitions
  inf.trans = inf.trans[inf.trans$p.boot==1,]
  
  drugs = cabbage.fit$feature.names
  old.data = as.matrix(sens_mat[sens_mat$species==species & sens_mat$date < 41500,
                                6:ncol(sens_mat)])
  old.data.r <- old.data[, match(drugs, colnames(old.data))]
  colnames(old.data.r) == drugs
  new.data = as.matrix(sens_mat[sens_mat$species==species & sens_mat$date >= 41500,
                                6:ncol(sens_mat)])
  new.data.r <- new.data[, match(drugs, colnames(new.data))]
  colnames(new.data.r) == drugs
  
  old.states.raw = apply(old.data.r, 1, BinToDec)
  new.states.obs = apply(new.data.r, 1, BinToDec)
  
  # different putative initial states
  #old.states = old.states.raw
  old.states = c(old.states.raw, rep(0, 50))
  #old.states = rep(0, 100)
  
  # assign probability weights to different patterns with numbers of steps from some putative initial state
  weights = matrix(0, nrow=5, ncol=2**length(drugs))
  for(old.state in unique(old.states)) {
    weights[1,old.state+1] = length(which(old.states==old.state))/length(old.states)
  }
  for(i in 2:5) {
    for(old.state in old.states) {
      for(new.state in inf.trans$To[inf.trans$From == old.state]) {
        weights[i, new.state+1] = weights[i, new.state+1] + 
          inf.trans$Probability[inf.trans$To == new.state & inf.trans$From == old.state]*
          weights[i-1, old.state+1]
      }    
    }
    old.states = which(weights[i,] != 0)-1
  }
  
  # produce df comparing predictions to observations across levels
  preds = colSums(weights)
  obs = rep(0, 2**length(drugs))
  for(new.state in unique(new.states.obs)) {
    obs[new.state + 1] = length(which(new.states.obs == new.state))/length(new.states.obs)
  }
  perf.df = data.frame(state = (1:2**length(drugs))-1, 
                       preds = preds,
                       obs = obs)
  perf.df$bin = sapply(perf.df$state, DecToBinS, len=length(drugs))
  perf.df$level = sapply(perf.df$state, DecToLevel, len=length(drugs))
  
}

####### TASK 4 -- embed new observations in previous model vis

embed.g = list()

for(speciesref in 1:3) {
  if(speciesref == 1) { 
    species = "Ec"
  } else {
    species = "Kp"
  }
  
  cabbage.fit = qfits[[speciesref]]
  
  old.data = as.matrix(sens_mat[sens_mat$species==species & sens_mat$date < 41500,
                                6:ncol(sens_mat)])
  old.data.r <- old.data[, match(drugs, colnames(old.data))]
  colnames(old.data.r) == drugs
  new.data = as.matrix(sens_mat[sens_mat$species==species & sens_mat$date >= 41500,
                                6:ncol(sens_mat)])
  new.data.r <- new.data[, match(drugs, colnames(new.data))]
  colnames(new.data.r) == drugs
  
  old.states.raw = apply(old.data.r, 1, BinToDec)
  new.states.obs = apply(new.data.r, 1, BinToDec)
  
  #pg = hyperinf::get_plot_graph(cabbage.fit, threshold = 5e-3)
  pg = hyperinf::get_plot_graph(cabbage.fit, threshold = 1e-10)
  pgg = pg$plot.graph
  E(pgg)$label3 = substr(E(pgg)$label, 2, 4)
  E(pgg)$label3[E(pgg)$Flux < 0.1] = ""
  
  V(pgg)$cat.new = V(pgg)$cat.old = 0
  V(pgg)$cat.new[names(V(pgg)) %in% new.states.obs] = 2
  V(pgg)$count.new = 0
  for(i in 1:length(names(V(pgg)))) {
    v = names(V(pgg))[i]
    if(v %in% new.states.obs) {
      V(pgg)$count.new[i] = length(which(new.states.obs == v))    
    }
  }
  V(pgg)$cat.old[names(V(pgg)) %in% old.states.raw] = 1
  V(pgg)$count.old = 0
  for(i in 1:length(names(V(pgg)))) {
    v = names(V(pgg))[i]
    if(v %in% old.states.raw) {
      V(pgg)$count.old[i] = length(which(old.states.raw == v))    
    }
  }
  
  nudge = 0.14
  embed.g[[speciesref]] =  ggraph(pgg, layout="sugiyama") + 
    geom_edge_link(color="#EEEEEE", 
                   aes(width=Flux, label=label3), 
                   label_size=2, label_colour = "#AAAAAA", check_overlap = FALSE) +
    geom_node_point(data = ~subset(.x, cat.old == 1),
                    shape = 19, alpha = 0.5, color = "#FF0000",
                    aes(size = count.old),
                    position = position_nudge(y = nudge)
    ) +
    
    geom_node_point(
      data = ~subset(.x, cat.new == 2),
      shape = 19, alpha = 0.5, color = "#0000BB",
      aes(size = count.new),
      position = position_nudge(y = -nudge)
    ) +
    geom_node_point(
      data = ~subset(.x, cat.new == 2 | cat.old == 1),
      shape = 3, alpha = 0.8, size=1, color = "#888888",
    ) +
    scale_edge_width(range = c(0, 5)) +
    theme_void() + 
    labs(colour="Period", size="# samples", width="Flux")
  
}

size.scale = scale_size_continuous(
  breaks = function(x) {
    rng <- range(x)
    pretty(seq(floor(rng[1]), ceiling(rng[2]), by = 1), n = 3)
  }
)

sf = 4
png("embed-nets.png", width=600*sf, height=250*sf, res=72*sf)
ggarrange(embed.g[[2]]+size.scale, embed.g[[1]]+size.scale, labels=c("Kp", "Ec"))
dev.off()

png("embed-nets-past.png", width=800*sf, height=250*sf, res=72*sf)
ggarrange(embed.g[[3]]+size.scale, embed.g[[2]]+size.scale, 
          embed.g[[1]]+size.scale, labels=c("Kp historic", "Kp", "Ec"),
          nrow=1)
dev.off()


#######

nudge = 0.1
ggraph(pgg, layout="sugiyama") + 
  geom_edge_link(color="#EEEEEE", 
                 aes(width=Flux)) +
  geom_node_point(data = ~subset(.x, cat.old == 1),
    shape = 19, alpha = 0.5, color = "#FF0000",
    aes(size = count.old),
    position = position_nudge(y = nudge)
  ) +
  
  geom_node_point(
    data = ~subset(.x, cat.new == 2),
    shape = 19, alpha = 0.5, color = "#0000BB",
    aes(size = count.new),
    position = position_nudge(y = -nudge)
  ) +
  geom_node_point(
    data = ~subset(.x, cat.new == 2 | cat.old == 1),
    shape = 3, alpha = 0.8, size=2, color = "#000000",
  ) +
  scale_edge_width(range = c(0, 5)) +
  theme_void() + 
  labs(colour="Period", size="# samples", width="Flux")

