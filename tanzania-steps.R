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
####### TASK 3 -- embed new observations in previous model vis
####### TASK 4 -- compare trained predictions to proportional null model

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

small_df = new_df[!is.na(new_df$species),c("ID", "MICROBESNG number", "admission date", "species", "drug", "drug_outcome")]
colnames(small_df)[2] = "SeqID"

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

write.table(sens_mat, file="new-phenotypes-clean.csv", row.names=FALSE, quote=FALSE)
drugs = c( "gentamicin", "trimethoprim-sulfamethoxazole", "ciprofloxacin",      
           "ceftazidime", "piperacillin-tazobactam", "meropenem", "amikacin"  )

sens_mat.r = sens_mat[,c(1:6, match(drugs, colnames(sens_mat)))]
write.table(sens_mat.r, file="new-phenotypes-clean-reduced.csv", row.names=FALSE, quote=FALSE)

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
  this.mat[[species]] = as.matrix(sens_mat[sens_mat$species == species,7:ncol(sens_mat)])
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

seqids = sens_mat$SeqID[!is.na(sens_mat$SeqID)]
lines <- readLines("fasta-files.txt")

# keep lines containing any pattern
matches <- lines[Reduce(`|`, lapply(seqids, function(p) grepl(p, lines)))]

# view or write out
matches
writeLines(matches, "filtered_lines.txt")

####### TASK 2 -- fits from CABBAGE data

drug.names = fits[[1]]$feature.names

qdf = read_parquet("phenotype.parquet")

qdf = qdf[qdf$species %in% c("Klebsiella pneumoniae", "Escherichia coli") &
            #as.numeric(qdf$collection_year) <= 2015 &
            qdf$antibiotic_name %in% drug.names, 
          c("BioSample_ID", "species", "country", "geographical_region", "collection_year", "antibiotic_name", "resistance_phenotype")]

qdf$species[qdf$species == "Klebsiella pneumoniae"] = "Kp"
qdf$species[qdf$species == "Escherichia coli"] = "Ec"

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

### now prune only to our focus set

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
  #filter(geographical_region == "Africa") %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = antibiotic_name,
    values_from = value,
    values_fn = max,   # collapse duplicates safely
    values_fill = 0
  )

qfits = qfits2 = qfits3 = qthis.mat = list()
for(species in unique(qsens_mat$species)) {
  qthis.mat[[species]] = as.matrix(qsens_mat[qsens_mat$species == species,7:ncol(qsens_mat)])
  if(run.inference == TRUE) {
    qfits[[species]] = hyperinf(qthis.mat[[species]], method="hyperhmm", boot.parallel=10)
    if(FALSE) {
      # run HyperTraPS analysis (not used in final version)
      set.seed(1)
      qfits2[[species]] = hyperinf(qthis.mat[[species]][sample(1:nrow(qthis.mat[[species]]), 200),], 
                                   method="hypertraps")
      qfits3[[species]] = hyperinf(qthis.mat[[species]][sample(1:nrow(qthis.mat[[species]]), 200),], 
                                   method="hypertraps", model=1)
    }
  }
}

if(run.inference == TRUE) {
  save(qfits, file = "qtanzania-fits-major.Rdata")
} else {
  load("qtanzania-fits-major.Rdata")
}

# use HyperHMM or HyperTraPS for predictions?
w.tech = "hyperhmm"
#w.tech = "hypertraps2"
#w.tech = "hypertraps1"
if(w.tech == "hyperhmm") {
  qfitsw = qfits
} else if(w.tech == "hypertraps2") {
  qfitsw = qfits2
} else {
  qfitsw = qfits3
  pdf = data.frame()
  for(i in 1:7) {
    for(j in 1:7) {
      pdf = rbind(pdf, data.frame(feature=i,
                                  step=j,
                                  count=length(which(qfitsw$Ec$routes[,1:j]==i-1))))
      propns = rep(0,7)
      for(i in 1:7) {
        propns[i] = mean(pdf$count[pdf$feature == i])/max(pdf$count)
      }
    }
  }
}


qplot.comp = plot_hyperinf_comparative(qfitsw[1:2], expt.names=names(qfitsw),
                                       feature.names = substr(qfits[[1]]$feature.names, 1, 3),
                                       style = "full")

qplot.comp

png("both-fits-major.png", width=800*sf, height=600*sf, res=72*sf)
ggarrange(plot.comp, qplot.comp)
dev.off()

####### TASK 3 -- embed new observations in previous model vis

embed.g = embed.g.label = list()

for(species in names(qfitsw)) {
  
  cabbage.fit = qfitsw[[species]]
  
  # old.data and new.data will both store Tanzanian observations
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
  
  # get transition graph from trained model
  pg = hyperinf::get_plot_graph(cabbage.fit, threshold = 1e-10)
  pgg = pg$plot.graph
  E(pgg)$label3 = substr(E(pgg)$label, 2, 4)
  E(pgg)$label3[E(pgg)$Flux < 0.1] = ""
  
  # label vertices corresponding to new observations
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
  
  # produce styled visualisation
  nudge = 0.14
  embed.g[[species]] =  ggraph(pgg, layout="sugiyama") + 
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
  
  embed.g.label[[species]] =  ggraph(pgg, layout="sugiyama") + 
    geom_edge_link(color="#EEEEEE", 
                   aes(width=Flux, label=label3), 
                   label_size=2, label_colour = "#AAAAAA", check_overlap = FALSE) +
    geom_node_point(data = ~subset(.x, cat.old == 1),
                    shape = 19, alpha = 0.5, color = "#FF0000",
                    aes(size = count.old),
                    position = position_nudge(y = nudge)
    ) +
    geom_node_text(data = ~subset(.x, cat.old == 1), 
                   color = "#FF0000", size=2,
                   aes(label=name),
                   position = position_nudge(y = nudge)
    ) +
    
    geom_node_point(
      data = ~subset(.x, cat.new == 2),
      shape = 19, alpha = 0.5, color = "#0000BB",
      aes(size = count.new),
      position = position_nudge(y = -nudge)
    ) +
    geom_node_text(data = ~subset(.x, cat.new == 2), 
                   color = "#0000BB", size=2,
                   aes(label=name),
                   position = position_nudge(y = nudge)
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
png("embed-nets-major-raw.png", width=600*sf, height=250*sf, res=72*sf)
ggarrange(plotlist=embed.g.label, labels=names(qfits))
dev.off()

ggarrange(plotlist=embed.g.label, labels=names(qfits))

sf = 4
png("embed-nets-major.png", width=600*sf, height=250*sf, res=72*sf)
ggarrange(embed.g[[2]]+size.scale, embed.g[[1]]+size.scale, labels=c("Kp", "Ec"))
dev.off()


####### TASK 4 -- compare trained predictions to proportional null model

# null: what is the probability that a walker with no interactions emits this signal?
# trained: what is the probability that a walker on the hypercube emits this signal?

# trained: state probability (= sum of in flux) * (probability of emission at this level)
# null: state probability * (probability of emission at this level)

# plackett-luce model for ordering likelihood of x given rates lambda
likelihood_binary <- function(x, lambda) {
  # x: binary vector (0/1)
  # lambda: positive rates, same length as x
  
  stopifnot(length(x) == length(lambda))
  L <- length(lambda)
  S <- which(x == 1)
  
  # memoisation environment
  memo <- new.env(parent = emptyenv())
  
  # helper: encode a set as a string key
  key_of <- function(set) {
    if (length(set) == 0) return("empty")
    paste(sort(set), collapse = ",")
  }
  
  # recursive function f(S)
  f <- function(S_current, remaining) {
    # S_current: indices still required to come before complement
    # remaining: indices still in the pool
    
    if (length(S_current) == 0) return(1)
    
    key <- paste0(key_of(S_current), "|", key_of(remaining))
    if (exists(key, envir = memo, inherits = FALSE)) {
      return(memo[[key]])
    }
    
    total_lambda <- sum(lambda[remaining])
    
    val <- 0
    for (i in S_current) {
      # probability i is chosen next
      p_i <- lambda[i] / total_lambda
      
      # recurse with i removed
      val <- val + p_i * f(
        S_current = setdiff(S_current, i),
        remaining = setdiff(remaining, i)
      )
    }
    
    memo[[key]] <- val
    return(val)
  }
  
  fS <- f(S_current = S, remaining = seq_len(L))
  
  # uniform over steps 0..L
  return(fS / (L + 1))
}

# what structure of null model are we going to use? 
null.choice = 1
null.choice.probs = 1

# initialise structures for comparison plots and data
comp.plots = comp.plots.label = list()
comp.data = list()
comp.lik = data.frame()

# loop over Ec and Kp
for(species in names(qfitsw)) {
  # pull CABBAGE data and new tanzanian observations, including unique sets and decimal/binary
  cabbage.fit = qfitsw[[species]]
  train.data = unique(cabbage.fit$data$obs)
  colnames(train.data) = cabbage.fit$feature.names
  
  tanz.obs = as.matrix(sens_mat[sens_mat$species==species, 6:ncol(sens_mat)])
  tanz.obs.r <- tanz.obs[, match(drugs, colnames(tanz.obs))]
  colnames(tanz.obs.r) == drugs
  
  comp.data[[species]] = tanz.obs.r
  
  tanz.obs.raw = apply(tanz.obs.r, 1, BinToDec)
  tanz.obs.r.uniq = unique(tanz.obs.r)
  tanz.obs.uniq = apply(tanz.obs.r.uniq, 1, BinToDec)
  
  # what are the feature probabilities under our chosen null model?
  if(null.choice == 1) {
    feature.probs = colMeans(train.data)
  } else if(null.choice == 2) {
    feature.probs = rep(0, ncol(train.data))
    for(i in 1:ncol(train.data)) {
      this.set = qdf[qdf$species == species &
                       qdf$antibiotic_name == colnames(train.data)[i] &
                       !is.na(qdf$resistance_phenotype), ]
      this.prop = nrow(this.set[this.set$resistance_phenotype == "resistant", ]) / 
        nrow(this.set)
      feature.probs[i] = this.prop
    }
    names(feature.probs) = colnames(train.data)
  }
  
  # loop through unique states in new observations
  for(i in 1:length(tanz.obs.uniq)) {
    test.state = tanz.obs.uniq[i]
    test.state.bin = tanz.obs.r.uniq[i,]
    
    # what is the model likelihood for this state?
    if(w.tech == "hyperhmm") {
      fluxes.in = cabbage.fit$transitions[#cabbage.fit$transitions$p.boot==1 & 
        cabbage.fit$transitions$To == test.state,]
    } else {
      fluxes.in = cabbage.fit$edges[cabbage.fit$edges$To == test.state,]
    }
    trained.p.state = (sum(fluxes.in$Flux)/(length(feature.probs)+1)) / length(unique(fluxes.in$p.boot))
    if(test.state == 0) { trained.p.state = 1/(length(feature.probs) + 1) }
    # and the null model prob?
    if(null.choice.probs == 1) {
      null.p.state = likelihood_binary(test.state.bin, feature.probs)
      if(null.p.state == 0) { null.p.state = -0.05 }
    } else if(null.choice.probs == 2) {
      null.ps = abs( (1-test.state.bin) - feature.probs ) 
      null.p.state = prod(null.ps)
    }
    # add this to the growing dataset
    if(new.states.obs[i] != 0 ) {
      comp.lik = rbind(comp.lik, data.frame(species = species,
                                            obs = test.state,
                                            count = length(which(tanz.obs.raw == test.state)),
                                            label = "",
                                            trained.p.state = trained.p.state,
                                            null.p.state = null.p.state))
    }
  }
  # add labelled and unlabelled plots
  comp.lik$label = comp.lik$obs
  comp.plots[[species]] = ggplot(comp.lik[comp.lik$species == species,], aes(x=null.p.state, y=trained.p.state, label=label)) + 
    geom_abline(color="#CCCCCC", linewidth=3) +
    geom_point(aes(size=count), alpha=0.5) + 
    #    geom_text(size=3) +
    labs(x= "P(state) from\nprevalence null", y="P(state)\nfrom model") +
    theme_minimal()
  comp.plots.label[[species]] = ggplot(comp.lik[comp.lik$species == species,], aes(x=null.p.state, y=trained.p.state, label=label)) + 
    geom_abline(color="#CCCCCC", linewidth=3) +
    geom_point(aes(size=count), alpha=0.5) + 
    geom_text(size=3) +
    labs(x= "P(state) from\nprevalence null", y="P(state)\nfrom model") +
    theme_minimal()
}

# visualisation of explicit likelihoods
if(FALSE) {
  comp.lik.plot = comp.lik
  comp.lik.plot$trained.p.state = comp.lik.plot$trained.p.state*8
  comp.lik.plot$trained.p.state[comp.lik.plot$trained.p.state < 1e-3] = 1e-3
  comp.lik.plot$null.p.state = comp.lik.plot$null.p.state*8
  comp.lik.plot$null.p.state[comp.lik.plot$null.p.state < 1e-3] = 1e-3
  library(ggbeeswarm)
  library(ggrepel)
  
  ggarrange(
    ggplot(comp.lik.plot, aes(x=species, y=trained.p.state)) + 
      geom_beeswarm(alpha = 0.5, aes(size=count)) + 
      geom_text_repel(aes(label=label), size = 2, nudge_x=0.2, segment.size = 0.1) +
      scale_y_log10(limits = c(0.0008,1.1), breaks = c(0.001, 0.01, 0.1, 1), labels=c("< 1e-3", 0.01, 0.1, 1)) + 
      coord_flip(),
    ggplot(comp.lik.plot, aes(x=species, y=null.p.state)) + 
      geom_beeswarm(alpha = 0.5, aes(size=count)) + 
      geom_text_repel(aes(label=label), size = 2, nudge_x=0.2, segment.size = 0.1) +
      scale_y_log10(limits = c(0.0008,1.1), breaks = c(0.001, 0.01, 0.1, 1), labels=c("< 1e-3", 0.01, 0.1, 1)) +
      coord_flip(),
    nrow = 2
  )
}

# which new states are on high probability pathways? 
length(which(comp.lik.plot$species == "Kp" & comp.lik.plot$trained.p.state > 0.01))
length(which(comp.lik.plot$species == "Kp" & comp.lik.plot$trained.p.state < 0.01))
length(which(comp.lik.plot$species == "Ec" & comp.lik.plot$trained.p.state > 0.01))
length(which(comp.lik.plot$species == "Ec" & comp.lik.plot$trained.p.state < 0.01))

# likelihood ratios trained vs null
comp.lik$ratio = comp.lik$trained.p.state/comp.lik$null.p.state
comp.lik$ratio.plot = comp.lik$ratio

# likelihood ratio without the states HyperHMM can't predict
prod(comp.lik$ratio[comp.lik$species=="Ec" & comp.lik$trained.p.state > 1e-4])
prod(comp.lik$ratio[comp.lik$species=="Kp" & comp.lik$trained.p.state > 1e-4])

# likelihood ratio without those states and the "strange" ones
prod(comp.lik$ratio[comp.lik$species=="Ec" & comp.lik$trained.p.state > 1e-4 &
                      comp.lik$obs != 108 & comp.lik$obs != 100])
prod(comp.lik$ratio[comp.lik$species=="Kp" & comp.lik$trained.p.state > 1e-4 &
                      comp.lik$obs != 108 & comp.lik$obs != 100])

# produce labelled and unlabelled versions of the likelihood comparison plots
sf = 4
png("embed-nets-probs.png", width=600*sf, height=350*sf, res=72*sf)
ggarrange(embed.g[[2]]+size.scale, embed.g[[1]]+size.scale, 
          comp.plots[[2]] + 
            labs(size = "# samples") + 
            #   scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1), 
            #                       labels=c("= 0", "0", "0.05", "0.1")) +
            scale_size_continuous(breaks = c(0,5,10)),
          comp.plots[[1]] + 
            labs(size = "# samples") #+ 
          #    scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1), 
          #                       labels=c("= 0", "0", "0.05", "0.1"))
          ,
          heights = c(2,1),
          labels=c("A i", "B i", "ii", "ii"))
dev.off()

sf = 4
png("embed-nets-label-probs.png", width=600*sf, height=650*sf, res=72*sf)
ggarrange(embed.g.label[[2]]+size.scale, embed.g.label[[1]]+size.scale, 
          comp.plots.label[[2]] + 
            labs(size = "# samples") + 
            #  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1), 
            #                     labels=c("= 0", "0", "0.05", "0.1")) +
            scale_size_continuous(breaks = c(0,5,10)),
          comp.plots.label[[1]] + 
            labs(size = "# samples") + 
            #scale_x_continuous(breaks = c(-0.05, 0, 0.02, 0.04), 
            #                  labels=c("= 0", "0", "0.02", "0.04")),
            heights = c(2,2),
          labels=c("A i", "B i", "ii", "ii"))
dev.off()



