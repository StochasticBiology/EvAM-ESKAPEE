library(arrow)
library(dplyr)
library(tidyr)
library(hyperinf)
library(ggrepel)
library(ggpubr)
library(hypermk2)
library(igraph)
library(ggraph)

# 5: Ec, Kp, Ab, Pa x 8
# 0: Ec, Kp, Sa, Ab, Ef, Es x 4
# 8: Ec, Kp x 10
# totals so far: 5/TRUE, 0/TRUE, 1/FALSE, 
# covariates so far: 5/TRUE/33 [geo region]; 4/TRUE/35 [age decade]; 7/TRUE/33 [TB, geo region]
expt = 0
run.hmk2 = TRUE
plot.hmk2 = TRUE
run.diagnostics = FALSE
cov.index = FALSE
sf = 2

if(cov.index != FALSE) {
  fname = paste0("eskapee-phylo-fits-", expt, "-cov-", cov.index, "-", run.hmk2, ".Rdata", collapse="")
} else {
  fname = paste0("eskapee-phylo-fits-", expt, "-", run.hmk2, ".Rdata", collapse="")
}

load(fname)

drug.labels = all.fits$fit.hmm[[1]]$feature.names
drug.labels
plot.drugs = ggtexttable(data.frame(Drug=drug.labels))
fit.hmm = all.fits$fit.hmm
fit.hmm.phy = all.fits$fit.hmm.phy
data.set = all.fits$data.set
tree.set = all.fits$tree.set

if(run.hmk2 == TRUE) {
  fit.mk2 = all.fits$fit.hmk2
}

################
if(run.diagnostics == TRUE) {
  fit.list = list()
  fit.names = c()
  j = 1
  for(i in 1:length(fit.mk2)) {
    fit.list[[j]] = c(fit.list, fit.mk2[[i]])
    fit.names = c(fit.names, names(fit.mk2)[[i]])
    j = j+1
    fit.list[[j]] = c(fit.list, fit.hmm.phy[[i]])
    fit.names = c(fit.names, names(fit.hmm.phy)[[i]])
    j = j+1
  }
  plot.d.om = plot_hyperinf_ordering_matrices(fit.list, type="relative", expt.names=fit.names)
  plot.d.om.1 = plot_hyperinf_ordering_matrices(fit.list, type="absolute", expt.names=fit.names)
  
  if(FALSE) {
    plot.d.om = plot_hyperinf_ordering_matrices(c(fit.mk2, fit.hmm.phy), 
                                              type = "relative",
                                  expt.names = c(names(fit.mk2),names(fit.hmm.phy)))                                                                            
  
  plot.d.om.1 = plot_hyperinf_ordering_matrices(c(fit.mk2, fit.hmm.phy), 
                                              type = "absolute",
                                              expt.names = c(names(fit.mk2),names(fit.hmm.phy)))                                                                            
  }
  
  if(FALSE) {
  plot.d.com.full = plot_hyperinf_comparative(c(fit.mk2, fit.hmm.phy), 
                                              expt.names = c(names(fit.mk2),names(fit.hmm.phy)),
                            style = "full",
                            feature.names = substr(drug.labels, start=1, stop=3)) 
  
  
    plot.d.com.lim = plot_hyperinf_comparative(c(fit.mk2, fit.hmm.phy), 
                                              expt.names = c(names(fit.mk2), names(fit.hmm.phy)),
                                              style = "limited",
                                              feature.names = substr(drug.labels, start=1, stop=3)) 
  }
  plot.d.com.full.mk2 = plot_hyperinf_comparative(c(fit.mk2), 
                                              expt.names = c(names(fit.mk2)),
                                              style = "full",
                                              feature.names = substr(drug.labels, start=1, stop=3)) 
  
  plot.d.com.full.hmm.phy = plot_hyperinf_comparative(c(fit.hmm.phy), 
                                                  expt.names = c(names(fit.hmm.phy)),
                                                  style = "full",
                                                  feature.names = substr(drug.labels, start=1, stop=3)) 
  
  fig.name = paste0("eskapee-phylo-diagnostics-", expt, "-trellis.png", collapse="")
  sf = 3
  png(fig.name, width=1200*sf, height=1000*sf, res=72*sf)
  print(ggarrange(plot.d.om, plot.d.om.1, 
                  plot.d.com.full.mk2, plot.d.com.full.hmm.phy,
                  nrow=2, ncol=2))
  dev.off()
}

if(plot.hmk2 == TRUE) {
  fit.plot = fit.mk2
} else {
  fit.plot = fit.hmm.phy
}

drug.labels = fit.plot[[1]]$feature.names
drug.labels
plot.drugs = ggtexttable(data.frame(Drug=drug.labels))

################

to.plot = 1:length(fit.plot)
#to.plot = 1:2
fname.index = to.plot[length(to.plot)]

if(cov.index != FALSE) {
  fig.name.1 = paste0("eskapee-phylo-comparison-mk2-", expt, "-cov-", cov.index, "-", plot.hmk2, "-", fname.index, ".png", collapse="")
  fig.name.2 = paste0("eskapee-phylo-comparison-simple-mk2-", expt, "-cov-", cov.index, "-", plot.hmk2, "-", fname.index, ".png", collapse="")
  fig.name.3 = paste0("eskapee-phylo-comparison-AICs-", expt, "-cov-", cov.index, "-", plot.hmk2, "-", fname.index, ".png", collapse="")
} else {
  fig.name.1 = paste0("eskapee-phylo-comparison-mk2-", expt, "-", plot.hmk2, "-", fname.index, ".png", collapse="")
  fig.name.2 = paste0("eskapee-phylo-comparison-simple-mk2-", expt, "-", plot.hmk2, "-", fname.index, ".png", collapse="")
  fig.name.3 = paste0("eskapee-phylo-comparison-AICs-", expt, "-", plot.hmk2, "-", fname.index, ".png", collapse="")
}

plot.om = plot_hyperinf_ordering_matrices(fit.plot[to.plot], 
                                          type = "relative",
                                          expt.names = names(fit.plot)[to.plot])  
plot.om.abs = plot_hyperinf_ordering_matrices(fit.plot[to.plot], 
                                          type = "absolute",
                                          expt.names = names(fit.plot)[to.plot])  
plot.net = plot_hyperinf_comparative(fit.plot[to.plot], 
                                     expt.names = names(fit.plot)[to.plot], 
                                     style= "full", threshold = 0.15,
                                     feature.names = substr(drug.labels, start=1, stop=3))

data.plots = list()
for(bug in names(data.set)[to.plot]) {
  data.plots[[bug]] = plot_hyperinf_data(data.set[[bug]], tree.set[[bug]],
                                         bmargin = 100,
                                         feature.names = substr(drug.labels, start=1, stop=3))
}

compare.plot = ggarrange( ggarrange(plotlist=data.plots, 
                                    labels=names(data.set)[to.plot], nrow=1),
                          ggarrange(ggarrange(plot.om, plot.om.abs, nrow=2),
                                    plot.net), nrow=2, heights=c(1,2))

sf = 2
png(fig.name.1, width=1000*sf, height=800*sf, res=72*sf)
print(compare.plot)
dev.off()

sf = 2
png(fig.name.2, width=800*sf, height=400*sf, res=72*sf)
print( ggarrange(plot.om, plot.net), nrow=1)
dev.off()


############## dAICs and interactions

nexpt = length(all.fits$data.set)

this.pc = this.null = this.null.AIC = this.hmk2.AIC = list()

for(i in 1:nexpt) {
  this.data = all.fits$data.set[[i]]
  this.tree = all.fits$tree.set[[i]]
  this.tree$edge.length = abs(this.tree$edge.length)
  #this.tree$edge.length <- this.tree$edge.length / max(this.tree$edge.length)
  this.mat = as.matrix(this.data[,2:ncol(this.data)])
  this.pc[[i]] = phylo_correlations(this.mat, this.tree)
  this.null[[i]] = hypermk2_independent(this.mat, this.tree)
  this.null.AIC[[i]] = this.null[[i]]$AIC
  this.hmk2.AIC[[i]] = all.fits$fit.hmk2[[i]]$fitted_mk$AIC
  c(this.null.AIC, this.hmk2.AIC)
} 


get_initials <- function(x) {
  sub("^([A-Za-z])[A-Za-z]*\\s+([A-Za-z])[A-Za-z]*.*$", "\\1\\2", x)
}

# delta AICs for null model comparison
daic.df = data.frame(bug=get_initials(names(all.fits$fit.hmk2)[rep(1:nexpt, 2)]),
                     AICtype = rep(c("Null", "HMk2"), each=nexpt),
                     vals=c(unlist(this.null.AIC), unlist(this.hmk2.AIC)))

plot.daic = ggplot(daic.df, aes(x=bug, y=vals, fill=AICtype)) + 
  geom_col(position="dodge") + scale_y_log10() +
  labs(x = "Species", y = "AIC", fill = "Model") +
  theme_minimal()


int.dfs = data.frame()
for(i in 1:nexpt) {
  test.daic = daic.df$vals[daic.df$bug == unique(daic.df$bug)[i] & daic.df$AICtype == "HMk2"] - 
    daic.df$vals[daic.df$bug == unique(daic.df$bug)[i] & daic.df$AICtype == "Null"] 
  this.ints = which(this.pc[[i]]$dAICs < test.daic/ncol(this.pc[[i]]$dAICs), arr.ind = TRUE)
  if(nrow(this.ints) > 0) {
  int.dfs = rbind(int.dfs, data.frame(expt=i, 
                                      from=this.ints[this.ints[,1]>this.ints[,2],1],
                                      to=this.ints[this.ints[,1]>this.ints[,2],2]))
  }
}
int.dfs
drug.names = colnames(this.pc[[i]]$dAICs)
drug.names.3 = substr(drug.names, 1, 3)


edges <- int.dfs %>%
  count(from, to, name = "weight")
edges$from = drug.names.3[edges$from]
edges$to = drug.names.3[edges$to]

# interaction graphs

g = graph_from_data_frame(edges)
plot.ints = ggraph(g) + geom_edge_link(aes(width = weight), alpha = 0.6) +
  geom_node_label(aes(label=name)) + theme_void()

png(fig.name.3, width=600*sf, height=200*sf, res=72*sf)
ggarrange(plot.daic, plot.ints, nrow=1, labels=c("A", "B"))
dev.off()

###########


test.daic = daic.df$vals[daic.df$bug == "Ec" & daic.df$AICtype == "HMk2"] - 
  daic.df$vals[daic.df$bug == "Ec" & daic.df$AICtype == "Null"] 

this.pc[[1]]$dAICs 

this.ints.df <- as.data.frame(this.pc[[1]]$dAICs) %>%
  mutate(drug1 = rownames(.)) %>%
  pivot_longer(-drug1, names_to = "drug2", values_to = "dAIC") %>%
  filter(!is.na(dAIC), dAIC < 0) %>%
  # remove duplicate pairs (since matrix is symmetric)
  rowwise() %>%
  mutate(pair = paste(sort(c(drug1, drug2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(dAIC)  # most negative first


# ^ issues with 1e100 values

this.null$by.feature[[1]]$fitted_mk$AIC
this.null$by.feature[[2]]$fitted_mk$AIC
this.null$by.feature[[3]]$fitted_mk$AIC
this.null$by.feature[[4]]$fitted_mk$AIC
this.null$by.feature[[5]]$fitted_mk$AIC
this.null$by.feature[[6]]$fitted_mk$AIC

i = 1
plot_hyperinf_data(all.fits$data.set[[i]], all.fits$tree.set[[i]])

##############

### particular cases
plot_hyperinf_ordering_matrices(fit.hmm.phy[to.plot])
plot_hyperinf_compare_orderings(fit.hmm.phy[[4]], fit.hmm.phy[[1]], 
                                thetastep = 3,
                                expt.names = names(data.set)[c(4,1)])

