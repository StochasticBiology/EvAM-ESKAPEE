library(arrow)
library(dplyr)
library(tidyr)
library(hyperinf)
library(ggrepel)
library(ggpubr)
library(hypermk2)

# covariates so far: 5/TRUE/33 [geo region]; 4/TRUE/35 [age decade]; 7/TRUE/33 [TB, geo region]
expt = 7
run.hmk2 = TRUE
run.diagnostics = TRUE
cov.index = 33

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
  plot.d.om = plot_hyperinf_ordering_matrices(c(fit.mk2, fit.hmm.phy), 
                                  expt.names = c(paste(names(fit.mk2), "Mk2"),
                                    paste(names(fit.hmm.phy), "HMM")))                                                                            
  
  plot.d.com.full = plot_hyperinf_comparative(c(fit.mk2, fit.hmm.phy), 
                            expt.names = c(paste(names(fit.mk2), "Mk2"),
                                           paste(names(fit.hmm.phy), "HMM")),
                            style = "full",
                            feature.names = substr(drug.labels, start=1, stop=3)) 
  
  plot.d.com.lim = plot_hyperinf_comparative(c(fit.mk2, fit.hmm.phy), 
                                              expt.names = c(paste(names(fit.mk2), "Mk2"),
                                                             paste(names(fit.hmm.phy), "HMM")),
                                              style = "limited",
                                              feature.names = substr(drug.labels, start=1, stop=3)) 
  
  
  fig.name = paste0("eskapee-phylo-diagnostics-", expt, "-trellis.png", collapse="")
  sf = 3
  png(fig.name, width=1200*sf, height=500*sf, res=72*sf)
  print(ggarrange(plot.d.om, plot.d.com.lim))
  dev.off()
}

drug.labels = all.fits$fit.hmm[[1]]$feature.names
drug.labels
plot.drugs = ggtexttable(data.frame(Drug=drug.labels))
fit.hmm = all.fits$fit.hmm
fit.hmm.phy = all.fits$fit.hmm.phy

################

to.plot = 1:length(fit.hmm.phy)
#to.plot = 1:2
fname.index = to.plot[length(to.plot)]

if(cov.index != FALSE) {
  fig.name.1 = paste0("eskapee-phylo-comparison-hmm-", expt, "-cov-", cov.index, "-", run.hmk2, "-", fname.index, ".png", collapse="")
  fig.name.2 = paste0("eskapee-phylo-comparison-simple-hmm-", expt, "-cov-", cov.index, "-", run.hmk2, "-", fname.index, ".png", collapse="")
} else {
  fig.name.1 = paste0("eskapee-phylo-comparison-hmm-", expt, "-", run.hmk2, "-", fname.index, ".png", collapse="")
  fig.name.2 = paste0("eskapee-phylo-comparison-simple-hmm-", expt, "-", run.hmk2, "-", fname.index, ".png", collapse="")
}

plot.om = plot_hyperinf_ordering_matrices(fit.hmm.phy[to.plot], 
                                          expt.names = names(fit.hmm.phy)[to.plot])  
plot.net = plot_hyperinf_comparative(fit.hmm.phy[to.plot], 
                                     expt.names = names(fit.hmm.phy)[to.plot], 
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
                          ggarrange(plot.om, plot.net), nrow=2)

sf = 2
png(fig.name.1, width=1000*sf, height=800*sf, res=72*sf)
print(compare.plot)
dev.off()

sf = 2
png(fig.name.2, width=800*sf, height=400*sf, res=72*sf)
print( ggarrange(plot.om, plot.net), nrow=1)
dev.off()
