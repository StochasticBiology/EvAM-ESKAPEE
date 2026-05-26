library(arrow)
library(dplyr)
library(tidyr)
library(hyperinf)
library(ggrepel)
library(ggpubr)
library(hypermk2)

# totals so far: 5/TRUE, 0/TRUE, 
# covariates so far: 5/TRUE/33 [geo region]; 4/TRUE/35 [age decade]; 7/TRUE/33 [TB, geo region]
expt = 0
run.hmk2 = TRUE
plot.hmk2 = FALSE
run.diagnostics = TRUE
cov.index = FALSE

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
} else {
  fig.name.1 = paste0("eskapee-phylo-comparison-mk2-", expt, "-", plot.hmk2, "-", fname.index, ".png", collapse="")
  fig.name.2 = paste0("eskapee-phylo-comparison-simple-mk2-", expt, "-", plot.hmk2, "-", fname.index, ".png", collapse="")
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


##############

##############

### particular cases
plot_hyperinf_ordering_matrices(fit.hmm.phy[to.plot])
plot_hyperinf_compare_orderings(fit.hmm.phy[[4]], fit.hmm.phy[[1]], 
                                thetastep = 3,
                                expt.names = names(data.set)[c(4,1)])

