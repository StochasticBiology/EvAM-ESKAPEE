library(arrow)
library(dplyr)
library(tidyr)
library(hyperinf)
library(ggrepel)
library(ggpubr)
library(hypermk2)

expt = 1
run.hmk2 = FALSE
run.diagnostics = FALSE
fname = paste0("eskapee-phylo-fits-", expt, "-", run.hmk2, ".Rdata", collapse="")

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
  plot_hyperinf_ordering_matrices(list(fit.mk2[[1]], fit.hmm[[1]], fit.hmm.phy[[1]]), expt.names = c("Mk2", "HMM", "HMM+Phy"))                                                                            
  
  plot_hyperinf_comparative(list(fit.mk2[[1]], fit.hmm[[1]], fit.hmm.phy[[1]]), expt.names = c("Mk2", "HMM", "HMM+Phy"),
                            feature.names = substr(drug.labels, start=1, stop=3)) 
  
  
  trellis.plot = ggarrange(plot_hyperinf_data(data.set[[1]], tree.set[[1]], bmargin=100, feature.names=gsub("-", "-\n", fit.mk2[[1]]$feature.names)),
                           plot_hyperinf_ordering_matrices(list(fit.mk2[[1]], fit.hmm[[1]], fit.hmm.phy[[1]]), expt.names = c("Mk2", "HMM", "HMM+Phy")) + theme(axis.text.x = element_text(angle=45, hjust=1)),
                           plot_hyperinf_data(data.set[[2]], tree.set[[2]], bmargin=100, feature.names=gsub("-", "-\n", fit.mk2[[1]]$feature.names)),
                           plot_hyperinf_ordering_matrices(list(fit.mk2[[2]], fit.hmm[[2]], fit.hmm.phy[[2]]), expt.names = c("Mk2", "HMM", "HMM+Phy")) + theme(axis.text.x = element_text(angle=45, hjust=1)), labels=c("Ec", "", "Kp", "")
  )
  
  fig.name = paste0("eskapee-phylo-diagnostics-", expt, "-trellis.png", collapse="")
  sf = 2
  png(fig.name, width=800*sf, height=800*sf, res=72*sf)
  print(trellis.plot)
  dev.off()
}

drug.labels = all.fits$fit.hmm[[1]]$feature.names
drug.labels
plot.drugs = ggtexttable(data.frame(Drug=drug.labels))
fit.hmm = all.fits$fit.hmm
fit.hmm.phy = all.fits$fit.hmm.phy

################

plot.om = plot_hyperinf_ordering_matrices(fit.hmm.phy, expt.names = ESKAPEE)  
plot.net = plot_hyperinf_comparative(fit.hmm.phy, expt.names = ESKAPEE, 
                                     style= "full", threshold = 0.15,
                                     feature.names = substr(drug.labels, start=1, stop=3))

data.plots = list()
for(bug in ESKAPEE) {
  data.plots[[bug]] = plot_hyperinf_data(data.set[[bug]], tree.set[[bug]],
                                         bmargin = 100,
                                         feature.names = substr(drug.labels, start=1, stop=3))
}

compare.plot = ggarrange( ggarrange(plotlist=data.plots, labels=ESKAPEE, nrow=1),
                          ggarrange(plot.om, plot.net), nrow=2)

fig.name = paste0("eskapee-phylo-comparison-hmm-", expt, "-", run.hmk2, ".png", collapse="")
sf = 2
png(fig.name, width=1000*sf, height=800*sf, res=72*sf)
print(compare.plot)
dev.off()
