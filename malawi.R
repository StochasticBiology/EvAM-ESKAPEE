library(hyperinf)
library(hypertrapsct)
library(ggpubr)
sf = 2

  name <- load("~/Dropbox/klebevo/kp-evolution-inference-curate-wip/kleborate-analysis/all_models.Rdata")
  country.list <- get(name)

ggarrange(  plotHypercube.bubbles(country.list$Malawi[[1]]) + labs(size="Probability"),
  plotHypercube.sampledgraph2(country.list$Malawi[[1]], node.labels = FALSE,
                              no.times = TRUE, truncate = 6, edge.label.angle = "along",
                              edge.label.size=3))
############
# pull from Kp paper

  plot_hyperinf(country.list$Malawi[[1]])
# here we see differences
fit.x = multiple_fits_to_booted_fit(country.list$Malawi)
  
  plot_hyperinf(fit.x$boots[[1]])
  
fit.y = multiple_fits_to_booted_fit(country.list$South_Korea)

#plot_hyperinf_bubbles(fit.x)
#plot_hyperinf_bubbles(fit.y, sqrt.trans=TRUE)

c.plot = plot_hyperinf_compare_orderings(fit.x, fit.y, sqrt.trans=TRUE, expt.names=c("Gambia", "South Korea"))

png("c-plot-1.png", width=600*sf, height=400*sf, res=72*sf)
print(c.plot)
dev.off()

# here there's lots of noise
fit.x = multiple_fits_to_booted_fit(country.list$Venezuela)
fit.y = multiple_fits_to_booted_fit(country.list$Nigeria)

c.plot = plot_hyperinf_compare_orderings(fit.x, fit.y, sqrt.trans=TRUE, expt.names=c("Venezuela", "Nigeria"))

png("c-plot-2.png", width=600*sf, height=400*sf, res=72*sf)
print(c.plot)
dev.off()
