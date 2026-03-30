# EvAM-ESKAPEE
Evolutionary accumulation modelling for AMR in ESKAPEE pathogens

Kp-test.R – compare countries from Olav’s Kp dataset that fall at extreme PC points

Bugs-covariates-boot.R – covariates (sex, country, host, year, latitude, region, age) with promises about number of levels and number of features, for Mb, Sa, Kp, Ec, Pa

Bugs-uncertain.R – similar to above, older, supports missing data with HyperLAU

Eksape-boot.R – compare sets of ESKAPEE pathogen dynamics for max set of available drugs. Also a collection of backs-and-forth HMM-TraPS

 Eskape-rev.R – HyperMk look at ESKAPEE pathogens

Proof-of-principle.R – simple datasets, est phylo vs CS, rev vs irrev

Genomes.R (uses mash-batch.sh) – pipeline downloading and phylogeneticising actual genome data
needs `pip install ncbi-genome-download` and `brew install mash`
