# EvAM-ESKAPEE
Evolutionary accumulation modelling for AMR in ESKAPEE pathogens

For old dev content, see `old-versions`.

Uses hyperinf https://github.com/StochasticBiology/hyperinf including HyperMk2 and HyperHMM. Also uses `ncbi-genome-download` and `mash`, and pulls data from the CABBAGE database.

Various evolutionary accumulation modelling approaches for comparing MDR evolution across different pathogens (and covariates). HyperMk2 and HyperHMM outputs are compared when questions are small-scale. HyperHMM alone taken forward to explore bigger questions. Phylogenies are constructed using batch-downloaded genomes mash'd, presence/absence features are from CABBAGE.

`eskapee-inference.R` does data curation and inference of evolutionary pathways for a given set of bugs and drugs. `eskapee-inference-covariate.R` also looks at behaviour across levels of a factor covariate (geographical region; age decade). `eskapee-plotting.R` produces plots summarising the dynamics in these different cases. Across these codes, an "experiment" (bug-drug combination) is labelled by an integer, a covariate is labelled by another, and whether HyperMk2 is used or not is given as a Boolean in output. The curated data inputs and fitted model outputs are written to an Rdata file by `eskapee-inference*.R`; this Rdata file is read in by `eskapee-plotting.R` and plot output files are produced.

For example, `eskapee-phylo-fits-4-cov-35-TRUE.Rdata` is covariate 35 (geographical region) for experiment 4 (the 5 top drugs in Ec and Kp), including a HyperMk2 run. Output files will be labelled with this same pattern, with a final appended integer saying how many elements of the full experiment were included (in case we just want to look at a subset of bugs, for example).

Data curation is through script `mash-batch-cli.sh`. The R code first identifies accessions for all genomes in the current experiment. The script then downloads genomes 400 at a time into folders labelled by species, experiment, and covariate (if applicable). These are sketched in batches by `mash` and the genome records discarded to retain disk space. When all genomes have been sketched, a distance matrix is produced for each species. This is read in to R and neighbour-joining used to estimate the tree for this bug-experiment.

*Next to do:* hypothesis testing for covariate influence.

References
-----

Dickens, E., Derelle, R., Beardmore, R., Suresh, A., Uplekar, S., Azov, A. G., Gurbich, T. A., Houdaigui, B. E., Keatley, J., Ochkalova, S., Koci, O., Rahman, N. M., Shivalikanjli, A., Winterbottom, A., Yordanova, G., Parkinson, H., Yates, A. D., Finn, R. D., Lees, J. A., & Chindelevitch, L. (2026). A comprehensive AMR genotype-phenotype database (CABBAGE) (p. 2025.11.12.688105). bioRxiv. https://doi.org/10.1101/2025.11.12.688105

Ondov, B. D., Treangen, T. J., Melsted, P., Mallonee, A. B., Bergman, N. H., Koren, S., & Phillippy, A. M. (2016). Mash: Fast genome and metagenome distance estimation using MinHash. Genome Biology, 17(1), 132. https://doi.org/10.1186/s13059-016-0997-x
