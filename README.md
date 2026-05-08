# EvAM-ESKAPEE
Evolutionary accumulation modelling for AMR in ESKAPEE pathogens

For old dev content, see `old-versions`.

Uses hyperinf https://github.com/StochasticBiology/hyperinf including HyperMk2 and HyperHMM. Also uses `ncbi-genome-download` and `mash`, and pulls data from the CABBAGE database.

Various evolutionary accumulation modelling approaches for comparing MDR evolution across different pathogens (and covariates). First HyperMk2 and HyperHMM outputs are compared for a small subset question. Then HyperHMM is taken forward to explore bigger questions. Phylogenies are constructed using downloaded genomes mash'd, presence/absence features are from CABBAGE.

References
-----

Dickens, E., Derelle, R., Beardmore, R., Suresh, A., Uplekar, S., Azov, A. G., Gurbich, T. A., Houdaigui, B. E., Keatley, J., Ochkalova, S., Koci, O., Rahman, N. M., Shivalikanjli, A., Winterbottom, A., Yordanova, G., Parkinson, H., Yates, A. D., Finn, R. D., Lees, J. A., & Chindelevitch, L. (2026). A comprehensive AMR genotype-phenotype database (CABBAGE) (p. 2025.11.12.688105). bioRxiv. https://doi.org/10.1101/2025.11.12.688105

Ondov, B. D., Treangen, T. J., Melsted, P., Mallonee, A. B., Bergman, N. H., Koren, S., & Phillippy, A. M. (2016). Mash: Fast genome and metagenome distance estimation using MinHash. Genome Biology, 17(1), 132. https://doi.org/10.1186/s13059-016-0997-x
