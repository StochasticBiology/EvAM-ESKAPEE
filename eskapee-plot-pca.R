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
# 4: Ec, Kp x 5
# 1: Ec, Kp x 12
# totals so far: 5/TRUE, 0/TRUE, 1/FALSE, 8/TRUE
# covariates so far: 5/TRUE/33 [geo region]; 4/TRUE/35 [age decade]
expt = 5
run.hmk2 = TRUE
plot.hmk2 = FALSE
cov.index.set = c(FALSE, 33)
sf = 4
cov.index = FALSE

summary.list = list()
summary.names = c()

for(cov.index in cov.index.set) {
if(cov.index != FALSE) {
  fname = paste0("eskapee-phylo-fits-", expt, "-cov-", cov.index, "-", run.hmk2, ".Rdata", collapse="")
} else {
  fname = paste0("eskapee-phylo-fits-", expt, "-", run.hmk2, ".Rdata", collapse="")
}

get_initials <- function(x) {
  sub("^([A-Za-z])[A-Za-z]*\\s+([A-Za-z])[A-Za-z]*.*$", "\\1\\2", x)
}

load(fname)

drug.labels = all.fits$fit.hmm[[1]]$feature.names

plot.drugs = ggtexttable(data.frame(Drug=drug.labels))
fit.hmm = all.fits$fit.hmm
fit.hmm.phy = all.fits$fit.hmm.phy
data.set = all.fits$data.set
tree.set = all.fits$tree.set

if(run.hmk2 == TRUE) {
  fit.mk2 = all.fits$fit.hmk2
}


if(plot.hmk2 == TRUE) {
  fit.plot = fit.mk2
} else {
  fit.plot = fit.hmm.phy
}

drug.labels = fit.plot[[1]]$feature.names
drug.labels
plot.drugs = ggtexttable(data.frame(Drug=drug.labels))



for(bug in names(fit.plot)) {
  for(boot in 1:length(fit.plot[[bug]]$boots)) {
    this.om = ordering_matrix(fit.plot[[bug]]$boots[[boot]])
    this.vec = as.vector(this.om)
    split.bug = strsplit(bug, "-")[[1]]
    if(length(split.bug) == 2) {
      this.label = paste0(get_initials(split.bug[1]), "-", substr(split.bug[2], 1, 2), boot)
    } else {
      this.label = paste0(get_initials(bug), "-", boot)
    }
    summary.list[[length(summary.list)+1]] = this.vec
    summary.names = c(summary.names, this.label)
  }
}
}

summary.mat = do.call(rbind, summary.list)
rownames(summary.mat) = summary.names

############

X = summary.mat
train_idx <- nchar(rownames(X)) < 6

X_train <- X[train_idx, , drop = FALSE]
X_other <- X[!train_idx, , drop = FALSE]

pca <- prcomp(X_train, center = TRUE)#, scale. = TRUE)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)

scores_train <- as.data.frame(pca$x)
scores_train$label <- rownames(X_train)
scores_train$group <- "train"

# project supplementary data (DO NOT refit PCA)
scores_other_mat <- predict(pca, newdata = X_other)

scores_other <- as.data.frame(scores_other_mat)
scores_other$label <- rownames(X_other)
scores_other$group <- "other"

# full plot
scores_all <- rbind(scores_train, scores_other)
scores_all$group_id <- gsub("[0-9]+$", "", scores_all$label)
scores_all$label2 = ""
scores_all$label2[grep("9", scores_all$label)] = scores_all$label[grep("9", scores_all$label)]
scores_all$label2 <- gsub("[-]*[0-9]+$", "", scores_all$label2)

scores_all_1 = scores_all
pca.plot.1 = ggplot(scores_all_1, aes(PC1, PC2)) +
  
  # points colored by group type (train vs other)
  geom_point(aes(color = factor(group, levels=c("train", "other"))), size = 2, alpha = 0.7) +
  
  # ellipses by biological/label group
  stat_ellipse(aes(group = group_id, 
                   fill = factor(group_id, levels = unique(scores_all_1$group_id))),
               geom = "polygon",
               alpha = 0.25,
               color = NA) +
  
  # labels (optional; can be noisy)
  geom_text_repel(aes(label = label2),
                  size = 3, force_pull = 0.01,
                  min.segment.length = 0,
                  max.overlaps = 20) +
  
  theme_minimal() +
  guides(fill = "none") +
  scale_colour_manual(values = c(rep("black", 1), rep("#FFAAAA", 1))) +
  scale_fill_manual(values = c(rep("blue", 4), 
                               rep("#FFCCCC", 2),
                               rep("#FFCC44", 4),
                               rep("#CCAA44", 1),
                               rep("#AA444488", 1)))  +
  xlab(sprintf("PC1 (%.1f%%)", 100 * var_explained[1])) +
  ylab(sprintf("PC2 (%.1f%%)", 100 * var_explained[2])) +
  theme(legend.position = "none")

#################################################################

# 5: Ec, Kp, Ab, Pa x 8
# 0: Ec, Kp, Sa, Ab, Ef, Es x 4
# 8: Ec, Kp x 10
# 4: Ec, Kp x 5
# 1: Ec, Kp x 12
# totals so far: 5/TRUE, 0/TRUE, 1/FALSE, 8/TRUE
# covariates so far: 5/TRUE/33 [geo region]; 4/TRUE/35 [age decade]
expt = 0
run.hmk2 = TRUE
plot.hmk2 = FALSE
cov.index.set = c(FALSE)
sf = 4
cov.index = FALSE

summary.list = list()
summary.names = c()

for(cov.index in cov.index.set) {
  if(cov.index != FALSE) {
    fname = paste0("eskapee-phylo-fits-", expt, "-cov-", cov.index, "-", run.hmk2, ".Rdata", collapse="")
  } else {
    fname = paste0("eskapee-phylo-fits-", expt, "-", run.hmk2, ".Rdata", collapse="")
  }
  
  get_initials <- function(x) {
    sub("^([A-Za-z])[A-Za-z]*\\s+([A-Za-z])[A-Za-z]*.*$", "\\1\\2", x)
  }
  
  load(fname)
  
  drug.labels = all.fits$fit.hmm[[1]]$feature.names
  
  plot.drugs = ggtexttable(data.frame(Drug=drug.labels))
  fit.hmm = all.fits$fit.hmm
  fit.hmm.phy = all.fits$fit.hmm.phy
  data.set = all.fits$data.set
  tree.set = all.fits$tree.set
  
  if(run.hmk2 == TRUE) {
    fit.mk2 = all.fits$fit.hmk2
  }
  
  
  if(plot.hmk2 == TRUE) {
    fit.plot = fit.mk2
  } else {
    fit.plot = fit.hmm.phy
  }
  
  drug.labels = fit.plot[[1]]$feature.names
  drug.labels
  plot.drugs = ggtexttable(data.frame(Drug=drug.labels))
  
  
  
  for(bug in names(fit.plot)) {
    for(boot in 1:length(fit.plot[[bug]]$boots)) {
      this.om = ordering_matrix(fit.plot[[bug]]$boots[[boot]])
      this.vec = as.vector(this.om)
      split.bug = strsplit(bug, "-")[[1]]
      if(length(split.bug) == 2) {
        this.label = paste0(get_initials(split.bug[1]), "-", substr(split.bug[2], 1, 2), boot)
      } else {
        this.label = paste0(get_initials(bug), "-", boot)
      }
      summary.list[[length(summary.list)+1]] = this.vec
      summary.names = c(summary.names, this.label)
    }
  }
}

summary.mat = do.call(rbind, summary.list)
rownames(summary.mat) = summary.names

############

X = summary.mat
train_idx <- nchar(rownames(X)) < 6

X_train <- X[train_idx, , drop = FALSE]
X_other <- X[!train_idx, , drop = FALSE]

pca <- prcomp(X_train, center = TRUE)#, scale. = TRUE)
var_explained_2 <- pca$sdev^2 / sum(pca$sdev^2)

scores_train <- as.data.frame(pca$x)
scores_train$label <- rownames(X_train)
scores_train$group <- "train"

# project supplementary data (DO NOT refit PCA)
scores_other_mat <- predict(pca, newdata = X_other)

scores_other <- as.data.frame(scores_other_mat)
scores_other$label <- rownames(X_other)
scores_other$group <- "other"

# full plot
scores_all <- rbind(scores_train, scores_other)
scores_all$group_id <- gsub("[0-9]+$", "", scores_all$label)
scores_all$label2 = ""
scores_all$label2[grep("9", scores_all$label)] = scores_all$label[grep("9", scores_all$label)]
scores_all$label2 <- gsub("[-]*[0-9]+$", "", scores_all$label2)

scores_all_2 = scores_all
pca.plot.2 = ggplot(scores_all_2, aes(PC1, PC2)) +
  
  # points colored by group type (train vs other)
  geom_point(aes(color = factor(group, levels=c("train", "other"))), size = 2, alpha = 0.7) +
  
  # ellipses by biological/label group
  stat_ellipse(aes(group = group_id, 
                   fill = factor(group_id, levels = unique(scores_all_2$group_id))),
               geom = "polygon",
               alpha = 0.25,
               color = NA) +
  
  # labels (optional; can be noisy)
  geom_text_repel(aes(label = label2),
                  size = 3, force_pull = 0.01,
                  min.segment.length = 0,
                  max.overlaps = 20) +
  
  theme_minimal() +
  guides(fill = "none") +
  scale_colour_manual(values = c(rep("black", 2))) + #, rep("#FFAAAA", 1))) +
  scale_fill_manual(values = c(rep("blue", 6)))   +
  xlab(sprintf("PC1 (%.1f%%)", 100 * var_explained_2[1])) +
  ylab(sprintf("PC2 (%.1f%%)", 100 * var_explained_2[2])) +
  theme(legend.position = "none")

#######

png("pca-plots.png", width=700*sf, height=180*sf, res=72*sf)
ggarrange(pca.plot.2, pca.plot.1, nrow=1)
dev.off()
