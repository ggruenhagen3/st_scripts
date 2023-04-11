# Input ========================================================================
# mz.dataset = "st"; mm.dataset = "oritzg"; isGE = F
# Read Input
args = commandArgs(trailingOnly=TRUE)
mz.dataset = as.character(args[1])
mm.dataset = as.character(args[2])
message(paste0("Running comparison on the following datasets: mz.dataset=", mz.dataset, ", mm.dataset=", mm.dataset))

# Load Libraries
message("Loading Libraries")
if (! "ggh4x" %in% installed.packages()) { stop(paste0("package ggh4x is not available. Ensure that the SeuratDisk conda environment is active.")) }
library("ggplot2")
library("scales")
library("viridisLite")
library("RColorBrewer")
library("ggh4x")
library("parallel")

# Read Data
message("Reading Data")
samc_folder = "~/scratch/bcs/samc/"
meta = as.matrix(read.csv(paste0(samc_folder, mz.dataset, "_", mm.dataset, ".csv"),   row.names = 1))

# Select dataset specific variables: metadata column and color palette
col.pal = viridis(100)
if (grepl("bb", mz.dataset)) { mz_col = "mz_good_names" } else { mz_col = "mz_struct" }
if (grepl("tasic", mm.dataset)) {
  mm_col = "mm_cluster"
  col.pal = rev(brewer.pal(11, "PRGn")[1:6])
} else if (grepl("saunders", mm.dataset)) {
  col.pal = rev(brewer.pal(11, "PiYG")[1:6])
} else if (grepl("oritz", mm.dataset)) {
  mm_col = "mm_g_parent"
  col.pal = brewer.pal(9, "Blues")
} else if (grepl("zeisel", mm.dataset)) {
  mm_col = "mm_ClusterName"
  col.pal = brewer.pal(9, "Greens")
} else if (grepl("turtle", mm.dataset)) {
  col.pal = brewer.pal(11, "BrBG")[6:11]
} else if (grepl("axolotl", mm.dataset)) {
  col.pal = magma(100)
} else if (grepl("bird", mm.dataset)) {
  col.pal = rev(brewer.pal(11, "PuOr")[1:6])
}

# Input Checks
if (! mz_col            %in% colnames(meta) ) { stop(paste0("mz_col not in metadata, mz_col=", mz_col)) }
if (! mm_col            %in% colnames(meta) ) { stop(paste0("mm_col not in metadata, mm_col=", mm_col)) }
if (! 'leiden_clusters' %in% colnames(meta) ) { stop(paste0("leiden_clusters not in metadata"))         }

# Find distribution of cichild clusters in combined clusters
mz_both = unclass(table(meta[,mz_col], meta[,'leiden_clusters']))
mz_both = mz_both[which(rownames(mz_both) != "unassigned"),]
mz_both = mz_both / rowSums(mz_both)

# Find distribution of mouse clusters in combined clusters
mm_both = unclass(table(meta[,mm_col], meta[,'leiden_clusters']))
mm_both = mm_both[which(rownames(mm_both) != "unassigned"),]
mm_both = mm_both / rowSums(mm_both)

# Find Correlation
message("Finding Correlations")
mzmm.cor = cor(t(mz_both), t(mm_both))
mzmm.melt = reshape2::melt(mzmm.cor)
colnames(mzmm.melt) = c("mz.cluster", "mm.cluster", "cor")

permCluster = function() {
  mz_perm = meta[,mz_col]
  mz_perm[which(mz_perm != "unassigned")] = sample(mz_perm[which(mz_perm != "unassigned")])
  mz_perm_both = unclass(table(mz_perm, meta[,'leiden_clusters']))
  mz_perm_both = mz_perm_both[which(rownames(mz_perm_both) != "unassigned"),]
  mz_perm_both = mz_perm_both / rowSums(mz_perm_both)
  
  mm_perm = meta[,mm_col]
  mm_perm[which(mm_perm != "unassigned")] = sample(mm_perm[which(mm_perm != "unassigned")])
  mm_perm_both = unclass(table(mm_perm, meta[,'leiden_clusters']))
  mm_perm_both = mm_perm_both[which(rownames(mm_perm_both) != "unassigned"),]
  mm_perm_both = mm_perm_both / rowSums(mm_perm_both)
  perm_cor = cor(t(mz_perm_both), t(mm_perm_both))
  return(perm_cor)
}

# Permutations
message("Performing Permutations")
nperm = 1000
perm_cor = mclapply(1:nperm, function(x) permCluster(), mc.cores = 20)
perm_cor_vector = unlist(perm_cor)
mzmm.melt$p   = sapply(1:nrow(mzmm.melt), function(x) 1 - (length(which(perm_cor_vector<mzmm.melt$cor[x])) / length(perm_cor_vector)) )
mzmm.melt$bh  = p.adjust(mzmm.melt$p, method = "BH")
mzmm.melt$bon = p.adjust(mzmm.melt$p, method = "bonferroni")
mzmm.melt$p_sig   = mzmm.melt$p   < 0.05
mzmm.melt$bh_sig  = mzmm.melt$bh  < 0.05
mzmm.melt$bon_sig = mzmm.melt$bon < 0.05
mzmm.melt$p0 = mzmm.melt$p == 0
bh_sig_mat = reshape2::acast(mzmm.melt, mz.cluster ~ mm.cluster, value.var = "bh_sig")

# Arrange clusters
arrange_mat = mzmm.cor
mz.order  = hclust(dist(arrange_mat), method = "complete")
mzmm.melt$mz.cluster = factor(as.vector(mzmm.melt$mz.cluster), levels = mz.order$labels[mz.order$order])
mouse.order = hclust(dist(t(arrange_mat)), method = "complete")
mzmm.melt$mm.cluster = factor(as.vector(mzmm.melt$mm.cluster), levels = mouse.order$labels[mouse.order$order])

# Save Data
message("Saving Supplementary Data")
write.csv(mzmm.melt, paste0(samc_folder, mz.dataset, "_", mm.dataset, "_sup.csv"))

# Plot
message("Plotting")
ggplot(mzmm.melt, aes(x = mm.cluster, y = mz.cluster, fill = cor)) + geom_raster() + scale_fill_gradientn(colors = col.pal, limits = c(0, 1), oob=squish) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + force_panelsizes(rows = unit(nrow(mzmm.cor)/8, "in"), cols = unit(ncol(mzmm.cor)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.4, color = "black") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1, color = "white")
ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, ".pdf"), width = (ncol(mzmm.cor)/5), height = (nrow(mzmm.cor)/5))
message("All Done")