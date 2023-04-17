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
if (grepl("bb", mz.dataset)) { 
  mz_col = "mz_good_names" 
} else if (grepl("st", mz.dataset)) { 
  mz_col = "mz_struct" 
} else if (grepl("zeisel", mz.dataset)) {
  mz_col = "ms_ClusterName"
} else {
  mz_col = "mz_good_names" 
}
if (grepl("tasic", mm.dataset)) {
  mm_col = "mm_cluster"
  col.pal = rev(brewer.pal(11, "PRGn")[1:6])
} else if ( mm.dataset == "saunders" ) {
  mm_col = "mm_tissue_cluster"
  col.pal = rev(brewer.pal(11, "PiYG")[1:6])
} else if ( mm.dataset == "saunders_sub" ) {
  mm_col = "mm_tissue_subcluster"
  col.pal = rev(brewer.pal(11, "PiYG")[1:6])
} else if (grepl("oritz", mm.dataset)) {
  mm_col = "mm_g_parent"
  col.pal = brewer.pal(9, "Blues")
} else if (grepl("zeisel", mm.dataset)) {
  mm_col = "mm_ClusterName"
  col.pal = brewer.pal(9, "Greens")
} else if (grepl("turtle", mm.dataset)) {
  mm_col = "cp_detail"
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
  return(as.vector(perm_cor))
}

# Permutations
message("Performing Permutations")
nperm = 1000
perm_cor = mclapply(1:nperm, function(x) permCluster(), mc.cores = 20)
perm_cor_vector = unlist(perm_cor)
# mzmm.melt$p   = sapply(1:nrow(mzmm.melt), function(x) 1 - (length(which(perm_cor_vector<mzmm.melt$cor[x])) / length(perm_cor_vector)) )
rank_df = rbind(data.frame(id = 1:nrow(mzmm.melt), value = mzmm.melt$cor, real = T),
                data.frame(id = "perm", value = perm_cor_vector, real = F))
rank_df = rank_df[order(rank_df$value, decreasing = T),]
rank_df$num_perm_greater = cumsum(!rank_df$real)
rank_df$p = rank_df$num_perm_greater / length(perm_cor_vector)
mzmm.melt$p = rank_df$p[match(1:nrow(mzmm.melt), rank_df$id)]

# Multiple Hypothesis Test Correction
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

# Proportion Plot
message("Plotting SAMap Cluster Proportion")
if (mz.dataset == "vert2") {
  mm_over_mm_cluster = unclass(table(meta[,"species"], meta[,"leiden_clusters"]))
  mm_over_mm_cluster = mm_over_mm_cluster / colSums(mm_over_mm_cluster)
  species_count = unclass(table(meta[,"species"]))
  species_count_adj = species_count / sum(species_count)
  print(mm_over_mm_cluster[1:5])
  print(species_count_adj)
  print(dim(mm_over_mm_cluster))
  relative_prop = t(t(mm_over_mm_cluster) / species_count_adj)
  relative_prop = t(t(relative_prop) / colSums(relative_prop))
  df_prop = reshape2::melt(relative_prop)
  print(head(df_prop))
  df_prop = rbind(data.frame(prop = relative_prop, species = 'mm', cluster = names(relative_prop)), 
                  data.frame(prop = 1-relative_prop, species = 'mz', cluster = names(relative_prop)))
  df_prop$cluster = factor(as.numeric(df_prop$cluster), levels = as.character(sort(unique(as.numeric(df_prop$cluster)))))
  df_prop$prop = df_prop$prop * 100
  df_prop$color = "goldenrod1"
  df_prop$color[which(df_prop$species == "mm")] =  colorRampPalette(col.pal)(100)[80]
  df_prop$color[which(df_prop$species == "cp")] =  colorRampPalette(brewer.pal(11, "BrBG")[6:11])(100)[80]
  df_prop$color[which(df_prop$species == "tg")] =  colorRampPalette(rev(brewer.pal(11, "PuOr")[1:6]))(100)[80]
  df_prop$color[which(df_prop$species == "am")] =  magma(100)[80]
  # print(ggplot(df_prop, aes(x = prop, y = cluster, fill = color)) + geom_bar(stat='identity') + scale_x_continuous(expand = c(0,0), name = "") + ylab("") + theme_classic() + theme(axis.text = element_text(size = 10)) + geom_vline(xintercept = 50, linetype = "dashed", color = "gray40") + geom_vline(xintercept = 25, linetype = "dashed", color = "gray60") + geom_vline(xintercept = 75, linetype = "dashed", color = "gray60") + scale_fill_identity()) 
  # ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, "_prop.pdf"), width = 2.5, height = length(relative_prop)*0.15, limitsize = F)
  print(ggplot(df_prop, aes(x = cluster, y = prop, fill = color)) + geom_bar(stat='identity') + scale_y_continuous(expand = c(0,0), name = "") + xlab("") + theme_classic() + theme(axis.text = element_text(size = 10)) + geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") + geom_hline(yintercept = 25, linetype = "dashed", color = "gray60") + geom_hline(yintercept = 75, linetype = "dashed", color = "gray60") + scale_fill_identity())
  ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, "_prop.pdf"), width = length(relative_prop)*0.225, height = 2.5, limitsize = F)
} else {
  mz_species = meta[which(meta[,mz_col] != "unassigned")[1], "species"]
  mm_species = meta[which(meta[,mm_col] != "unassigned")[1], "species"]
  mm_over_mm_cluster = unclass(table(meta[,"species"], meta[,"leiden_clusters"]))
  mm_over_mm_cluster = mm_over_mm_cluster[mm_species,] / mm_over_mm_cluster[mz_species,]
  mm_over_mz = length(which(meta[,mm_col] != "unassigned")) / length(which(meta[,mz_col] != "unassigned"))
  relative_prop = (mm_over_mm_cluster / mm_over_mz) / ((mm_over_mm_cluster / mm_over_mz) + 1)
  df_prop = rbind(data.frame(prop = relative_prop, species = 'mm', cluster = names(relative_prop)), 
                  data.frame(prop = 1-relative_prop, species = 'mz', cluster = names(relative_prop)))
  df_prop$cluster = factor(as.numeric(df_prop$cluster), levels = as.character(sort(unique(as.numeric(df_prop$cluster)))))
  df_prop$prop = df_prop$prop * 100
  df_prop$color = "goldenrod1"
  df_prop$color[which(df_prop$species == "mm")] =  colorRampPalette(col.pal)(100)[80]
  # print(ggplot(df_prop, aes(x = prop, y = cluster, fill = color)) + geom_bar(stat='identity') + scale_x_continuous(expand = c(0,0), name = "") + ylab("") + theme_classic() + theme(axis.text = element_text(size = 10)) + geom_vline(xintercept = 50, linetype = "dashed", color = "gray40") + geom_vline(xintercept = 25, linetype = "dashed", color = "gray60") + geom_vline(xintercept = 75, linetype = "dashed", color = "gray60") + scale_fill_identity()) 
  # ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, "_prop.pdf"), width = 2.5, height = length(relative_prop)*0.15, limitsize = F)
  print(ggplot(df_prop, aes(x = cluster, y = prop, fill = color)) + geom_bar(stat='identity') + scale_y_continuous(expand = c(0,0), name = "") + xlab("") + theme_classic() + theme(axis.text = element_text(size = 10)) + geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") + geom_hline(yintercept = 25, linetype = "dashed", color = "gray60") + geom_hline(yintercept = 75, linetype = "dashed", color = "gray60") + scale_fill_identity())
  ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, "_prop.pdf"), width = length(relative_prop)*0.225, height = 2.5, limitsize = F)
}

# Plot
message("Plotting Celltype Mapping")
ggplot(mzmm.melt, aes(x = mm.cluster, y = mz.cluster, fill = cor)) + geom_raster() + scale_fill_gradientn(colors = col.pal, limits = c(0, 1), oob=squish) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + force_panelsizes(rows = unit(nrow(mzmm.cor)/8, "in"), cols = unit(ncol(mzmm.cor)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.4, color = "black") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1, color = "white")
ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, ".pdf"), width = (ncol(mzmm.cor)/5)+2, height = (nrow(mzmm.cor)/5)+2, limitsize = F)
message("All Done")