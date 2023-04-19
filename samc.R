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
  mm_col = "am_cluster"
  col.pal = magma(100)
} else if (grepl("bird", mm.dataset)) {
  mm_col = "tg_cluster_orig2"
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
mm_over_mm_cluster = unclass(table(meta[,"species"], meta[,"leiden_clusters"]))
mm_over_mm_cluster = mm_over_mm_cluster / colSums(mm_over_mm_cluster)
species_count = as.vector(unclass(table(meta[,"species"])))
species_count_adj = species_count / sum(species_count)
print(mm_over_mm_cluster[1:5])
print(species_count_adj)
print(dim(mm_over_mm_cluster))
relative_prop = mm_over_mm_cluster / species_count_adj
relative_prop = t(t(relative_prop) / colSums(relative_prop))
df_prop = reshape2::melt(relative_prop)
colnames(df_prop) = c("species", "cluster", "prop")
df_prop$cluster = factor(as.numeric(df_prop$cluster), levels = as.character(sort(unique(as.numeric(df_prop$cluster)))))
df_prop$prop = df_prop$prop * 100
df_prop$color = "goldenrod1"
df_prop$color[which(df_prop$species == "mm")] =  colorRampPalette(brewer.pal(9, "Greens"))(100)[80]
df_prop$color[which(df_prop$species == "cp")] =  colorRampPalette(brewer.pal(11, "BrBG")[6:11])(100)[80]
df_prop$color[which(df_prop$species == "tg")] =  colorRampPalette(rev(brewer.pal(11, "PuOr")[1:6]))(100)[80]
df_prop$color[which(df_prop$species == "am")] =  magma(100)[80]  num_mat=unclass(table(meta[,"species"], meta[,"leiden_clusters"]))
num_df=reshape2::melt(num_mat)
df_prop$num = num_df[,3]
num_mat_notSpecies = abs(sweep(num_mat, 2, colSums(num_mat)))
num_df_notSpecies = reshape2::melt(num_mat_notSpecies)
df_prop$num2 = num_df_notSpecies[,3]
num_mat_notCluster = rowSums(num_mat) - num_mat
num_df_notCluster = reshape2::melt(num_mat_notCluster)
df_prop$num3 = num_df_notCluster[,3]
tot_num = sum(rowSums(num_mat))
df_prop$num4 = tot_num - df_prop$num - df_prop$num2 - df_prop$num3
df_prop$p = unlist(mclapply(1:nrow(df_prop), function(x) chisq.test(matrix(c(df_prop$num[x], df_prop$num2[x], df_prop$num3[x], df_prop$num4[x]), ncol=2))$p.value, mc.cores=20))
df_prop$v = unlist(mclapply(1:nrow(df_prop), function(x) vcd::assocstats(matrix(c(df_prop$num[x], df_prop$num2[x], df_prop$num3[x], df_prop$num4[x]), ncol=2))$cramer, mc.cores=20))
print(ggplot(df_prop, aes(x = cluster, y = prop, fill = color)) + geom_bar(stat='identity') + scale_y_continuous(expand = c(0,0), name = "") + xlab("") + theme_classic() + theme(axis.text = element_text(size = 10)) + geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") + geom_hline(yintercept = 25, linetype = "dashed", color = "gray60") + geom_hline(yintercept = 75, linetype = "dashed", color = "gray60") + scale_fill_identity())
ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, "_prop.pdf"), width = ncol(relative_prop)*0.225, height = 2.5, limitsize = F)
if (mz.dataset == "vert2") {
  # cp_col = "cp_detail"
  # tg_col = "tg_cluster_orig2"
  # am_col = "am_cluster"
  # 
  # prop_thresh  = 10
  # prop_thresh2 = .75
  # df_prop_specific = df_prop[which(df_prop$prop < prop_thresh),]
  # meta = cbind(meta, cpSpecificJoint = F)
  # meta[which(meta[,"species"] == "cp" & meta[,"leiden_clusters"] %in% as.vector(df_prop_specific$cluster[which(df_prop_specific$species == "cp")]) ), "cpSpecificJoint"]=T
  # cp_dist = unclass(table(meta[,cp_col], meta[,"cpSpecificJoint"]))
  # cp_dist = cp_dist / rowSums(cp_dist)
  # 
  # meta = cbind(meta, mzSpecificJoint = F)
  # meta[which(meta[,"species"] == "mz" & meta[,"leiden_clusters"] %in% as.vector(df_prop_specific$cluster[which(df_prop_specific$species == "mz")]) ), "mzSpecificJoint"]=T
  # mz_dist = unclass(table(meta[,mz_col], meta[,"mzSpecificJoint"]))
  # mz_dist = mz_dist / rowSums(mz_dist)
  # 
  # meta = cbind(meta, mmSpecificJoint = F)
  # meta[which(meta[,"species"] == "mm" & meta[,"leiden_clusters"] %in% as.vector(df_prop_specific$cluster[which(df_prop_specific$species == "mm")]) ), "mmSpecificJoint"]=T
  # mm_dist = unclass(table(meta[,mm_col], meta[,"mmSpecificJoint"]))
  # mm_dist = mm_dist / rowSums(mm_dist)
  # 
  # meta = cbind(meta, tgSpecificJoint = F)
  # meta[which(meta[,"species"] == "tg" & meta[,"leiden_clusters"] %in% as.vector(df_prop_specific$cluster[which(df_prop_specific$species == "tg")]) ), "tgSpecificJoint"]=T
  # tg_dist = unclass(table(meta[,tg_col], meta[,"tgSpecificJoint"]))
  # tg_dist = tg_dist / rowSums(tg_dist)
  # 
  # meta = cbind(meta, amSpecificJoint = F)
  # meta[which(meta[,"species"] == "am" & meta[,"leiden_clusters"] %in% as.vector(df_prop_specific$cluster[which(df_prop_specific$species == "am")]) ), "amSpecificJoint"]=T
  # am_dist = unclass(table(meta[,am_col], meta[,"amSpecificJoint"]))
  # am_dist = am_dist / rowSums(am_dist)
  # 
  # species_specific = data.frame()
  # if (ncol(cp_dist) > 1 && length(which(cp_dist[,"TRUE"] > prop_thresh2))) { species_specific = rbind(species_specific, data.frame(cluster = rownames(cp_dist)[which(cp_dist[,"TRUE"] > prop_thresh2)], species = "cp")) }
  # if (ncol(mz_dist) > 1 && length(which(mz_dist[,"TRUE"] > prop_thresh2))) { species_specific = rbind(species_specific, data.frame(cluster = rownames(mz_dist)[which(mz_dist[,"TRUE"] > prop_thresh2)], species = "mz")) }
  # if (ncol(tg_dist) > 1 && length(which(tg_dist[,"TRUE"] > prop_thresh2))) { species_specific = rbind(species_specific, data.frame(cluster = rownames(tg_dist)[which(tg_dist[,"TRUE"] > prop_thresh2)], species = "tg")) }
  # if (ncol(am_dist) > 1 && length(which(am_dist[,"TRUE"] > prop_thresh2))) { species_specific = rbind(species_specific, data.frame(cluster = rownames(am_dist)[which(am_dist[,"TRUE"] > prop_thresh2)], species = "am")) }
  # if (ncol(mm_dist) > 1 && length(which(mm_dist[,"TRUE"] > prop_thresh2))) { species_specific = rbind(species_specific, data.frame(cluster = rownames(mm_dist)[which(mm_dist[,"TRUE"] > prop_thresh2)], species = "mm")) }
} else {
  # print(ggplot(df_prop, aes(x = prop, y = cluster, fill = color)) + geom_bar(stat='identity') + scale_x_continuous(expand = c(0,0), name = "") + ylab("") + theme_classic() + theme(axis.text = element_text(size = 10)) + geom_vline(xintercept = 50, linetype = "dashed", color = "gray40") + geom_vline(xintercept = 25, linetype = "dashed", color = "gray60") + geom_vline(xintercept = 75, linetype = "dashed", color = "gray60") + scale_fill_identity()) 
  # ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, "_prop.pdf"), width = 2.5, height = length(relative_prop)*0.15, limitsize = F)
  # print(ggplot(df_prop, aes(x = cluster, y = prop, fill = color)) + geom_bar(stat='identity') + scale_y_continuous(expand = c(0,0), name = "") + xlab("") + theme_classic() + theme(axis.text = element_text(size = 10)) + geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") + geom_hline(yintercept = 25, linetype = "dashed", color = "gray60") + geom_hline(yintercept = 75, linetype = "dashed", color = "gray60") + scale_fill_identity())
  # ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, "_prop.pdf"), width = length(relative_prop)*0.225, height = 2.5, limitsize = F)
  v_thresh1 = 0.1
  df_prop_specific = df_prop[which(df_prop$v > v_thresh1),]
  not_mz_species = unique(meta[,"species"])
  not_mz_species = not_mz_species[which(not_mz_species != "mz")]
  species_specific = data.frame()
  
  mz_specific_joint_cluster = as.vector(df_prop_specific$cluster[which(df_prop_specific$species == "mz" & df_prop_specific$prop > 50)])
  if (length(mz_specific_joint_cluster) > 0) {
    meta = cbind(meta, mzSpecificJoint = F)
    meta[which(meta[,"species"] == "mz" & meta[,"leiden_clusters"] %in% mz_specific_joint_cluster ), "mzSpecificJoint"]=T
    mz_dist = unclass(table(meta[,mz_col], meta[,"mzSpecificJoint"]))
    mz_dist = mz_dist[which(rownames(mz_dist) != "unassigned"),]
    mz_dist_df = data.frame(num = mz_dist[,"TRUE"], num2 = sum(mz_dist[,"TRUE"])-mz_dist[,"TRUE"], num3 = mz_dist[,"FALSE"])
    mz_dist_df$num4 = length(which(meta[,"species"] == "mz")) - mz_dist_df$num - mz_dist_df$num2 - mz_dist_df$num3
    mz_dist_df$v = unlist(mclapply(1:nrow(mz_dist_df), function(x) vcd::assocstats(matrix(c(mz_dist_df$num[x], mz_dist_df$num2[x], mz_dist_df$num3[x], mz_dist_df$num4[x]), ncol=2))$cramer, mc.cores=20))
    mz_dist_df$prop1 = mz_dist_df$num  / mz_dist_df$num3
    mz_dist_df$prop2 = mz_dist_df$num2 / mz_dist_df$num4
    mz_specific = rownames(mz_dist_df)[which(mz_dist_df$v > v_thresh1 & mz_dist_df$prop1 > mz_dist_df$prop2)]
    if (length(mz_specific) > 0) {
      species_specific = rbind(species_specific, data.frame(cluster = mz_specific, species = "mz"))
    }
  }
  
  mm_specific_joint_cluster = as.vector(df_prop_specific$cluster[which(df_prop_specific$species == not_mz_species & df_prop_specific$prop > 50)])
  if (length(mm_specific_joint_cluster) > 0) {
    meta = cbind(meta, mmSpecificJoint = F)
    meta[which(meta[,"species"] == not_mz_species & meta[,"leiden_clusters"] %in% mm_specific_joint_cluster ), "mmSpecificJoint"]=T
    mm_dist = unclass(table(meta[,mm_col], meta[,"mmSpecificJoint"]))
    mm_dist = mm_dist[which(rownames(mm_dist) != "unassigned"),]
    mm_dist_df = data.frame(num = mm_dist[,"TRUE"], num2 = sum(mm_dist[,"TRUE"])-mm_dist[,"TRUE"], num3 = mm_dist[,"FALSE"])
    mm_dist_df$num4 = length(which(meta[,"species"] == not_mz_species)) - mm_dist_df$num - mm_dist_df$num2 - mm_dist_df$num3
    mm_dist_df$v = unlist(mclapply(1:nrow(mm_dist_df), function(x) vcd::assocstats(matrix(c(mm_dist_df$num[x], mm_dist_df$num2[x], mm_dist_df$num3[x], mm_dist_df$num4[x]), ncol=2))$cramer, mc.cores=20))
    mm_dist_df$prop1 = mm_dist_df$num  / mm_dist_df$num3
    mm_dist_df$prop2 = mm_dist_df$num2 / mm_dist_df$num4
    mm_specific = rownames(mm_dist_df)[which(mm_dist_df$v > v_thresh1 & mm_dist_df$prop1 > mm_dist_df$prop2)]
    if (length(mm_specific) > 0) {
      species_specific = rbind(species_specific, data.frame(cluster = mm_specific, species = not_mz_species))
    }
  }
  write.csv(species_specific, paste0(samc_folder, mz.dataset, "_", mm.dataset, "_species_specific.csv"))
}

# Plot
message("Plotting Celltype Mapping")
ggplot(mzmm.melt, aes(x = mm.cluster, y = mz.cluster, fill = cor)) + geom_raster() + scale_fill_gradientn(colors = col.pal, limits = c(0, 1), oob=squish) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + force_panelsizes(rows = unit(nrow(mzmm.cor)/8, "in"), cols = unit(ncol(mzmm.cor)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.4, color = "black") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1, color = "white")
ggsave(paste0(samc_folder, mz.dataset, "_", mm.dataset, ".pdf"), width = (ncol(mzmm.cor)/5)+2, height = (nrow(mzmm.cor)/5)+2, limitsize = F)
message("All Done")