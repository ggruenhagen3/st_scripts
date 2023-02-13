# Input ========================================================================
# Load Libraries
source("~/scratch/st/st_scripts/st_f.R")

# Mouse Object
mouse.dataset = "zeisel"
rds.path = list()
rds.path[["oritz"]]    = "~/scratch/bcs/data/mst_norm.rds"
rds.path[["saunders"]] = "~/scratch/bcs/data/mouse_w_pc_down_norm.rds"
rds.path[["tasic"]]       = "~/scratch/bcs/data/tasic_norm_020823.rds"
rds.path[["zeisel"]]      = "~/scratch/bcs/data/l5_tel_norm.rds"
mouse.path = rds.path[[mouse.dataset]]
mouse = readRDS(mouse.path)

# Set Identity of Mouse Object
mouse.ident = list()
mouse.ident[["oritiz"]] = mouse$cluster_name
mouse.ident[["saunders"]] = paste0(mouse$region, "_", mouse$subcluster)
mouse.ident[["tasic"]] = mouse$subclass
mouse.ident[["zeisel"]] = mouse$ClusterName
Idents(mouse) = mouse.ident[[mouse.dataset]]

# Mouse DEG
mouse.deg.path = list()
mouse.deg.path[["oritz"]]    = read.csv("~/scratch/bcs/results/mst_cluster_markers_020723.csv")
mouse.deg.path[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_cluster_region_subcluster_020723.csv")
mouse.deg.path[["tasic"]]    = read.csv("~/scratch/bcs/results/tasic_subclass_020823.csv")
mouse.deg.path[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_020723.csv")
# mouse.deg.path[["zeisel"]]      = read.csv("~/scratch/bcs/results/l5_cluster_markers_tax_region_020723.csv")
mouse.deg = mouse.deg.path[[mouse.dataset]]

# Mouse Glut
mouse.glut.cells = list()
mouse.glut.cells[["saunders"]] = colnames(mouse)[which( (grepl("Slc17a6", mouse$class_marker) | grepl("Slc17a7", mouse$class_marker)) & !(grepl("Gad1", mouse$class_marker) | grepl("Gad2", mouse$class_marker)) & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) > 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) == 0 )]
mouse.glut.cells[["tasic"]]    = colnames(mouse)[which( colSums(l5@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) > 0 & colSums(l5@assays$RNA@counts[c("Gad1","Gad2"),]) == 0 )]
mouse.glut.cells[["zeisel"]]   = colnames(mouse)[which( startsWith(mouse$ClusterName, "TEGLU") & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) > 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) == 0 )]
mouse.deg.path.glut = list()
mouse.deg.path.glut[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_cluster_region_subcluster_glut_020923.csv")
mouse.deg.path.glut[["tasic"]]    = read.csv("~/scratch/bcs/results/tasic_subclass_glut_020823.csv")
mouse.deg.path.glut[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_glut_020823.csv")

# Mouse GABA
mouse.gaba.cells = list()
mouse.gaba.cells[["saunders"]] = colnames(mouse)[which( !(grepl("Slc17a6", mouse$class_marker) | grepl("Slc17a7", mouse$class_marker)) & (grepl("Gad1", mouse$class_marker) | grepl("Gad2", mouse$class_marker)) & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) == 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) > 0 )]
mouse.gaba.cells[["zeisel"]]   = colnames(mouse)[which( startsWith(mouse$ClusterName, "TEINH") | startsWith(mouse$ClusterName, "MSN") & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),])==0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),])>0 )]
mouse.deg.path.gaba = list()
mouse.deg.path.gaba[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_cluster_region_subcluster_gaba_020923.csv")
mouse.deg.path.gaba[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_gaba_020923.csv")

isGlut = F
if (isGlut) {
  message("Using mouse glutamatergic clusters")
  mouse = subset(mouse, cells=mouse.glut.cells[[mouse.dataset]])
  mouse.deg = mouse.deg.path.glut[[mouse.dataset]]
}

isGABA = T
if (isGABA) {
  message("Using mouse inhibitory clusters")
  mouse = subset(mouse, cells=mouse.gaba.cells[[mouse.dataset]])
  mouse.deg = mouse.deg.path.gaba[[mouse.dataset]]
}

# Cichlid Object and DEG
isBB = T
gene_info = read.csv("~/scratch/m_zebra_ref/gene_info_3.csv")
convert15 = read.csv("~/scratch/st/data/convert15.csv")
convert53 = read.csv("~/scratch/st/data/convert53.csv")
if (isBB) { 
  if (isGlut) {
    message("Using cichlid glutamatergic clusters")
    mz = readRDS("~/scratch/st/data/bb_glut.rds")
    mz$good_names53 = factor(convert53$new[match(mz$seurat_clusters, convert53$old)], levels = rev(convert53$new))
    Idents(mz) = "good_names53"
    mz.deg = read.csv("~/scratch/st/data/bb53_glut_deg.csv")
  } else if (isGABA) {
    message("Using cichlid glutamatergic clusters")
    mz = readRDS("~/scratch/st/data/bb_gaba.rds")
    mz$good_names53 = factor(convert53$new[match(mz$seurat_clusters, convert53$old)], levels = rev(convert53$new))
    Idents(mz) = "good_names53"
    mz.deg = read.csv("~/scratch/st/data/bb53_gaba_deg.csv")
  } else {
    mz = readRDS("~/scratch/brain/data/bb_demux_102021.rds")
    mz$good_names53 = factor(convert53$new[match(mz$seurat_clusters, convert53$old)], levels = rev(convert53$new))
    Idents(mz) = "good_names53"
    mz.deg = read.csv("~/scratch/brain/results/bb53_deg_w_one_012323.csv")
  }
  mz.deg$cluster = stringr::str_replace(mz.deg$cluster, "Astro", "RG")
} else { 
  mz = readRDS("~/scratch/st/data/st_c2b2_hi_012723.rds") 
  Idents(mz) = "structure"
  mz.deg = read.csv("~/scratch/st/data/c2b2_structure_deg.csv")
}
mz.deg$one_to_one_human = gene_info$one_to_one_human[match(mz.deg$gene, gene_info$seurat_name)]

# Find Correlations ============================================================
# Find Common Gene Set
common.gene.set = sort(unique(toupper(mouse.deg$gene)))
common.gene.set = common.gene.set[which(common.gene.set %in% mz.deg$one_to_one_human)]
mz.common.gene.set = mz.deg$gene[match(common.gene.set, mz.deg$one_to_one_human)]
mouse.common.gene.set = stringr::str_to_title(common.gene.set)

# Cichlid Average Expression
mz.assay.to.use = ifelse("SCT" %in% names(mz@assays), "SCT", "RNA")
mz.avg.exp = AverageExpression(mz, features = mz.common.gene.set, assays = mz.assay.to.use, slot = "data")[[1]]
mz.avg.exp.norm = log(mz.avg.exp+1) + 0.1
mz.avg.exp.norm = mz.avg.exp.norm / rowMeans(mz.avg.exp.norm)

# Mouse Average Expression
mouse.avg.exp = AverageExpression(mouse, features = mouse.common.gene.set, assays = "SCT", slot = "data")[[1]]
mouse.avg.exp.norm = log(mouse.avg.exp+1) + 0.1
mouse.avg.exp.norm = mouse.avg.exp.norm / rowMeans(mouse.avg.exp.norm)

# Correlation
mz.mouse.cor = cor(mz.avg.exp.norm, mouse.avg.exp.norm, method = "spearman")

# Permutation Function =========================================================
permCor = function(old.mat) {
  new.mat.list = lapply(1:nrow(old.mat), function(x) sample(old.mat[x,]))
  new.mat = do.call('rbind', new.mat.list)
  perm.cor = cor(new.mat, mouse.avg.exp.norm, method = "spearman")
  perm.cor.melt = reshape2::melt(perm.cor)
  return(perm.cor.melt[,3])
}

# Permutations =================================================================
n.perms = 1000
mz.mouse.cor.melt = reshape2::melt(mz.mouse.cor)
orig.gene.labels = unname(rownames(mz.avg.exp.norm))
perm.cor.list = mclapply(1:n.perms, function(x) permCor(mz.avg.exp.norm), mc.cores = 20)
perm.cor.mat = do.call('cbind', perm.cor.list)
# perm.supp = cbind(mz.mouse.cor.melt, perm.cor.mat)
# colnames(perm.supp) = c("mz.cluster", "mouse.cluster", "cor", paste0("perm.", 1:n.perms))
mz.mouse.cor.melt$num_perm_greater = unlist(lapply(1:nrow(mz.mouse.cor.melt), function(x) length(which(perm.cor.mat[x,] > mz.mouse.cor.melt$value[x]))))
mz.mouse.cor.melt$p.perm = (mz.mouse.cor.melt$num_perm_greater) / n.perms
mz.mouse.cor.melt$bon.perm = p.adjust(mz.mouse.cor.melt$p.perm, method = "bonferroni")

# Spearman Correlation Test ====================================================
mz.mouse.cor.p = matrix(NA, nrow = nrow(mz.mouse.cor), ncol = ncol(mz.mouse.cor), dimnames = list(rownames(mz.mouse.cor), colnames(mz.mouse.cor)))
for (mz.clust in colnames(mz.avg.exp.norm)) {
  for (mouse.clust in colnames(mouse.avg.exp.norm)) {
    mz.mouse.cor.p[mz.clust, mouse.clust] = cor.test(mz.avg.exp.norm[,mz.clust], mouse.avg.exp.norm[,mouse.clust], method = "spearman", alternative = "greater")$p.value
  }
}
mz.mouse.cor.bon = matrix(p.adjust(mz.mouse.cor.p, method = "bonferroni"), nrow = nrow(mz.mouse.cor), ncol = ncol(mz.mouse.cor), dimnames = list(rownames(mz.mouse.cor), colnames(mz.mouse.cor)))

# Output: Plots and CSV ========================================================
mz.mouse.cor.maxed.out.melt = reshape2::melt(mz.mouse.cor)
colnames(mz.mouse.cor.maxed.out.melt) = c("mz.cluster", "mouse.cluster", "cor")
mz.mouse.cor.maxed.out.melt[, c("mouse.region", "mouse.celltype")] = reshape2::colsplit(mz.mouse.cor.maxed.out.melt$mouse.cluster, "_", c('1', '2'))
mz.mouse.cor.maxed.out.melt = mz.mouse.cor.maxed.out.melt[,c("mz.cluster", "mouse.cluster", "mouse.region", "mouse.celltype", "cor")]
maxed.num = 0.35
mz.mouse.cor.maxed.out.melt$cor.maxed = mz.mouse.cor.maxed.out.melt$cor
mz.mouse.cor.maxed.out.melt$cor.maxed[which(mz.mouse.cor.maxed.out.melt$cor >  maxed.num)] =  maxed.num
mz.mouse.cor.maxed.out.melt$cor.maxed[which(mz.mouse.cor.maxed.out.melt$cor < -maxed.num)] = -maxed.num
mz.mouse.cor.maxed.out.melt$cor.test.p = reshape2::melt(mz.mouse.cor.p)[,3]
mz.mouse.cor.maxed.out.melt$cor.test.bon = reshape2::melt(mz.mouse.cor.bon)[,3]
mz.mouse.cor.maxed.out.melt$perm.test.p = mz.mouse.cor.melt$p.perm
mz.mouse.cor.maxed.out.melt$perm.test.bon = mz.mouse.cor.melt$bon.perm
mz.mouse.cor.maxed.out.melt$all.sig = mz.mouse.cor.maxed.out.melt$cor.test.bon < 0.05 & mz.mouse.cor.maxed.out.melt$perm.test.bon < 0.05
mz.order  = hclust(dist(mz.mouse.cor), method = "complete")
mz.mouse.cor.maxed.out.melt$mz.cluster = factor(mz.mouse.cor.maxed.out.melt$mz.cluster, levels = mz.order$labels[mz.order$order])
mouse.order = hclust(dist(t(mz.mouse.cor)), method = "complete")
mz.mouse.cor.maxed.out.melt$mouse.cluster = factor(mz.mouse.cor.maxed.out.melt$mouse.cluster, levels = mouse.order$labels[mouse.order$order])
ggplot(mz.mouse.cor.maxed.out.melt, aes(x = mouse.cluster, y = mz.cluster, fill = cor.maxed)) + geom_raster() + geom_point(data = mz.mouse.cor.maxed.out.melt[which(mz.mouse.cor.maxed.out.melt$all.sig),], size = 0.8) + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), n.breaks = 6, limits = c(-maxed.num, maxed.num)) + coord_fixed() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# ggplot(mz.mouse.cor.maxed.out.melt, aes(x = mz.cluster, y = mouse.cluster, fill = cor.maxed)) + geom_raster() + geom_point(data = mz.mouse.cor.maxed.out.melt[which(mz.mouse.cor.maxed.out.melt$all.sig),], size = 0.8) + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), n.breaks = 6, limits = c(-maxed.num, maxed.num)) + coord_fixed() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
cichlid_str = ifelse(isBB, "bb_", "mz_")
glut_str    = ifelse(isGlut, "glut_", "")
gaba_str    = ifelse(isGABA, "gaba_", "")
out_name = paste0("~/scratch/bcs/results/", cichlid_str, mouse.dataset, "_cluster_", glut_str, gaba_str, "cor")
ggsave(paste0(out_name, ".pdf"), width = 6, height = 6)

write.csv(mz.mouse.cor.maxed.out.melt, paste0(out_name, ".csv"))
message(paste0("rclone copy ", paste0(out_name, ".pdf"), " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/bcs/", mouse.dataset))
message(paste0("rclone copy ", paste0(out_name, ".csv"), " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/bcs/", mouse.dataset))
