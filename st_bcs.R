# Input ========================================================================
# Read Input
args = commandArgs(trailingOnly=TRUE)
mouse.dataset = args[1]
isGlut = as.logical(args[2])
isGABA = as.logical(args[3])
isNN   = as.logical(args[4])
isSub1 = as.logical(args[5])
isBB   = as.logical(args[6])
message(paste0("Running correlation comparison using the following arguments: mouse.dataset= ", mouse.dataset, ", isGlut=", isGlut, ", isGABA=", isGABA, ", isNN=", isNN, ", isSub1=", isSub1, ", isBB=", isBB))

out_name_overide = ""
max_overide = 0
if (length(args) > 5) { 
  if(is.na(as.numeric(args[7]))) { out_name_overide=args[7] } else { max_overide=args[7] }
  message(paste0("Overide Argument = ", args[7]))
}
# mouse.dataset = "zeisel"; isGlut = F; isGABA = F; isNN = F; isSub1 = F; isBB = T

# Load Libraries
message("Loading Libraries")
suppressMessages(source("~/scratch/st/st_scripts/st_f.R"))
suppressMessages(source("~/scratch/bcs/bcs_scripts/bcs_f.R"))
library("ggh4x")

# Mouse Object
message("Loading mouse object")
rds.path = list()
rds.path[["oritz"]]    = "~/scratch/bcs/data/mst_norm.rds"
rds.path[["oritzb"]]    = "~/scratch/bcs/data/oritz_b_raw.rds"
rds.path[["saunders"]] = "~/scratch/bcs/data/mouse_w_pc_down_norm.rds"
rds.path[["tasic"]]    = "~/scratch/bcs/data/tasic_norm_020823.rds"
rds.path[["tran"]]     = "~/scratch/bcs/data/tran_norm.rds"
rds.path[["zeisel"]]   = "~/scratch/bcs/data/l5_tel_norm.rds"
rds.path[["zei_yu"]]   = "~/scratch/bcs/data/zei_yu_norm.rds"
mouse.path = rds.path[[mouse.dataset]]
mouse = readRDS(mouse.path)

# Set Identity of Mouse Object
message("Setting mouse identity")
mouse.ident = list()
mouse.ident = switch(mouse.dataset,
                     "oritz" = mouse$ABA_parent,
                     "oritzb" = mouse$b_parent,
                     "saunders" = paste0(mouse$region, "_", mouse$subcluster),
                     "tasic" = mouse$subclass,
                     "tran" = Idents(mouse),
                     "zeisel" = mouse$ClusterName,
                     "zei_yu" = mouse$my.specific)
Idents(mouse) = mouse.ident

# Mouse DEG
message("Loading mouse DEGs")
mouse.deg.path = list()
# mouse.deg.path[["oritz"]]    = read.csv("~/scratch/bcs/results/mst_cluster_markers_020723.csv")
mouse.deg.path[["oritz"]]    = read.csv("~/scratch/bcs/results/mst_aba_parent_markers_020823.csv")
mouse.deg.path[["oritzb"]]    = read.csv("~/scratch/bcs/results/mst_b_parent_022423.csv")
mouse.deg.path[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_cluster_region_subcluster_020723.csv")
mouse.deg.path[["tasic"]]    = read.csv("~/scratch/bcs/results/tasic_subclass_020823.csv")
mouse.deg.path[["tran"]]     = read.csv("~/scratch/bcs/results/tran_broad_021523.csv")
mouse.deg.path[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_020723.csv")
# mouse.deg.path[["zeisel"]]      = read.csv("~/scratch/bcs/results/l5_cluster_markers_tax_region_020723.csv")
mouse.deg.path[["zei_yu"]]   = read.csv("~/scratch/bcs/results/zei_yu_broad_022023.csv")
mouse.deg = mouse.deg.path[[mouse.dataset]]

if (isGlut) {
  mouse.glut.cells = switch(mouse.dataset,
                            "saunders" = colnames(mouse)[which( (grepl("Slc17a6", mouse$class_marker) | grepl("Slc17a7", mouse$class_marker)) & !(grepl("Gad1", mouse$class_marker) | grepl("Gad2", mouse$class_marker)) & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) > 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) == 0 )],
                            "tasic"    = colnames(mouse)[which( colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) > 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) == 0 )],
                            "tran"     = colnames(mouse)[which( grepl("Excit", mouse$broad) & colSums(mouse@assays$RNA@counts[c("SLC17A6","SLC17A7"),]) > 0 & colSums(mouse@assays$RNA@counts[c("GAD1","GAD2"),]) == 0 )],
                            "zeisel"   = colnames(mouse)[which( startsWith(mouse$ClusterName, "TEGLU") & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) > 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) == 0 )],
                            "zei_yu"   = colnames(mouse)[which( (grepl("TEGLU", mouse$my.specific) | grepl("ExN", mouse$my.specific)) & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) > 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) == 0 )])
  mouse.deg.path.glut = list()
  mouse.deg.path.glut[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_subcluster_glut_022423.csv")
  mouse.deg.path.glut[["tasic"]]    = read.csv("~/scratch/bcs/results/tasic_subclass_glut_020823.csv")
  mouse.deg.path.glut[["tran"]]    = read.csv("~/scratch/bcs/results/tran_broad_glut_021623.csv")
  mouse.deg.path.glut[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_glut_020823.csv")
  mouse.deg.path.glut[["zei_yu"]]   = read.csv("~/scratch/bcs/results/zei_yu_specific_glut_022023.csv")
  
  message("Using mouse glutamatergic clusters")
  mouse = subset(mouse, cells=mouse.glut.cells)
  mouse.deg = mouse.deg.path.glut[[mouse.dataset]]
}
if (isGABA) {
  mouse.gaba.cells = switch(mouse.dataset,
                            "saunders" = colnames(mouse)[which( !(grepl("Slc17a6", mouse$class_marker) | grepl("Slc17a7", mouse$class_marker)) & (grepl("Gad1", mouse$class_marker) | grepl("Gad2", mouse$class_marker)) & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) == 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) > 0 )],
                            "tran"     = colnames(mouse)[which( grepl("Inhib", mouse$broad) & colSums(mouse@assays$RNA@counts[c("SLC17A6","SLC17A7"),]) == 0 & colSums(mouse@assays$RNA@counts[c("GAD1","GAD2"),]) > 0 )],
                            "zeisel"   = colnames(mouse)[which( startsWith(mouse$ClusterName, "TEINH") | startsWith(mouse$ClusterName, "MSN") | startsWith(mouse$ClusterName, "OBINH") & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),])==0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),])>0 )],
                            "zei_yu"   = colnames(mouse)[which( (startsWith(mouse$my.specific, "TEINH") | startsWith(mouse$my.specific, "MSN") | startsWith(mouse$my.specific, "OBINH") | grepl("InN", mouse$my.specific)) & colSums(mouse@assays$RNA@counts[c("Slc17a6","Slc17a7"),]) == 0 & colSums(mouse@assays$RNA@counts[c("Gad1","Gad2"),]) > 0 )])
  mouse.deg.path.gaba = list()
  mouse.deg.path.gaba[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_subcluster_gaba_022423.csv")
  mouse.deg.path.gaba[["tran"]]     = read.csv("~/scratch/bcs/results/tran_broad_gaba_021623.csv")
  mouse.deg.path.gaba[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_gaba_020923.csv")
  mouse.deg.path.gaba[["zei_yu"]]   = read.csv("~/scratch/bcs/results/zei_yu_specific_gaba_022023.csv")
  
  message("Using mouse inhibitory clusters")
  mouse = subset(mouse, cells=mouse.gaba.cells)
  mouse.deg = mouse.deg.path.gaba[[mouse.dataset]]
}
if (isNN) {
  mouse.nn.cells = switch(mouse.dataset,
                          "saunders" = colnames(mouse)[which( paste0(mouse$region, "_", mouse$cluster) %in% c("FC_10", "FC_12", "FC_13", "FC_14", "FC_8", "FC_9", "PC_10", "PC_11", "PC_12", "PC_13", "PC_14", "PC_8", "PC_9", "GP_10", "GP_11", "GP_4", "GP_5", "GP_6", "GP_7", "GP_8", "GP_9", "HC_10", "HC_12", "HC_13", "HC_15", "HC_16", "HC_17", "HC_7", "HC_8", "HC_9", "STR_1", "STR_2", "STR_3", "STR_4", "STR_5", "STR_6", "STR_7", "STR_8", "STR_9") | paste0(mouse$region, "_", mouse$subcluster) %in% c("FC_11-1", "FC_11-3", "FC_11-4") )],
                          "tran"     = colnames(mouse)[which( !grepl("Excit", mouse$broad) & !grepl("Inhib", mouse$broad) )],
                          "zeisel"   = colnames(mouse)[which( mouse$ClusterName %in% c("SZNBL", "OEC", "RGDG", "RGSZ", "ACOB", "ACTE1", "ACTE2", "DGNBL1", "DGNBL2", "OBNBL1", "OBNBL2", "OBNBL3", "OBNBL4", "OBNBL5") )])
  mouse.deg.path.nn = list()
  mouse.deg.path.nn[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_subcluster_nn_022423.csv")
  mouse.deg.path.nn[["tran"]]     = read.csv("~/scratch/bcs/results/tran_broad_nn_021623.csv")
  mouse.deg.path.nn[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_nn_021623.csv")
  
  message("Using mouse non-neuronal clusters")
  mouse = subset(mouse, cells=mouse.nn.cells)
  mouse.deg = mouse.deg.path.nn[[mouse.dataset]]
}
if (isSub1) {
  message("Using mouse sub1 clusters")
  # mouse = subset(mouse, cells=colnames(mouse)[which(mouse$ABA_parent %in% c( "Isocortex", "Olfactory areas", "Cortical subplate", "Retrohippocampal region", "Hippocampal region" ))])
  # mouse.deg = read.csv("~/scratch/bcs/results/mst_sub1_aba_parent_markers_022123.csv")
  mouse = subset(mouse, cells = colnames(mouse)[which(mouse$b_parent %in% c( "Isocortex", "CA1", "DG", "CA2", "CA3", "RHP", "HIP-other", "PA", "TR", "U_CTX", "PAA", "LA", "TT", "DP", "PIR", "EP", "BLA", "CLA", "BMA", "AON" ))])
  mouse.deg = read.csv("~/scratch/bcs/results/oritzb_sub1_022723.csv")
}

# Cichlid Object and DEG
message("Loading cichlid data")
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
    message("Using cichlid gabaergic clusters")
    mz = readRDS("~/scratch/st/data/bb_gaba.rds")
    mz$good_names53 = factor(convert53$new[match(mz$seurat_clusters, convert53$old)], levels = rev(convert53$new))
    Idents(mz) = "good_names53"
    mz.deg = read.csv("~/scratch/st/data/bb53_gaba_deg.csv")
  } else if (isNN) {
    message("Using cichlid non-neuronal clusters")
    mz = readRDS("~/scratch/st/data/bb_nn.rds")
    mz.deg = read.csv("~/scratch/st/data/bb53_nn_deg.csv")
  } else if (isSub1) {
    message("Using bb sub1 glutamatergic clusters")
    mz = readRDS("~/scratch/st/data/bb_glut.rds")
    mz$good_names53 = factor(convert53$new[match(mz$seurat_clusters, convert53$old)], levels = rev(convert53$new))
    # mz = subset(mz, cells=colnames(mz)[which(!mz$good_names53 %in% c("15.7_Glut", "15.6_Glut", "14_Glut", "9.5_Glut", "9.8_Glut", "8.7_Glut"))])
    # mz.deg = read.csv("~/scratch/st/data/bb53_glut_sub1_deg.csv")
    mz = subset(mz, cells=colnames(mz)[which(!mz$good_names53 %in% c("15.7_Glut", "15.6_Glut", "14_Glut", "9.5_Glut"))])
  }else {
    mz = readRDS("~/scratch/brain/data/bb_demux_102021.rds")
    mz$good_names53 = factor(convert53$new[match(mz$seurat_clusters, convert53$old)], levels = rev(convert53$new))
    Idents(mz) = "good_names53"
    mz.deg = read.csv("~/scratch/brain/results/bb53_deg_w_one_012323.csv")
  }
  mz.deg$cluster = stringr::str_replace(mz.deg$cluster, "Astro", "RG")
} else { 
  mz = readRDS("~/scratch/st/data/st_c2b2_hi_022023.rds") 
  Idents(mz) = "structure"
  mz.deg = read.csv("~/scratch/bcs/results/c2b2_structure_deg_022123.csv")
  if (isSub1) {
    message("Using cichlid spatial sub1 clusters")
    # mz = subset(mz, cells = colnames(mz)[which(mz$structure %in% c( "Dm-2r", "Dc-3", "Dd", "Dm-3", "Dc-1/2", "NT", "Dp", "Dl-g", "Dl-v" ))])
    mz = subset(mz, cells = colnames(mz)[which(mz$structure %in% c( "Dm-2r", "Dc-3", "Dd", "Dm-3", "Dc-1/2", "NT", "Dp", "Dl-g", "Dl-v", "Dc-4", "Dm-2c", "Dc-5", "Dm-1" ))])
    mz = subset(mz, cells = colnames(mz)[which(mz$structure %in% c( "Dm-2r", "Dc-3", "Dd", "Dm-3", "Dc-1/2", "Dp", "Dl-g", "Dl-v", "Dc-4", "Dm-2c", "Dc-5", "Dm-1" ))])
    mz.deg = read.csv("~/scratch/bcs/results/c2b2_sub3_deg_022123.csv")
    # mz.deg = read.csv("~/scratch/bcs/results/c2b2_sub2_deg_022123.csv")
    # mz.deg = read.csv("~/scratch/bcs/results/c2b2_sub1_deg_022123.csv")
  }
}
mz.deg$one_to_one_human = gene_info$one_to_one_human[match(mz.deg$gene, gene_info$seurat_name)]

# Find Correlations ============================================================
# Find Common Gene Set
message("Finding correlations")
common.gene.set = sort(unique(toupper(mouse.deg$gene)))
common.gene.set = common.gene.set[which(common.gene.set %in% mz.deg$one_to_one_human)]
mz.common.gene.set = mz.deg$gene[match(common.gene.set, mz.deg$one_to_one_human)]
mouse.common.gene.set = common.gene.set
if (mouse.dataset != "tran") { mouse.common.gene.set = stringr::str_to_title(common.gene.set) }

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
# permCor = function(old.mat) {
#   new.mat = old.mat[sample(rownames(old.mat)),]
#   perm.cor = cor(new.mat, mouse.avg.exp.norm, method = "spearman")
#   perm.cor.melt = reshape2::melt(perm.cor)
#   return(perm.cor.melt[,3])
# }

# Permutations =================================================================
message("Performing permutations")
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
message("Performing spearman correlation test")
mz.mouse.cor.p = matrix(NA, nrow = nrow(mz.mouse.cor), ncol = ncol(mz.mouse.cor), dimnames = list(rownames(mz.mouse.cor), colnames(mz.mouse.cor)))
for (mz.clust in colnames(mz.avg.exp.norm)) {
  for (mouse.clust in colnames(mouse.avg.exp.norm)) {
    mz.mouse.cor.p[mz.clust, mouse.clust] = cor.test(mz.avg.exp.norm[,mz.clust], mouse.avg.exp.norm[,mouse.clust], method = "spearman", alternative = "greater")$p.value
  }
}
mz.mouse.cor.bon = matrix(p.adjust(mz.mouse.cor.p, method = "bonferroni"), nrow = nrow(mz.mouse.cor), ncol = ncol(mz.mouse.cor), dimnames = list(rownames(mz.mouse.cor), colnames(mz.mouse.cor)))

# Output: Plots and CSV ========================================================
message("Plotting")
mz.mouse.cor.maxed.out.melt = reshape2::melt(mz.mouse.cor)
colnames(mz.mouse.cor.maxed.out.melt) = c("mz.cluster", "mouse.cluster", "cor")
mz.mouse.cor.maxed.out.melt[, c("mouse.region", "mouse.celltype")] = reshape2::colsplit(mz.mouse.cor.maxed.out.melt$mouse.cluster, "_", c('1', '2'))
mz.mouse.cor.maxed.out.melt = mz.mouse.cor.maxed.out.melt[,c("mz.cluster", "mouse.cluster", "mouse.region", "mouse.celltype", "cor")]
if (max_overide != 0) {
  maxed.num = max_overide
} else {
  maxed.num = plyr::round_any(as.numeric(quantile(mz.mouse.cor.maxed.out.melt$cor, 0.99)), .05, f=ceiling)
  maxed.num = ifelse(maxed.num > max(mz.mouse.cor.maxed.out.melt$cor), max.num - 0.5, maxed.num)
}
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
ggplot(mz.mouse.cor.maxed.out.melt, aes(x = mouse.cluster, y = mz.cluster, fill = cor.maxed)) + geom_raster() + geom_point(data = mz.mouse.cor.maxed.out.melt[which(mz.mouse.cor.maxed.out.melt$all.sig),], size = 1.2, color = "gray60") + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), n.breaks = 6, limits = c(-maxed.num, maxed.num)) + coord_fixed() + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 10), axis.line=element_blank()) + force_panelsizes(rows = unit(nrow(mz.mouse.cor)/8, "in"), cols = unit(ncol(mz.mouse.cor)/8, "in"))
# ggplot(mz.mouse.cor.maxed.out.melt, aes(x = mz.cluster, y = mouse.cluster, fill = cor.maxed)) + geom_raster() + geom_point(data = mz.mouse.cor.maxed.out.melt[which(mz.mouse.cor.maxed.out.melt$all.sig),], size = 0.8) + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), n.breaks = 6, limits = c(-maxed.num, maxed.num)) + coord_fixed() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
cichlid_str = ifelse(isBB, "bb_", "mz_")
glut_str    = ifelse(isGlut, "glut_", "")
gaba_str    = ifelse(isGABA, "gaba_", "")
nn_str      = ifelse(isNN, "nn_", "")
sub1_str      = ifelse(isSub1, "sub1_", "")
out_name = paste0("~/scratch/bcs/results/", cichlid_str, mouse.dataset, "_cluster_", glut_str, gaba_str, nn_str, sub1_str, "cor")
if (out_name_overide != "") { out_name = out_name_overide }
col.width = 2
ggsave(paste0(out_name, ".pdf"), width = (ncol(mz.mouse.cor)/5) + col.width, height = nrow(mz.mouse.cor)/5)
write.csv(mz.mouse.cor.maxed.out.melt, paste0(out_name, ".csv"))
message(paste0("rclone copy ", paste0(out_name, ".pdf"), " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/bcs/", mouse.dataset))
system(paste0("rclone copy ", paste0(out_name, ".pdf"), " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/bcs/", mouse.dataset))
# message(paste0("rclone copy ", paste0(out_name, ".csv"), " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/bcs/", mouse.dataset))
message("Done")
