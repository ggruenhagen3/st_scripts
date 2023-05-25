processOrtholog = function(ortho) {
  ortho = ortho[order(ortho$value, decreasing = T),]
  ortho$id = paste0(ortho$mz_gene, "_", ortho$mm_gene)
  mz_top = ortho[which(!duplicated(ortho$mz_gene)),]
  mm_top = ortho[which(!duplicated(ortho$mm_gene)),]
  ortho_reciprocal_top_hit = ortho[which(ortho$id %in% mz_top$id & ortho$id %in% mm_top$id),]
  
  ortho_reciprocal_top_hit$X = NULL
  colnames(ortho_reciprocal_top_hit) = c("value", "gene1", "gene2")
  ortho_reciprocal_top_hit$gene1 = stringr::str_sub(ortho_reciprocal_top_hit$gene1, 4, 50)
  ortho_reciprocal_top_hit$gene2 = stringr::str_sub(ortho_reciprocal_top_hit$gene2, 4, 50)
  return(ortho_reciprocal_top_hit)
}

geneOverlap = function(df1, df2, cluster1, cluster2, gene_converter, strict = T, return_genes = F) {
  if (!"gene1" %in% colnames(gene_converter) | !"gene2" %in% colnames(gene_converter)) { stop("Columns gene1 and gene2 are required in gene_converter, but not found.") }
  if (strict) {
    df1$pct.dif = df1$pct.1 - df1$pct.2
    df1$strict = df1$pct.2 < 0.15 & df1$pct.dif > 0.15 & df1$p_val_adj < 1e-5
    df1 = df1[which(df1$strict),]
    df2$pct.dif = df2$pct.1 - df2$pct.2
    df2$strict = df2$pct.2 < 0.15 & df2$pct.dif > 0.15 & df2$p_val_adj < 1e-5
    df2 = df2[which(df2$strict),]
  }
  df1 = df1[which(df1$cluster == cluster1),]
  df2 = df2[which(df2$cluster == cluster2),]
  
  df1_genes = df1$gene
  df2_genes = df2$gene
  
  df1_genes2 = gene_converter$gene2[which(gene_converter$gene1 %in% df1_genes)]
  ovlp = unique(df1_genes2[which(df1_genes2 %in% df2_genes)])
  smallest_cluster = min(length(df1_genes), length(df2_genes))
  pct = length(ovlp) / smallest_cluster
  pct = ifelse(length(ovlp) == 0, 0, pct)
  
  if (return_genes) {
    df_ortho = merge(df1, gene_converter, by.x = "gene", by.y = "gene1")
    df_ortho = merge(df_ortho, df2, by.x = "gene2", by.y = "gene", suffixes = c("_1", "_2"))
    df_ortho = df_ortho[order(df_ortho$p_val_adj_1, df_ortho$p_val_adj_2),]
    df_ortho = df_ortho[which(!duplicated(df_ortho$gene2)),]
    df_ortho = df_ortho[order(df_ortho$p_val_adj_2, df_ortho$p_val_adj_1),]
    df_ortho = df_ortho[which(!duplicated(df_ortho$gene)),]
    df_ortho = df_ortho[,which(!colnames(df_ortho) %in% c("X_1", "X_2"))]
    colnames(df_ortho)[which(colnames(df_ortho) == "gene")] = "gene1"
    return(df_ortho)
  }
  return(pct)
}


geneOverlap2 = function(df1, df2, cluster1, cluster2, hgnc1, hgnc2, strict = T, return_genes = F) {
  if (strict) {
    df1$pct.dif = df1$pct.1 - df1$pct.2
    df1$strict = df1$pct.2 < 0.15 & df1$pct.dif > 0.15 & df1$p_val_adj < 1e-5
    df1 = df1[which(df1$strict),]
    df2$pct.dif = df2$pct.1 - df2$pct.2
    df2$strict = df2$pct.2 < 0.15 & df2$pct.dif > 0.15 & df2$p_val_adj < 1e-5
    df2 = df2[which(df2$strict),]
  }
  df1 = df1[which(df1$cluster == cluster1),]
  df2 = df2[which(df2$cluster == cluster2),]
  
  df1$hgnc = hgnc1$hgnc[match(df1$gene, hgnc1$gene)]
  df2$hgnc = hgnc2$hgnc[match(df2$gene, hgnc2$gene)]
  
  df1_genes = df1$hgnc[which(!is.na(df1$hgnc))]
  df2_genes = df2$hgnc[which(!is.na(df2$hgnc))]
  
  ovlp = unique(df1_genes[which(df1_genes %in% df2_genes)])
  smallest_cluster = min(length(df1_genes), length(df2_genes))
  pct = length(ovlp) / smallest_cluster
  pct = ifelse(length(ovlp) == 0, 0, pct)
  
  if (return_genes) {
    df_ortho = merge(df1[which(!is.na(df1$hgnc)),], df2[which(!is.na(df1$hgnc)),], by = "hgnc", suffixes = c("_1", "_2"))
    df_ortho = df_ortho[order(df_ortho$p_val_adj_1, df_ortho$p_val_adj_2),]
    df_ortho = df_ortho[which(!duplicated(df_ortho$gene_2)),]
    df_ortho = df_ortho[order(df_ortho$p_val_adj_2, df_ortho$p_val_adj_1),]
    df_ortho = df_ortho[which(!duplicated(df_ortho$gene_1)),]
    df_ortho = df_ortho[,which(!colnames(df_ortho) %in% c("X_1", "X_2"))]
    # colnames(df_ortho)[which(colnames(df_ortho) == "gene")] = "gene1"
    return(df_ortho)
  }
  return(pct)
}
mz.dataset = "st"
mm.dataset = "turtle"
mouse.dataset = mm.dataset

# mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/", mz.dataset, "_", mm.dataset, "_gene_ortholog2.csv"))
# mz_mm_gene_map = processOrtholog(mz_mm_gene_map)
# mz.hgnc = read.csv(paste0("~/scratch/bcs/data/cichlid_egg_many.csv"))
# if (mm.dataset == "zeisel" || mm.dataset == "oritz") {
#   mm.hgnc = read.csv(paste0("~/scratch/bcs/data/", mm.dataset, "_human_many.csv"))
#   mm.hgnc = data.frame(X = 1:nrow(mm.hgnc), gene = mm.hgnc$mouse, hgnc = mm.hgnc$human, isOne = !mm.hgnc$multiple.mouse & !mm.hgnc$multiple.human)
# } else {
#   mm.hgnc = read.csv(paste0("~/scratch/bcs/data/", mm.dataset, "_egg_many.csv"))
# }
mz.hgnc = ortho[["cichlid"]]
mm.hgnc = ortho[[mm.dataset]]
mz_deg = read.csv("~/scratch/brain/results/bb53_deg_012323.csv")
if (mz.dataset == "st") {
  mz_deg = read.csv("~/scratch/bcs/results/c2b2_brain_structure_deg_042023.csv")
}
mz_deg$cluster = stringr::str_replace(as.character(as.vector(mz_deg$cluster)), "Astro", "RG")
mz_deg$cluster = stringr::str_replace(as.character(as.vector(mz_deg$cluster)), "OB gc", "OB-gc")
mz_deg$cluster = stringr::str_replace(as.character(as.vector(mz_deg$cluster)), "OB gml", "OB-gml")

message("Loading mouse DEGs")
mouse.deg.path = list()
mouse.deg.path[["axolotl"]]    = read.csv("~/scratch/bcs/results/axolotl_cluster_markers_022123.csv")
mouse.deg.path[["bird"]]    = read.csv("~/scratch/bcs/results/bird_markers_041923.csv")
# mouse.deg.path[["oritz"]]    = read.csv("~/scratch/bcs/results/mst_cluster_markers_020723.csv")
mouse.deg.path[["oritz"]]    = read.csv("~/scratch/bcs/results/oritzg_markers_042023.csv")
mouse.deg.path[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_cluster_region_subcluster_020723.csv")
mouse.deg.path[["tasic"]]    = read.csv("~/scratch/bcs/results/tasic_subclass_020823.csv")
mouse.deg.path[["tran"]]     = read.csv("~/scratch/bcs/results/tran_broad_021523.csv")
# mouse.deg.path[["turtle"]]   = read.csv("~/scratch/bcs/results/turtle_cluster_markers_041923.csv")
mouse.deg.path[["turtle"]]   = read.csv("~/scratch/bcs/results/turtle_cp_pallial_area_markers_042623.csv")
# mouse.deg.path[["turtle_neurons"]]   = read.csv("~/scratch/bcs/results/turtle_neurons_area_markers_042023.csv")
mouse.deg.path[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_020723.csv")
# mouse.deg.path[["zeisel"]]      = read.csv("~/scratch/bcs/results/l5_cluster_markers_tax_region_020723.csv")
mouse.deg.path[["zei_yu"]]   = read.csv("~/scratch/bcs/results/zei_yu_broad_022023.csv")
mouse.deg = mouse.deg.path[[mouse.dataset]]
mm_deg = mouse.deg

mm_clusters = unique(mm_deg$cluster)
mz_clusters = unique(mz_deg$cluster)
ovlp_df = expand.grid(mm_clusters, mz_clusters)
colnames(ovlp_df) = c("mm_cluster", "mz_cluster")
ovlp_df$mz_cluster = stringr::str_replace(as.character(as.vector(ovlp_df$mz_cluster)), "Astro", "RG")
ovlp_df$id = paste0(ovlp_df$mz_cluster, "_", ovlp_df$mm_cluster)
ovlp_df$pct_loose  = unlist(mclapply(1:nrow(ovlp_df), function(x) geneOverlap2(mz_deg, mm_deg, as.vector(ovlp_df$mz_cluster[x]), as.vector(ovlp_df$mm_cluster[x]), mz.hgnc, mm.hgnc, strict = F, return_genes = F), mc.cores = 20))
ovlp_df$pct_strict = unlist(mclapply(1:nrow(ovlp_df), function(x) geneOverlap2(mz_deg, mm_deg, as.vector(ovlp_df$mz_cluster[x]), as.vector(ovlp_df$mm_cluster[x]), mz.hgnc, mm.hgnc, strict = T, return_genes = F), mc.cores = 20))

hit_sup = read.csv(paste0("~/scratch/bcs/results/", mz.dataset, "_", mm.dataset, "_sig_hits.csv"))
hit_sup$X = NULL
# hit_sup$id = paste0(hit_sup$mz.cluster, "_", hit_sup$mm.cluster)
hit_sup$id = paste0(hit_sup$mz_name, "_", hit_sup$mm_name)
# ovlp_df[,colnames(hit_sup)[3:ncol(hit_sup)]] = hit_sup[match(ovlp_df$id, hit_sup$id), colnames(hit_sup)[3:ncol(hit_sup)]]
ovlp_df[,"p0"] = hit_sup[match(ovlp_df$id, hit_sup$id), "p0"]
ovlp_df[which(is.na(ovlp_df[,"p0"])),"p0"] = FALSE

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
  # col.pal = magma(100)
  col.pal = rev(c(magma(100)[40:100], magma(100)[80:100]))
} else if (grepl("bird", mm.dataset)) {
  mm_col = "tg_cluster_orig2"
  col.pal = rev(brewer.pal(11, "PuOr")[1:6])
}

ovlp_df$pct_loose = ovlp_df$pct_loose * 100
ggplot(ovlp_df, aes(x = pct_loose, fill = p0, color = p0)) + geom_density(alpha = 0.75) + scale_color_manual(values = c(colorRampPalette(colors = col.pal)(100)[50], colorRampPalette(colors = col.pal)(100)[100])) + scale_fill_manual(values = c(colorRampPalette(colors = col.pal)(100)[30], colorRampPalette(colors = col.pal)(100)[80])) + theme_classic() + scale_x_continuous(expand = c(0,0), name = "") + scale_y_continuous(expand = expansion(mult = c(0, .05)), name = "") + theme(axis.text = element_text(size = 10))
ggsave(paste0("~/scratch/bcs/results/", mz.dataset, "_", mm.dataset, "_pct_ovlp", ".pdf"), width = 5, height = 4)# ovlp_df$sig_level = plyr::revalue(as.character(ovlp_df$bh_sig), c("TRUE" = "sig", "FALSE" = "none"))
system(paste0("rclone copy ~/scratch/bcs/results/", mz.dataset, "_", mm.dataset, "_pct_ovlp", ".pdf dropbox:BioSci-Streelman/George/Brain/spatial/analysis/samap/"))
# ovlp_df$sig_level[which(ovlp_df$p0)] = "very sig"
# ggplot(ovlp_df, aes(x = pct_loose, fill = sig_level, color = sig_level)) + geom_density(alpha = 0.6) + scale_color_manual(values = c(colorRampPalette(colors = col.pal)(100)[50], colorRampPalette(colors = col.pal)(100)[75], colorRampPalette(colors = col.pal)(100)[100])) + scale_fill_manual(values = c(colorRampPalette(colors = col.pal)(100)[30], colorRampPalette(colors = col.pal)(100)[55], colorRampPalette(colors = col.pal)(100)[80])) + theme_classic() + scale_x_continuous(expand = c(0,0), name = "") + scale_y_continuous(expand = expansion(mult = c(0, .05)), name = "") + theme(axis.text = element_text(size = 10))

# col.pal2 = RColorBrewer::brewer.pal(9, "YlOrBr")[1:7]
# ovlp_df$pct_loose = ovlp_df$pct_loose * 100
# ovlp_df_p0_p = ovlp_df
# ovlp_df_p0_p$isGABA = grepl("GABA", ovlp_df_p0_p$mz_cluster)
# ovlp_df_p0_p$isGlut = grepl("Glut", ovlp_df_p0_p$mz_cluster)
# # ovlp_df_p0_p$isNN   = ! ( ovlp_df_p0_p$isGABA | ovlp_df_p0_p$isGlut )
# # ovlp_df_p0_p = ovlp_df_p0[which( (ovlp_df_p0$isGABA | ovlp_df_p0$isGlut) & !(ovlp_df_p0$isGABA & ovlp_df_p0$isGlut) ),]
# ovlp_df_p0_p$class1 = "nn"
# ovlp_df_p0_p$class1[which(ovlp_df_p0_p$isGABA)] = "GABA"
# ovlp_df_p0_p$class1[which(ovlp_df_p0_p$isGlut)] = "Glut"
# ovlp_df_p0_p$mmIsGABA = grepl("INH", ovlp_df_p0_p$mm_cluster) | grepl("MSN", ovlp_df_p0_p$mm_cluster) | grepl("OBDOP", ovlp_df_p0_p$mm_cluster) | grepl("SZNBL", ovlp_df_p0_p$mm_cluster) | ovlp_df_p0_p$mm_cluster %in% c("OBNBL3", "OBNBL4")
# ovlp_df_p0_p$mmIsGlut = grepl("TEGLU", ovlp_df_p0_p$mm_cluster) | grepl("CR", ovlp_df_p0_p$mm_cluster) | ovlp_df_p0_p$mm_cluster %in% c("OBNBL1", "OBNBL2") | grepl("DGNBL", ovlp_df_p0_p$mm_cluster) | grepl("SEPNBL", ovlp_df_p0_p$mm_cluster) | grepl("TECHO", ovlp_df_p0_p$mm_cluster)
# # ovlp_df_p0_p$mmIsNN   = ! ( ovlp_df_p0_p$mmIsGABA | ovlp_df_p0_p$mmIsGlut )
# ovlp_df_p0_p$class2 = "nn"
# ovlp_df_p0_p$class2[which(ovlp_df_p0_p$mmIsGABA)] = "GABA"
# ovlp_df_p0_p$class2[which(ovlp_df_p0_p$mmIsGlut)] = "Glut"
# ovlp_df_p0_p$class = paste0(ovlp_df_p0_p$class1, "_", ovlp_df_p0_p$class2)
# ovlp_df_p0_p = ovlp_df_p0_p[which(ovlp_df_p0_p$class %in% c("GABA_GABA", "Glut_Glut", "nn_nn")),]
# ovlp_df_p0_p = ovlp_df_p0_p[which(ovlp_df_p0_p$p0),]
# ggplot(ovlp_df_p0_p, aes(x = pct_loose, fill = class, color = class)) + geom_density(alpha = 0.75) + scale_color_manual(values = c(colorRampPalette(colors = col.pal)(100)[20], colorRampPalette(colors = col.pal)(100)[50], colorRampPalette(colors = col.pal)(100)[100])) + scale_fill_manual(values = c(colorRampPalette(colors = col.pal)(100)[10], colorRampPalette(colors = col.pal)(100)[30], colorRampPalette(colors = col.pal)(100)[80])) + theme_classic() + scale_x_continuous(expand = c(0,0), name = "") + scale_y_continuous(expand = expansion(mult = c(0, .05)), name = "") + theme(axis.text = element_text(size = 10))
# ggsave(paste0("~/scratch/bcs/results/bb_zeisel_pct_ovlp_class", ".pdf"), width = 5, height = 2.5)

ovlp_df_axolotl$species = "axolotl"
ovlp_df_axolotl$col = "#F96D98FF"
ovlp_df_axolotl$col[which(ovlp_df_axolotl$p0)] = "#CD4071FF"
ovlp_df_axolotl$col  = factor(ovlp_df_axolotl$col, levels = c("#F96D98FF", "#CD4071FF"))
ovlp_df_turtle$species = "turtle"
ovlp_df_turtle$col = "#8CB1A6"
ovlp_df_turtle$col[which(ovlp_df_turtle$p0)] = "#004034"
ovlp_df_turtle$col = factor(ovlp_df_turtle$col, levels = c("#8CB1A6", "#004034"))
ovlp_df_bird$species = "bird"
ovlp_df_bird$col = "#FFB084"
ovlp_df_bird$col[which(ovlp_df_bird$p0)] = "#CD700E"
ovlp_df_bird$col = factor(ovlp_df_bird$col, levels = c("#FFB084", "#CD700E"))
ovlp_df_vert = rbind(ovlp_df_axolotl, ovlp_df_turtle, ovlp_df_bird)
ovlp_df_vert$species = factor(ovlp_df_vert$species, c("axolotl", "turtle", "bird"))
ggplot(ovlp_df_vert, aes(x = pct_loose, fill = col, color = col)) + geom_density(alpha = 0.75) + scale_color_identity() + scale_fill_identity() + theme_classic() + scale_x_continuous(expand = c(0,0), name = "") + scale_y_continuous(expand = expansion(mult = c(0, .01)), name = "") + theme(axis.text = element_text(size = 10)) + facet_wrap(~ species, ncol = 1)
ggsave("~/scratch/bcs/results/vert_pct_ovlp_dark.pdf", width = 6, height = 10)# ovlp_df$sig_level = plyr::revalue(as.character(ovlp_df$bh_sig), c("TRUE" = "sig", "FALSE" = "none"))
system("rclone copy ~/scratch/bcs/results/vert_pct_ovlp_dark.pdf dropbox:BioSci-Streelman/George/Brain/spatial/analysis/samap/")

ovlp_df_turtle_neurons$species = "turtle"
ovlp_df_turtle_neurons$col = "#8CB1A6"
ovlp_df_turtle_neurons$col[which(ovlp_df_turtle_neurons$p0)] = "#004034"
ovlp_df_turtle_neurons$col = factor(ovlp_df_turtle_neurons$col, levels = c("#8CB1A6", "#004034"))
ovlp_df_oritz$species = "oritz"
ovlp_df_oritz$col = "#FDA284"
ovlp_df_oritz$col[which(ovlp_df_oritz$p0)] = "#8F0813"
ovlp_df_oritz$col = factor(ovlp_df_oritz$col, levels = c("#FDA284", "#8F0813"))
ovlp_df_vert = rbind(ovlp_df_oritz, ovlp_df_turtle_neurons)
ovlp_df_vert$species = factor(ovlp_df_vert$species, levels = c("turtle", "oritz"))
ggplot(ovlp_df_vert, aes(x = pct_loose, fill = col, color = col)) + geom_density(alpha = 0.75) + scale_color_identity() + scale_fill_identity() + theme_classic() + scale_x_continuous(expand = c(0,0), name = "") + scale_y_continuous(expand = expansion(mult = c(0, .01)), name = "") + theme(axis.text = element_text(size = 10)) + facet_wrap(~ species, ncol = 1)
ggsave(paste0("~/scratch/bcs/results/st_pct_ovlp_dark.pdf"), width = 3.5, height = 5)# ovlp_df$sig_level = plyr::revalue(as.character(ovlp_df$bh_sig), c("TRUE" = "sig", "FALSE" = "none"))
system("rclone copy ~/scratch/bcs/results/st_pct_ovlp_dark.pdf dropbox:BioSci-Streelman/George/Brain/spatial/analysis/samap/")

se <- function(x) sd(x)/sqrt(length(x))
ovlp_df = ovlp_df_oritz
this_p0 = ovlp_df$pct_loose[which(ovlp_df$p0)]
this_not_p0 = ovlp_df$pct_loose[which(!ovlp_df$p0)]
t.test(this_p0, this_not_p0)$p.value
t.test(this_p0, this_not_p0)$statistic
mean(this_p0)
se(this_p0)
mean(this_not_p0)
se(this_not_p0)

# Gene Heatmap 
message("Loading mouse object")
rds.path = list()
rds.path[["oritz"]]    = "~/scratch/bcs/data/mst_norm.rds"
rds.path[["oritzb"]]   = "~/scratch/bcs/data/oritz_b_raw.rds"
rds.path[["saunders"]] = "~/scratch/bcs/data/mouse_w_pc_down_norm.rds"
rds.path[["tasic"]]    = "~/scratch/bcs/data/tasic_norm_020823.rds"
rds.path[["tran"]]     = "~/scratch/bcs/data/tran_norm.rds"
rds.path[["turtle"]]   = "/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/turtle_neurons_norm.rds"
rds.path[["zeisel"]]   = "~/scratch/bcs/data/l5_tel_norm.rds"
rds.path[["zei_down"]] = "~/scratch/bcs/data/zei_down_norm_022823.rds"
rds.path[["zei_yu"]]   = "~/scratch/bcs/data/zei_yu_norm.rds"
mouse.path = rds.path[[mm.dataset]]
mouse = readRDS(mouse.path)

message("Setting mouse identity")
mouse.ident = list()
mouse.ident = switch(mm.dataset,
                     "oritz" = mouse$ABA_parent,
                     "oritzb" = mouse$b_parent,
                     "saunders" = paste0(mouse$region, "_", mouse$subcluster),
                     "tasic" = mouse$subclass,
                     "tran" = Idents(mouse),
                     "turtle" = mouse$clusters,
                     "zeisel" = mouse$ClusterName,
                     "zei_down" = mouse$ClusterName,
                     "zei_yu" = mouse$my.specific)
Idents(mouse) = mouse.ident

message("Loading cichlid data")
gene_info = read.csv("~/scratch/m_zebra_ref/gene_info_3.csv")
convert15 = read.csv("~/scratch/st/data/convert15.csv")
convert53 = read.csv("~/scratch/st/data/convert53.csv")
mz = readRDS("~/scratch/brain/data/bb_demux_102021.rds")
mz$good_names53 = factor(convert53$new[match(mz$seurat_clusters, convert53$old)], levels = rev(convert53$new))
Idents(mz) = mz$good_names53

ovlp_df = ovlp_df_turtle_neurons
ovlp_df_p0 = ovlp_df[which(ovlp_df$p0),]
# ovlp_df_p0 = ovlp_df_p0[order(ovlp_df_p0$cor, decreasing = T),]
# ovlp_df_p0 = ovlp_df_p0[which(!duplicated(ovlp_df_p0$mz_cluster)),]
# ovlp_df_p0 = ovlp_df_p0[which(!duplicated(ovlp_df_p0$mm_cluster)),]
res = mclapply(1:nrow(ovlp_df_p0), function(x) geneOverlap2(mz_deg, mm_deg, as.vector(ovlp_df_p0$mz_cluster[x]), as.vector(ovlp_df_p0$mm_cluster[x]), mz.hgnc, mm.hgnc, strict = F, return_genes = T), mc.cores = 20)
df_ortho = do.call('rbind', res)
df_ortho_write = df_ortho
colnames(df_ortho_write) = c("hgnc", "")
ovlp_df_p0$mz_cluster = stringr::str_replace(as.character(as.vector(ovlp_df_p0$mz_cluster)), "Astro", "RG")

df_ortho$id = paste0(df_ortho$cluster_1, "_", df_ortho$cluster_2)
df_ortho_backup = df_ortho
# df_ortho = df_ortho[which(df_ortho$id %in% c("2.2_Oligo_OEC", "4.1_GABA_MSN1", "7_GABA_MSN5", "15.1_GABA/Glut_TEINH1", "15.3_GABA_TEINH21", "8.1_Glut_TEGLU6", "8.9_Glut_TEGLU23", "9.5_Glut_DGNBL1", "15.6_Glut_OBNBL2")),]
df_ortho = df_ortho[which(df_ortho$id %in% c("OB-gc_MOB", "Vv_LSX", "Vd-r_ACB", "Dl-vv_CA3", "Dl-g_VIS")),]
# df_ortho = df_ortho[which(df_ortho$id %in% c("Dl-vv_DMC", "Dp_aLC", "Dm-2r_pDVR")),]

mz = ScaleData(mz, features = unique(df_ortho$gene1))
mz.assay.to.use = "RNA"
mz.exp = AverageExpression(mz, features = unique(df_ortho$gene1), assays = mz.assay.to.use, slot = "scale.data")[[1]]
mz.exp.order = mz.exp[df_ortho$gene1, unique(df_ortho$cluster_1)]
mz.exp.order.id = mz.exp.order
rownames(mz.exp.order.id) = 1:nrow(mz.exp.order.id)
mz.exp.order.df = reshape2::melt(mz.exp.order.id)
mz.exp.order.df$Var2 = factor(mz.exp.order.df$Var2, levels = unique(df_ortho$cluster_1))
# mz.exp.order.df$Var1 = factor(rev(mz.exp.order.df$Var1))

mouse = ScaleData(mouse, features = unique(df_ortho$gene2))
mouse.assay.to.use = "SCT"
mouse.exp = AverageExpression(mouse, features = unique(df_ortho$gene2), assays = mouse.assay.to.use, slot = "scale.data")[[1]]
mouse.exp.order = mouse.exp[df_ortho$gene2, unique(df_ortho$cluster_2)]
mouse.exp.order.id = mouse.exp.order
rownames(mouse.exp.order.id) = 1:nrow(mouse.exp.order.id)
mouse.exp.order.df = reshape2::melt(mouse.exp.order.id)
mouse.exp.order.df$Var2 = factor(mouse.exp.order.df$Var2, levels = unique(df_ortho$cluster_2))
# mouse.exp.order.df$Var1 = factor(rev(mouse.exp.order.df$Var1))

library(scales)
col.pal2 = RColorBrewer::brewer.pal(9, "YlOrBr")[1:7]
pdf('~/scratch/bcs/results/bb_zeisel_heatmap_exp.pdf', width = nrow(mz.exp.order   )/60, height = (ncol(mz.exp.order)+ncol(mouse.exp.order))/4)
p1 = ggplot(mz.exp.order.df,    aes(x = Var1, y = Var2, fill = value)) + geom_raster() + scale_fill_gradientn(colors = col.pal2, limits = c(0, 2.5), oob=squish) + theme_classic() + scale_y_discrete(name = "", expand=c(0,0)) + scale_x_discrete(name = "", expand=c(0,0)) + force_panelsizes(rows = unit(ncol(mz.exp.order   )/6, "in"), cols = unit(nrow(mz.exp.order   )/110, "in")) + theme(axis.line=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
p2 = ggplot(mouse.exp.order.df, aes(x = Var1, y = Var2, fill = value)) + geom_raster() + scale_fill_gradientn(colors = col.pal,  limits = c(0, 2.5), oob=squish) + theme_classic() + scale_y_discrete(name = "", expand=c(0,0)) + scale_x_discrete(name = "", expand=c(0,0)) + force_panelsizes(rows = unit(ncol(mouse.exp.order)/6, "in"), cols = unit(nrow(mouse.exp.order)/110, "in")) + theme(axis.line=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
cowplot::plot_grid(plotlist = list(p1, p2), ncol = 1)
dev.off()

col.pal2 = RColorBrewer::brewer.pal(9, "YlOrBr")[1:7]
pdf('~/scratch/bcs/results/bb_zeisel_heatmap_exp_t.pdf', width = 2+(ncol(mz.exp.order)+ncol(mouse.exp.order))/4, height = nrow(mz.exp.order   )/80)
p1 = ggplot(mz.exp.order.df,    aes(x = Var2, y = Var1, fill = value)) + geom_raster() + scale_fill_gradientn(colors = col.pal2, limits = c(0, 2.5), oob=squish) + theme_classic() + scale_y_discrete(name = "", expand=c(0,0)) + scale_x_discrete(name = "", expand=c(0,0), position = "top") + force_panelsizes(rows = unit(nrow(mz.exp.order   )/110, "in"), cols = unit(ncol(mz.exp.order   )/6, "in")) + theme(axis.line=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(vjust = -1, hjust = 0, angle = 45, size = 10))
p2 = ggplot(mouse.exp.order.df, aes(x = Var2, y = Var1, fill = value)) + geom_raster() + scale_fill_gradientn(colors = col.pal,  limits = c(0, 2.5), oob=squish) + theme_classic() + scale_y_discrete(name = "", expand=c(0,0)) + scale_x_discrete(name = "", expand=c(0,0), position = "top") + force_panelsizes(rows = unit(nrow(mouse.exp.order)/110, "in"), cols = unit(ncol(mouse.exp.order)/6, "in")) + theme(axis.line=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(vjust = -1, hjust = 0, angle = 45, size = 10))
cowplot::plot_grid(plotlist = list(p1, p2), ncol = 2)
dev.off()

# mz.avg.exp = AverageExpression(mz, features = mz.common.gene.set, assays = mz.assay.to.use, slot = "data")[[1]]
# mz.avg.exp.norm = log(mz.avg.exp+1) + 0.1
# mz.avg.exp.norm = mz.avg.exp.norm / rowMeans(mz.avg.exp.norm)
# 
# mouse.avg.exp = AverageExpression(mouse, features = mouse.common.gene.set, assays = "SCT", slot = "data")[[1]]
# mouse.avg.exp.norm = log(mouse.avg.exp+1) + 0.1
# mouse.avg.exp.norm = mouse.avg.exp.norm / rowMeans(mouse.avg.exp.norm)

# Plotting
mz.dataset = "bb"
mm.dataset = "zeisel"
mz_col = "mz_good_names"
mm_col = "mm_g_parent"
samc_folder = "~/scratch/bcs/samc/"
meta = as.matrix(read.csv(paste0(samc_folder, mz.dataset, "_", mm.dataset, ".csv"),   row.names = 1))
obj_fpath = paste0("~/scratch/bcs/data/mz_mm_", mz.dataset, "_oritiz_b")
obj_fpath = paste0("~/scratch/bcs/data/bb_zeisel")
Convert(paste0(obj_fpath, "_3.h5ad"), dest = "h5seurat", overwrite = TRUE)
merged = LoadH5Seurat(paste0(obj_fpath, "_3.h5seurat"), meta.data = FALSE, misc = FALSE)
merged@meta.data = cbind(merged@meta.data, meta)
merged@meta.data[which(merged@meta.data[,mz_col] == "unassigned"), mz_col] = NA
merged@meta.data[which(merged@meta.data[,mm_col] == "unassigned"), mm_col] = NA

# col.pal = RColorBrewer::brewer.pal(9, "Greens")
col.pal = RColorBrewer::brewer.pal(9, "Blues")
yellows = c("#E06C00", "#ff8800", "#ffa200", "#ffaa00", "#ffb700", "#ffc300", "#FFDD00", "#ffea00")
mz_merged = subset(merged, cells = colnames(merged)[which(merged$species == "mz")])
mm_merged = subset(merged, cells = colnames(merged)[which(merged$species == "mm")])
Idents(mm_merged) = mm_merged$leiden_clusters
Idents(mz_merged) = mz_merged$leiden_clusters
num_clusters = max(as.numeric(merged$leiden_clusters))+1
png("~/scratch/bcs/results/st_oritzg_split.png", width = 6400, height = 3200, res = 600)
p1 = DimPlot(mz_merged, label = T, repel = F, label.size = 3, cols = colorRampPalette(colors = yellows)(num_clusters)) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend()
p2 = DimPlot(mm_merged, label = T, repel = F, label.size = 3, cols = colorRampPalette(colors = col.pal)(num_clusters)) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend()
print(cowplot::plot_grid(plotlist = list(p1, p2)))
dev.off()
# Cairo::Cairo("~/scratch/bcs/results/mz_mm_zei_split.png", width = 3200, height = 1600, res = 300)
# print(DimPlot(merged, split.by = "species", label = F, cols = c("goldenrod1", colorRampPalette(col.pal)(100)[80])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend())
# dev.off()

Idents(merged) = merged$leiden_clusters
pdf("~/scratch/bcs/results/st_oritzg_umap_clusters.pdf", width = 7, height = 7)
print(DimPlot(merged, order = T, label = T) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
dev.off()

Idents(merged) = merged@meta.data[,mz_col]
pdf("~/scratch/bcs/results/st_oritzg_umap_mz.pdf", width = 7, height = 7)
print(DimPlot(merged, order = T, label = T) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
dev.off()

Idents(merged) = merged@meta.data[,mm_col]
pdf("~/scratch/bcs/results/st_oritzg_umap_mm.pdf", width = 7, height = 7)
print(DimPlot(merged, order = T, label = T) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
dev.off()

# Calculate mouse-specific clusters
mz_species = meta[which(meta[,mz_col] != "unassigned")[1], "species"]
mm_species = meta[which(meta[,mm_col] != "unassigned")[1], "species"]
mm_over_mm_cluster = unclass(table(meta[,"species"], meta[,"leiden_clusters"]))
mm_over_mm_cluster = mm_over_mm_cluster[mm_species,] / mm_over_mm_cluster[mz_species,]
mm_over_mz = length(which(meta[,mm_col] != "unassigned")) / length(which(meta[,mz_col] != "unassigned"))
relative_prop = (mm_over_mm_cluster / mm_over_mz) / ((mm_over_mm_cluster / mm_over_mz) + 1)
df_prop = rbind(data.frame(prop = relative_prop, species = 'mm', cluster = names(relative_prop)), 
                data.frame(prop = 1-relative_prop, species = 'mz', cluster = names(relative_prop)))
prop_thresh  = .75
prop_thresh2 = .75
mm_specific_joint = names(relative_prop)[which(relative_prop > prop_thresh)]
meta = cbind(meta, mmSpecificJoint=F)
meta[which(meta[,"species"] == 'mm' & meta[,"leiden_clusters"] %in% mm_specific_joint), "mmSpecificJoint"] = T
pct_mm_in_mm_specific_joint = unclass(table(meta[,mm_col], meta[,"mmSpecificJoint"]))
pct_mm_in_mm_specific_joint = pct_mm_in_mm_specific_joint[,"TRUE"] / rowSums(pct_mm_in_mm_specific_joint)
pct_mm_in_mm_specific_joint = pct_mm_in_mm_specific_joint[which(pct_mm_in_mm_specific_joint > prop_thresh2)]
mm_specific_mm = names(pct_mm_in_mm_specific_joint)

mz_specific_joint = names(relative_prop)[which(relative_prop < (1-prop_thresh))]
meta = cbind(meta, mzSpecificJoint=F)
meta[which(meta[,"species"] == 'mz' & meta[,"leiden_clusters"] %in% mz_specific_joint), "mzSpecificJoint"] = T
pct_mz_in_mz_specific_joint = unclass(table(meta[,mz_col], meta[,"mzSpecificJoint"]))
pct_mz_in_mz_specific_joint = pct_mz_in_mz_specific_joint[,"TRUE"] / rowSums(pct_mz_in_mz_specific_joint)
pct_mz_in_mz_specific_joint = pct_mz_in_mz_specific_joint[which(pct_mz_in_mz_specific_joint > prop_thresh2)]
mz_specific_mz = names(pct_mz_in_mz_specific_joint)

if (length(mz_specific_mz) > 0) {
  species_specific_df_mz = data.frame(cluster = mz_specific_mz, species = "mz")
  write.csv(species_specific_df_mz, paste0("~/scratch/bcs/samc/", mz.dataset, "_", mm.dataset, "_species_specific.csv"))
}
if (length(mm_specific_mm) > 0) {
  species_specific_df_mm = data.frame(cluster = mm_specific_mm, species = "mm")
  write.csv(species_specific_df_mm, paste0("~/scratch/bcs/samc/", mz.dataset, "_", mm.dataset, "_species_specific.csv"))
}
if (length(mz_specific_mz) > 0 && length(mm_specific_mm) > 0) {
  species_specific_df = rbind(species_specific_df_mz, species_specific_df_mm)
  write.csv(species_specific_df, paste0("~/scratch/bcs/samc/", mz.dataset, "_", mm.dataset, "_species_specific.csv"))
}

# Core Conserved Celltype Markers ==============================================
vert = sort(c("mouse", "axolotl", "bird", "turtle", "cichlid"))
deg = list()
deg[["mouse"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_020723.csv")
deg[["axolotl"]] = read.csv("~/scratch/bcs/results/axolotl_cluster_markers_022123.csv")
deg[["turtle"]]  = read.csv("~/scratch/bcs/results/turtle_cluster_markers_041923.csv")
deg[["bird"]]    = read.csv("~/scratch/bcs/results/bird_markers_042823.csv")
deg[["cichlid"]] = read.csv("~/scratch/brain/results/bb53_deg_012323.csv")
deg[["cichlid"]]$cluster = stringr::str_replace(as.character(as.vector(deg[["cichlid"]]$cluster)), "Astro", "RG")

ortho = list()
for (v1 in vert[1:length(vert)]) {
  if (v1 == "mouse") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/mouse_human_many.csv"))
    mz_mm_gene_map = mz_mm_gene_map[,c("gene", "hgnc")]
  } else if (v1 == "turtle") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/turtle_ortho_final.csv"))
    mz_mm_gene_map$hgnc = mz_mm_gene_map$hgnc.many
  } else if (v1 == "cichlid") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/gene_info_4.csv"))
    mz_mm_gene_map$hgnc = mz_mm_gene_map$human_reasonable
    mz_mm_gene_map$gene = mz_mm_gene_map$seurat_name
  } else if (v1 == "bird") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/bird_ortho_final.csv"))
    mz_mm_gene_map$hgnc = mz_mm_gene_map$hgnc.many
  } else if (v1 == "axolotl") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/", v1, "_egg_many.csv"))  
  }
  mz_mm_gene_map = mz_mm_gene_map[,c("gene", "hgnc")]
  mz_mm_gene_map$hgnc[which(mz_mm_gene_map$hgnc == "")] = NA
  ortho[[v1]] = mz_mm_gene_map
}
ortho[["zeisel"]] = ortho[["mouse"]]
ortho[["saunders"]] = ortho[["mouse"]]
ortho[["oritz"]] = ortho[["mouse"]]
length(which(ortho[["mouse"]]$hgnc %in% ortho[["axolotl"]]$hgnc & ortho[["mouse"]]$hgnc %in% ortho[["turtle"]]$hgnc & ortho[["mouse"]]$hgnc %in% ortho[["bird"]]$hgnc & ortho[["mouse"]]$hgnc %in% ortho[["cichlid"]]$hgnc))

ortho2 = list()
for (v1 in vert) {
  if (v1 == "cichlid") {
    ortho2[[v1]] = ortho[[v1]][which(ortho[[v1]]$gene %in% rownames(v[[v1]]@assays$RNA@counts)),]
  } else {
    ortho2[[v1]] = ortho[[v1]][which(ortho[[v1]]$gene %in% rownames(v[[v1]]@assays$SCT@counts)),]
  }
  ortho2[[v1]] = ortho2[[v1]][which(!duplicated(paste0(ortho2[[v1]]$gene, "_", ortho2[[v1]]$hgnc))),]
  ortho2[[v1]]$dup = duplicated(ortho2[[v1]]$gene) | duplicated(ortho2[[v1]]$gene, fromLast = TRUE)
  ortho2[[v1]]$dup[which(ortho2[[v1]]$dup & ortho2[[v1]]$hgnc == toupper(ortho2[[v1]]$gene))] = FALSE
  ortho2[[v1]] = ortho2[[v1]][which(!ortho2[[v1]]$dup),]
}
common_gene = ortho2[["mouse"]]$hgnc[which(ortho2[["mouse"]]$hgnc %in% ortho2[["axolotl"]]$hgnc & ortho2[["mouse"]]$hgnc %in% ortho2[["turtle"]]$hgnc & ortho2[["mouse"]]$hgnc %in% ortho2[["bird"]]$hgnc & ortho2[["mouse"]]$hgnc %in% ortho2[["cichlid"]]$hgnc)]
length(common_gene)

# common_name = c("OLIG", "OB", "MGE1", "MGE2", "DGNBL", "CA3")
# # common_name = c("OLIG", "OB", "MSN", "MGE1", "MGE2", "DGNBL", "CA3")
# # common[[common_name[1]]] = data.frame(name = common_name[1], celltype = c("ACTE1", "ACTE2", "EPEN0", "EPEN2", "EPEN6", "EPEN8", "EPEN9", "EPEN10"),
# #                                       species = c("mouse", "mouse", "axolotl", "turtle", "bird", "cichlid"))
# common = list()
# common[[common_name[1]]] = data.frame(name = common_name[1], celltype = c("OEC", "OLIG15", "tsOlig", "Oligo", "2.2_Oligo"),
#                                       species = c("mouse",  "axolotl", "turtle", "bird", "cichlid"))
# common[[common_name[2]]] = data.frame(name = common_name[2], celltype = c("OBDOP2", "GABA1", "GABA3", "GABA10", "i01", "GABA-1-1", "5.1_GABA", "5.2_GABA"),
#                                       species = c("mouse",  "axolotl",  "axolotl",  "axolotl", "turtle", "bird", "cichlid", "cichlid"))
# # common[[common_name[3]]] = data.frame(name = common_name[3], celltype = c("MSN1", "MSN2", "GABA11", "i05", "i06", "MSN3", "4.1_GABA", "4.2_GABA", "4.4_GABA", "4.7_GABA"),
# #                                       species = c("mouse",  "mouse",  "axolotl",  "turtle", "turtle", "bird", "cichlid", "cichlid", "cichlid", "cichlid"))
# common[[common_name[3]]] = data.frame(name = common_name[3], celltype = c("TEINH17", "TEINH18", "GABA2", "GABA6", "i07", "i11", "i12", "i13", "GABA-3", "GABA-4", "6_GABA"),
#                                       species = c("mouse",  "mouse",  "axolotl",  "axolotl",  "turtle", "turtle", "turtle", "turtle", "bird", "bird", "cichlid"))
# common[[common_name[4]]] = data.frame(name = common_name[4], celltype = c("TEINH19", "TEINH21", "GABA17", "i08", "i09", "i10", "GABA-2", "15.3_GABA"),
#                                       species = c("mouse",  "mouse",  "axolotl", "turtle", "turtle", "turtle", "bird", "cichlid"))
# common[[common_name[5]]] = data.frame(name = common_name[5], celltype = c("DGNBL1", "DGNBL2", "NB4", "NB7", "NB8", "NB1", "NB3", "NB11", "tsNPCs", "Pre-2", "9.5_Glut"),
#                                       species = c("mouse",  "mouse",  "axolotl",  "axolotl",  "axolotl",  "axolotl",  "axolotl",  "axolotl", "turtle", "bird", "cichlid"))
# common[[common_name[6]]] = data.frame(name = common_name[6], celltype = c("TEGLU23", "GLUT7", "e34", "HVC_Glut-3", "8.5_Glut", "8.9_Glut"),
#                                       species = c("mouse", "axolotl",  "turtle", "bird", "cichlid", "cichlid"))

common_name = c("OLIG", "OB", "MGE1", "MGE2", "DGNBL", "CA3")
common = list()
common[[common_name[1]]] = data.frame(name = common_name[1], celltype = c("OEC", "OLIG15", "tsOlig", "Oligo", "2.2_Oligo"),
                                      species = c("mouse",  "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[2]]] = data.frame(name = common_name[2], celltype = c("OBDOP2", "GABA3", "i01", "GABA-1-1", "5.2_GABA"),
                                      species = c("mouse",  "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[3]]] = data.frame(name = common_name[3], celltype = c("TEINH17", "GABA2", "i07", "GABA-3", "6_GABA"),
                                      species = c("mouse", "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[4]]] = data.frame(name = common_name[4], celltype = c("TEINH21", "GABA17", "i08", "GABA-2", "15.3_GABA"),
                                      species = c("mouse", "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[5]]] = data.frame(name = common_name[5], celltype = c("DGNBL1", "NB1", "tsNPCs", "Pre-2", "9.5_Glut"),
                                      species = c("mouse",  "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[6]]] = data.frame(name = common_name[6], celltype = c("TEGLU23", "GLUT7", "e34", "HVC_Glut-3", "8.9_Glut"),
                                      species = c("mouse", "axolotl",  "turtle", "bird", "cichlid"))
common_celltype = do.call('rbind', common)
# Maybe MSN's should just 4.1 and 4.2_GABA

common_celltype_genes = data.frame()
common_celltype_genes4 = data.frame()
for (name_i in common_name) {
  cluster1 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == "cichlid")]
  for (v1 in vert[which(vert != "cichlid")]) {
    print(v1)
    cluster2 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == v1)]
    
    df1 = deg[["cichlid"]]
    df2 = deg[[v1]]
    df1 = df1[which(df1$cluster %in% cluster1),]
    df2 = df2[which(df2$cluster %in% cluster2),]
    
    df1$hgnc = ortho[["cichlid"]]$hgnc[match(df1$gene, ortho[["cichlid"]]$gene)]
    df2$hgnc = ortho[[v1]]$hgnc[match(df2$gene, ortho[[v1]]$gene)]
    
    ovlp = unique(df1$hgnc[which(df1$hgnc %in% df2$hgnc & !is.na(df1$hgnc) & df1$hgnc != "")])
    if (length(ovlp) > 0) {
      common_celltype_genes = rbind(common_celltype_genes, data.frame(name = name_i, species = v1, gene = ovlp))
    }
  }
  gene_rep = table(common_celltype_genes$gene[which(common_celltype_genes$name == name_i)])
  print(table(gene_rep))
  if (length(which(gene_rep == 4)) > 0) { common_celltype_genes4 = rbind(common_celltype_genes4, data.frame(name = name_i, gene = names(gene_rep)[which(gene_rep == 4)])) }
}
write.csv(common_celltype_genes4, "~/scratch/bcs/results/vert_cons_genes.csv")

common_pairwise = data.frame()
for (name_i in common_name[which(common_name != "MSN")]) {
  for (v1 in vert) {
    cluster1 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == v1)]
    for (v2 in vert[1:(which(vert == v1)-1)]) {
      cluster2 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == v2)]
      
      df1 = deg[[v1]]
      df2 = deg[[v2]]
      df1 = df1[which(df1$cluster %in% cluster1),]
      df2 = df2[which(df2$cluster %in% cluster2),]
      
      df1$hgnc = ortho[[v1]]$hgnc[match(df1$gene, ortho[[v1]]$gene)]
      df2$hgnc = ortho[[v2]]$hgnc[match(df2$gene, ortho[[v2]]$gene)]
      
      ovlp = unique(df1$hgnc[which(df1$hgnc %in% df2$hgnc & !is.na(df1$hgnc))])
      min_num = min(length(df1$hgnc[which(!is.na(df1$hgnc))]), length(df2$hgnc[which(!is.na(df2$hgnc))]))
      pct = length(ovlp) / min_num
      common_pairwise = rbind(common_pairwise, data.frame(v1 = v1, v2 = v2, name = name_i, pct = pct, num = length(ovlp)))
    }
  }
}
library("scales")
library("viridis")
common_pairwise = common_pairwise[which(common_pairwise$v1 != common_pairwise$v2),]
pdf(paste0("~/scratch/bcs/results/vert_ovlp_pct2.pdf"), width = 8, height = 8)
# ggplot(common_pairwise, aes(x = v1, y = v2, size = pct, color = num)) + geom_point() + facet_wrap(~ name, ncol = 2) + theme_classic() + theme(axis.text.x = element_text(size = 10)) + coord_fixed() + scale_color_viridis(limits=c(0,250), breaks = c(0, 250), oob=squish) + scale_size_continuous(limits=c(0,1))
# ggplot(common_pairwise, aes(x = v1, y = v2, size = pmax(pmin(num, 250),0), color = pct)) + geom_point() + facet_wrap(~ name, ncol = 2) + theme_classic() + theme(axis.text.x = element_text(size = 10)) + coord_fixed() + scale_color_viridis(limits=c(0,0.5), breaks = c(0, 0.5), oob=squish) + scale_size_continuous(limits=c(0,250), breaks = c(0,250))
ggplot(common_pairwise, aes(x = v1, y = v2, size = pmax(pmin(num, 250),0), color = pct)) + geom_point() + facet_wrap(~ name, ncol = 2) + theme_classic() + theme(axis.text.x = element_text(size = 10)) + coord_fixed() + scale_color_viridis() + xlab("") + ylab("")
dev.off()

common_celltype_genes_stats = data.frame()
for (name_i in common_name[which(common_name != "MSN")]) {
  for (v1 in vert) {
    cluster1 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == v)]
    
    df1 = deg[[v1]]
    df1 = df1[which(df1$cluster %in% cluster1),]
    df1$hgnc = ortho[[v1]]$hgnc[match(df1$gene, ortho[[v1]]$gene)]
    
    high_cons_genes = common_celltype_genes4$gene[which(common_celltype_genes4$name == name_i)]
    df1 = df1[which(df1$hgnc %in% high_cons_genes),]
    df1 = df1[order(df1$p_val_adj, decreasing = F),]
    df1 = df1[which(!duplicated(df1$hgnc)),]
    common_celltype_genes_stats = rbind(common_celltype_genes_stats, data.frame(name = name_i, species = v1, gene = df1$hgnc, p_val_adj = df1$p_val_adj, log2FC = df1$avg_log2FC))
  }
}
neg_log_bon_thresh = 300
common_celltype_genes_stats$neg_log_bon = -log10(common_celltype_genes_stats$p_val_adj)
common_celltype_genes_stats$neg_log_bon[which(common_celltype_genes_stats$neg_log_bon > neg_log_bon_thresh)] = neg_log_bon_thresh
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "cichlid")] = "goldenrod1"
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "mouse")] = "#006d2c"
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "axolotl")] = "#bb3665"
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "turtle")] = "#004034"
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "bird")] = "#CD700E"
pdf(paste0("~/scratch/bcs/results/vert_ovlp_genes.pdf"), width = 8, height = 8)
ggplot(common_celltype_genes_stats, aes(x = neg_log_bon, y = gene, color = species, size = log2FC)) + geom_point() + scale_color_identity() + theme_bw() + facet_grid(name ~ ., scales = "free_y", space = "free_y") + theme(strip.text.x = element_blank()) + scale_x_continuous(limits=c(0, neg_log_bon_thresh+10))
dev.off()

# Common Gene Patterns =========================================================
v = list()
v[["cichlid"]] = readRDS("~/scratch/brain/data/bb_demux_102021.rds")
v[["mouse"]]   = readRDS("~/scratch/bcs/data/l5_tel_norm.rds")
v[["axolotl"]] = readRDS("~/scratch/bcs/data/axolotl_norm_w_cluster.rds")
v[["turtle"]]  = readRDS("~/scratch/bcs/data/turtle_norm_cp_detail.rds")
v[["bird"]]    = readRDS("~/scratch/bcs/data/bird_smallx_033023_norm.rds")

convert53 = read.csv("~/scratch/st/data/convert53.csv")
v[["cichlid"]]$good_names53 = factor(convert53$new[match(v[["cichlid"]]$seurat_clusters, convert53$old)], levels = rev(convert53$new))
# turtle_meta = read.csv("~/scratch/bcs/data/turtle_meta.csv")
# v[["turtle"]]$cp_pallial_area = turtle_meta$cp_pallial_area

Idents(v[["cichlid"]]) = "good_names53"
Idents(v[["mouse"]])   = v[["mouse"]]$ClusterName
Idents(v[["axolotl"]]) = v[["axolotl"]]$cluster
Idents(v[["turtle"]])  = v[["turtle"]]$cp_detail
Idents(v[["bird"]])    = v[["bird"]]$cluster_orig2

findSpeciesFC = function(i) {
  c1 = common_celltype$name[i]
  v1 = common_celltype$species[i]
  this_ident = common_celltype$celltype[i]
  print(paste0('common celltype = ', c1, ', species = ', v1, ", species celltype = ", this_ident))
  
  this_genes = ortho2[[v1]]$gene[which(ortho2[[v1]]$hgnc %in% common_gene)]
  this_fc = FoldChange(v[[v1]], features = this_genes, ident.1 = this_ident)
  this_fc$gene = rownames(this_fc)
  this_fc$hgnc = ortho2[[v1]]$hgnc[match(this_fc$gene, ortho2[[v1]]$gene)]
  this_fc$name = c1
  this_fc$species = v1
  this_fc$ident = this_ident
  print(paste0('*DONE* common celltype = ', c1, ', species = ', v1, ", species celltype = ", this_ident))
  return(this_fc)
}

v_common_fc = parallel::mclapply(1:nrow(common_celltype), function(x) findSpeciesFC(x), mc.cores = 1)
fc_df = do.call('rbind', v_common_fc)

hgnc_counts = fc_df[which(fc_df$name == "OLIG"),]
hgnc_counts = table(hgnc_counts$hgnc, hgnc_counts$species)
one_hgnc   = rownames(hgnc_counts)[which(rowSums(hgnc_counts) == 5)] 
multi_hgnc = rownames(hgnc_counts)[which(rowSums(hgnc_counts) != 5 & rowSums(hgnc_counts > 0) == 5)]
one_common = data.frame()
for (c1 in common_name) {
  this_one = data.frame(hgnc = c(one_hgnc, multi_hgnc), isOne = c(rep("TRUE", length(one_hgnc)), rep("FALSE", length(multi_hgnc))))
  for (v1 in vert) {
    this_df = fc_df[which(fc_df$name == c1 & fc_df$species == v1 & fc_df$hgnc %in% c(one_hgnc, multi_hgnc)),]
    multi_df = this_df[order(this_df$avg_log2FC, decreasing = T),]
    multi_df = multi_df[which(multi_df$hgnc %in% c(one_hgnc, multi_hgnc)),]
    multi_df = multi_df[which(!duplicated(multi_df$hgnc)),]
    multi_df = multi_df[match(c(one_hgnc, multi_hgnc), multi_df$hgnc),]
    this_one = cbind(this_one, multi_df$avg_log2FC)
  }
  colnames(this_one) = c("hgnc", "isOne", vert)
  this_one$name = c1
  one_common = rbind(one_common, this_one)
}
one_common$same_sign = rowSums(one_common[,vert] > 0)
table(one_common$same_sign, one_common$name)
one_common2 = one_common[which(one_common$same_sign == 5 & rowSums(one_common[,vert] > 0.1) == 5),]

all_hit = one_common2
rownames(all_hit) = paste0("G", 1:nrow(all_hit))
this_other = -all_hit[,vert]
colnames(this_other) = paste0(colnames(this_other), "_other")
all_hit = cbind(all_hit[,vert], this_other)
all_hit = as.matrix(all_hit)
fc_thresh = 1
all_hit[which(all_hit >  fc_thresh)] =  fc_thresh
all_hit[which(all_hit < -fc_thresh)] = -fc_thresh
gene_annot = data.frame(annot = one_common2$name, row.names = rownames(all_hit))
pheatmap::pheatmap(all_hit, annotation_row = gene_annot, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, filename = paste0("~/scratch/st/results/loose_cons_genes_hit.pdf"))

hit_genes = one_common2$hgnc
hit_combos = paste0(one_common2$hgnc, "_", one_common2$name)
for (v1 in vert) {
  this_df = one_common[which(one_common$hgnc %in% hit_genes), c("hgnc", "name", v1)]
  this_mat = reshape2::acast(this_df, hgnc ~ name, value.var = v1)
  this_mat2 = this_mat[hit_genes,]
  rownames(this_mat2) = 1:nrow(this_mat2)
  this_df2 = reshape2::melt(this_mat2)
  this_df2$Var1 = factor(this_df2$Var1, levels = 1:nrow(this_mat2))
  print(ggplot(this_df2, aes(y = Var1, x = Var2, fill = value)) + geom_raster() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-2, 2), oob = scales::squish ) + scale_y_discrete(labels = hit_genes, name = "", expand = c(0,0)) + scale_x_discrete(name = "", expand = c(0,0)) + theme_classic() + theme(axis.line=element_blank(), axis.text.y = element_text(face = ifelse(one_common2$hgnc %in% common_celltype_genes4$gene, "bold", "plain"))))
  ggsave(paste0("~/scratch/st/results/loose_cons_genes_hit_", v1, ".pdf"), width = 5, height = 12)
  
  # pheatmap::pheatmap(this_mat2, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, filename = paste0("~/scratch/st/results/loose_cons_genes_hit_", v1, ".pdf"))
  # this_df$comobo = paste0(this_df$hgnc, "_", this_df$name)
  # this_df$hgnc_id = 
  # this_df$hgnc = factor(this_df)
  # print(ggplot(this_df, aes_string(x = "name", y = "hgnc", fill = v1)) + geom_tile() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-1, 1), oob = scales::squish ))
  # ggsave(paste0("~/scratch/st/results/loose_cons_genes_hit_", v1, ".pdf"), width = 5, height = 6)
}
