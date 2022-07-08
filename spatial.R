#*******************************************************************************
# Load Libraries ===============================================================
#*******************************************************************************
wdstr = substr(getwd(), 1, 12)
switch(wdstr,
       "C:/Users/mil" = { main_path = "C:/Users/miles/Downloads/";        },
       "/home/george" = { main_path = "~/research/"                       },
       "/storage/scr" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/hom" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/cod" = { main_path = "/storage/scratch1/6/ggruenhagen3/" })
brain_dir = paste0(main_path, "brain/")
data_dir  = paste0(main_path, "st/data/")
out_dir   = paste0(main_path, "st/results/")
if (main_path == "/storage/scratch1/6/ggruenhagen3/") { data_dir = "/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/data/st/data/" }
source(paste0(brain_dir, "/brain_scripts/all_f.R"))
setwd(out_dir)

#*******************************************************************************
# Load Objects =================================================================
#*******************************************************************************
all_merge = qs::qread(paste0(data_dir, "st_070822.qs"))
spo = qs::qread(paste0(data_dir, "st_obj_list_070822.qs"))

#*******************************************************************************
# Integration with BB ==========================================================
#*******************************************************************************
bb = readRDS(paste0(brain_dir, "data/bb_sct_070522.rds"))
bb_convert15 = data.frame(old = 0:14, new = c("8_Glut", "9_Glut", "4_GABA", "15_GABA/Glut", "1_RGC/MG", "10_Glut", "5_GABA", "11_Glut", "6_GABA", "2_OPC/Oligo", "12_Glut", "13_Glut", "14_Glut", "3_Peri", "7_GABA"))
bb_convert53 = data.frame(old = 0:52, new = c("4.1_GABA", "10.1_Glut", "15.1_GABA/Glut", "9.1_Glut", "8.1_Glut", "1.1_RGC", "6_GABA", "5.1_GABA", "9.2_Glut", "8.2_Glut", "15.2_GABA", "11.1_Glut", "8.3_Glut", "8.4_Glut", "9.3_Glut", "4.2_GABA", "8.5_Glut", "5.2_GABA", "8.6_Glut", "8.7_Glut", "1.2_RGC", "4.3_GABA", "4.4_GABA", "9.4_Glut", "9.5_Glut", "8.8_Glut", "9.6_Glut", "4.5_GABA", "12_Glut", "8.9_Glut", "10.2_Glut", "2.1_OPC", "15.3_GABA", "11.2_Glut", "15.4_GABA", "4.6_GABA", "9.7_Glut", "13_Glut", "14_Glut", "4.7_GABA", "11.3_Glut", "9.8_Glut", "8-9_Glut", "15.5_GABA/Glut", "4.8_GABA", "1.3_MG", "2.2_Oligo", "15.6_Glut", "8.10_Glut", "8.11_Glut", "3_Peri", "15.7_Glut", "7_GABA"))
bb$names15 = bb_convert15$new[match(bb$seuratclusters15, bb_convert15$old)]
bb$names53 = bb_convert53$new[match(bb$seuratclusters53, bb_convert53$old)]

# ================ #
# Primary Clusters #
# ================ #
anchors = FindTransferAnchors(reference = bb, query = all_merge, normalization.method = "SCT", npcs = 50)
predictions15 = TransferData(anchorset = anchors, refdata = bb$names15, prediction.assay = TRUE, weight.reduction = all_merge[["pca"]], dims = 1:50)
all_merge[["predictions15"]] = predictions15

DefaultAssay(all_merge) = "predictions15"
SpatialFeaturePlot(all_merge, features = c("8-Glut"), pt.size.factor = 3) + plot_layout(ncol = 4)
all_merge$bb15 = GetTransferPredictions(all_merge, assay = "predictions15")
all_merge$bb15 = str_replace(all_merge$bb15, "-", "_")
all_merge$bb15[which(all_merge$bb15 == "Unassigned")] = NA
all_merge$bb15 = factor(all_merge$bb15, levels = c("15_GABA/Glut", "14_Glut", "13_Glut", "12_Glut", "11_Glut", "10_Glut", "9_Glut", "8_Glut", "7_GABA", "6_GABA", "5_GABA", "4_GABA", "3_Peri", "2_OPC/Oligo", "1_RGC/MG"))
# all_merge$bb15name = factor(bb_convert15$new[match(all_merge$bb15, bb_convert15$old)], levels = c("15_GABA/Glut", "14_Glut", "13_Glut", "12_Glut", "11_Glut", "10_Glut", "9_Glut", "8_Glut", "7_GABA", "6_GABA", "5_GABA", "4_GABA", "3_Peri", "2_OPC/Oligo", "1_RGC/MG"))
Idents(all_merge) = "bb15"

# Paint the "Scores" for each bb cluster on the spatial plots for all samples
for(bb_clust in unique(GetTransferPredictions(all_merge))) {
  if (bb_clust != "Unassigned") {
    Cairo::Cairo(file = paste0(out_dir, "all_bb15_integration_", bb_clust, ".png"), width = 1800, height = 1800, res = 150)
    print(SpatialFeaturePlot(all_merge, features = bb_clust, pt.size.factor = 3) + plot_layout(ncol = 4))
    dev.off()
  }
}

# DimPlot of spots painted by their predicted bb cluster
Cairo::Cairo(file = paste0(out_dir, "all_bb15_integration.png"), width = 2000, height = 1800, res = 150)
print(DimPlot(all_merge, label = T, pt.size = 2, label.size = 6, label.box = F) + theme_void())
dev.off()

# ================== #
# Secondary Clusters #
# ================== #
anchors = FindTransferAnchors(reference = bb, query = all_merge, normalization.method = "SCT", npcs = 50)
predictions53 = TransferData(anchorset = anchors, refdata = bb$names53, prediction.assay = TRUE, weight.reduction = all_merge[["pca"]], dims = 1:50)
all_merge[["predictions53"]] = predictions53

DefaultAssay(all_merge) = "predictions53"
all_merge$bb53 = GetTransferPredictions(all_merge, assay = "predictions53")
all_merge$bb53 = str_replace(all_merge$bb53, "-", "_")
all_merge$bb53[which(all_merge$bb53 == "Unassigned")] = NA
all_merge$bb53 = factor(all_merge$bb53, levels = c("15.7_Glut", "15.6_Glut", "15.5_GABA/Glut", "15.4_GABA", "15.3_GABA", "15.2_GABA", "15.1_GABA/Glut", "14_Glut", "13_Glut", "12_Glut", "11.3_Glut", "11.2_Glut", "11.1_Glut", "10.2_Glut", "10.1_Glut", "9.8_Glut", "9.7_Glut", "9.6_Glut", "9.5_Glut", "9.4_Glut", "9.3_Glut", "9.2_Glut", "9.1_Glut", "8-9_Glut", "8.11_Glut", "8.10_Glut", "8.9_Glut", "8.8_Glut", "8.7_Glut", "8.6_Glut", "8.5_Glut", "8.4_Glut", "8.3_Glut", "8.2_Glut", "8.1_Glut", "7_GABA", "6_GABA", "5.2_GABA", "5.1_GABA", "4.8_GABA", "4.7_GABA", "4.6_GABA", "4.5_GABA", "4.4_GABA", "4.3_GABA", "4.2_GABA", "4.1_GABA", "3_Peri", "2.2_Oligo", "2.1_OPC", "1.3_MG", "1.2_RGC", "1.1_RGC"))
# all_merge$bb53name = factor(bb_convert15$new[match(all_merge$bb53, bb_convert15$old)], levels = c("15_GABA/Glut", "14_Glut", "13_Glut", "12_Glut", "11_Glut", "10_Glut", "9_Glut", "8_Glut", "7_GABA", "6_GABA", "5_GABA", "4_GABA", "3_Peri", "2_OPC/Oligo", "1_RGC/MG"))
Idents(all_merge) = "bb53"

for(bb_clust in unique(GetTransferPredictions(all_merge))) {
  if (bb_clust != "Unassigned") {
    Cairo::Cairo(file = paste0(out_dir, "all_bb53_integration_", bb_clust, ".png"), width = 1800, height = 1800, res = 150)
    print(SpatialFeaturePlot(all_merge, features = bb_clust, pt.size.factor = 3) + plot_layout(ncol = 4))
    dev.off()
  }
}

# DimPlot of spots painted by their predicted bb cluster
Cairo::Cairo(file = paste0(out_dir, "all_bb53_integration.png"), width = 2000, height = 1800, res = 150)
print(DimPlot(all_merge, label = T, pt.size = 2, label.size = 6, label.box = F) + theme_void())
dev.off()

# ================== #
# Overlap of Markers #
# ================== #
st.umap.2 = read.csv("~/Downloads/all_merge_umap2_loose_degs_raw.csv")
st.umap.2 = st.umap.2[which(st.umap.2$p_val_adj < 0.05 & abs(st.umap.2$avg_log2FC) > 0.1 & st.umap.2$pct.1 > 0.05),]
bb15.deg = read.csv("~/research/brain/data/bb_all_cluster_15_degs.csv")
bb15.deg = bb15.deg[which(bb15.deg$p_val_adj < 0.05 & abs(bb15.deg$avg_logFC) > 0.1 & bb15.deg$pct.1 > 0.05),]
bb15.deg$old = bb15.deg$cluster
bb15.deg$new = factor(bb_convert15$new[match(bb15.deg$cluster, bb_convert15$old)], levels = convert15$new.full)
bb15.deg$cluster = bb15.deg$new
bb53.deg = read.csv("~/research/brain/data/bb_all_cluster_53_degs.csv")
bb53.deg = bb53.deg[which(bb53.deg$p_val_adj < 0.05 & abs(bb53.deg$avg_logFC) > 0.1 & bb53.deg$pct.1 > 0.05),]
bb53.deg$old = bb53.deg$cluster
bb53.deg$new = factor(bb_convert53$new[match(bb53.deg$cluster, bb_convert53$old)], levels = rev(convert53$new))
bb53.deg$cluster = bb53.deg$new

# Number of Markers per cluster
ggplot(as.data.frame(table(st.umap.2$cluster)), aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cluster") + ylab("Number of Markers") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + NoLegend() + ggtitle("ST - UMAP2") 
ggplot(as.data.frame(table(bb15.deg$cluster)), aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cluster") + ylab("Number of Markers") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + NoLegend() + ggtitle("BB - Primary")
ggplot(as.data.frame(table(bb53.deg$cluster)), aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cluster") + ylab("Number of Markers") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + NoLegend() + ggtitle("BB - Secondary")

st.umap.2$hgnc = st.umap.2$gene
bb15.deg$hgnc = bb15.deg$gene
bb15.deg$cluster = as.vector(bb15.deg$cluster)
bb53.deg$hgnc = bb53.deg$gene
bb53.deg$cluster = as.vector(bb53.deg$cluster)
my.list = list(st.umap.2, bb15.deg)
names(my.list) = c("ST", "BB15")
bh_out15 = bigHeatmap(my.list, pdf.name = "st_w_bb15.pdf", single.org.1 = "ST", single.org.2 = "BB15", cluster.includes.org = T)
bh_out15_1 = bh_out15[[1]]
bh_out15_1$org.cluster.2 = factor(bh_out15_1$org.cluster.2, levels = rev(convert15$new.full))
Cairo::Cairo(file = paste0(out_dir, "all_bb15_marker_overlap_heatmap_num.png"), width = 2000, height = 1000, res = 200)
ggplot(bh_out15_1, aes(org.cluster.2, org.cluster.1, fill=num)) + geom_tile() + scale_fill_viridis_c(name = "") + guides(color = 'none') + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ggtitle("Number of Overlapping Markers")
dev.off()
Cairo::Cairo(file = paste0(out_dir, "all_bb15_marker_overlap_heatmap_pct.png"), width = 2000, height = 1000, res = 200)
ggplot(bh_out15_1, aes(org.cluster.2, org.cluster.1, fill=pct.both)) + geom_tile() + scale_fill_viridis_c(name = "") + guides(color = 'none') + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ggtitle("Percent of Overlapping Markers")
dev.off()

my.list = list(st.umap.2, bb53.deg)
names(my.list) = c("ST", "BB53")
bh_out53 = bigHeatmap(my.list, pdf.name = "st_w_bb53.pdf", single.org.1 = "ST", single.org.2 = "BB53", cluster.includes.org = T)
bh_out53_1 = bh_out53[[1]]
bh_out53_1$org.cluster.2 = factor(bh_out53_1$org.cluster.2, levels = convert53$new)

Cairo::Cairo(file = paste0(out_dir, "all_bb53_marker_overlap_heatmap_num.png"), width = 2000, height = 1000, res = 200)
ggplot(bh_out53_1, aes(org.cluster.2, org.cluster.1, fill=num)) + geom_tile() + scale_fill_viridis_c(name = "") + guides(color = 'none') + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ggtitle("Number of Overlapping Markers")
dev.off()
Cairo::Cairo(file = paste0(out_dir, "all_bb53_marker_overlap_heatmap_pct.png"), width = 2000, height = 1000, res = 200)
ggplot(bh_out53_1, aes(org.cluster.2, org.cluster.1, fill=pct.both)) + geom_tile() + scale_fill_viridis_c(name = "") + guides(color = 'none') + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ggtitle("Percent of Overlapping Markers")
dev.off()

#*******************************************************************************
# BHVE vs CTRL DEGs ============================================================
#*******************************************************************************

#*******************************************************************************
# Initial clustering ===========================================================
#*******************************************************************************

# Load Data
dir_of_sr_dirs = "~/Downloads/sp_data/" # Folder where all the individual samples are kept
spo = list()
spo[["c2a"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_295_A1_fr/outs/"))
spo[["c2b"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_295_B1_fr/outs/"))
spo[["c2c"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_295_C1_fr/outs/"))
spo[["c2d"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_295_D1_fr/outs/"))
spo[["b2a"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_296_A1_fr/outs/"))
spo[["b2b"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_296_B1_fr/outs/"))
spo[["b2c"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_296_C1_fr/outs/"))
spo[["b2d"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_296_D1_fr/outs/"))
spo[["c1a"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_297_A1_fr/outs/"))
spo[["c1b"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_297_B1_fr/outs/"))
spo[["c1c"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_297_C1_fr/outs/"))
spo[["c1d"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_297_D1_fr/outs/"))
spo[["b1c"]] = Load10X_Spatial(paste0(dir_of_sr_dirs, "/JTS12_293_C1_fr/outs/"))

# Add MetaData
spo[["c2a"]]$fish = "c2"; spo[["c2a"]]$fnum = 295; spo[["c2a"]]$cond = "CTRL"; spo[["c2a"]]$area = "a"; spo[["c2a"]]$sample = "c2a"; 
spo[["c2b"]]$fish = "c2"; spo[["c2b"]]$fnum = 295; spo[["c2b"]]$cond = "CTRL"; spo[["c2b"]]$area = "b"; spo[["c2b"]]$sample = "c2b";
spo[["c2c"]]$fish = "c2"; spo[["c2c"]]$fnum = 295; spo[["c2c"]]$cond = "CTRL"; spo[["c2c"]]$area = "c"; spo[["c2c"]]$sample = "c2c";
spo[["c2d"]]$fish = "c2"; spo[["c2d"]]$fnum = 295; spo[["c2d"]]$cond = "CTRL"; spo[["c2d"]]$area = "d"; spo[["c2d"]]$sample = "c2d";
spo[["b2a"]]$fish = "b2"; spo[["b2a"]]$fnum = 296; spo[["b2a"]]$cond = "BHVE"; spo[["b2a"]]$area = "a"; spo[["b2a"]]$sample = "b2a";
spo[["b2b"]]$fish = "b2"; spo[["b2b"]]$fnum = 296; spo[["b2b"]]$cond = "BHVE"; spo[["b2b"]]$area = "b"; spo[["b2b"]]$sample = "b2b";
spo[["b2c"]]$fish = "b2"; spo[["b2c"]]$fnum = 296; spo[["b2c"]]$cond = "BHVE"; spo[["b2c"]]$area = "c"; spo[["b2c"]]$sample = "b2c";
spo[["b2d"]]$fish = "b2"; spo[["b2d"]]$fnum = 296; spo[["b2d"]]$cond = "BHVE"; spo[["b2d"]]$area = "d"; spo[["b2d"]]$sample = "b2d";
spo[["c1a"]]$fish = "c1"; spo[["c1a"]]$fnum = 297; spo[["c1a"]]$cond = "CTRL"; spo[["c1a"]]$area = "a"; spo[["c1a"]]$sample = "c1a";
spo[["c1b"]]$fish = "c1"; spo[["c1b"]]$fnum = 297; spo[["c1b"]]$cond = "CTRL"; spo[["c1b"]]$area = "b"; spo[["c1b"]]$sample = "c1b";
spo[["c1c"]]$fish = "c1"; spo[["c1c"]]$fnum = 297; spo[["c1c"]]$cond = "CTRL"; spo[["c1c"]]$area = "c"; spo[["c1c"]]$sample = "c1c";
spo[["c1d"]]$fish = "c1"; spo[["c1d"]]$fnum = 297; spo[["c1d"]]$cond = "CTRL"; spo[["c1d"]]$area = "d"; spo[["c1d"]]$sample = "c1d";
spo[["b1c"]]$fish = "b1"; spo[["b1c"]]$fnum = 293; spo[["b1c"]]$cond = "BHVE"; spo[["b1c"]]$area = "c"; spo[["b1c"]]$sample = "b1c"; 

# Rename Spots: Add the sample name to the spot names
for (s in names(spo)) {
  print(s)
  spo[[s]] = RenameCells(spo[[s]], paste0(s))
  spo[[s]]@images$slice1@key = s
}

# Extraneous cells
remove_cells = data.frame(cell = c("c2a_ACCTCCGTTATTCACC-1", "c2a_GAAATATGCTTGAATG-1", "c2a_ACGCTTAGTGTCTCTC-1", "c2a_AATGACTGTCAGCCGG-1", "c2b_GCTTAGGGAAGCGGTA-1", "c2b_GCTTAGGGAAGCGGTA-1", "c2b_AGGCAATACGGAGGAC-1", "c2b_TGCATGAGTAGATTCG-1", "c2b_TTAGGTGTGACTGGTC-1", "c2c_AGGACTTATAGGAGAA-1", "c2c_GACGCATACCCGTCGG-1", "c2c_CCCTGACTAACAAATT-1", "c2d_GTGGCCTAATATCATT-1", "c2d_TTAAACAGAGTCCCGC-1", "c2d_AACAGCTGTGTGGCAA-1", "c2d_CACTCCTATGTAAGAT-1", "c2d_ACGTTAGATTTGCCCG-1", "c2d_TGCTTCCCAAGCAGTA-1", "b2a_CTGTTGGCTCTTCTGA-1", "b2a_AATTACGAGACCCATC-1", "b2a_ACCGAAGAGTCTGGTT-1", "b2a_CCGTTCCGAATCTCGG-1", "b2a_CGAACGGCCGGACAAC-1", "b2a_GGTCGTAAGCTCGCAC-1", "b2a_TCCATCAATACTAATC-1", "b2a_TCGCGCGTTTACATGA-1", "b2a_TGAGGCATGTACTGTG-1", "b2a_TGATTCAGGTCCCGCG-1", "b2a_CACCTAATCAGTTTAC-1", "b2a_ATTCTTCGTACTTATG-1", "b2b_GGACTCGTGAGTGGTC-1", "b2b_TCGTAAGCTCCGAGGA-1", "b2d_GACAGCCAGACCTGAC-1", "b2d_GATAACTCGCACTGTG-1", "c1a_TCTCCACAAGTTGAAT-1", "c1a_TGGTATCGCATCCCAA-1", "c1a_CAGTCTGTATACTGGG-1", "c1c_AGTGGCGTCTGAAGGT-1", "c1d_ACAATCCATTTAAACC-1", "c1d_TACGAGAACTTCACGT-1", "c1d_TCTTACGGCATCCGAC-1", "c1d_TTACCATTGATTACCC-1", "c1d_GTACTCCTGGGTATGC-1", "c1d_CGGTTATCCAACAGTG-1", "c1d_GGGTCACCGTGACGGT-1", "c1d_TGTTGTCAAGAAGTCT-1", "b1c_ACAGGTGTGTTGTTGC-1", "b1c_AAACAATCTACTAGCA-1"))
remove_cells$sample = reshape2::colsplit(remove_cells$cell, "_", c('1', '2'))[,1]
for (s in unique(remove_cells$sample)) {
  spo[[s]] = spo[[s]][, which(! colnames(spo[[s]]) %in% remove_cells$cell )]
}

# Cluster each sample separate. Try multiple clustering methods.
for (s in names(spo)) {
  print(s)
  
  # Default Method of Clustering
  print(ncol(spo[[s]]))
  spo[[s]] = subset(spo[[s]], subset = nCount_Spatial > 0)
  print(ncol(spo[[s]]))
  spo[[s]] = SCTransform(spo[[s]], assay = "Spatial", verbose = FALSE)
  spo[[s]] = RunPCA(spo[[s]], assay = "SCT", verbose = FALSE)
  spo[[s]] = FindNeighbors(spo[[s]], reduction = "pca", dims = 1:30)
  spo[[s]] = FindClusters(spo[[s]], verbose = FALSE)
  spo[[s]] = RunUMAP(spo[[s]], reduction = "pca", dims = 1:30)
  spo[[s]]$default_cluster = spo[[s]]$seurat_clusters
  p1 = DimPlot(spo[[s]], reduction = "umap", label = TRUE) + ggtitle("PCA Clustering - Default Res")
  p2 = SpatialDimPlot(spo[[s]], label = TRUE, label.size = 3) + NoLegend()
  
  # Clustering on UMAP Dimensions w/ default resolution (0.8)
  spo[[s]] = FindNeighbors(spo[[s]], reduction = "umap", dims = 1:2)
  spo[[s]] = FindClusters(spo[[s]], verbose = FALSE)
  spo[[s]]$umap_cluster = spo[[s]]$seurat_clusters
  p3 = DimPlot(spo[[s]], reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering - Default Res")
  p4 = SpatialDimPlot(spo[[s]], label = TRUE, label.size = 3) + NoLegend()
  
  # Clustering on UMAP Dimensions w/ low resolution (0.3)
  spo[[s]] = FindNeighbors(spo[[s]], reduction = "umap", dims = 1:2)
  spo[[s]] = FindClusters(spo[[s]], verbose = FALSE, resolution = 0.3)
  spo[[s]]$umap_cluster2 = spo[[s]]$seurat_clusters
  p5 = DimPlot(spo[[s]], reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering - Lower Res")
  p6 = SpatialDimPlot(spo[[s]], label = TRUE, label.size = 3) + NoLegend()
  
  # Save the plots
  # Cairo::Cairo(file = paste0(out_dir, s, "_clustering.png"), width = 1800, height = 1200, res = 150)
  # print(wrap_plots(p1, p3, p5, p2, p4, p6, ncol = 3))
  # dev.off()
}

# Plot Quality Metrics for Each Spot: nCount (Number of UMIs) and nFeature (Number of Genes)
out_dir = "~/research/brain/sp_results/"
sp.fp.list = list()
sp.fp.nFeature.list = list()
for (s in names(spo)) {
  Idents(spo[[s]]) = 1
  # nCount Plots
  plot1 = VlnPlot(spo[[s]], features = "nCount_Spatial", pt.size = 0.1) + NoLegend()  + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  plot2 = SpatialFeaturePlot(spo[[s]], features = "nCount_Spatial", stroke = 0, pt.size.factor = 3) + theme(legend.position = "right") + ggtitle(paste0(s, " nCount"))
  plot3 = ggplot(data.frame(nCount = spo[[s]]$nCount_Spatial), aes(x = nCount)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nCount")
  sp.fp.list[[s]] = SpatialFeaturePlot(spo[[s]], features = "nCount_Spatial", stroke = 0) + NoLegend()
  
  # nFeature Plots
  plot4 = VlnPlot(spo[[s]], features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()  + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  plot5 = SpatialFeaturePlot(spo[[s]], features = "nFeature_Spatial", stroke = 0, pt.size.factor = 3) + theme(legend.position = "right") + ggtitle(paste0(s, " nFeature"))
  plot6 = ggplot(data.frame(nFeature = spo[[s]]$nFeature_Spatial), aes(x = nFeature)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nFeature")
  sp.fp.nFeature.list[[s]] = SpatialFeaturePlot(spo[[s]], features = "nFeature_Spatial", stroke = 0) + NoLegend()
  
  # Save the plots
  Cairo::Cairo(file = paste0(out_dir, s, "_nCount_nFeature.png"), width = 1800, height = 1200, res = 150)
  print(wrap_plots(plot1, plot3, plot2, plot4, plot6, plot5, ncol = 3))
  dev.off()
}

# Plot of nCount for all samples
Cairo::Cairo(file = paste0(out_dir, "all_nCount.png"), width = 2400, height = 2400, res = 150)
print(wrap_plots(sp.fp.list, ncol = 4))
dev.off()

# Plot of nFeature for all samples
Cairo::Cairo(file = paste0(out_dir, "all_nFeature.png"), width = 2400, height = 2400, res = 150)
print(wrap_plots(sp.fp.nFeature.list, ncol = 4))
dev.off()

# Merge all the samples
all_merge = merge(spo[[names(spo)[1]]], spo[[names(spo)[2]]])
for (s in names(spo)[3:length(spo)]) {
  print(s)
  all_merge = merge(all_merge, spo[[s]])
}

# Merged object clustering
all_merge = subset(all_merge, subset = nCount_Spatial > 0)
all_merge = SCTransform(all_merge, assay = "Spatial", verbose = FALSE)
all_merge = RunPCA(all_merge, assay = "SCT", verbose = FALSE)
all_merge = RunUMAP(all_merge, reduction = "pca", dims = 1:30)
all_merge = FindNeighbors(all_merge, reduction = "umap", dims = 1:2)
all_merge = FindClusters(all_merge, verbose = FALSE, resolution = 0.55)
all_merge$cluster = all_merge$seurat_clusters
all_merge$all_cluster = all_merge$seurat_clusters
p1 = DimPlot(all_merge, reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering - Brianna Res")

# Default Method of Clustering
# all_merge = FindNeighbors(all_merge, reduction = "pca", dims = 1:30)
# all_merge = FindClusters(all_merge, verbose = FALSE)
# all_merge = RunUMAP(all_merge, reduction = "pca", dims = 1:30)
# all_merge$all_default_cluster = all_merge$seurat_clusters
# p1 = DimPlot(all_merge, reduction = "umap", label = TRUE) + ggtitle("PCA Clustering - Default Res")

# Clustering on UMAP Dimensions w/ default resolution (0.8)
all_merge = FindNeighbors(all_merge, reduction = "umap", dims = 1:2)
all_merge = FindClusters(all_merge, verbose = FALSE)
all_merge$all_cluster_res_80 = all_merge$seurat_clusters
p2 = DimPlot(all_merge, reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering - Default Res")

# Clustering on UMAP Dimensions w/ low resolution (0.3)
all_merge = FindNeighbors(all_merge, reduction = "umap", dims = 1:2)
all_merge = FindClusters(all_merge, verbose = FALSE, resolution = 0.3)
all_merge$all_cluster_res_30 = all_merge$seurat_clusters
p3 = DimPlot(all_merge, reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering - Lower Res")

# Reset to the clustering that Brianna likes best
all_merge$seurat_clusters = all_merge$cluster
Idents(all_merge) = all_merge$cluster

# Save the Clustering Plots
# Cairo::Cairo(file = paste0(out_dir, "all_clustering.png"), width = 2000, height = 600, res = 150)
# print(wrap_plots(p1, p2, p3, ncol = 3))
# dev.off()

# Add the clusters found from the merged object into the separate objects and plot
all_on_separate  = list()
for (s in names(spo)) {
  spo[[s]]$cluster              = all_merge$cluster[match(colnames(spo[[s]]), colnames(all_merge))]
  spo[[s]]$all_cluster_res_80   = all_merge$all_cluster_res_80[match(colnames(spo[[s]]),    colnames(all_merge))]
  spo[[s]]$all_cluster_res_30   = all_merge$all_cluster_res_30[match(colnames(spo[[s]]),   colnames(all_merge))]
  
  Idents(spo[[s]]) = spo[[s]]$cluster
  p1 = SpatialDimPlot(spo[[s]], label = TRUE, label.size = 3) + NoLegend()
  
  all_on_separate[[s]] = p1
}

# Clusters from merged object back on separate objects
Cairo::Cairo(file = paste0(out_dir, "all_cluster_on_separate.png"), width = 2400, height = 2400, res = 150)
print(wrap_plots(all_on_separate_umap2, ncol = 4))
dev.off()

# Show the clusters by sample
all_merge$sample = factor(all_merge$sample, levels = unique(all_merge$sample))
Cairo::Cairo(file = paste0(out_dir, "all_clustering_split_by_sample_color_by_umap2.png"), width = 2400, height = 2400, res = 150)
print(DimPlot(all_merge, reduction = "umap", label = T, split.by = "sample", ncol = 4))
dev.off()

names(all_merge@images) = levels(all_merge$sample)

# Before Saving, Go to the Integration w/ BB section to add that data into the object
Idents(all_merge) = all_merge$cluster
DefaultAssay(all_merge) = "SCT"
all_merge@active.assay = "SCT"
for (s in names(spo)) {
  spo[[s]]$bb15 = all_merge$bb15[match(colnames(spo[[s]]), colnames(all_merge))]
  spo[[s]]$bb53 = all_merge$bb53[match(colnames(spo[[s]]), colnames(all_merge))]
}

# Save the merged object
saveRDS(all_merge,   paste0(out_dir, "st_070822.rds"))
qs::qsave(all_merge, paste0(out_dir, "st_070822.qs"))

# Save the list of samples
saveRDS(spo,   paste0(out_dir, "st_obj_list_070822.rds"))
qs::qsave(spo, paste0(out_dir, "st_obj_list_070822.qs"))

# Save each sample separately
for (s in names(spo)) {
  saveRDS(spo[[s]], paste0(out_dir, "/sample_objs/", s, ".rds"))
}

#*******************************************************************************
# Trash Can ====================================================================
#*******************************************************************************

# # Look at clustering at multiple resolutions
# for (this.res in seq(0.3, 0.8, by = 0.05)) {
#   all_merge = FindNeighbors(all_merge, reduction = "umap", dims = 1:2)
#   all_merge = FindClusters(all_merge, verbose = FALSE, resolution = this.res)
#   all_merge$all_umap_cluster2 = all_merge$seurat_clusters
#   p2 = DimPlot(all_merge, reduction = "umap", label = TRUE) + ggtitle(paste0("UMAP Clustering - ", this.res))
#   
#   Cairo::Cairo(file = paste0(out_dir, "all_cluster_multi_res/all_umap_res_", this.res, ".png"), width = 1200, height = 1200, res = 150)
#   print(p2)
#   dev.off()
# } # resolution = 0.55 looks the best

# # Manual Rotation of Images
# coords2 = coords
# coords2$px = coords2$imagecol
# coords2$py = -coords2$imagerow
# coords3 = coords2[, c("px", "py")]
# alpha = 45 # positive rotates counterclockwise.
# rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
# coords.test = coords2[, c("imagerow", "imagecol")]
# coords.test <- as.data.frame(t(rotm %*% (t(coords.test))))
# coords3 <- as.data.frame(t(rotm %*% (t(coords3))))
# coords2$rot.px = coords3[,1]
# coords2$rot.py = coords3[,2]
# 
# big.img.raster.real = img.grob$raster
# big.img.raster.dummy = big.img.raster.real
# big.img.raster.dummy[which(big.img.raster.dummy != "#000000ff")] = "#000000ff"
# big.img.raster.dummy[my.x.min:my.x.max, my.y.min:my.y.max] = "#ff0000"
# m.img.dummy = magick::image_read(big.img.raster.dummy)
# m.img.dummy = magick::image_rotate(m.img.dummy, degrees = -alpha)
# m.img.dummy.grob = rasterGrob(m.img.dummy)
# m.img.dummy.grob$raster[which(!m.img.dummy.grob$raster %in% bad.dummy.cols)] = "#ff0000"
# 
# good.idx = as.data.frame(which(m.img.dummy.grob$raster == "#ff0000ff", arr.ind = T))
# new.x.min = min(good.idx[,1])
# new.x.max = max(good.idx[,1])
# new.y.min = min(good.idx[,2])
# new.y.max = max(good.idx[,2])
# m.img.real = magick::image_read(big.img.raster.real)
# m.img.real = magick::image_rotate(m.img.real, degrees = -alpha)
# m.img.real.grob = rasterGrob(m.img.real)
# m.img.real.grob.raster = m.img.real.grob$raster
# m.img.real.grob.raster = m.img.real.grob.raster[new.x.min:new.x.max, new.y.min:new.y.max]
# m.img.real.grob$raster = m.img.real.grob.raster
# m.img.real.grob$width  = unit(x = 1, units = "npc")
# m.img.real.grob$height = unit(x = 1, units = "npc")
# 
# big.img.raster.df = melt(as.matrix(img.grob$raster))
# colnames(big.img.raster.df)[1:2] = c("orig.x", "orig.y")
# big.img.raster.df$contains_spot = F
# big.img.raster.df$contains_spot[which(big.img.raster.df$orig.x > my.x.min & big.img.raster.df$orig.x < my.x.max & big.img.raster.df$orig.y > my.y.min & big.img.raster.df$orig.y < my.y.max)] = T
# big.img.raster.df[, c("rot.x", "rot.y")] <- as.data.frame(t(rotm %*% (t( big.img.raster.df[,1:2] ))))
# big.img.raster.df$m.x = 0
# big.img.raster.df$m.y = 0