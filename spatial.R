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
gene_info = read.table(paste0(main_path, "/all_research/gene_info_2.txt"), header = T, stringsAsFactors = F)
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

#------------------------------------#
# Cell2location Secondary - Separate #
#------------------------------------#
all_merge_hi$class2_81_84 = "#FDE72500"
my.thresh = 0.4
this.samples = names(spo)
for (s in this.samples) {
  s_mean = read.csv(paste0(out_dir, "cell2location/bb_secondary/", s, "/means.csv"))
  rownames(s_mean) = s_mean$X
  s_mean = s_mean[,2:ncol(s_mean)]
  colnames(s_mean) = str_replace(colnames(s_mean), "meanscell_abundance_w_sf_", "")
  col_order = convert53$new[which(convert53$new %in% colnames(s_mean))]
  s_mean = s_mean[, c( col_order, paste0("RGC", 0:10) )]
  df84_81_2 = data.frame(class2 = s_mean[, "RGC2"], class84 = s_mean[, "8.4_Glut"], class81 = s_mean[, "8.1_Glut"], row.names = rownames(s_mean))
  
  min_max_value = min(c(max(df84_81_2$class81), max(df84_81_2$class84)))
  
  df84_81_2$class84[which(df84_81_2$class84 > min_max_value)] = min_max_value
  df84_81_2$class81[which(df84_81_2$class81 > min_max_value)] = min_max_value
  
  df84_81_2$dif84_81 = df84_81_2$class84 - df84_81_2$class81
  df84_81_2$dif84_81_col = viridis(100)[cut(df84_81_2$dif84_81,100)] # scaled across all samples
  df84_81_2$dif84_81_col[which(df84_81_2$dif84_81 == 0)] = "#FDE72500"
  df84_81_2$dif84_81_col[which(df84_81_2$class2 >= my.thresh)] = "#FF9E24"
  all_merge_hi$class2_81_84[rownames(df84_81_2)] = df84_81_2$dif84_81_col
}

pdf(paste0(out_dir, paste0("cell2loc_sep_class2_81_84_", my.thresh, "_2.pdf")), width = 12, height = 12, onefile = F)
print(myMultiSFP(all_merge_hi, feature = "class2_81_84", pt.size.multiplier = 1.3, pal = colorRampPalette(viridis(100)), rm.zero = T, col.ident = T, high.res = T ))
dev.off()


#-------------------------#
# Cell2location Secondary #
#-------------------------#
c2l_mean = read.csv("cell2location_spatial_output_means.csv")
rownames(c2l_mean) = c2l_mean$X
colnames(c2l_mean)[2:ncol(c2l_mean)] = as.character(0:52)
c2l_mean$top1 = colnames(c2l_mean[,2:ncol(c2l_mean)])[apply(c2l_mean[,2:ncol(c2l_mean)],1,which.max)]
c2l_mean$top1 = convert53$new[match(c2l_mean$top1, convert53$old)]
c2l_mean$top1 = factor(c2l_mean$top1, levels = convert53$new)
all_merge$bb53_c2l = c2l_mean$top1 
Idents(all_merge) = all_merge$bb53_c2l
DimPlot(all_merge, order = T, pt.size = 0.8, label = T) + coord_fixed()
all_merge$bb53_8_1_Glut = c2l_mean[, "4"]
all_merge$bb53_8_4_Glut = c2l_mean[, "13"]

# Paint 8.1 and 8.4
Cairo::Cairo(file = paste0(out_dir, "bb53_c2l_8_1_glut.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "bb53_8_1_Glut", pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5)
dev.off()
Cairo::Cairo(file = paste0(out_dir, "bb53_c2l_8_4_glut.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "bb53_8_4_Glut", pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5)
dev.off()

# Visualize 8.1 vs 8.4 differences
all_merge$diff_84_81 = c2l_mean[, "13"] - c2l_mean[, "4"]
Cairo::Cairo(file = paste0(out_dir, "bb53_c2l_8_4_glut_diff.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "diff_84_81", pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5)
dev.off()
Cairo::Cairo(file = paste0(out_dir, "bb53_c2l_8_4_glut_diff_diffscale.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "diff_84_81", pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5, same.col.scale = F)
dev.off()
scaled_diff = unlist(lapply(levels(all_merge$sample), function(x) scale(all_merge$diff_84_81[which(all_merge$sample == x)]) ))
names(scaled_diff) = unlist(lapply(levels(all_merge$sample), function(x) colnames(all_merge)[which(all_merge$sample == x)] ))
all_merge$scaled_diff = scaled_diff[colnames(all_merge)]
Cairo::Cairo(file = paste0(out_dir, "bb53_c2l_8_4_glut_diff_scaled.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "scaled_diff", pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5)
dev.off()
Cairo::Cairo(file = paste0(out_dir, "bb53_c2l_8_4_glut_diff_scaled_diffscale.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "scaled_diff", pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5, same.col.scale = F)
dev.off()

# Paint every cluster
for (i in 0:52) {
  good_name = convert53$new[match(i, convert53$old)]
  good_name = str_replace(good_name, "/", "_")
  all_merge@meta.data[,paste0("bb53_", i)] = c2l_mean[, as.character(i)]
  Cairo::Cairo(file = paste0(out_dir, "/cell2location/", good_name, ".png"), width = 1800, height = 1800, res = 150)
  print(allSamplesSFP(all_merge, paste0("bb53_", i), pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5))
  dev.off()
}

Cairo::Cairo(file = paste0(out_dir, "number_of_cells.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "num_cell", pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5)
dev.off()
Cairo::Cairo(file = paste0(out_dir, "number_of_cells_diffscale.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "num_cell", pal = colorRampPalette(viridis(100)), pt.size.multiplier = 1.5, same.col.scale = F)
dev.off()

col_df = data.frame(cluster = sort(unique(all_merge$bb53_c2l)), cols = hue_pal()(length(unique(all_merge$bb53_c2l))))
all_merge$bb53_c2l_col = col_df$cols[match(all_merge$bb53_c2l, col_df$cluster)]
Cairo::Cairo(file = paste0(out_dir, "top_cluster.png"), width = 1800, height = 1800, res = 150)
allSamplesSFP(all_merge, "bb53_c2l_col", pal = colorRampPalette(viridis(100)), col.ident = T, pt.size.multiplier = 1.5)
dev.off()

c2l_mean_raw = read.csv("cell2location_spatial_output_means.csv")
rownames(c2l_mean_raw) = c2l_mean_raw$X
c2l_mean_raw$X = NULL
colnames(c2l_mean_raw) = as.character(0:52)

c2l_q05 = read.csv("cell2location_spatial_output_q05.csv")
rownames(c2l_q05) = c2l_q05$X
c2l_q05 = c2l_q05[,2:ncol(c2l_q05)]
colnames(c2l_q05) = as.character(0:52)

c2l_q95 = read.csv("cell2location_spatial_output_q95.csv")
rownames(c2l_q95) = c2l_q95$X
c2l_q95 = c2l_q95[,2:ncol(c2l_q95)]
colnames(c2l_q95) = as.character(0:52)

c2l_std = read.csv("cell2location_spatial_output_stds.csv")
rownames(c2l_std) = c2l_std$X
c2l_std = c2l_std[,2:ncol(c2l_std)]
colnames(c2l_std) = as.character(0:52)
c2l_spot_which_max_num = apply(c2l_mean_raw,1,which.max)
std.df = data.frame(max_num = apply(c2l_mean_raw,1,max), max_num_idx = c2l_spot_which_max_num, std = sapply(1:nrow(c2l_std), function(x) c2l_std[x, c2l_spot_which_max_num[x]]))
ggplot(std.df, aes(x = max_num, y = std)) + geom_point(alpha = 0.25) + theme_bw() + scale_x_continuous(expand = c(0,0), "Maximum # of Cells") + scale_y_continuous(expand = c(0,0), name = "STD") + ggtitle("STD for Cell Type with Max # of Cells in a Spot")

# Percentage
c2l_mean_pct = read.csv("cell2location_spatial_output_means.csv")
rownames(c2l_mean_pct) = c2l_mean_pct$X
c2l_mean_pct = c2l_mean_pct[,2:ncol(c2l_mean_pct)]
c2l_spot_sum = rowSums(c2l_mean_pct)
c2l_spot_max_num = apply(c2l_mean_pct,1,max)
c2l_mean_pct = c2l_mean_pct / c2l_spot_sum
colnames(c2l_mean_pct) = as.character(0:52)

all_merge$num_cell = c2l_spot_sum
all_merge$max_num = c2l_spot_max_num
all_merge$max_pct = apply(c2l_mean_pct,1,max)
# all_merge$num_cell[which(all_merge$num_cell >= 40)] = 40
all_merge$max_num[which(all_merge$max_num >= 25)] = 25
Idents(all_merge) = all_merge$cluster
FeaturePlot(all_merge, "num_cell", order = T, pt.size = 0.8, label = F) + theme_void() + coord_fixed() + scale_color_viridis()
FeaturePlot(all_merge, "max_num", order = T, pt.size = 0.8, label = F) + theme_void() + coord_fixed() + scale_color_viridis()
FeaturePlot(all_merge, "max_pct", order = T, pt.size = 0.8, label = F) + theme_void() + coord_fixed() + scale_color_viridis()

stat.df = data.frame(max_num = all_merge$max_num, max_pct = all_merge$max_pct, combo = all_merge$max_num * all_merge$max_pct)
ggplot(stat.df, aes(x = max_num, y = max_pct)) + geom_point(alpha = 0.25) + theme_bw() + xlab("Max # Cells from a Cell Type") + ylab("Max % from a Cell Type")

# Number of Cell Types with >10% in a spot
my.pct = 0.1
count.my.pct = sapply(1:nrow(c2l_mean_pct), function(x) length(which(c2l_mean_pct[x,] > my.pct)))
count.my.pct.df = as.data.frame(table(count.my.pct))
ggplot(count.my.pct.df, aes(x = count.my.pct, y = Freq)) + geom_bar(stat = "identity") + theme_classic() + xlab("# of Cell Types > 10% in a Spot") + scale_y_continuous(expand = c(0,0), name = "# of Spots")

# c2l_mean_pct2 = c2l_mean_pct
# colnames(c2l_mean_pct2) = convert53$new[match(colnames(c2l_mean_pct2), convert53$old)]
# c2l_mean_pct2 = c2l_mean_pct2[, convert53$new]
# c2l_mean_pct2$other = sapply(1:nrow(c2l_mean_pct), function(x) sum(c2l_mean_pct))

c2l_mean_num2 = c2l_mean[, as.character(0:52)]
colnames(c2l_mean_num2) = convert53$new[match(colnames(c2l_mean_num2), convert53$old)]
c2l_mean_num2 = c2l_mean_num2[, convert53$new]
# c2l_mean_num2_neglible = sapply(1:nrow(c2l_mean_num2), function(x) sum(c2l_mean_num2[x,which(c2l_mean_num2[x,] < 0.5)]))
# c2l_mean_num2 = as.matrix(c2l_mean_num2)
# c2l_mean_num2[which(c2l_mean_num2 < 0.5)] = 0
# c2l_mean_num2 = as.data.frame(c2l_mean_num2)
# c2l_mean_num2$neglible = c2l_mean_num2_neglible
# c2l_mean_num2[1:5, 50:54]

# Plot percentages as pi chart
library("scatterpie")
coords = GetTissueCoordinates(object = all_merge_hi, image = "c2c")
# coords = GetTissueCoordinates(object = spo[["c2c"]])
coords$spot = 1:nrow(coords)
coords[, colnames(c2l_mean_num2)] = c2l_mean_num2[rownames(coords),]
rownames(coords) = NULL
pdf("c2c_pi_chart.pdf", width = 12, height = 12, onefile = F)
# print(ggplot() + geom_scatterpie(data = coords, aes(x=imagecol, y=-imagerow, group = spot), cols = colnames(c2l_mean_num2), color=NA, alpha = 0.8) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void()) + coord_fixed() + scale_fill_manual(values = c(hue_pal()(48), 'gray60'), name = "cell type")
print(ggplot() + geom_scatterpie(data = coords, aes(x=imagecol, y=-imagerow, group = spot), cols = colnames(c2l_mean_num2), color=NA, alpha = 0.8) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed())
dev.off()

c2l_mean_num2_round = round(c2l_mean_num2)
coords_multi = data.frame()
for (i in 1:nrow(c2l_mean_num2_round)) {
  this_non_zero_cell_type = colnames(c2l_mean_num2_round)[which(c2l_mean_num2_round[i,] != 0)]
  for (this_ct in this_non_zero_cell_type) {
    coords_multi = rbind(coords_multi, 
                         data.frame(spot = rep(rownames(c2l_mean_num2_round)[i], c2l_mean_num2_round[i, this_ct]), 
                                    imagerow = rep(coords$imagerow[i], c2l_mean_num2_round[i, this_ct]), 
                                    imagecol = rep(coords$imagecol[i], c2l_mean_num2_round[i, this_ct]), 
                                    ct = rep(this_ct, c2l_mean_num2_round[i, this_ct])) ) 
  }
  
}
coords_multi$value = all_merge$cluster[coords_multi$spot]
coords_multi$ct = factor(coords_multi$ct, levels = convert53$new)

img.grob = GetImage(all_merge_hi, image = "c2c")
img.grob.test = img.grob$raster[min(coords$imagerow):max(coords$imagerow), min(coords$imagecol):max(coords$imagecol)]
img.grob.test.grob = img.grob
img.grob.test.grob$raster = img.grob.test

pdf("c2c_with_mini_cells_no_circle.pdf", width = 12, height = 12, onefile = F)
# print(ggplot(coords_multi, aes(x=imagecol, y=-imagerow)) + annotation_custom(img.grob.test.grob) + ggforce::geom_mark_hull(aes(fill = value, group = spot), expand = unit(4, "mm"), radius = unit(4, "mm"), size = 0, color = NA) + geom_point(size = 0.75, position = position_jitter(width = 2.5, height = 2.5), aes(color = ct)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed() + guides(fill="none"))
# print(ggplot(coords_multi, aes(x=imagecol, y=-imagerow)) + annotation_custom(img.grob.test.grob) + ggforce::geom_mark_hull(aes(fill = value, group = spot), expand = unit(4, "mm"), radius = unit(4, "mm"), size = 0, color = NA, alpha = 0.5) + geom_point(size = 0.75, position = position_jitter(width = 8, height = 8), aes(color = ct)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed() + guides(fill="none"))
print(ggplot(coords_multi, aes(x=imagecol, y=-imagerow)) + geom_point(size = 0.75, position = position_jitter(width = 8, height = 8), aes(color = ct)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed() + guides(fill="none"))
dev.off()

# ---------------- #
# Primary Clusters #
# ---------------- #
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
for(bb_clust in bb_convert15$new) {
  if (bb_clust != "Unassigned") {
    print(bb_clust)
    bb_clust = str_replace(bb_clust, "_", "-")
    fname = bb_clust; fname = str_replace(fname, "2-OPC/Oligo", "2"); fname = str_replace(fname, "15-GABA\\/Glut", "15"); fname = str_replace(fname, "1-RGC\\/MG", "1")
    print(bb_clust)
    Cairo::Cairo(file = paste0(out_dir, "all_bb15_integration_", fname, ".png"), width = 1800, height = 1800, res = 150)
    p1 = SpatialFeaturePlot(all_merge, features = bb_clust, pt.size.factor = 3) + plot_layout(ncol = 4)
    p2 = FeaturePlot(all_merge, features = bb_clust, pt.size = 0.7, order = T) + theme_void() + NoLegend()
    print(p1 + p2)
    dev.off()
  }
}

# DimPlot of spots painted by their predicted bb cluster
Cairo::Cairo(file = paste0(out_dir, "all_bb15_integration.png"), width = 2000, height = 1800, res = 150)
print(DimPlot(all_merge, label = T, pt.size = 2, label.size = 6, label.box = F) + theme_void())
dev.off()

# ------------------ #
# Secondary Clusters #
# ------------------ #
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

for(bb_clust in "15.5_GABA/Glut") {
  if (bb_clust != "Unassigned") {
    print(bb_clust)
    bb_clust = str_replace(bb_clust, "_", "-")
    fname = bb_clust; fname = str_replace(fname, "15.1-GABA\\/Glut", "15.1"); fname = str_replace(fname, "15.5-GABA\\/Glut", "15.5");
    print(bb_clust)
    Cairo::Cairo(file = paste0(out_dir, "all_bb53_integration_", fname, ".png"), width = 1800, height = 1800, res = 150)
    p1 = SpatialFeaturePlot(all_merge, features = bb_clust, pt.size.factor = 3) + plot_layout(ncol = 4)
    p2 = FeaturePlot(all_merge, features = bb_clust, pt.size = 0.7, order = T) + theme_void() + NoLegend()
    print(p1 + p2)
    dev.off()
  }
}

# DimPlot of spots painted by their predicted bb cluster
Cairo::Cairo(file = paste0(out_dir, "all_bb53_integration.png"), width = 2000, height = 1800, res = 150)
print(DimPlot(all_merge, label = T, pt.size = 2, label.size = 6, label.box = F) + theme_void())
dev.off()

# --------------- #
# RGC Subclusters #
# --------------- #
# rgc_sub_sct = SCTransform(rgc_sub)
bb$rgc_sub = bb$names53
bb$rgc_sub[colnames(rgc_sub)] = paste0("RGC", rgc_sub$seurat_clusters)
anchors = FindTransferAnchors(reference = bb, query = all_merge, normalization.method = "SCT", npcs = 50)
predictionsRGC = TransferData(anchorset = anchors, refdata = bb$rgc_sub, prediction.assay = TRUE, weight.reduction = all_merge[["pca"]], dims = 1:50)
all_merge[["predictionsRGC"]] = predictionsRGC

DefaultAssay(all_merge) = "predictionsRGC"
all_merge$bbrgc = GetTransferPredictions(all_merge, assay = "predictionsRGC")

for(bb_clust in paste0("RGC", as.character(0:11))) {
  if (bb_clust != "Unassigned") {
    print(bb_clust)
    Cairo::Cairo(file = paste0(out_dir, "all_bbrgc_integration_test_", bb_clust, ".png"), width = 1800, height = 1800, res = 150)
    p1 = SpatialFeaturePlot(all_merge, features = bb_clust, pt.size.factor = 3) + plot_layout(ncol = 4)
    p2 = FeaturePlot(all_merge, features = bb_clust, pt.size = 0.7, order = T) + theme_void() + NoLegend()
    print(p1 + p2)
    dev.off()
  }
}

all_merge$bb.best = apply(all_merge@assays$predictionsRGC[, 1:ncol(all_merge)], 2, which.max)
all_merge$bb.best = rownames(all_merge@assays$predictionsRGC)[all_merge$bb.best]
all_merge$bb.best = factor(all_merge$bb.best, levels = c(str_replace(convert53$new[3:52], "_", "-"), "RGC2", "RGC5", "RGC6", "RGC8") )
Idents(all_merge) = all_merge$bb.best
DimPlot(all_merge, label = T, pt.size = 0.8)
SpatialDimPlot(all_merge, label = F, pt.size.factor = 3, cells.highlight = colnames(all_merge)[which(all_merge$bb.best == "RGC2")]) + NoLegend() + ggtitle("RGC2") + plot_layout(ncol = 4)

rgc1.deg = bbrgc.deg[which(bbrgc.deg$cluster == 1),]
bb$calc_fc = 0
bb$calc_fc[colnames(rgc_sub)[which(rgc_sub$seurat_clusters == 1)]] = 1
Idents(bb) = bb$calc_fc
rgc1.stats = FoldChange(bb, ident.1 = 1, ident.2 = 0)
rgc1.deg = rgc1.deg[which(rgc1.deg$gene %in% rownames(rgc1.stats)[which(rgc1.stats$pct.2 < 0.20)]),]

rgc2.deg = bbrgc.deg[which(bbrgc.deg$cluster == 1),]
bb$calc_fc = 0
bb$calc_fc[colnames(rgc_sub)[which(rgc_sub$seurat_clusters == 2)]] = 1
Idents(bb) = bb$calc_fc
rgc2.stats = FoldChange(bb, ident.1 = 1, ident.2 = 0)
rgc2.deg = rgc2.deg[which(rgc2.deg$gene %in% rownames(rgc2.stats)[which(rgc2.stats$pct.2 < 0.20)]),]

genes.in.both = rgc1.deg$gene[which(rgc1.deg$gene %in% rgc2.deg$gene)]
rgc1.deg = rgc1.deg[which(!rgc1.deg$gene %in% genes.in.both),]
rgc2.deg = rgc2.deg[which(!rgc2.deg$gene %in% genes.in.both),]

# ------------------ #
# Overlap of Markers #
# ------------------ #
st.umap.2 = read.csv("~/research/st/results/all_merge_umap2_loose_degs_raw.csv")
st.umap.2 = st.umap.2[which(st.umap.2$p_val_adj < 0.05 & abs(st.umap.2$avg_log2FC) > 0.1 & st.umap.2$pct.1 > 0.05),]
bb15.deg = read.csv("~/research/brain/data/bb_all_cluster_15_degs.csv")
bb15.deg = bb15.deg[which(bb15.deg$p_val_adj < 0.05 & abs(bb15.deg$avg_logFC) > 0.1 & bb15.deg$avg_logFC > 0 & bb15.deg$pct.1 > 0.05),]
bb15.deg$old = bb15.deg$cluster
bb15.deg$new = factor(bb_convert15$new[match(bb15.deg$cluster, bb_convert15$old)], levels = convert15$new.full)
bb15.deg$cluster = bb15.deg$new
bb53.deg = read.csv("~/research/brain/data/bb_all_cluster_53_degs.csv")
bb53.deg = bb53.deg[which(bb53.deg$p_val_adj < 0.05 & abs(bb53.deg$avg_logFC) > 0.1 & bb53.deg$avg_logFC > 0 & bb53.deg$pct.1 > 0.05),]
bb53.deg$old = bb53.deg$cluster
bb53.deg$new = factor(bb_convert53$new[match(bb53.deg$cluster, bb_convert53$old)], levels = rev(convert53$new))
bb53.deg$cluster = bb53.deg$new
bbrgc.deg = read.csv("~/research/st/results/rgc_sub_markers_for_brianna.csv")
bbrgc.deg = bbrgc.deg[which(bbrgc.deg$p_val_adj < 0.05 & abs(bbrgc.deg$avg_log2FC) > 0.1 & bbrgc.deg$avg_log2FC > 0 & bbrgc.deg$pct.1 > 0.05),]


# Number of Markers per cluster
ggplot(as.data.frame(table(st.umap.2$cluster)), aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cluster") + ylab("Number of Markers") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + NoLegend() + ggtitle("ST - UMAP2") 
ggplot(as.data.frame(table(bb15.deg$cluster)), aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cluster") + ylab("Number of Markers") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + NoLegend() + ggtitle("BB - Primary")
ggplot(as.data.frame(table(bb53.deg$cluster)), aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cluster") + ylab("Number of Markers") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + NoLegend() + ggtitle("BB - Secondary")

st.umap.2$hgnc = st.umap.2$gene
bb15.deg$hgnc = bb15.deg$gene
bb15.deg$cluster = as.vector(bb15.deg$cluster)
bb53.deg$hgnc = bb53.deg$gene
bb53.deg$cluster = as.vector(bb53.deg$cluster)
bbrgc.deg$hgnc = bbrgc.deg$gene
bbrgc.deg$cluster = as.vector(bbrgc.deg$cluster)
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

st.umap.2$cluster = paste0("ST.", st.umap.2$cluster)
bbrgc.deg$cluster = paste0("BBRGC.", bbrgc.deg$cluster)
my.list = list(st.umap.2, bbrgc.deg)
names(my.list) = c("ST", "BBRGC")
bh_outrgc = bigHeatmap(my.list, pdf.name = "st_w_bbrgc.pdf", single.org.1 = "ST", single.org.2 = "BBRGC", cluster.includes.org = T)
bh_outrgc_1 = bh_outrgc[[1]]
bh_outrgc_1$org.cluster.2 = factor(bh_outrgc_1$org.cluster.2, levels = convertrgc$new)

Cairo::Cairo(file = paste0(out_dir, "all_bbrgc_marker_overlap_heatmap_num.png"), width = 2000, height = 1000, res = 200)
ggplot(bh_outrgc_1, aes(org.cluster.2, org.cluster.1, fill=num)) + geom_tile() + scale_fill_viridis_c(name = "") + guides(color = 'none') + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ggtitle("Number of Overlapping Markers")
dev.off()
Cairo::Cairo(file = paste0(out_dir, "all_bbrgc_marker_overlap_heatmap_pct.png"), width = 2000, height = 1000, res = 200)
ggplot(bh_outrgc_1, aes(org.cluster.2, org.cluster.1, fill=pct.both)) + geom_tile() + scale_fill_viridis_c(name = "") + guides(color = 'none') + theme_classic() + coord_fixed() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ggtitle("Percent of Overlapping Markers")
dev.off()

# 8.1_Glut, 8.4_Glut, RGC2 Figure
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

this.pal = c("#ed2828", "#3F4788", "#FDE725")
my.thresh = 50
df84_81_2 = data.frame(class1 = as.vector(range01(all_merge@assays$predictionsRGC["RGC1",]))*100, class2 = as.vector(range01(all_merge@assays$predictionsRGC["RGC2",]))*100, class81 = as.vector(range01(all_merge@assays$predictionsRGC["8.1-Glut",]))*100, class84 = as.vector(range01(all_merge@assays$predictionsRGC["8.4-Glut",]))*100)
# df84_81_2$class = colnames(df84_81_2)[max.col(df84_81_2, ties.method="first")]
# df_col = data.frame(class = c("class2", "class81", "class84"), col = this.pal)
# df84_81_2$class_col = df_col$col[match(df84_81_2$class, df_col$class)]
# df84_81_2$alpha = as.character(round(apply(df84_81_2[, c("class2", "class81", "class84")], 1, max)))
# df84_81_2$alpha[which(nchar(df84_81_2$alpha) == 1)] = paste0("0", df84_81_2$alpha[which(nchar(df84_81_2$alpha) == 1)])
# df84_81_2$alpha[which(df84_81_2$alpha == "100")] = "ff"
# df84_81_2$col = paste0(df84_81_2$class_col, df84_81_2$alpha)
df84_81_2$dif2_1 = df84_81_2$class2 - df84_81_2$class1
# df84_81_2$dif2_1_col = colorRampPalette(brewer.pal(11, "Spectral"))(100)[cut(df84_81_2$dif2_1,100)]
df84_81_2$dif2_1_col = magma(100)[cut(df84_81_2$dif2_1,100)]
df84_81_2$dif84_81 = df84_81_2$class84 - df84_81_2$class81
df84_81_2$dif84_81_col = viridis(100)[cut(df84_81_2$dif84_81,100)]
df84_81_2$dif84_81_col[which(df84_81_2$dif84_81 == 0)] = "#FDE72500"
df84_81_2$dif84_81_col[which(df84_81_2$class2 >= my.thresh)] = "#FF9E24"
df84_81_2$scale_col = df84_81_2$dif84_81_col
df84_81_2$scale_col[which(df84_81_2$class2 >= my.thresh | df84_81_2$class1 >= my.thresh)] = df84_81_2$dif2_1_col[which(df84_81_2$class2 >= my.thresh | df84_81_2$class1 >= my.thresh)]
all_merge_hi$class2_81_84 = df84_81_2$dif84_81_col

pdf(paste0(out_dir, paste0("class2_81_84_4_", my.thresh, ".pdf")), width = 12, height = 12)
print(allSamplesSFP(all_merge_hi, "class2_81_84", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100)), rm.zero = T, col.ident = T, high.res = T ))
dev.off()

# ---------------- #
# Zack Rabbit Hole #
# -----------------#
# Trying to distinguish b/w RGC1 and RGC2
bb$tmp = "other"
bb$tmp[colnames(rgc_sub)[which(rgc_sub$seurat_clusters == 1)]] = "1"
bb$tmp[colnames(rgc_sub)[which(rgc_sub$seurat_clusters == 2)]] = "2"
deg1 = FindAllMarkers(bb, min.pct = 0.05, logfc.threshold = 0.1)
Idents(bb) = bb$tmp
deg1 = FindAllMarkers(bb, min.pct = 0.05, logfc.threshold = 0.1, only.pos = T)
deg1.1 = deg1[which(deg1$p_val_adj < 0.05 & deg1$avg_log2FC > 0 & deg1$cluster == "1"),]
deg1.2 = deg1[which(deg1$p_val_adj < 0.05 & deg1$avg_log2FC > 0 & deg1$cluster == "2"),]
common.genes = deg1.1$gene[which(deg1.1$gene %in% deg1.2$gene)]
deg1.1.no2 = deg1.1[which(!deg1.1$gene %in% common.genes),]
deg1.2.no1 = deg1.2[which(!deg1.2$gene %in% common.genes),]

all_merge$rgc1.top = colSums(all_merge@assays$Spatial@counts[deg1.1.no2$gene[1:75],] > 0)
all_merge$rgc2.top = colSums(all_merge@assays$Spatial@counts[deg1.2.no1$gene[1:75],] > 0)
allSamplesSFP(all_merge, "rgc1.top")

bb$rgc1.top = colSums(bb@assays$RNA@counts[deg1.1.no2$gene[1:75],] > 0)
bb$rgc2.top = colSums(bb@assays$RNA@counts[deg1.2.no1$gene[1:75],] > 0)
FeaturePlot(bb, "rgc1.top", order = T, pt.size = 0.8)

bb$tmp1 = bb$tmp
bb$tmp1[which(bb$tmp == 2)] = "other"
bb$tmp2 = bb$tmp
bb$tmp2[which(bb$tmp == 1)] = "other"
Idents(bb) = bb$tmp2
rgc2.stats.for.rgc1deg = FoldChange(bb, ident.1 = "2", ident.2 = "other", features = deg1.1.no2$gene)
Idents(bb) = bb$tmp1
rgc1.stats.for.rgc2deg = FoldChange(bb, ident.1 = "1", ident.2 = "other", features = deg1.2.no1$gene)

deg1.1.no2[, paste0("rgc2.", colnames(rgc2.stats.for.rgc1deg))] = rgc2.stats.for.rgc1deg
deg1.2.no1[, paste0("rgc1.", colnames(rgc1.stats.for.rgc2deg))] = rgc1.stats.for.rgc2deg

zrgc = read.csv("~/Downloads/rgc_subclusters123_for_brianna_062922.csv")
zrgc = zrgc[order(zrgc$zack_marker_score, decreasing = T),]
for (clust in 1:2) {
  for (gene in zrgc$gene[which(zrgc$rgc_subcluster == clust)]) {
    p1 = FeaturePlot(bb, gene, order = T, pt.size = 0.1) + theme_void() + NoLegend() + ggtitle(paste0(clust , ": ", zrgc$zack_marker_score[which(zrgc$gene == gene)])) 
    p2 = FeaturePlot(rgc_sub, gene, order = T, pt.size = 0.75) + theme_void() + NoLegend() + theme(title = element_blank()) 
    Cairo::Cairo(file = paste0("zrgc_", clust, "_", gene, "_spatial_peripheral.png"), width = 1200, height = 1200, res = 150)
    print(allSamplesSFP(all_merge_rgc, gene, pt.size.multiplier = 0.7))
    dev.off()
    # Cairo::Cairo(file = paste0("zrgc_", clust, "_", gene, "_sn.png"), width = 600, height = 1200, res = 150)
    # print(wrap_plots(p1, p2, ncol = 1))
    # dev.off()
  }
}

bb15.deg4 = bb15.deg[which(bb15.deg$cluster == 4),]
all_merge$rgc_score = colSums(all_merge@assays$Spatial@counts[bb15.deg4$gene[order(bb15.deg4$pct.dif, decreasing = T)[1:20]],] > 0)
all_merge$rgc_score = colSums(all_merge@assays$Spatial@counts[bb15.deg4$gene[1:20],] > 0)
allSamplesSFP(all_merge, "rgc_score")

all_merge_rgc = subset(all_merge, cells = colnames(all_merge)[which(all_merge$rgc_score > 6)])
SpatialFeaturePlot(all_merge_rgc, "nCount_Spatial", pt.size.factor = 3) + plot_layout(ncol = 4)

bb$zdef = bb$good_names
bb$zdef[which(bb$zdef %in% c("1.1_Astro", "1.2_Astro"))] = "RGC"
bb$zdef[colnames(rgc_sub)[which(rgc_sub$seurat_clusters %in% c(1, 2))]] = "RGC12"
Idents(bb) = bb$zdef
DimPlot(bb, label = T) + NoLegend()

all_merge_rgc12 = subset(all_merge, cells = colnames(all_merge)[which(all_merge@assays$predictionsZDEF["RGC12",] > 0.4)])
all_merge_rgc12$rgc1.top.znew = colSums(all_merge@assays$Spatial@counts[c("LOC101464051", "cdon", "LOC101465090", "LOC101478557", "cep112"),] > 0)
all_merge_rgc12$rgc2.top.znew = colSums(all_merge@assays$Spatial@counts[c("LOC101467223", "LOC101482938", "LOC101466608", "nr2f2", "LOC101475144"),] > 0)
allSamplesSFP(all_merge_rgc12, "LOC101475144", pt.size.multiplier = 0.7)
allSamplesSFP(all_merge_rgc12, "rgc1.top.znew", assay = "predictionsZDEF", pt.size.multiplier = 0.7)

rgc_sub12 = rgc_sub[,which(rgc_sub$seurat_clusters %in% c(1,2))]
rgc_sub12 = SCTransform(rgc_sub12)
all_merge_rgc12@active.assay = "SCT"
anchors = FindTransferAnchors(reference = rgc_sub12, query = all_merge_rgc12, normalization.method = "SCT", npcs = 50)
predictions12 = TransferData(anchorset = anchors, refdata = rgc_sub12$seurat_clusters, prediction.assay = TRUE, weight.reduction = all_merge_rgc12[["pca"]], dims = 1:50)
all_merge_rgc12[["predictions12"]] = predictions12
allSamplesSFP(all_merge_rgc12, "1", assay = "predictions12", pt.size.multiplier = 0.7)

anchor.fc1 = anchor.fc1[order(anchor.fc1$avg_log2FC, decreasing = T),]
all_merge$tmp = colSums(all_merge@assays$Spatial@counts[head(rownames(anchor.fc1), 10),] > 0)
all_merge$tmp2 = colSums(all_merge@assays$Spatial@counts[tail(rownames(anchor.fc1), 10),] > 0)
allSamplesSFP(all_merge, "tmp", pt.size.multiplier = 0.7)
allSamplesSFP(all_merge, "tmp2", pt.size.multiplier = 0.7)
rgc_sub$tmp = colSums(rgc_sub@assays$RNA@counts[head(rownames(anchor.fc1), 10),] > 0)
rgc_sub$tmp2 = colSums(rgc_sub@assays$RNA@counts[tail(rownames(anchor.fc1), 10),] > 0)
FeaturePlot(rgc_sub, "tmp", order = T, pt.size = 1.5) + scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))
FeaturePlot(rgc_sub, "tmp2", order = T, pt.size = 1.5) + scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))
bb$tmp = colSums(bb@assays$RNA@counts[head(rownames(anchor.fc1), 10),] > 0)
bb$tmp2 = colSums(bb@assays$RNA@counts[tail(rownames(anchor.fc1), 10),] > 0)
FeaturePlot(bb, "tmp", order = T, pt.size = 0.8) + scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))
FeaturePlot(bb, "tmp2", order = T, pt.size = 0.8) + scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))

# Hippocampus kernel density
slide.seq = qs::qread("~/research/st/data/sshippo_processed.qs")
zdeg = read.csv("~/Downloads/rgc12_8_1_and_8_4_Glut_markers_071622_minpct20.csv")
zdeg$hgnc = gene_info$human[match(zdeg$gene, gene_info$mzebra)]
zdeg$hgnc[which(is.na(zdeg$hgnc))] = toupper(zdeg$gene[which(is.na(zdeg$hgnc))])
glut81_hgnc = zdeg$hgnc[which(zdeg$cluster == 4 & zdeg$hgnc %in% rownames(slide.seq))]
glut84_hgnc = zdeg$hgnc[which(zdeg$cluster == 13 & zdeg$hgnc %in% rownames(slide.seq))]
slide.seq$glut81 = colSums(slide.seq@assays$Spatial@counts[glut81_hgnc,] > 0) >= 2
slide.seq$glut84 = colSums(slide.seq@assays$Spatial@counts[glut84_hgnc,] > 0) >= 2

# gene = "COL12A1"
coords = GetTissueCoordinates(object = slide.seq)
coords$value = slide.seq$glut81
# coords$value = slide.seq@assays$Spatial@counts[gene,] > 0
coords = coords[order(coords$value, decreasing = F),]
ggplot(coords, aes(x=y, y=-x, color = value, alpha = value > 0)) + geom_point(size = 0.5) + scale_alpha_manual(values = c(0.03, 1), guide = 'none') + scale_color_manual(values = c("gray40", "black")) + theme_void() + coord_fixed()
coords$value = as.numeric(coords$value)
coords = coords[,c("x", "y", "value")]

this_grid_length = 75
num_neighbors = 4
this_interpolate = F
plot_alpha = 0.9
this_alpha_hull = 5

# KNN Smoothing
this_grid = expand.grid(seq(min(coords$x), max(coords$x), length = this_grid_length), seq(min(coords$y), max(coords$y), length = this_grid_length))
colnames(this_grid) = c("x", "y")
test_knn <- kknn::kknn(value ~ ., train = coords, test = this_grid, kernel = "gaussian", k = num_neighbors)
# for (row in 1:nrow(this_grid)) {
#   for (col in 1:ncol(this_grid)) {
#     
#   }
# }
this_grid$value = test_knn$fitted.values

c1 = c(mean(coords$x), mean(coords$y))
coords_max_dist = max(apply(coords[, c("x", "y")],1,function(x,c1) {(sqrt((x[1] - c1[1])^2+(x[2]-c1[2])^2))},c1))
c1 = c(mean(this_grid$x), mean(this_grid$y))
this_grid$dist = apply(this_grid[, c("x", "y")],1,function(x,c1) {(sqrt((x[1] - c1[1])^2+(x[2]-c1[2])^2))},c1)
this_grid_relevant = this_grid[which(this_grid$dist <= coords_max_dist),]

# this_grid_relevant = this_grid
# temp = rev(brewer.pal(11,"Spectral"))
temp = viridis(100)
pal = colorRampPalette(temp)
ggplot(this_grid_relevant, aes(x=y, y=-x, color = value, fill = value)) + geom_raster(alpha = plot_alpha, interpolate = this_interpolate) + scale_fill_gradientn(colors=pal(50)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed()

# Hippocampus sections
slide.seq <- SeuratData::LoadData("ssHippo")
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)

brain <- SeuratData::LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", ncells = 3000, verbose = FALSE)
brain <- RunPCA(brain)
brain <- RunUMAP(brain, dims = 1:30)
brain <- FindNeighbors(brain, dims = 1:30)
brain <- FindClusters(brain, resolution = 0.3, verbose = FALSE)

post <- SeuratData::LoadData("stxBrain", type = "posterior1")
post <- SCTransform(post, assay = "Spatial", ncells = 3000, verbose = FALSE)
post <- RunPCA(post)
post <- RunUMAP(post, dims = 1:30)
post <- FindNeighbors(post, dims = 1:30)
post <- FindClusters(post, resolution = 0.3, verbose = FALSE)

slide.seq = qs::qread("~/scratch/brain/data/sshippo_processed.qs")
bb.mouse = qs::qread("/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/data/bcs/data/cichlid_human.qs")
bb.mouse = SCTransform(bb.mouse)
anchors <- FindTransferAnchors(reference = bb.mouse, query = slide.seq, normalization.method = "SCT", npcs = 50)
predictions53 <- TransferData(anchorset = anchors, refdata = bb.mouse$names53, prediction.assay = TRUE, weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions53"]] <- predictions53

anchors <- FindTransferAnchors(reference = bb.mouse, query = brain, normalization.method = "SCT", npcs = 50)
predictions53 <- TransferData(anchorset = anchors, refdata = bb.mouse$names53, prediction.assay = TRUE, weight.reduction = brain[["pca"]], dims = 1:2)
brain[["predictions53"]] <- predictions53

anchors <- FindTransferAnchors(reference = bb.mouse, query = post, normalization.method = "SCT", npcs = 50)
predictions15 <- TransferData(anchorset = anchors, refdata = bb.mouse$seuratclusters15, prediction.assay = TRUE, weight.reduction = post[["pca"]], dims = 1:2)
post[["predictions15"]] <- predictions15


all_merge_hgnc = convertSeuratObject(all_merge, ortho[[vstr]], assay = "Spatial")
for (i in names(all_merge@images)) { all_merge_hgnc[[i]] = all_merge[[i]] }
all_merge_hgnc = qs::qread("~/scratch/brain/sp_results/st_hgnc_070822.qs")
all_merge_hgnc = SCTransform(all_merge_hgnc, assay = "Spatial")
all_merge_hgnc = RunPCA(all_merge_hgnc, assay = "SCT", verbose = FALSE)
all_merge_hgnc = RunUMAP(all_merge_hgnc, reduction = "pca", dims = 1:30)
all_merge_hgnc = FindNeighbors(all_merge_hgnc, reduction = "umap", dims = 1:2)
all_merge_hgnc = FindClusters(all_merge_hgnc, verbose = FALSE, resolution = 0.55)
anchors <- FindTransferAnchors(reference = bb.mouse, query = all_merge_hgnc, normalization.method = "SCT", npcs = 50)
predictions53 <- TransferData(anchorset = anchors, refdata = bb.mouse$names53, prediction.assay = TRUE, weight.reduction = all_merge_hgnc[["pca"]], dims = 1:50)
all_merge_hgnc[["predictions53"]] <- predictions53

DefaultAssay(slide.seq) = "predictions53"
# for (i in c("8.1_Glut", "8.4_Glut", "8.6_Glut")) {
for (i in unique(bb.mouse$seuratclusters15)) {
  print(i)
  i = str_replace(i, "_", "-")
  # slide.seq$highlight = slide.seq@assays$predictions53@data[as.character(i),] > 0.05
  # coords = GetTissueCoordinates(object = slide.seq)
  # coords$value = factor(as.character(slide.seq$highlight), levels = c("TRUE", "FALSE"))
  Cairo::Cairo(file = paste0("~/scratch/brain/sp_results/post/sshippo_", i, ".png"), width = 1600, height = 1600, res = 150)
  print(SpatialFeaturePlot(post, features = as.character(i), alpha = c(0.3, 0.9)) + scale_fill_gradientn(colors = viridis(100)))
  # coords = GetTissueCoordinates(object = slide.seq)
  # coords$value = colSums(slide.seq@assays$predictions15[as.character(i),])
  # coords$value = slide.seq@assays$predictions15[as.character(i),] > 0
  # coords = coords[order(coords$value, decreasing = F),]
  # print(ggplot(coords[which(coords$value > 0),], aes(x=y, y=-x, color = value)) + geom_point(size = 0.5) + scale_color_viridis() + theme_void() + coord_fixed())
  # print(ggplot(coords, aes(x=y, y=-x, color = value, alpha = value)) + geom_point(size = 0.5) + scale_color_manual(values = c("black", "white")) + scale_alpha_manual(values = c(1, 0)) + theme_void() + coord_fixed() + NoLegend())
  dev.off()
}

bb15.deg = read.csv("~/scratch/brain/data/bb_all_markers_15clusters_102820_more_info.csv")
bb53.deg = read.csv("~/scratch/brain/data/bb_all_markers_53clusters_102720_more_info.csv")

bb15.deg$mouse = str_to_title(bb15.deg$human)
bb53.deg$mouse = str_to_title(bb53.deg$human)

deg4  = bb53.deg[which(bb53.deg$cluster == 4  & !is.na(bb53.deg$human) & bb53.deg$human %in% rownames(slide.seq) & bb53.deg$p_val_adj < 0.05 & bb53.deg$pct.1 > 0.4 & bb53.deg$pct.2 < 0.2),]
deg13 = bb53.deg[which(bb53.deg$cluster == 13 & !is.na(bb53.deg$human) & bb53.deg$human %in% rownames(slide.seq) & bb53.deg$p_val_adj < 0.05 & bb53.deg$pct.1 > 0.4 & bb53.deg$pct.2 < 0.2),]

common_degs = deg4$gene[which(deg4$gene %in% deg13$gene)]
deg4 = deg4[which(!deg4$gene %in% common_degs),]
deg13 = deg13[which(!deg13$gene %in% common_degs),]

bb$ident4  = bb$seuratclusters53 == 4
bb$ident13 = bb$seuratclusters53 == 13
Idents(bb) = bb$ident13
deg4_fc_in_13 = FoldChange(bb, ident.1 = "TRUE", ident.2 = "FALSE", features = deg4$gene)
Idents(bb) = bb$ident4
deg13_fc_in_4 = FoldChange(bb, ident.1 = "TRUE", ident.2 = "FALSE", features = deg13$gene)
deg4[,paste0("in_13_", colnames(deg4_fc_in_13))] = deg4_fc_in_13
deg13[,paste0("in_4_", colnames(deg13_fc_in_4))] = deg13_fc_in_4

deg4 = deg4[which(deg4$in_13_avg_log2FC < 0   & deg4$in_13_pct.1 < 0.2),]
deg13 = deg13[which(deg13$in_4_avg_log2FC < 0 & deg13$in_4_pct.1 < 0.2),]

bb$deg4.norm = colSums(bb@assays$RNA@counts[deg84_81_84$gene[1:10],] > 0)
bb$deg4 = colSums(bb@assays$RNA@counts[deg4$gene[1:5],])
FeaturePlot(bb, "deg4.norm", order = T) + theme_void() + theme(title=element_blank()) + coord_fixed()

bb$deg13.norm = colSums(bb@assays$RNA@counts[deg13$gene[1:5],] > 0)
bb$deg13 = colSums(bb@assays$RNA@counts[deg84_81_84$gene[1:5],])
FeaturePlot(bb, "deg13.norm", order = T) + theme_void() + theme(title=element_blank()) + coord_fixed()

gene = "COL12A1"
coords = GetTissueCoordinates(object = slide.seq)
coords$value = slide.seq@assays$SCT@data[gene,]
coords = coords[order(coords$value, decreasing = F),]
ggplot(coords[which(coords$value > 0),], aes(x=y, y=-x, color = value)) + geom_point(size = 0.5) + scale_color_viridis() + theme_void() + coord_fixed() + NoLegend() + ggtitle(gene)

coords = GetTissueCoordinates(object = slide.seq)
# coords$value = colSums(slide.seq@assays$Spatial@counts[deg84_81_81$human[1:5],] > 0)
coords$value = colSums(slide.seq@assays$Spatial@data[deg84_81_84$human[1:5],])
coords = coords[order(coords$value, decreasing = F),]
ggplot(coords, aes(x=y, y=-x, color = value, alpha = value > 0)) + geom_point(size = 0.5) + scale_alpha_manual(values = c(0.03, 1), guide = 'none') + scale_color_viridis() + theme_void() + coord_fixed()

brain$value = colSums(brain@assays$Spatial@counts[deg84_81_84$mouse[1:10],] > 0)
SpatialFeaturePlot(brain, "value") + scale_fill_gradientn(colors = viridis(100)) + theme(legend.position="right")

post$value = colSums(post@assays$Spatial@counts[deg13$mouse[1:5],] > 0)
SpatialFeaturePlot(post, "value") + scale_fill_gradientn(colors = viridis(100)) + theme(legend.position="right")

bb$glut84_81 = "none"
bb$glut84_81[which(bb$seuratclusters53 == 4)] = "81_Glut"
bb$glut84_81[which(bb$seuratclusters53 == 13)] = "84_Glut"
Idents(bb) = bb$glut84_81
deg84_81 = FindMarkers(bb, ident.1 = "81_Glut", ident.2 = "84_Glut", min.pct = 1e-100, logfc.threshold = 0, only.pos = F)
bb$glut84_and_81 = "none"
bb$glut84_and_81[which(bb$seuratclusters53 %in% c(4, 13))] = "84_81"
Idents(bb) = bb$glut84_and_81
deg84_81$gene = rownames(deg84_81)
deg84_81_in_other = FoldChange(bb, ident.1 = "none", ident.2 = "84_81", features = deg84_81$gene)
deg84_81[,paste0("other", colnames(deg84_81_in_other))] = deg84_81_in_other
deg84_81$human = gene_info$human[match(deg84_81$gene, gene_info$mzebra)]
deg84_81$mouse = str_to_title(deg84_81$human)

deg84_81_81 = deg84_81[which(deg84_81$p_val_adj < 0.05 & !is.na(deg84_81$human) & deg84_81$human %in% rownames(slide.seq) & deg84_81$pct.1 > 0.4 & deg84_81$pct.2 < 0.3 & deg84_81$otherpct.1 < 0.2),]
deg84_81_84 = deg84_81[which(deg84_81$p_val_adj < 0.05 & !is.na(deg84_81$human) & deg84_81$human %in% rownames(slide.seq) & deg84_81$pct.2 > 0.4 & deg84_81$pct.1 < 0.3 & deg84_81$otherpct.1 < 0.2),]

#*******************************************************************************
# BHVE vs CTRL DEGs ============================================================
#*******************************************************************************
b2.samples = c("b2a", "b2b", "b2c", "b2d")
c2.samples = c("c2a", "c2b", "c2c", "c2d")
all_merge_b2_c2 = subset(all_merge, cells = colnames(all_merge)[which(all_merge$sample %in% c(b2.samples, c2.samples))])

# --------------------------- #
# Real B2 vs C2 Cluster Level #
# --------------------------- #
all_merge_b2_c2$cluster_cond = paste0(all_merge_b2_c2$cluster, "_", all_merge_b2_c2$cond)
Idents(all_merge_b2_c2) = all_merge_b2_c2$cluster_cond
b2c2.deg.cluster = data.frame()
for (i in 0:31) {
  if (length(which(Idents(all_merge_b2_c2) == paste0(i, "_BHVE"))) >= 3 && length(which(Idents(all_merge_b2_c2) == paste0(i, "_CTRL"))) >= 3) {
    print(i)
    this.deg = FindMarkers(all_merge_b2_c2, ident.1 = paste0(i, "_BHVE"), ident.2 = paste0(i, "_CTRL"), only.pos = F, min.pct = 1e-100, logfc.threshold = 0)
    this.deg$gene = rownames(this.deg)
    this.deg$cluster = i
    b2c2.deg.cluster = rbind(b2c2.deg.cluster, this.deg)
  }
}

b2c2.deg.cluster.strict = b2c2.deg.cluster[which(b2c2.deg.cluster$p_val_adj < 0.05 & abs(b2c2.deg.cluster$avg_log2FC) > 0.1 & b2c2.deg.cluster$pct.1 > 0.1),]

# ---------------------------------------------------------------- #
# Compare Real B2 vs C2 Differences to Permutations of the Samples #
# ---------------------------------------------------------------- #

# Real B2 v2 C2 Differences
viridis.purple = "#440154"; viridis.yellow = "#FDE725";
realb2c2 = FindMarkers(all_merge_b2_c2, ident.1 = "BHVE", ident.2 = "CTRL", only.pos = F, min.pct = 1e-100, logfc.threshold = 0)
realb2c2$perm = "real"
realb2c2$gene = rownames(realb2c2)
realb2c2$neg_log_bon = -log10(realb2c2$p_val_adj)
realb2c2$pct.dif = realb2c2$pct.1 - realb2c2$pct.2
realb2c2$pct.ratio = realb2c2$pct.1 / realb2c2$pct.2
realb2c2$abs.avg_log2FC = abs(realb2c2$avg_log2FC)
realb2c2$col = "gray80"
realb2c2$col[which(realb2c2$p_val_adj < 0.05 & realb2c2$avg_log2FC > 0.1 & realb2c2$pct.1 > 0.05)] = viridis.yellow
realb2c2$col[which(realb2c2$p_val_adj < 0.05 & realb2c2$avg_log2FC < 0.1 & realb2c2$pct.2 > 0.05)] = viridis.purple
realb2c2 = realb2c2[order(realb2c2$col, decreasing = T),]
realb2c2$hgnc = gene_info$human[match(realb2c2$gene, gene_info$mzebra)]
realb2c2$ens = gene_info$ens[match(realb2c2$gene, gene_info$mzebra)]
realb2c2$label = tolower(realb2c2$hgnc)
realb2c2$label[which( is.na(realb2c2$label) )] = realb2c2$gene[which( is.na(realb2c2$label) )]
x.max = max(realb2c2$abs.avg_log2FC)

Cairo::Cairo(file = paste0("b2_c2_real_volcano.png"), width = 900, height = 600, res = 150)
print(ggplot(realb2c2, aes(x = avg_log2FC, y = neg_log_bon, color = col, size = abs(pct.dif))) + geom_point(aes(alpha = neg_log_bon)) + scale_color_identity() + scale_size_continuous(range = c(1, 5), name = "|% Difference|") + scale_alpha_continuous(range = c(0.4, 1), guide = 'none') + xlab(expression("Average -"*Log["2"]*" Fold Change")) + ylab(expression("-"*Log["10"]*"(adjusted p)")) + theme_bw() + scale_y_sqrt(expand = c(0,0)) + scale_x_continuous(limits = c(-x.max, x.max)) + geom_label_repel(data = realb2c2[which(realb2c2$neg_log_bon > 90),], aes(label = label, fill = paste0(col, "40")), color="black", size = 3.25, label.padding = 0.2, nudge_y = 0.07, box.padding = 0.3) + scale_fill_identity())
dev.off()

# Make a heatmap
realb2c2.sig = realb2c2[which(realb2c2$col %in% c(viridis.purple, viridis.yellow)),]
realb2c2.sig = realb2c2.sig[order(realb2c2.sig$col, realb2c2.sig$neg_log_bon, decreasing = T),]
fc.mat = matrix(0L, nrow = nrow(realb2c2.sig), ncol = 8, dimnames = list(realb2c2.sig$gene, c(c2.samples, b2.samples)))
for (s in c(c2.samples, b2.samples)) {
  all_merge$this.s = all_merge$cond
  all_merge$this.s[which(all_merge$sample == s)] = "this.s"
  Idents(all_merge) = all_merge$this.s
  if (startsWith(s, "c")) {
    fc.mat[,s] = FoldChange(all_merge, ident.1 = "this.s", ident.2 = "BHVE", features = realb2c2.sig$gene)[,1]
  } else {
    fc.mat[,s] = FoldChange(all_merge, ident.1 = "this.s", ident.2 = "CTRL", features = realb2c2.sig$gene)[,1]
  }
}

fc.mat2 = fc.mat
fc.mat2[which(fc.mat2 > 0.5)] = 0.5
fc.mat2[which(fc.mat2 < -0.5)] = -0.5
pheatmap::pheatmap(fc.mat2, color = viridis(100), show_rownames = F, angle_col = 0)

# Permutations of Samples
even_perm_deg_df = data.frame()
even_perm_deg_df_num = data.frame()
i = 1
for (bi in 1:4) {
  b2.b1 = b2.samples[bi]
  for (bj in (bi+1):4) {
    if (bj <= 4) {
      
      b2.b2 = b2.samples[bj]
      
      for (ci in 1:4) {
        c2.b3 = c2.samples[ci]
        for (cj in (ci+1):4) {
          c2.b4 = c2.samples[cj]
          if (cj <= 4) {
            print(i)
            this.id = i
            # this.id = paste0("(", bi, ", ",  bj, ")", ", ", "(", ci, ", ",  cj, ")")
            # print(this.id)
            new_b2 = c(b2.b1, b2.b2, c2.b3, c2.b4)
            
            all_merge_b2_c2$tmp = "C"
            all_merge_b2_c2$tmp[which(all_merge$sample %in% new_b2)] = "B"
            Idents(all_merge_b2_c2) = all_merge_b2_c2$tmp
            this.deg = FindMarkers(all_merge_b2_c2, ident.1 = "B", ident.2 = "C", only.pos = T, min.pct = 1e-100, logfc.threshold = 0)
            this.deg$perm = this.id
            this.deg$gene = rownames(this.deg)
            even_perm_deg_df = rbind(even_perm_deg_df, this.deg)
            
            even_perm_deg_df_num = rbind(even_perm_deg_df_num, data.frame(b2.b1, b2.b2, c2.b3, c2.b4, length(which(this.deg$p_val_adj < 0.05 & this.deg$pct.1 >= 0.05 & this.deg$avg_log2FC >= 0.25))))
            
            i = i + 1
          }
        }
      }
      
    }
  }
}

perm_df = read.csv("b2_c2_perm_with_neg_071422.csv")
perm_df$X = NULL
perm_df$perm = paste0("perm: ", perm_df$perm)
perm_df$neg_log_bon = -log10(perm_df$p_val_adj)
perm_and_real = rbind(realb2c2[, colnames(perm_df)], perm_df)
perm_and_real$perm = factor(perm_and_real$perm, levels = c("real", paste0("perm: ", 1:49)))
perm_and_real$col = "gray80"
perm_and_real$col[which(perm_and_real$p_val_adj < 0.05 & perm_and_real$avg_log2FC > 0.1 & perm_and_real$pct.1 > 0.05)] = viridis.yellow
perm_and_real$col[which(perm_and_real$p_val_adj < 0.05 & perm_and_real$avg_log2FC < 0.1 & perm_and_real$pct.2 > 0.05)] = viridis.purple
ggplot(perm_and_real, aes(x = avg_log2FC, y = neg_log_bon, color = col)) + geom_point(size = 0.6) + scale_color_identity() + xlab(expression("Average -"*Log["2"]*" Fold Change")) + ylab(expression("-"*Log["10"]*"(adjusted p)")) + theme_classic() + facet_wrap(~ perm)

#*******************************************************************************
# Cluster markers ==============================================================
#*******************************************************************************
st.clust.deg = read.csv("all_merge_cluster_markers_loose_070822.csv")
st.clust.deg = st.clust.deg[which(st.clust.deg$p_val_adj < 0.05 & abs(st.clust.deg$avg_log2FC) > 0.25 & st.clust.deg$pct.1 > 0.1),]
st.clust.deg$hgnc = gene_info$human[match(st.clust.deg$gene, gene_info$mzebra)]
st.clust.deg$ens = gene_info$ens[match(st.clust.deg$gene, gene_info$mzebra)]
st.clust.deg$avg_logFC = st.clust.deg$avg_log2FC
st.clust.deg$abs.avg_logFC = abs(st.clust.deg$avg_logFC)
st.clust.deg$pct.dif = st.clust.deg$pct.1 - st.clust.deg$pct.2
st.clust.deg$abs.pct.dif = abs(st.clust.deg$pct.dif)
st.clust.deg$cluster_gene = paste0(st.clust.deg$cluster, "_", st.clust.deg$gene)
st.clust.deg = st.clust.deg[which(!duplicated(st.clust.deg$cluster_gene)),]
n_gene_appears = data.frame(table(st.clust.deg$gene))
st.clust.deg$n_gene_appears = n_gene_appears[match(st.clust.deg$gene, n_gene_appears[,1]),2]
appears = data.frame()
for (gene in unique(st.clust.deg$gene)) {
  rows = st.clust.deg[which(st.clust.deg$gene == gene),]
  appears = rbind(appears, t(c(gene, paste(rows$cluster, collapse = ", ") )))
}
st.clust.deg$DEG_in = appears$V2[match(st.clust.deg$gene, appears$V1)]
st.clust.deg$mz_description = gene_info$mzebra_description[match(st.clust.deg$gene, gene_info$mzebra)]
write.csv(st.clust.deg, "all_merge_cluster_markers_standard_070822.csv")

st.clust.deg.sub = st.clust.deg[which(st.clust.deg$n_gene_appears <= 3),]
write.csv(st.clust.deg.sub, "all_merge_cluster_markers_3n_070822.csv")
st.clust.deg.sub = st.clust.deg[which(st.clust.deg$pct.2 < 0.3 & st.clust.deg$pct.1 > 0.5),]
write.csv(st.clust.deg.sub, "all_merge_cluster_markers_subsetpct2_070822.csv")

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

# High Res
all_merge_hi = all_merge
for (s in levels(all_merge_hi$sample)) {
  real.slice = all_merge_hi@images[[s]]
  id = paste0(all_merge_hi$fnum[which(all_merge_hi$sample == s)][1], "_", toupper(all_merge_hi$area[which(all_merge_hi$sample == s)][1]), "1")
  print(paste0(s, ": ", id))
  hires = Read10X_Image(paste0("~/Downloads/sp_data/spatial_modified/", id, "/"), image.name = "tissue_hires_image.png")
  rownames(hires@coordinates) = paste0(s, "_", rownames(hires@coordinates))
  hires@coordinates = hires@coordinates[rownames(real.slice@coordinates),]
  hires@scale.factors$lowres = hires@scale.factors$hires
  hires@assay = real.slice@assay
  hires@key = real.slice@key
  all_merge_hi@images[[s]]= hires
  # SpatialFeaturePlot(b2d, feature = "egr1", images = "slice1", crop = T)
}

pdf(paste0(out_dir, "high_res_cluster.pdf"), width = 12, height = 12, onefile = F)
print(allSamplesSFP( all_merge_hi, "cluster.num", pal = scales::hue_pal(), rm.zero = F, high.res = T, points.as.text = T, pt.size.multiplier = 1.4 ))
dev.off()

#*******************************************************************************
# Gene Info 2 ==================================================================
#*******************************************************************************
pat = read.table("~/research/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", sep = "\t", header = FALSE, fill = TRUE)
gtf = read.table("~/research/all_research/gtf_processed.gtf", sep = "\t", header = FALSE, fill = TRUE)
ncbi_data = data.table::fread("~/research/bcs/data/gene_orthology/cichlid_ncbi_datasets.tsv", data.table = F)
gene_info = data.table::fread("~/research/all_research/gene_info.txt", data.table = F)
gene_info = gene_info[which(! duplicated(gene_info$mzebra) ),]
ncbi_data$LOC = paste0("LOC", ncbi_data$`NCBI GeneID`)
ncbi_data$name = gtf$V10[match(ncbi_data$LOC, gtf$V11)]

# Genes in the genome that weren't included in the first gene info
forgotten_genes = as.data.frame(matrix(NA, nrow = length(which(!gtf$V10 %in% gene_info$mzebra)), ncol = ncol(gene_info)))
colnames(forgotten_genes) = colnames(gene_info)
forgotten_genes$mzebra = gtf$V10[which(!gtf$V10 %in% gene_info$mzebra)]

gene_info2 = rbind(gene_info, forgotten_genes)
gene_info2$mzebra = factor(gene_info2$mzebra, levels = gtf$V10)
gene_info2 = gene_info2[order(gene_info2$mzebra),]
gene_info2$mzebra = as.vector(gene_info2$mzebra)
gene_info2 = gene_info2[, c(1:4, 7:10, 5, 6)]

gene_info2$loc_id = gtf$V11
gene_info2$nd_symbol = ncbi_data$Symbol[match(gene_info2$mzebra, ncbi_data$name)]
gene_info2$nd_description = ncbi_data$Description[match(gene_info2$mzebra, ncbi_data$name)]
gene_info2$nd_symbol_minus1 = substr(gene_info2$nd_symbol,1,nchar(gene_info2$nd_symbol)-1)
gene_info2$ens_clean = str_replace(gene_info2$ens, " (1 of many)", "")
gene_info2$ens_minus1 = substr(gene_info2$ens_clean,1,nchar(gene_info2$ens_clean)-1)

gene_info2$mzebra_desc_clean = reshape2::colsplit(gene_info2$mzebra_description, " \\[", c('1','2'))[,1]
gene_info2$human_desc_clean = reshape2::colsplit(gene_info2$human_description, " \\[", c('1','2'))[,1]
gene_info2$mz_human_desc_match = gene_info2$mzebra_desc_clean != "" & (gene_info2$mzebra_desc_clean == gene_info2$human_desc_clean)
# gene_info2$mz_in_human_desc = sapply(1:nrow(gene_info2), function(x) gene_info2$mzebra_desc_clean[x] != "" & grepl(gene_info2$mzebra_desc_clean[x], gene_info2$human_desc_clean[x], fixed = TRUE))
# gene_info2$human_in_mz_desc = sapply(1:nrow(gene_info2), function(x) gene_info2$human_desc_clean[x] != "" & grepl(gene_info2$human_desc_clean[x], gene_info2$mzebra_desc_clean[x], fixed = TRUE))

gene_info2$human_agrees = gene_info2$human_mart == gene_info2$human_pat
gene_info2$something_agrees = gene_info2$human_agrees
gene_info2$something_agrees[which(!gene_info2$something_agrees & (tolower(gene_info2$human_mart) == tolower(gene_info2$ens) | tolower(gene_info2$human_mart) == tolower(gene_info2$ens_minus1) | tolower(gene_info2$human_mart) == tolower(gene_info2$nd_symbol) | tolower(gene_info2$human_mart) == tolower(gene_info2$nd_symbol_minus1) ) )] = T
gene_info2$something_agrees[which(!gene_info2$something_agrees & (tolower(gene_info2$human_pat) == tolower(gene_info2$ens) | tolower(gene_info2$human_pat) == tolower(gene_info2$ens_minus1) | tolower(gene_info2$human_pat) == tolower(gene_info2$nd_symbol) | tolower(gene_info2$human_pat) == tolower(gene_info2$nd_symbol_minus1) ) )] = T
gene_info2$something_agrees[which(!gene_info2$something_agrees & gene_info2$mz_human_desc_match)] = T

gene_info2$label = gene_info2$nd_symbol
gene_info2$label[which( startsWith(gene_info2$label, "LOC") & !startsWith(gene_info2$mzebra, "LOC") )] = gene_info2$mzebra[which( startsWith(gene_info2$label, "LOC") & !startsWith(gene_info2$mzebra, "LOC") )]
gene_info2$label[which( startsWith(gene_info2$label, "zgc:") | startsWith(gene_info2$label, "si:") )] = gene_info2$mzebra[which( startsWith(gene_info2$label, "zgc:") | startsWith(gene_info2$label, "si:") )]
gene_info2$label[which( startsWith(gene_info2$label, "LOC") & gene_info2$something_agrees)] = tolower(gene_info2$human[which( startsWith(gene_info2$label, "LOC") & gene_info2$something_agrees)])
gene_info2 = gene_info2[, c("mzebra", "label", "nd_symbol", "ens", "ens_mart", "ens_me", "human", "human_mart", "human_pat", "mzebra_description", "nd_description", "human_description")]
write.table(gene_info2, "~/research/all_research/gene_info_2.txt", sep = "\t")

# test = gene_info2[which(gene_info2$something_agrees & startsWith(gene_info2$label, "LOC")),]

#*******************************************************************************
# Sample Stats =================================================================
#*******************************************************************************
num_spots_fish = as.data.frame(table(all_merge$fish))
num_spots_fish[,1] = toupper(num_spots_fish[,1])
fish_stats = xlsx::read.xlsx2(paste0(data_dir, "st sample info (official).xlsx"), sheetIndex = 1, startRow = 2)
fish_stats$num_spots = num_spots_fish[match(fish_stats$Fish.ID, num_spots_fish[,1]),2]
fish_stats_melt = reshape2::melt(fish_stats, id.var = "Fish.ID")
fish_stats_melt = fish_stats_melt[which(fish_stats_melt$variable %in% c("BM..g.", "Std..Length..cm.", "Whole.Brain..g.", "Tel..g.", "num_spots")),]
fish_stats_melt$value = as.numeric(as.vector(fish_stats_melt$value))
fish_stats_melt$my.facet = "Whole Fish"
fish_stats_melt$my.facet[which(fish_stats_melt$variable %in% c("Whole.Brain..g.", "Tel..g."))] = "Brain"
fish_stats_melt$variable = plyr::revalue(fish_stats_melt$variable, c("BM..g." = "Body Mass (g)", "Std..Length..cm." = "Std. Length (cm)", "Whole.Brain..g." = "Whole Brain (g)", "Tel..g." = "Tel (g)", "num_spots" = "# of Spots"))
ggplot(fish_stats_melt, aes(x = Fish.ID, y = value, color = variable, group = variable)) + geom_point(size = 4) + geom_smooth() + theme_classic() + xlab("Fish") + ylab("Value") + facet_wrap(~my.facet, scales = "free")
ggplot(fish_stats_melt, aes(x = Fish.ID, y = value, color = variable, group = variable)) + geom_point(size = 4) + geom_smooth() + theme_classic() + xlab("Fish") + ylab("Value") + facet_wrap(~variable, ncol = 1, scales = "free")

#*******************************************************************************
# BB Figure ====================================================================
#*******************************************************************************

# Prep BB Data
bb = readRDS(paste0(brain_dir, "data/bb_sct_070522.rds"))
rgc_sub = readRDS("~/research/brain/data/rgc_subclusters_reclustered_q_c_nb_scores.rds")
bb_convert15 = data.frame(old = 0:14, new = c("8_Glut", "9_Glut", "4_GABA", "15_GABA/Glut", "1_RGC/MG", "10_Glut", "5_GABA", "11_Glut", "6_GABA", "2_OPC/Oligo", "12_Glut", "13_Glut", "14_Glut", "3_Peri", "7_GABA"))
bb_convert53 = data.frame(old = 0:52, new = c("4.1_GABA", "10.1_Glut", "15.1_GABA/Glut", "9.1_Glut", "8.1_Glut", "1.1_RGC", "6_GABA", "5.1_GABA", "9.2_Glut", "8.2_Glut", "15.2_GABA", "11.1_Glut", "8.3_Glut", "8.4_Glut", "9.3_Glut", "4.2_GABA", "8.5_Glut", "5.2_GABA", "8.6_Glut", "8.7_Glut", "1.2_RGC", "4.3_GABA", "4.4_GABA", "9.4_Glut", "9.5_Glut", "8.8_Glut", "9.6_Glut", "4.5_GABA", "12_Glut", "8.9_Glut", "10.2_Glut", "2.1_OPC", "15.3_GABA", "11.2_Glut", "15.4_GABA", "4.6_GABA", "9.7_Glut", "13_Glut", "14_Glut", "4.7_GABA", "11.3_Glut", "9.8_Glut", "8-9_Glut", "15.5_GABA/Glut", "4.8_GABA", "1.3_MG", "2.2_Oligo", "15.6_Glut", "8.10_Glut", "8.11_Glut", "3_Peri", "15.7_Glut", "7_GABA"))
bb$names15 = bb_convert15$new[match(bb$seuratclusters15, bb_convert15$old)]
bb$names53 = bb_convert53$new[match(bb$seuratclusters53, bb_convert53$old)]
bb$rgc_sub = bb$names53
bb$rgc_sub[colnames(rgc_sub)] = paste0("RGC", rgc_sub$seurat_clusters)

# Function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Integrate Separately
spo = qs::qread(paste0(data_dir, "st_obj_list_070822.qs"))
this.pal = c("#ed2828", "#3F4788", "#FDE725")
my.thresh = 40
spo_hi = list()
for (s in names(spo)) {
  message(s)
  message("Integrating with BB Secondary with RGC")
  names(spo[[s]]@images) = s
  anchors = FindTransferAnchors(reference = bb, query = spo[[s]], normalization.method = "SCT", npcs = 50, verbose = F)
  predictionsRGC = TransferData(anchorset = anchors, refdata = bb$rgc_sub, prediction.assay = TRUE, weight.reduction = spo[[s]][["pca"]], dims = 1:50, verbose = F)
  spo[[s]][["predictionsRGC"]] = predictionsRGC
  
  message("Integrating with BB Primary")
  anchors = FindTransferAnchors(reference = bb, query = spo[[s]], normalization.method = "SCT", npcs = 50, verbose = F)
  predictions15 = TransferData(anchorset = anchors, refdata = bb$names15, prediction.assay = TRUE, weight.reduction = spo[[s]][["pca"]], dims = 1:50, verbose = F)
  spo[[s]][["predictions15"]] = predictions15

  message("Making High Res Image")
  spo_hi[[s]] = spo[[s]]
  real.slice = spo_hi[[s]]@images[[s]]
  id = paste0(spo_hi[[s]]$fnum[which(spo_hi[[s]]$sample == s)][1], "_", toupper(spo_hi[[s]]$area[which(spo_hi[[s]]$sample == s)][1]), "1")
  print(paste0(s, ": ", id))
  hires = Read10X_Image(paste0("~/Downloads/sp_data/spatial_modified/", id, "/"), image.name = "tissue_hires_image.png")
  rownames(hires@coordinates) = paste0(s, "_", rownames(hires@coordinates))
  hires@coordinates = hires@coordinates[rownames(real.slice@coordinates),]
  hires@scale.factors$lowres = hires@scale.factors$hires
  hires@assay = real.slice@assay
  hires@key = real.slice@key
  spo_hi[[s]]@images[[s]]= hires
  
  df84_81_2 = data.frame(class1 = as.vector(range01(spo[[s]]@assays$predictionsRGC["RGC1",]))*100, class2 = as.vector(range01(spo[[s]]@assays$predictionsRGC["RGC2",]))*100, class9 = as.vector(range01(spo[[s]]@assays$predictions15["9-Glut",]))*100, class81 = as.vector(range01(spo[[s]]@assays$predictionsRGC["8.1-Glut",]))*100, class84 = as.vector(range01(spo[[s]]@assays$predictionsRGC["8.4-Glut",]))*100)
  df84_81_2$dif84_81 = df84_81_2$class84 - df84_81_2$class81
  if (!all(is.na(df84_81_2$dif84_81))) {
    df84_81_2$dif84_81_col = viridis(100)[cut(df84_81_2$dif84_81,100)]
    df84_81_2$dif84_81_col[which(df84_81_2$dif84_81 == 0)] = "#FDE72500"
    df84_81_2$dif84_81_col[which(df84_81_2$class2 >= my.thresh)] = "#FF9E24"
    df84_81_2$dif84_81_col[which(df84_81_2$class9 >= my.thresh)] = "#5b35bd"
    spo_hi[[s]]$class2_81_84 = df84_81_2$dif84_81_col
    # mySingleSFP(spo_hi[[s]], feature = "class2_81_84", assay = "Spatial", slot = "data", col.ident = T, my.pt.size = 3.5)
  } else {
    message("Setting all to 0.")
    spo_hi[[s]]$class2_81_84 = 0
  }
}

# for (s in names(spo_hi)) {
#   df84_81_2 = data.frame(class1 = as.vector(range01(spo[[s]]@assays$predictionsRGC["RGC1",]))*100, class2 = as.vector(range01(spo[[s]]@assays$predictionsRGC["RGC2",]))*100, class81 = as.vector(range01(spo[[s]]@assays$predictionsRGC["8.1-Glut",]))*100, class84 = as.vector(range01(spo[[s]]@assays$predictionsRGC["8.4-Glut",]))*100)
#   df84_81_2$dif2_1 = df84_81_2$class2 - df84_81_2$class1
#   df84_81_2$dif2_1_col = magma(100)[cut(df84_81_2$dif2_1,100)]
#   df84_81_2$dif84_81 = df84_81_2$class84 - df84_81_2$class81
#   if (!all(is.na(df84_81_2$dif84_81))) {
#     df84_81_2$dif84_81_col = viridis(100)[cut(df84_81_2$dif84_81,100)]
#     df84_81_2$dif84_81_col[which(df84_81_2$dif84_81 == 0)] = "#FDE72500"
#     df84_81_2$dif84_81_col[which(df84_81_2$class2 >= my.thresh)] = "#FF9E24"
#     df84_81_2$scale_col = df84_81_2$dif84_81_col
#     df84_81_2$scale_col[which(df84_81_2$class2 >= my.thresh | df84_81_2$class1 >= my.thresh)] = df84_81_2$dif2_1_col[which(df84_81_2$class2 >= my.thresh | df84_81_2$class1 >= my.thresh)]
#     spo_hi[[s]]$class2_81_84 = df84_81_2$dif84_81_col
#     # mySingleSFP(spo_hi[[s]], feature = "class2_81_84", assay = "Spatial", slot = "data", col.ident = T, my.pt.size = 3.5)
#   } else {
#     message("Setting all to 0.")
#     spo_hi[[s]]$class2_81_84 = 0
#   }
# }

# dif84_81_abs_max = max(abs(df84_81_2$dif84_81))
# my_thresh_raw = 0.25
# for (s in names(spo_hi)) {
#   df84_81_2 = data.frame(class1 = as.vector(range01(spo[[s]]@assays$predictionsRGC["RGC1",]))*100, class2 = as.vector(spo[[s]]@assays$predictionsRGC["RGC2",]), class81 = as.vector(range01(spo[[s]]@assays$predictionsRGC["8.1-Glut",]))*100, class84 = as.vector(range01(spo[[s]]@assays$predictionsRGC["8.4-Glut",]))*100)
#   df84_81_2$dif2_1 = df84_81_2$class2 - df84_81_2$class1
#   df84_81_2$dif2_1_col = magma(100)[cut(df84_81_2$dif2_1,100)]
#   df84_81_2$dif84_81 = df84_81_2$class84 - df84_81_2$class81
#   if (!all(is.na(df84_81_2$dif84_81))) {
#     # df84_81_2$dif84_81_col = viridis(100)[cut(df84_81_2$dif84_81,100)]
#     df84_81_2$dif84_81_col = "#FDE72500"
#     df84_81_2$dif84_81_col[which(df84_81_2$class2 >= my_thresh_raw)] = "#FF9E24"
#     df84_81_2$dif84_81_col[which(df84_81_2$class84 >= 20 & df84_81_2$class81 <= 20)] = "#FDE725FF"
#     df84_81_2$dif84_81_col[which(df84_81_2$class81 >= 20 & df84_81_2$class84 <= 20)] = "#440154FF"
#     spo_hi[[s]]$class2_81_84 = df84_81_2$dif84_81_col
#   } else {
#     message("Setting all to 0.")
#     spo_hi[[s]]$class2_81_84 = 0
#   }
# }
my.thresh2 = 30
for (s in names(spo)) {
  df84_81_2 = data.frame(class1 = as.vector(range01(spo[[s]]@assays$predictionsRGC["RGC1",]))*100, class2 = as.vector(range01(spo[[s]]@assays$predictionsRGC["RGC2",]))*100, class9 = as.vector(range01(spo[[s]]@assays$predictions15["9-Glut",]))*100, class81 = as.vector(range01(spo[[s]]@assays$predictionsRGC["8.1-Glut",]))*100, class84 = as.vector(range01(spo[[s]]@assays$predictionsRGC["8.4-Glut",]))*100)
  df84_81_2$dif84_81 = df84_81_2$class84 - df84_81_2$class81
  if (!all(is.na(df84_81_2$dif84_81))) {
    df84_81_2$dif84_81_col = viridis(100)[cut(df84_81_2$dif84_81,100)]
    df84_81_2$dif84_81_col[which(df84_81_2$dif84_81 == 0)] = "#FDE72500"
    df84_81_2$dif84_81_col[which(df84_81_2$class9 >= my.thresh)] = "#5b35bd"
    df84_81_2$dif84_81_col[which(df84_81_2$class2 >= my.thresh2)] = "#FF9E24"
    spo_hi[[s]]$class2_81_84 = df84_81_2$dif84_81_col
    # mySingleSFP(spo_hi[[s]], feature = "class2_81_84", assay = "Spatial", slot = "data", col.ident = T, my.pt.size = 3.5)
  } else {
    message("Setting all to 0.")
    spo_hi[[s]]$class2_81_84 = 0
  }
}


tmp = merge(spo_hi[[names(spo_hi)[1]]], spo_hi[[names(spo_hi)[2]]])
for (s in names(spo)[3:length(spo_hi)]) {
  print(s)
  tmp = merge(tmp, spo_hi[[s]])
}

pdf(paste0(out_dir, paste0("class2_81_84_11_", my.thresh, ".pdf")), width = 12, height = 12, onefile = F)
print(myMultiSFP(tmp, feature = "class2_81_84", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100)), rm.zero = T, col.ident = T, high.res = T ))
dev.off()

pdf(paste0(out_dir, paste0("9_Glut.pdf")), width = 12, height = 12, onefile = F)
print(myMultiSFP(all_merge_hi, assay = "predictions15", feature = "9-Glut", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100)), rm.zero = T, high.res = T ))
dev.off()

# 8.1_Glut, 8.4_Glut, RGC2 Figure
this.pal = c("#ed2828", "#3F4788", "#FDE725")
my.thresh = 40
df84_81_2 = data.frame(class1 = as.vector(range01(all_merge_for_bb@assays$predictionsRGC["RGC1",]))*100, class2 = as.vector(range01(all_merge_for_bb@assays$predictionsRGC["RGC2",]))*100, class81 = as.vector(range01(all_merge_for_bb@assays$predictionsRGC["8.1-Glut",]))*100, class84 = as.vector(range01(all_merge_for_bb@assays$predictionsRGC["8.4-Glut",]))*100)
df84_81_2$dif2_1 = df84_81_2$class2 - df84_81_2$class1
df84_81_2$dif2_1_col = magma(100)[cut(df84_81_2$dif2_1,100)]
df84_81_2$dif84_81 = df84_81_2$class84 - df84_81_2$class81
# df84_81_2$dif84_81_col = viridis(100)[cut(df84_81_2$dif84_81,100)] # scaled across all samples
col_values_by_sample = unlist(lapply(c("c2c", "b2b", "b1c"), function(x) viridis(100)[cut(df84_81_2$dif84_81[which(all_merge_for_bb$sample == x)],100)])) # scaled by sample
names(col_values_by_sample) = unlist(lapply(c("c2c", "b2b", "b1c"), function(x) colnames(all_merge_for_bb)[which(all_merge_for_bb$sample == x)])) # scaled by sample
df84_81_2$dif84_81_col = col_values_by_sample[colnames(all_merge_for_bb)]
df84_81_2$dif84_81_col[which(df84_81_2$dif84_81 == 0)] = "#FDE72500"
df84_81_2$dif84_81_col[which(df84_81_2$class2 >= my.thresh)] = "#FF9E24"
df84_81_2$scale_col = df84_81_2$dif84_81_col
df84_81_2$scale_col[which(df84_81_2$class2 >= my.thresh | df84_81_2$class1 >= my.thresh)] = df84_81_2$dif2_1_col[which(df84_81_2$class2 >= my.thresh | df84_81_2$class1 >= my.thresh)]
all_merge_for_bb_hi$class2_81_84 = df84_81_2$dif84_81_col

pdf(paste0(out_dir, paste0("class2_81_84_6_", my.thresh, "_only3.pdf")), width = 12, height = 12)
print(myMultiSFP(all_merge_for_bb_hi, samples = c("c2c", "b2b", "b1c"), "class2_81_84", pt.size.multiplier = 1.3, pal = colorRampPalette(viridis(100)), rm.zero = T, col.ident = T, high.res = T ))
dev.off()

#*******************************************************************************
# Trash Can ====================================================================
#*******************************************************************************

means = read.csv(paste0(out_dir, "cell2location/bb_secondary/cell2location_spatial_output_means.csv"))
rownames(means) = means$X
means$X = NULL
colnames(means) = str_replace(colnames(means), "meanscell_abundance_w_sf_", "")
means_round = round(means)
all_merge$sum = rowSums(means)

for (i in as.character(0)) {
  all_merge$tmp = means_round[,i]
  # pdf(paste0(out_dir, "cell2location/bb_secondary/results/", i, ".pdf"), width = 8, height = 8, onefile = F)
  Cairo::Cairo(file = paste0(out_dir, "cell2location/bb_secondary/results/", i, ".png"), width = 1800, height = 1800, res = 150)
  print(myMultiSFP(all_merge, feature = "tmp", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100))))
  dev.off()
  
  pct_sample = unlist(lapply(levels(all_merge$sample), function(x) range01(all_merge$tmp[which(all_merge$sample ==x)])*100))
  names(pct_sample) = unlist(lapply(levels(all_merge$sample), function(x) colnames(all_merge)[which(all_merge$sample ==x)]))
  all_merge$tmp2 = 0
  all_merge$tmp2[names(pct_sample)] = pct_sample
  # pdf(paste0(out_dir, "cell2location/bb_secondary/results/", i, "_max.pdf"), width = 8, height = 8, onefile = F)
  Cairo::Cairo(file = paste0(out_dir, "cell2location/bb_secondary/results/", i, "_max.png"), width = 1800, height = 1800, res = 150)
  print(myMultiSFP(all_merge, feature = "tmp2", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100))))
  dev.off()
  
  all_merge$tmp = (means_round[,i] / all_merge$sum) * 100
  # pdf(paste0(out_dir, "cell2location/bb_secondary/results/", i, "_pct.pdf"), width = 8, height = 8, onefile = F)
  Cairo::Cairo(file = paste0(out_dir, "cell2location/bb_secondary/results/", i, "_pct.png"), width = 1800, height = 1800, res = 150)
  print(myMultiSFP(all_merge, feature = "tmp", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100))))
  dev.off()
  
  pct_sample = unlist(lapply(levels(all_merge$sample), function(x) range01(all_merge$tmp[which(all_merge$sample ==x)])*100))
  names(pct_sample) = unlist(lapply(levels(all_merge$sample), function(x) colnames(all_merge)[which(all_merge$sample ==x)]))
  all_merge$tmp2 = 0
  all_merge$tmp2[names(pct_sample)] = pct_sample
  # pdf(paste0(out_dir, "cell2location/bb_secondary/results/", i, "_pct_max.pdf"), width = 8, height = 8, onefile = F)
  Cairo::Cairo(file = paste0(out_dir, "cell2location/bb_secondary/results/", i, "_pct_max.png"), width = 1800, height = 1800, res = 150)
  print(myMultiSFP(all_merge, feature = "tmp2", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100))))
  dev.off()
}

# Label Tissue Halves
angle.df = as.data.frame(c("c2a" = 195, "c2b" = 208, "c2c" = -115, "c2d" = 210,
                           "b2a" =  -90, "b2b" =  -95, "b2c" =  -100, "b2d" =  -85,
                           "c1a" =  -120, "c1b" =  -95, "c1c" =   -118, "c1d" =  -90,
                           "b1a" =   0, "b1b" =   0, "b1c" = 120, "b1d" =   0))
colnames(angle.df) = "angle"

s = "b1c"
spo[[s]]$half = "neither"

coords = GetTissueCoordinates(object = spo[[s]])
coords$spot = rownames(coords)
colnames(coords)[1:2] = c("UMAP_1", "UMAP_2")
coords$UMAP_1 = - coords$UMAP_1

M = as.matrix(coords[, c("UMAP_1", "UMAP_2")])
my.degree = angle.df[s, "angle"]
my.radians = my.degree * pi/180
alpha = my.radians
#rotation matrix
rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
#shift, rotate, shift back
M2 <- t(rotm %*% (
  t(M)-c(M[1,1],M[1,2])
)+c(M[1,1],M[1,2]))
coords$UMAP_1 = M2[,1]
coords$UMAP_2 = M2[,2]
p = ggplot(coords, aes(x=UMAP_2, y=UMAP_1)) + geom_point(size = 4, stroke = 0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + NoLegend()
a = CellSelector(p)

spo[[s]]$half[a] = "left"
spo[[s]]$half[which(spo[[s]]$half == "neither")] = "right"
Idents(spo[[s]]) = spo[[s]]$half
SpatialDimPlot(spo[[s]])

half_df = data.frame(cell = unlist(lapply(names(spo), function(x) colnames(spo[[x]]))), half = unlist(lapply(names(spo), function(x) spo[[x]]$half)))
write.csv(half_df, "~/research/st/data/st_halves_labelled.csv")
all_merge$half = half_df$half[match(colnames(all_merge), half_df$cell)]
all_merge$half_binary = as.numeric(as.vector(plyr::revalue(all_merge$half, c("left" = "0", "right" = "1", "neither" = "3"))))

# GGFORCE
library(ggforce)
mySingleSFP(c2c, "cluster", assay = "SCT", slot = "data", my.pt.size = 4, discrete = F)
coords = GetTissueCoordinates(object = c2c)
coords$value = c2c$cluster
coords$spot = rownames(coords)
ggplot(coords[which(coords$value != 29),], aes(x=imagecol, y=-imagerow, color = value)) + geom_point(size = 4, stroke = 0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + NoLegend() + ggforce::geom_mark_hull(aes(fill = value), expand = unit(3, "mm"))
ggplot(coords[which(coords$value != 29),], aes(x=imagecol, y=-imagerow, color = value)) + geom_point(size = 4, stroke = 0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + NoLegend() + ggforce::geom_mark_hull(aes(fill = value, group = spot), expand = unit(3, "mm"))
ggplot(coords[which(coords$value != 29),], aes(x=imagecol, y=-imagerow, color = value)) + geom_point(size = 4, stroke = 0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + NoLegend() + ggforce::geom_mark_hull(aes(fill = value), expand = unit(3, "mm")) + ggforce::geom_voronoi_segment(aes(group = 1))
ggplot(coords[which(coords$value != 29),], aes(x=imagecol, y=-imagerow, color = value)) + geom_point(size = 4, stroke = 0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + NoLegend() + ggforce::geom_voronoi_tile(aes(fill = value, group = 1), alpha = 0.2) + ggforce::geom_voronoi_segment(aes(group = 1)) + geom_point(color = "black")

ggplot(coords[which(coords$value != 29),], aes(x=imagecol, y=-imagerow, color = value)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + NoLegend() + ggforce::geom_voronoi_segment(aes(group = 1))


# Find the composition of clusters by samples and samples by clusters
sample.df  = as.data.frame(table(all_merge$sample))
cluster.df = as.data.frame(table(all_merge$cluster))
sample.cluster.df = as.data.frame(table(paste0(all_merge$sample, "_", all_merge$cluster)))
colnames(sample.cluster.df)[1] = "sample_cluster"
sample.cluster.df[, c("sample", "cluster")] = reshape2::colsplit(sample.cluster.df$sample_cluster, "_", c('1', '2'))
sample.cluster.df$sample_n  = sample.df[match(sample.cluster.df$sample, sample.df[,1]), 2]
sample.cluster.df$cluster_n = cluster.df[match(sample.cluster.df$sample, sample.df[,1]), 2]
sample.cluster.df$pct.of.cluster = sample.cluster.df$Freq / sample.cluster.df$cluster_n
sample.cluster.df$pct.of.sample  = sample.cluster.df$Freq / sample.cluster.df$sample_n
ggplot(sample.cluster.df, aes(x = sample, y = pct.of.cluster, fill = sample))  + geom_bar(stat = "identity") + xlab("Sample")  + ylab("% of Cluster") + scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + facet_wrap(~ cluster)
sample.cluster.df$cluster = factor(sample.cluster.df$cluster, levels = sort(unique(sample.cluster.df$cluster)))
ggplot(sample.cluster.df, aes(x = cluster, y = pct.of.sample, fill = cluster)) + geom_bar(stat = "identity") + xlab("Cluster") + ylab("% of Sample")  + scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + facet_wrap(~ sample)

markers8184 = read.csv('~/Downloads/rgc12_8_1_and_8_4_Glut_markers_071622_minpct20.csv')
markers8184$hgnc = gene_info$human[match(markers8184$mzebra, gene_info$mzebra)]
markers8184$mouse = str_to_title(markers8184$hgnc)
allSamplesSFP(all_merge, "LOC101477092", pal = colorRampPalette(colors = viridis(100)))


Idents(all_merge) = all_merge$sample
pw.dif.df = data.frame()
for (i in 1:length(names(spo))) {
  s1 = names(spo)[i]
  print(s1)
  for ( j in (i+1):length(names(spo)) ) {
    s2 = names(spo)[j]
    this.deg = FindMarkers(all_merge, ident.1 = s1, ident.2 = s2, only.pos = T, min.pct = 1e-100, logfc.threshold = 0)
    this.deg$gene = rownames(this.deg)
    this.deg$s1 = s1
    this.deg$s2 = s2
    pw.dif.df = rbind(pw.dif.df, this.deg)
  }
}
pw.dif.df$s1_s2 = paste0(pw.dif.df$s1, "_", pw.dif.df$s2)
write.csv(pw.dif.df, "pairwise_sample_degs_loose_071122.csv")

pw.dif.df.p = as.data.frame(table( pw.dif.df$s1_s2[which(pw.dif.df$p_val_adj < 0.05 & pw.dif.df$pct.1 > 0.05 & pw.dif.df$avg_log2FC > 0.25)] ))
pw.dif.df.p[, c("s1", "s2")] = reshape2::colsplit(pw.dif.df.p$Var1, pattern = "_", c("1", "2"))
pw.dif.df.p$s1 = factor(pw.dif.df.p$s1, levels = names(spo))
pw.dif.df.p$s2 = factor(pw.dif.df.p$s2, levels = names(spo))
ggplot(pw.dif.df.p, aes(x = s2, y = s1, fill = Freq)) + geom_tile() + scale_fill_gradientn(colours = viridis(100)) + xlab("Sample 1") + ylab("Sample 2") + theme_classic() +  ggtitle("Number of DEGs b/w Samples") + coord_fixed()

pw.dif.df.p = pw.dif.df.p[which( (startsWith(as.vector(pw.dif.df.p$s1), "c2") | startsWith(as.vector(pw.dif.df.p$s1), "b2")) & (startsWith(as.vector(pw.dif.df.p$s2), "c2") | startsWith(as.vector(pw.dif.df.p$s2), "b2"))  ),]
ggplot(pw.dif.df.p, aes(x = s2, y = s1, fill = Freq)) + geom_tile() + scale_fill_gradientn(colours = viridis(100)) + xlab("Sample 1") + ylab("Sample 2") + theme_classic() +  ggtitle("Number of DEGs b/w Samples") + coord_fixed()

pw.dif.df.p = pw.dif.df.p[which( !(startsWith(as.vector(pw.dif.df.p$s1), "c2d") | startsWith(as.vector(pw.dif.df.p$s1), "c1c")) & !(startsWith(as.vector(pw.dif.df.p$s2), "c2d") | startsWith(as.vector(pw.dif.df.p$s2), "c1c")) ),]
ggplot(pw.dif.df.p, aes(x = s2, y = s1, fill = Freq)) + geom_tile() + scale_fill_gradientn(colours = viridis(100)) + xlab("Sample 1") + ylab("Sample 2") + theme_classic() +  ggtitle("Number of DEGs b/w Samples") + coord_fixed()

# Get peripheral points
spots.rgc = c()
p_list = list()
start.r = 1000
for (s in levels(all_merge$sample)) {
  this.coords = spo[[s]]@images[["slice1"]]@coordinates
  this.coords$spot = rownames(this.coords)
  this.coords$angle = atan2(this.coords$imagerow - mean(this.coords$imagerow), this.coords$imagecol - mean(this.coords$imagecol))  # find angle of points from center
  c1 = c(mean(this.coords$imagerow), mean(this.coords$imagecol))
  this.coords$dist = apply(this.coords[, c("imagerow", "imagecol")],1,function(x,c1) {(sqrt((x[1] - c1[1])^2+(x[2]-c1[2])^2))},c1)
  # this.r = (nrow(this.coords)/500) * start.r
  this.r = (max(this.coords$dist)/1410) * start.r
  this.coords$peripheral = this.coords$dist > this.r
  spots.rgc = c(spots.rgc, this.coords$spot[which(this.coords$peripheral)])
  # ggplot(this.coords, aes(imagecol, -imagerow, color = angle, size = dist)) + geom_point()
  p_list[[s]] = ggplot(this.coords, aes(imagecol, -imagerow, color = peripheral)) + geom_point(size = 1.5) + theme_void() + ggtitle(s)
}
wrap_plots(p_list, ncol = 4)

all_merge_rgc = subset(all_merge, cells = spots.rgc)
allSamplesSFP(all_merge_rgc, "rgc1.top")


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