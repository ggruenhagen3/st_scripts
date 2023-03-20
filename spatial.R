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
gene_info = read.csv(paste0(main_path, "/all_research/gene_info_3.csv"))
all_merge = readRDS(paste0(data_dir, "st_c2b2_hi_012723.rds"))
convert15 = read.csv("~/research/brain/results/convert15.csv")
convert53 = read.csv("~/research/brain/results/convert53.csv")
# all_merge = qs::qread(paste0(data_dir, "st_070822.qs"))
all_merge_hi = qs::qread(paste0(data_dir, "all_merge_hi.qs"))
spo = qs::qread(paste0(data_dir, "st_obj_list_070822.qs"))

#*******************************************************************************
# Cell2location ================================================================
#*******************************************************************************

#-----------------------------------#
# Cell2location Secondary - Rounded #
#-----------------------------------#
# Load cell2location results
c2l_mean = read.csv(paste0(out_dir, "cell2location_c2b2_spatial53_output_means.csv")); rownames(c2l_mean) = c2l_mean$X; c2l_mean$X = NULL; colnames(c2l_mean) = as.character(0:52)
# convert the cluster numbers into cell type names
colnames(c2l_mean) = convert53$new[match(colnames(c2l_mean), convert53$old)]; c2l_mean = c2l_mean[, convert53$new]
# Prevent spots from having 0 cells after rounding
zero.cell.st = which(rowSums(round(c2l_mean)) == 0); for (i in zero.cell.st) { c2l_mean[i,which.max(c2l_mean[i,])] = 1 }
c2l_mean = round(c2l_mean)
c2l_mean_pct = c2l_mean / rowSums(c2l_mean)

# Get the celltype with the most cells per spot
st.celltype.char = unlist(lapply(1:nrow(c2l_mean), function(x) colnames(c2l_mean)[which.max(c2l_mean[x,])]))

c2l.to.cluster = expand.grid(cluster = levels(all_merge$cluster), celltype = colnames(c2l_mean))
sum.cells.per.cluster = unlist(lapply(levels(all_merge$cluster), function(x) sum(c2l_mean[colnames(all_merge)[which(all_merge$cluster == x)],]) ))
names(sum.cells.per.cluster) = levels(all_merge$cluster)
c2l.to.cluster$cluster.sum.num.cells = sum.cells.per.cluster[c2l.to.cluster$cluster]
c2l.to.cluster$num = unlist(lapply(1:nrow(c2l.to.cluster), function(x) sum(c2l_mean[colnames(all_merge)[which(all_merge$cluster == c2l.to.cluster$cluster[x])], c2l.to.cluster$celltype[x]]) ))
c2l.to.cluster$pct.of.cluster = c2l.to.cluster$num / c2l.to.cluster$cluster.sum.num.cells
pdf(paste0( "cluster_to_c2l.pdf"), width = 8, height = 5)
ggplot(c2l.to.cluster, aes(x = celltype, y = cluster, fill = pct.of.cluster)) + geom_raster() + coord_fixed() + scale_x_discrete(expand = c(0,0), name = "") + scale_y_discrete(expand = c(0,0), name = "") + theme_classic() + scale_fill_viridis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank())
dev.off()

all_merge$celltype = factor(st.celltype.char[match(colnames(all_merge), rownames(c2l_mean))], levels = convert53$new)
pdf(paste0( "b2_painted_by_best_celltype.pdf"), width = 16, height = 2)
myB2SFP(all_merge, "celltype", pal = hue_pal()(53), showLegend = T)
dev.off()

all_merge$num.cells = rowSums(c2l_mean)
all_merge$num.cells[which(all_merge$num.cells > 30)] = 30
pdf(paste0( "c2b2_num_cells.pdf"), width = 3, height = 3)
print(FeaturePlot(all_merge, "num.cells2", order = T) + scale_color_viridis() + theme_void() + NoLegend() + ggtitle("") + coord_fixed())
dev.off()
pdf(paste0( "b2_painted_by_best_celltype.pdf"), width = 3, height = 3)
Idents(all_merge) = all_merge$celltype
print(DimPlot(all_merge) + theme_void() + NoLegend() + ggtitle("") + coord_fixed())
dev.off()
pdf(paste0(out_dir, "b2_mini_cells.pdf"), width = 25, height = 3, onefile = F)
print(myB2SFP( all_merge, "ct", stsc.list = list(stsc.mat, stsc.meta), pal = scales::hue_pal()(length(levels(all_merge$ct))), showLegend = T, pt.size.multiplier = 1.1 ))
dev.off()
pdf(paste0( "c2b2_mean_num_cell_per_cluster.pdf"), width = 2.5, height = 5)
ggplot(st.num.cells.per.cluster, aes(x = seurat_clusters, y = num.cells, fill = seurat_clusters)) + geom_bar(stat = 'identity') + theme_classic() + scale_y_continuous(expand = c(0,0), name = "Mean # of Cells") + xlab("") + NoLegend() + coord_flip()
dev.off()
pdf(paste0( "c2b2_mean_num_cell_per_ct.pdf"), width = 3, height = 6)
ggplot(st.num.cells.per.ct, aes(x = colSums.c2l_mean., y = ct, fill = ct)) + geom_bar(stat = "identity") + theme_classic() + theme(panel.grid.major.x = element_line(color = "gray75", linetype = "dashed")) + scale_x_continuous(expand = c(0,0), name = "# of Cells") + NoLegend() + ylab("")
dev.off()

all_merge$num.cells.max = unlist(lapply(1:nrow(c2l_mean), function(x) max(c2l_mean[x,])))[match(colnames(all_merge), rownames(c2l_mean))]
pct.of.spot.bin = floor((all_merge$num.cells.max/all_merge$num.cells) / 0.2)
pct.of.spot.bin[which(pct.of.spot.bin %% 0.2 == 0)] = pct.of.spot.bin[which(pct.of.spot.bin %% 0.2 == 0)]-1
# pct.of.spot.bin[which(pct.of.spot.bin == 5)] = 4
df = as.data.frame(table(pct.of.spot.bin))
df$pct = df$Freq / sum(df$Freq) * 100
pdf(paste0(out_dir, "spot_comp.pdf"), width = 6, height = 2)
print(ggplot(df, aes(x = 1, y = pct, fill = pct.of.spot.bin)) + geom_bar(stat = "identity") + coord_flip() + scale_x_continuous(expand = c(0.8,0.8), name = "") + scale_y_continuous(expand = c(0,0), name = "") + theme_classic())
dev.off()

all_merge$num.celltypes = unlist(lapply(1:nrow(c2l_mean), function(x) length(which(c2l_mean[x,] != 0)) ))[match(colnames(all_merge), rownames(c2l_mean))]
num.celltypes.bin = floor(all_merge$num.celltypes / 3)
num.celltypes.bin[which(all_merge$num.celltypes  %% 3 == 0)] = num.celltypes.bin[which(all_merge$num.celltypes  %% 3 == 0)]-1
df = as.data.frame(table(num.celltypes.bin))
df$pct = df$Freq / sum(df$Freq) * 100
pdf(paste0(out_dir, "spot_comp2.pdf"), width = 6, height = 2)
print(ggplot(df, aes(x = 1, y = pct, fill = pct.of.spot.bin)) + geom_bar(stat = "identity") + coord_flip() + scale_x_continuous(expand = c(0.8,0.8), name = "") + scale_y_continuous(expand = c(0,0), name = "") + theme_classic())
dev.off()

# Jaccard Index of celltypes that overlap in spots
jcr.df = data.frame()
cor.df = data.frame()
for (i in 1:(ncol(c2l_mean)-1)) {
  for (j in (i+1):ncol(c2l_mean)) {
    this.int = length(which(c2l_mean[,i] > 0 & c2l_mean[,j] > 0))
    if ( this.int > 0 ) { this.int = sum(c2l_mean[which(c2l_mean[,i] > 0 & c2l_mean[,j] > 0), c(i,j)]) }
    this.union = sum(c2l_mean[,i]) + sum(c2l_mean[,j])
    jcr.df = rbind(jcr.df, data.frame(clust1 = colnames(c2l_mean)[i], clust2 = colnames(c2l_mean)[j], int = this.int, union = this.union, jcr = this.int/this.union))
    cor.df = rbind(cor.df, data.frame(clust1 = colnames(c2l_mean)[i], clust2 = colnames(c2l_mean)[j], cor = cor(c2l_mean[,i], c2l_mean[,j])))
  }
}

# Do celltypes show a bias on the anterior-posterior axis?
# TODO add a z-score
c2l_mean_tmp = c2l_mean[colnames(all_merge),]
c2l_mean_tmp = c2l_mean_pct[colnames(all_merge),]
c2l_mean_tmp$spot = rownames(c2l_mean_tmp)
c2l_mean_tmp$num.cells = all_merge$num.cells
c2l_mean_tmp$sh = all_merge$sh
c2l_mean_tmp_melt = reshape2::melt(c2l_mean_tmp, id.vars = c("spot", "num.cells", "sh"))
ct.sh.df = c2l_mean_tmp_melt %>%  group_by(variable, sh) %>% summarise(num = sum(value))
sh.df = c2l_mean_tmp_melt %>%  group_by(sh) %>% summarise(num = sum(value))
ct.sh.df$sh.num = sh.df$num[match(ct.sh.df$sh, sh.df$sh)]
ct.sh.df$pct = ct.sh.df$num / ct.sh.df$sh.num
ct.sh.df$sh.ordered = match(ct.sh.df$sh, c("c2dr", "c2dl", "c2cr", "c2cl", "c2br", "c2bl", "b2ar", "b2al", "b2br", "b2cr", "b2dr", "b2bl", "b2cl", "b2dl"))
ct.sh.df$variable = as.character(as.vector(ct.sh.df$variable))
# ct.sh.df = as.data.frame(ct.sh.df)
ct.sh.df = ct.sh.df %>% group_by(variable) %>% mutate(z.score=scale(pct))
# ct.sh.df$z.score.pos = ct.sh.df$z.score > 0
ct.sh.df$bias = F
for (i in unique(ct.sh.df$variable)) { 
  pos.slide.nums = ct.sh.df$sh.ordered[which(ct.sh.df$variable == i & ct.sh.df$z.score > 0)]
  diff.pos.slide.nums = diff(sort(pos.slide.nums))
  if (all(diff.pos.slide.nums == 1)) { ct.sh.df$bias[which(ct.sh.df$variable == i)] = T }
}

lm.res = as.data.frame(ct.sh.df %>% group_by(variable) %>% do(p = summary(lm(pct ~ sh.ordered, data = .))$coefficients[2,4], co = summary(lm(pct ~ sh.ordered, data = .))$coefficients[2,1]))
ggplot(ct.sh.df[which(ct.sh.df$variable %in% c("8.5_Glut")),], aes(x = sh.ordered, y = num)) + geom_point() + stat_smooth(method = "lm", col = "red") + theme_bw()
lm.res = data.frame(celltype = lm.res[,1], p = unlist(lm.res[,2]), co = unlist(lm.res[,3]))
lm.res$bh  = p.adjust(lm.res$p, method = "BH")
lm.res$bon = p.adjust(lm.res$p, method = "bonferroni")

lm.res$celltype[which(lm.res$p < 0.05)]
ct.sh.df.tmp = ct.sh.df[which(ct.sh.df$variable %in% lm.res$celltype[which(lm.res$p < 0.05)]),]
ct.sh.df.tmp$sign = sign(lm.res$co[match(ct.sh.df.tmp$variable, lm.res$celltype)])
ggplot(ct.sh.df.tmp, aes(x = sh.ordered, y = num, color = variable)) + geom_point() + stat_smooth(method = "lm", se = F) + theme_bw() + facet_wrap(~ sign)

ct.sh.df.mat = reshape2::acast(ct.sh.df, variable ~ sh.ordered, value.var = "pct")
ct.sh.df.mat.scale = scale(ct.sh.df.mat)
ct.sh.df.mat.scale.melt = reshape2::melt(ct.sh.df.mat.scale)
ct.sh.df.mat.scale.melt.hclust  = hclust(dist(ct.sh.df.mat.scale), method = "complete")
ct.sh.df.mat.scale.melt.order = ct.sh.df.mat.scale.melt.hclust$labels[ct.sh.df.mat.scale.melt.hclust$order]
ct.sh.df.mat.scale.melt$Var1 = factor(ct.sh.df.mat.scale.melt$Var1, levels = ct.sh.df.mat.scale.melt.order)
ggplot(ct.sh.df.mat.scale.melt, aes(x = Var2, y = Var1, fill = value)) + geom_raster() + scale_fill_viridis() + coord_fixed() + scale_x_continuous(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))

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
pdf("percent_and_count_scatter.pdf", width = 5, height = 5)
print(ggplot(stat.df, aes(x = max_num, y = max_pct)) + geom_point(alpha = 0.25) + theme_bw() + xlab("Max # Cells from a Cell Type") + ylab("Max % from a Cell Type"))
dev.off()

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
# pdf("c2c_pi_chart.pdf", width = 12, height = 12, onefile = F)
Cairo::Cairo("c2c_pi_chart.png", width = 2400, height = 1600, res = 240)
# print(ggplot() + geom_scatterpie(data = coords, aes(x=imagecol, y=-imagerow, group = spot), cols = colnames(c2l_mean_num2), color=NA, alpha = 0.8) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void()) + coord_fixed() + scale_fill_manual(values = c(hue_pal()(48), 'gray60'), name = "cell type")
print(ggplot() + geom_scatterpie(data = coords, aes(x=imagecol, y=-imagerow, group = spot), cols = colnames(c2l_mean_num2), color=NA, alpha = 0.8) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed())
dev.off()

c2l_mean_num2_round = round(c2l_mean_num2)
coords_multi = data.frame()
for (i in 1:nrow(c2l_mean_num2_round)) {
  if (startsWith(rownames(c2l_mean_num2_round)[i], "c2c_")) {
    this_non_zero_cell_type = colnames(c2l_mean_num2_round)[which(c2l_mean_num2_round[i,] != 0)]
    for (this_ct in this_non_zero_cell_type) {
      coords_multi = rbind(coords_multi, 
                           data.frame(spot = rep(rownames(c2l_mean_num2_round)[i], c2l_mean_num2_round[i, this_ct]), 
                                      imagerow = rep(coords[rownames(c2l_mean_num2_round)[i], "imagerow"], c2l_mean_num2_round[i, this_ct]), 
                                      imagecol = rep(coords[rownames(c2l_mean_num2_round)[i], "imagecol"], c2l_mean_num2_round[i, this_ct]), 
                                      ct = rep(this_ct, c2l_mean_num2_round[i, this_ct])) ) 
    }
  }
}
coords_multi$value = all_merge$cluster[coords_multi$spot]
coords_multi$best_ct = all_merge$bb53_c2l[coords_multi$spot]
coords_multi$ct = factor(coords_multi$ct, levels = convert53$new)

img.grob = GetImage(all_merge_hi, image = "c2c")
img.grob.test = img.grob$raster[min(coords$imagerow):max(coords$imagerow), min(coords$imagecol):max(coords$imagecol)]
img.grob.test.grob = img.grob
img.grob.test.grob$raster = img.grob.test

# pdf("c2c_with_mini_cells_no_circle.pdf", width = 12, height = 12, onefile = F)
# pdf("c2c_with_mini_cells_w_circle.pdf", width = 12, height = 8, onefile = F)
Cairo::CairoPNG("c2c_with_mini_cells_w_circle.png", width = 3000, height = 1950, res = 240)
# print(ggplot(coords_multi, aes(x=imagecol, y=-imagerow)) + annotation_custom(img.grob.test.grob) + ggforce::geom_mark_hull(aes(fill = value, group = spot), expand = unit(4, "mm"), radius = unit(4, "mm"), size = 0, color = NA) + geom_point(size = 0.75, position = position_jitter(width = 2.5, height = 2.5), aes(color = ct)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed() + guides(fill="none"))
# print(ggplot(coords_multi, aes(x=imagecol, y=-imagerow)) + annotation_custom(img.grob.test.grob) + ggforce::geom_mark_hull(aes(fill = value, group = spot), expand = unit(4, "mm"), radius = unit(4, "mm"), size = 0, color = NA, alpha = 0.5) + geom_point(size = 0.75, position = position_jitter(width = 8, height = 8), aes(color = ct)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed() + guides(fill="none"))
# print(ggplot(coords_multi, aes(x=imagecol, y=-imagerow)) + geom_point(size = 0.75, position = position_jitter(width = 8, height = 8), aes(color = ct)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_void() + coord_fixed() + guides(fill="none"))
print(ggplot(coords_multi, aes(x=imagecol, y=-imagerow)) + ggforce::geom_mark_hull(aes(fill = best_ct, group = spot), expand = unit(3.3, "mm"), radius = unit(3.3, "mm"), size = 0, color = NA, alpha = 0.5) + geom_point(size = 0.75, position = position_jitter(width = 7, height = 7), aes(color = ct)) + scale_x_continuous(expand=c(0.1,0.1)) + scale_y_continuous(expand=c(0.1,0.1)) + theme_void() + coord_fixed() + guides(fill="none"))
dev.off()

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

# Correlation w/ manual annotation
library(ggpmisc)
this.num = read.csv("~/Downloads/c2bl_num_cells.csv")
colnames(this.num) = c("Barcode", "num.cells")
this.num$Barcode = paste0("c2b_", this.num$Barcode)
this.num$cell2loc = all_merge$num.cells[match(this.num$Barcode, colnames(all_merge))]
cor(this.num[,2], this.num$cell2loc)
cor(this.num[,2], this.num$cell2loc, method = "spearman")
ggplot(this.num, aes(x = num.cells, y = cell2loc)) + geom_point() + stat_poly_line() + stat_poly_eq()

#*******************************************************************************
# BB Figure of Spatial =========================================================
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
bb$names15_4gaba = bb$names15
bb$names15_4gaba[which(bb$names15 == "4_GABA" & bb@assays$RNA@counts["htr1d",] > 0)] = "4_GABA htr1d"
bb$names15_4gaba[which(bb$names15 == "4_GABA" & bb@assays$RNA@counts["vipr2",] > 0)] = "4_GABA vipr2"
anchors = FindTransferAnchors(reference = bb, query = all_merge, normalization.method = "SCT", npcs = 50)
predictionsRGC = TransferData(anchorset = anchors, refdata = bb$rgc_sub, prediction.assay = TRUE, weight.reduction = all_merge[["pca"]], dims = 1:50)
all_merge[["predictionsRGC"]] = predictionsRGC
predictions15 = TransferData(anchorset = anchors, refdata = bb$names15_4gaba, prediction.assay = TRUE, weight.reduction = all_merge[["pca"]], dims = 1:50)
all_merge[["predictions15"]] = predictions15

# Function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

this.pal = c("#04D9FF", "#3F4788", "#FDE725")
my.thresh = 50
my.thesh2 = 60
# df84_81_2 = data.frame(class1 = as.vector(range01(all_merge@assays$predictionsRGC["RGC1",]))*100, class2 = as.vector(range01(all_merge@assays$predictionsRGC["RGC2",]))*100, class4 = as.vector(range01(all_merge@assays$predictions15["4-GABA",]))*100, class81 = as.vector(range01(all_merge@assays$predictionsRGC["8.1-Glut",]))*100, class84 = as.vector(range01(all_merge@assays$predictionsRGC["8.4-Glut",]))*100, class9 = as.vector(range01(all_merge@assays$predictions15["9-Glut",]))*100)
df84_81_2 = data.frame(class1 = as.vector(range01(all_merge@assays$predictionsRGC["RGC1",]))*100, class2 = as.vector(range01(all_merge@assays$predictionsRGC["RGC2",]))*100, class3 = as.vector(range01(all_merge@assays$predictionsRGC["RGC3",]))*100, class4_htr1d = as.vector(range01(all_merge@assays$predictions15["4-GABA htr1d",]))*100, class4_vipr2 = as.vector(range01(all_merge@assays$predictions15["4-GABA vipr2",]))*100, class81 = as.vector(range01(all_merge@assays$predictionsRGC["8.1-Glut",]))*100, class84 = as.vector(range01(all_merge@assays$predictionsRGC["8.4-Glut",]))*100, class9 = as.vector(range01(all_merge@assays$predictions15["9-Glut",]))*100)
df84_81_2$dif2_1 = df84_81_2$class2 - df84_81_2$class1
df84_81_2$dif2_1_col = magma(100)[cut(df84_81_2$dif2_1,100)]
df84_81_2$dif84_81 = df84_81_2$class84 - df84_81_2$class81
df84_81_2$dif84_81_col = viridis(100)[cut(df84_81_2$dif84_81,100)]
df84_81_2$dif84_81_col[which(df84_81_2$dif84_81 == 0)] = "#FDE72500"
df84_81_2$dif84_81_col[which(df84_81_2$class9 >= my.thesh2)] = "#d600ff"
df84_81_2$dif84_81_col[which(df84_81_2$class3 >= 15)] = "#0df705"
# df84_81_2$dif84_81_col[which(df84_81_2$class4 >= my.thresh)] = "#04D9FF"
df84_81_2$dif84_81_col[which(df84_81_2$class4_htr1d >= my.thresh)] = "#04D9FF"
df84_81_2$dif84_81_col[which(df84_81_2$class4_vipr2 >= my.thresh)] = "#04D9FF50"
df84_81_2$dif84_81_col[which(df84_81_2$class2 >= my.thresh)] = "#FF9E24"
all_merge_hi$class2_81_84 = df84_81_2$dif84_81_col

# pdf(paste0(out_dir, paste0("class2_81_84_4_9_4_", my.thresh, ".pdf")), width = 12, height = 12, onefile = F)
pdf(paste0(out_dir, paste0("test", my.thresh, ".pdf")), width = 12, height = 12, onefile = F)
print(myMultiSFP(all_merge_hi, "class2_81_84", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100)), rm.zero = T, col.ident = T, high.res = T ))
dev.off()

#*******************************************************************************
# Seurat Integration with BB ===================================================
#*******************************************************************************
bb = readRDS(paste0(brain_dir, "data/bb_sct_070522.rds"))
bb_convert15 = data.frame(old = 0:14, new = c("8_Glut", "9_Glut", "4_GABA", "15_GABA/Glut", "1_RGC/MG", "10_Glut", "5_GABA", "11_Glut", "6_GABA", "2_OPC/Oligo", "12_Glut", "13_Glut", "14_Glut", "3_Peri", "7_GABA"))
bb_convert53 = data.frame(old = 0:52, new = c("4.1_GABA", "10.1_Glut", "15.1_GABA/Glut", "9.1_Glut", "8.1_Glut", "1.1_RGC", "6_GABA", "5.1_GABA", "9.2_Glut", "8.2_Glut", "15.2_GABA", "11.1_Glut", "8.3_Glut", "8.4_Glut", "9.3_Glut", "4.2_GABA", "8.5_Glut", "5.2_GABA", "8.6_Glut", "8.7_Glut", "1.2_RGC", "4.3_GABA", "4.4_GABA", "9.4_Glut", "9.5_Glut", "8.8_Glut", "9.6_Glut", "4.5_GABA", "12_Glut", "8.9_Glut", "10.2_Glut", "2.1_OPC", "15.3_GABA", "11.2_Glut", "15.4_GABA", "4.6_GABA", "9.7_Glut", "13_Glut", "14_Glut", "4.7_GABA", "11.3_Glut", "9.8_Glut", "8-9_Glut", "15.5_GABA/Glut", "4.8_GABA", "1.3_MG", "2.2_Oligo", "15.6_Glut", "8.10_Glut", "8.11_Glut", "3_Peri", "15.7_Glut", "7_GABA"))
bb$names15 = bb_convert15$new[match(bb$seuratclusters15, bb_convert15$old)]
bb$names53 = bb_convert53$new[match(bb$seuratclusters53, bb_convert53$old)]

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

my.col = "#27e5fa"
all_merge_hi$glut11 = as.vector(range01(all_merge@assays$predictions15["11-Glut",])) * 100
all_merge_hi$glut11 = all_merge_hi$glut11 > 10
all_merge_hi$glut11 = plyr::revalue(as.character(all_merge_hi$glut11), c("TRUE" = my.col, "FALSE" = "#27e5fa00"))
pdf(paste0(out_dir, paste0("glu11", my.thresh, ".pdf")), width = 12, height = 12, onefile = F)
print(myMultiSFP(all_merge_hi, "glut11", pt.size.multiplier = 1.5, pal = colorRampPalette(viridis(100)), rm.zero = T, col.ident = T, high.res = T ))
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
# Correlations =================================================================
#*******************************************************************************
library("rhdf5")
h5f = H5Fopen("~/scratch/st/results/sct_data_cor.h5")
mat = h5f$name
h5closeAll()

# p-values
h5f = H5Fopen("~/scratch/st/results/sct_data_cor_one_p.h5")
mat_p = h5f$name
h5closeAll()

# mat_p_adj = do.call("rbind", lapply(1:nrow(mat_p), function(x) p.adjust(mat_p[x,], method = "bonferroni") ))

genes = data.frame(gene = rownames(all_merge@assays$SCT@data), sum = rowSums(all_merge@assays$SCT@data), raw_sum = rowSums(all_merge@assays$Spatial@counts[rownames(all_merge@assays$SCT@data),]))
rownames(mat) = colnames(mat) = genes$gene
rownames(mat_p_adj) = colnames(mat_p_adj) = genes$gene
# mat = mat[genes$gene[which(genes$raw_sum > 20)], genes$gene[which(genes$raw_sum > 20)]]
mat[upper.tri(mat, diag = T)] = NA
genes$cor_max = NA
genes$cor_max = unlist(lapply(rownames(mat), function(x) max(mat[x,], na.rm = T)))
genes$cor_max_p = c(NA, unlist(lapply(1:nrow(mat), function(x) mat_p[x, which.max(mat[x,])] )))
# genes$cor_max_p_adj = c(NA, unlist(lapply(1:nrow(mat), function(x) mat_p_adj[x, which.max(mat[x,])] )))
genes$cor_max_gene = c(NA, unlist(lapply(rownames(mat), function(x) colnames(mat)[which.max(mat[x,])] )))
genes$cor_max_raw_sum = genes$raw_sum[match(genes$cor_max_gene, genes$gene)]

no_na_ind = which(!is.na(mat), arr.ind = T)
mat_v = as.vector(mat_p[no_na_ind])
mat_v_adj = p.adjust(mat_v, method = "bonferroni")
mat_p_adj2 = matrix(NA, nrow = nrow(mat), ncol = ncol(mat), dimnames = list(genes$gene, genes$gene))
mat_p_adj2[no_na_ind] = mat_v_adj

cor_to_write = as.data.frame(which(mat_p_adj2 < 0.05, arr.ind=T))
cor_to_write$gene1 = rownames(mat)[cor_to_write$row]
cor_to_write$gene2 = rownames(mat)[cor_to_write$col]
cor_to_write$r = mat[cbind(cor_to_write$row, cor_to_write$col)]
cor_to_write$p = mat_p[cbind(cor_to_write$row, cor_to_write$col)]
cor_to_write$p_adj = mat_p_adj2[cbind(cor_to_write$row, cor_to_write$col)]
# cor_to_write = cor_to_write[which(! duplicated( paste0(cor_to_write$row, "_", cor_to_write$col) )),]
cor_to_write$row = cor_to_write$col = NULL
cor_to_write = cor_to_write[order(cor_to_write$r, decreasing = T),]
data.table::fwrite(cor_to_write, "~/scratch/st/results/spatail_two_tail_p_adj_sig.csv")

library("WGCNA")
this_rowSums = rowSums(all_merge@assays$Spatial@counts)
good_genes = names(this_rowSums)[which(this_rowSums > 0 & names(this_rowSums) %in% rownames(all_merge@assays$SCT@data))]
data_mat_c = t(all_merge@assays$SCT@data[good_genes,])
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft_c = pickSoftThreshold(data_mat_c, powerVector = powers, verbose = 5)
adjacency = adjacency(data_mat_c, type = "signed", power = 1)
TOM = adjacency
dissTOM = 1-TOM

df = data.frame(mod1 = 0, mod2 = 0, mod3 = 0, mod4 = 0, gene = colnames(data_mat_c), row.names = colnames(data_mat_c))
this.m = "ward.D"
geneTree = hclust(as.dist(dissTOM), method = this.m)
pred1 = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 1, pamRespectsDendro = FALSE, minClusterSize = 30);
pred2 = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30);
pred3 = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 3, pamRespectsDendro = FALSE, minClusterSize = 30);
pred4 = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = 30);
di1 = clValid::dunn(as.dist(dissTOM), clusters = pred1); di2 = clValid::dunn(as.dist(dissTOM), clusters = pred2); di3 = clValid::dunn(as.dist(dissTOM), clusters = pred3); di4 = clValid::dunn(as.dist(dissTOM), clusters = pred4);
df$mod1 = pred1; df$mod2 = pred2; df$mod3 = pred3; df$mod4 = pred4;
write.csv(df, "~/scratch/st/results/wgcna_power1_ds_1_3_minclustsize30_082922.csv")



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
st.clust.deg$label = gene_info$label[match(st.clust.deg$gene, gene_info$mzebra)]
st.clust.deg$nd = gene_info$nd_symbol[match(st.clust.deg$gene, gene_info$mzebra)]
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

# Exclusive markers
st.clust.deg = st.clust.deg[which(st.clust.deg$avg_log2FC > 0 & st.clust.deg$pct.dif > 0 & -log10(st.clust.deg$p_val_adj) > 4),]
st.clust.deg$exclusive_score = st.clust.deg$avg_log2FC * st.clust.deg$pct.dif * (1-st.clust.deg$pct.2)
st.clust.deg.exc = st.clust.deg[order(st.clust.deg$exclusive_score, decreasing = T),]
st.clust.deg.exc = st.clust.deg.exc %>% group_by(cluster) %>% slice(1:50)
write.csv(st.clust.deg.exc, "all_merge_cluster_exclusive_markers_standard_111522.csv")

st.clust.deg.sub = st.clust.deg[which(st.clust.deg$n_gene_appears <= 3),]
write.csv(st.clust.deg.sub, "all_merge_cluster_markers_3n_070822.csv")
st.clust.deg.sub = st.clust.deg[which(st.clust.deg$pct.2 < 0.3 & st.clust.deg$pct.1 > 0.5),]
write.csv(st.clust.deg.sub, "all_merge_cluster_markers_subsetpct2_070822.csv")

bri.clust.mark = read.csv("~/Downloads/ST_fig1_B2topgenes_111622.csv")
bri.clust.mark = bri.clust.mark[order(bri.clust.mark$cluster),]
bri.clust.mark[which(bri.clust.mark$label == "gabrp"),] = data.frame(cluster = 8, seurat_name = "hspg2", label = "hspg2")
bri.clust.mark[which(bri.clust.mark$label == "ucn"),] = data.frame(cluster = 16, seurat_name = "aqp9", label = "aqp9")
bri.clust.mark[which(bri.clust.mark$label == "trh"),] = data.frame(cluster = 18, seurat_name = "pax6", label = "pax6")
bri.clust.mark[which(bri.clust.mark$label == "id3"),] = data.frame(cluster = 24, seurat_name = "wdr49", label = "wdr49")
bri.clust.mark$seurat_name = factor(bri.clust.mark$seurat_name, levels = bri.clust.mark$seurat_name)
pdf(paste0(out_dir, "b2_cluster_markers_bri.pdf"), width = 12, height = 6, onefile = F)
DotPlot(all_merge, features = bri.clust.mark$seurat_name) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_color_viridis() + scale_x_discrete(labels = bri.clust.mark$label)
dev.off()

#*******************************************************************************
# B2 ===========================================================================
#*******************************************************************************
my.b2 = all_merge
for (s in c("b2al", "b2ar", "b2bl", "b2br", "b2cl", "b2cr", "b2dl", "b2dr")) {
  real.slice = all_merge@images[[substr(s, 1, 3)]]
  hires = Read10X_Image(paste0("~/Downloads/sp_data/b2_modified/296_", toupper(substr(s, 3, 3)), "1_", toupper(substr(s, 4, 4)), "/"), image.name = "tissue_hires_image.png")
  rownames(hires@coordinates) = paste0(substr(s, 1, 3), "_", rownames(hires@coordinates))
  hires@coordinates = hires@coordinates[rownames(real.slice@coordinates),]
  hires@scale.factors$lowres = hires@scale.factors$hires
  hires@assay = real.slice@assay
  hires@key = real.slice@key
  my.b2@images[[s]]= hires
  # SpatialFeaturePlot(b2d, feature = "egr1", images = "slice1", crop = T)
}
half_by_spot = read.csv(paste0(out_dir, "half.csv"))
my.b2$half = half_by_spot$half[match(colnames(my.b2), half_by_spot$cell)]
my.b2$sample.half = my.b2$sh = paste0(my.b2$sample, substr(my.b2$half, 1, 1))
my.b2$sh = my.b2$sample.half

all_merge$struct2 = factor(all_merge$struct, levels = c("Dm", "Dc-3", "Dd", "Dd-d", "Dl-g", "Dd-v", "Dl-v", "NT", "Dp", "Dc-5", "Dc-4", "Vs", "Vs-m", "Vs-l", "Vp", "Vd-c", "Vc", "Vv", "nPPa", "Vm", "Vl", "Vi", "E", "OB gml", "OB gc", "tract"))
this.pal = c("#1B9CA6", "#077278", "#12b0b3ff", "#00929aff", "#a7ebd4ff", "#089e90ff", "#0f9783ff", "#6ccab8ff", "#1BA668ff", "#6aa943ff", "#7fce4bff", "#f24839ff", "#fb717fff", "#fb8082ff", "#fea7f8ff", "#f367f1ff", "#fd91ffff", "#ff8da2ff", "#ED256E", "#FBCBC7", "#6BBD34", "#96e89bff", "#2fb824ff", "#C95EE6", "#D582EC", "#8135DE")
pdf(paste0(out_dir, "b2_brain_structure.pdf"), width = 16, height = 2, onefile = F)
print(myB2SFP(all_merge, "struct2", pal = this.pal, showLegend = T ))
dev.off()

# all_merge$struct2 = factor(all_merge$struct, levels = c("Dm-3", "Dc-1/2", "Dc-3", "Dc-4", "Dc-5", "Dd", "Dl-d", "Dl-g", "Dl-v", "Dm-1", "Dm-2c", "Dm-2r", "Dp", "NT", "OB gc", "OB gml", "SP-u", "Vc", "Vd-c", "Vd-r", "Vi", "Vl", "Vs", "Vv", "vVZ", "VZ"))
this.pal = c("#7dbdd4", "#077278", "#12b0b3ff", "#00929aff", "#6ccab8ff", "#00929aff", "#0f9783ff", "#00929aff", "#1BA668ff", "#2fb824ff", "#7fce4bff", "#4b9539ff", "#4db3a0", "#a7ebd4ff", "#fea7f8ff", "#f367f1ff", "#FBCBC7", "white", "#C95EE6", "#fd91ffff", "#f24839ff", "#fb8082ff", "#D582EC", "#ED256E", "#fb717fff", "#8135DE", "#434aab")
print(data.frame(struct2 = levels(all_merge$struct), colors = this.pal))
pdf(paste0(out_dir, "c2b2_brain_structure2.pdf"), width = 16, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "struct", pt.size.multiplier = 1.05, pal = this.pal, showLegend = T))
dev.off()

# all_merge$struct2 = factor(all_merge$struct, levels = c("Dm", "Dc-3", "Dd", "Dd-d", "Dl-g", "Dd-v", "Dl-v", "NT", "Dp", "Dc-5", "Dc-4", "Vs", "Vs-m", "Vs-l", "Vp", "Vd-c", "Vc", "Vv", "nPPa", "Vm", "Vl", "Vi", "E", "OB gml", "OB gc", "tract"))
# this.pal = c("#1ba598ff", "#08858bff", "#00a2abff", "#00929aff", "#12b0b3ff", "#089e90ff", "#0f9783ff", "#6ccab8ff", "#1BA668ff", "#6aa943ff", "#7fce4bff", "#f24839ff", "#fb717fff", "#fb8082ff", "#ff8da2ff", "#f367f1ff", "#fd91ffff", "#fea7f8ff", "#96e89bff", "#a7ebd4ff", "#6BBD34", "#88D158", "#2fb824ff", "#C95EE6", "#D177EA", "#8135DE")
# pdf(paste0(out_dir, "b2_test2.pdf"), width = 16, height = 2, onefile = F)
# # print(myB2SFP(all_merge, "struct2", pal = scales::hue_pal()(length(unique(all_merge$struct2))), showLegend = T ))
# print(myB2SFP(all_merge, "struct2", pal = this.pal, showLegend = T ))
# dev.off()

pdf(paste0(out_dir, "b2_cluster_half.pdf"), width = 16, height = 2, onefile = F)
print(myB2SFP(all_merge, "cluster.num", pal = scales::hue_pal()))
dev.off()

pdf(paste0(out_dir, "b2_cluster_markers.pdf"), width = 12, height = 6, onefile = F)
my.labels = gene_info$label[match(top.deg$gene, gene_info$seurat_name)]
print(DotPlot(my.b2, features = top.deg$gene) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_color_viridis() + scale_x_discrete(labels = my.labels))
dev.off()

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
names(all_merge@images) = names(spo)

# Merged object clustering
all_merge = subset(all_merge, subset = nCount_Spatial > 0)
all_merge = SCTransform(all_merge, assay = "Spatial", verbose = FALSE)
all_merge = RunPCA(all_merge, assay = "SCT", verbose = FALSE)
all_merge = RunUMAP(all_merge, reduction = "pca", dims = 1:30)
all_merge = FindNeighbors(all_merge, reduction = "umap", dims = 1:2)
this.res = 0.55
for (this.res in c(0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.7, 0.75, 0.8)) {
  all_merge = FindClusters(all_merge, verbose = FALSE, resolution = this.res)
  all_merge$cluster = all_merge$seurat_clusters
  all_merge$all_cluster = all_merge$seurat_clusters
  p1 = DimPlot(all_merge, reduction = "umap", label = TRUE) + ggtitle(this.res)
  p2 = DimPlot(all_merge, reduction = "umap", label = TRUE, pt.size = 1.25) + ggtitle(this.res)
  print(p1)
  
  my.b2$cluster = all_merge$cluster
  my.b2$cluster.num = as.numeric(as.vector(my.b2$cluster))
  # pdf(paste0(out_dir, "c2b2_cluster_half_", this.res, ".pdf"), width = 16, height = 4, onefile = F)
  Cairo::CairoPNG(paste0(out_dir, "c2b2_cluster_half_", this.res, ".png"), width = 3500, height = 1000, res = 240)
  print(myC2B2SFP(my.b2, "cluster.num", pal = scales::hue_pal()))
  dev.off()
  
  pdf(paste0(out_dir, "c2b2_cluster_half_rot_", this.res, ".pdf"), width = 14, height = 4, onefile = F)
  print(myC2B2SFP(my.b2, "cluster.num", pal = scales::hue_pal(), points.as.text = T, rot.text = T))
  dev.off()
}

my.b2 = all_merge
for (s in c("c2dr", "c2dl", "c2cr", "c2cl", "c2br", "c2bl", "b2al", "b2ar", "b2bl", "b2br", "b2cl", "b2cr", "b2dl", "b2dr")) {
  img.path = ifelse(substr(s,1,1) == "c", 'c2_modified/295_', 'b2_modified/296_')
  real.slice = all_merge@images[[substr(s, 1, 3)]]
  hires = Read10X_Image(paste0("~/Downloads/sp_data/", img.path, toupper(substr(s, 3, 3)), "1_", toupper(substr(s, 4, 4)), "/"), image.name = "tissue_hires_image.png")
  rownames(hires@coordinates) = paste0(substr(s, 1, 3), "_", rownames(hires@coordinates))
  hires@coordinates = hires@coordinates[rownames(real.slice@coordinates),]
  hires@scale.factors$lowres = hires@scale.factors$hires
  hires@assay = real.slice@assay
  hires@key = real.slice@key
  my.b2@images[[s]]= hires
  # SpatialFeaturePlot(b2d, feature = "egr1", images = "slice1", crop = T)
}
half_by_spot = read.csv(paste0(out_dir, "half.csv"))
half_by_spot$half[which(half_by_spot$half == "neither")] = "left"
my.b2$half = half_by_spot$half[match(colnames(my.b2), half_by_spot$cell)]
my.b2$sample.half = my.b2$sh = paste0(my.b2$sample, substr(my.b2$half, 1, 1))
my.b2$sh = my.b2$sample.half

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
# print( myMultiSFP( all_merge_hi, "cluster.num", pal = scales::hue_pal(), rm.zero = F, high.res = T, points.as.text = T, pt.size.multiplier = 1.4 ))
print( myMultiSFP( all_merge_hi, "cluster.num", pal = scales::hue_pal(), rm.zero = F, high.res = T, points.as.text = F, pt.size.multiplier = 1 ))
dev.off()

#*******************************************************************************
# Gene Info 2 ==================================================================
#*******************************************************************************
# TODO: remove "wu:", "sb:", "im:", "zmp:"
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
# 01/23/23
gene_info$human_agree = gene_info$human_mart == gene_info$human_pat
gene_info$human_agree[which(is.na(gene_info$human_agree))] = FALSE
gene_info$human_pat[which(gene_info$human_pat == "CantFind")] = NA
gene_info$human_pat[which(gene_info$human_pat == "")] = NA
gene_info$human_mart[which(gene_info$human_mart == "")] = NA
gene_info$seurat_name_has_human_pat = unlist(lapply(1:nrow(gene_info), function(x) grepl( tolower(gene_info$human_pat[x]), tolower(gene_info$seurat_name[x]), fixed = TRUE)  )) 
gene_info$ens_has_human_pat = unlist(lapply(1:nrow(gene_info), function(x) grepl( tolower(gene_info$human_pat[x]), tolower(gene_info$ens[x]), fixed = TRUE)  )) 
gene_info$nd_has_human_pat = unlist(lapply(1:nrow(gene_info), function(x) grepl( tolower(gene_info$human_pat[x]), tolower(gene_info$nd_symbol[x]), fixed = TRUE)  )) 
gene_info$seurat_name_has_human_mart = unlist(lapply(1:nrow(gene_info), function(x) grepl( tolower(gene_info$human_mart[x]), tolower(gene_info$seurat_name[x]), fixed = TRUE)  )) 
gene_info$ens_has_human_mart = unlist(lapply(1:nrow(gene_info), function(x) grepl( tolower(gene_info$human_mart[x]), tolower(gene_info$ens[x]), fixed = TRUE)  )) 
gene_info$nd_has_human_mart = unlist(lapply(1:nrow(gene_info), function(x) grepl( tolower(gene_info$human_mart[x]), tolower(gene_info$nd_symbol[x]), fixed = TRUE)  )) 
length(which(gene_info$seurat_name_has_human_pat | gene_info$ens_has_human_pat | gene_info$nd_has_human_pat | gene_info$seurat_name_has_human_mart | gene_info$ens_has_human_mart | gene_info$nd_has_human_mart))

gene_info$n_pat = rowSums(gene_info[,c("seurat_name_has_human_pat", "ens_has_human_pat", "nd_has_human_pat")])
gene_info$n_pat[which(gene_info$n_pat > 0)] = gene_info$n_pat[which(gene_info$n_pat > 0)] + 1
gene_info$n_pat[which(is.na(gene_info$n_pat))] = 0
gene_info$n_mart = rowSums(gene_info[,c("seurat_name_has_human_mart", "ens_has_human_mart", "nd_has_human_mart")])
gene_info$n_mart[which(gene_info$n_mart > 0)] = gene_info$n_mart[which(gene_info$n_mart > 0)] + 1
gene_info$n_mart[which(is.na(gene_info$n_mart))] = 0
gene_info$n_max_pat_mart = unlist(lapply(1:nrow(gene_info), function(x) max(gene_info$n_pat[x], gene_info$n_mart[x])))
gene_info$n_agree_any_human = gene_info$n_max_pat_mart + gene_info$human_agree
gene_info$n_pat_n_mart_dif = gene_info$n_pat - gene_info$n_mart
gene_info$greatest_human_n = unlist(lapply(1:nrow(gene_info), function(x) ifelse(gene_info$n_pat_n_mart_dif[x] > 0, gene_info$human_pat[x], gene_info$human_mart[x]) ))

gene_df = data.frame(gene = rownames(bb), bb_counts = rowSums(bb@assays$RNA@counts), st_counts = rowSums(all_merge@assays$Spatial@counts))
gene_df$bb_rank = rank(gene_df$bb_counts)
gene_df$st_rank = rank(gene_df$st_counts)
gene_df$mean_rank = rowMeans(gene_df[,c("bb_rank", "st_rank")])

gene_info$bb_st_counts_rank = gene_df$mean_rank[match(gene_info$seurat_name, gene_df$gene)]
gene_info_one_to_one = gene_info[which(gene_info$n_agree_any_human > 0),] %>% group_by(greatest_human_n) %>% arrange(desc(n_agree_any_human), bb_st_counts_rank) %>% slice(1)
gene_info$one_to_one_human = gene_info_one_to_one$greatest_human_n[match(gene_info$seurat_name, gene_info_one_to_one$seurat_name)]

st.clust.deg[, c("ens", "nd", "human_mart", "human_pat", "one_to_one_human", "n_sources", "bb_st_counts_rank")] = gene_info[match(st.clust.deg$gene, gene_info$seurat_name), c("ens", "nd_symbol", "human_mart", "human_pat", "one_to_one_human", "n_agree_any_human", "bb_st_counts_rank")]

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
# CDG ==========================================================================
#*******************************************************************************
cdg = c("cobl", "ddr1", "fhod3", "grik5", "LOC101476914", "LOC101477204", "LOC101479283", "LOC105941351", "nbeal2", "plekhf2", "plekhg4b", "wdr73")
all_merge$cdg_score = colSums(all_merge@assays$Spatial@counts[cdg,] > 0)
myMultiSFP(all_merge, "cdg_score", pt.size.multiplier = 1.3, pal = colorRampPalette(viridis(100)), rm.zero = T, high.res = F)

#*******************************************************************************
# Brain Structures =============================================================
#*******************************************************************************
# Take our annotations of Brain Strcutures and put them in the Seurat object
all_merge$structure = "unclassified"
for (s in levels(all_merge$sample)) {
  this.fish = substr(s, 1, 2)
  if (this.fish == "c2") { 
    this.struct = read.csv(paste0(out_dir, "structure_295_", toupper(substr(s, 3, 3)), "1_f1.csv"))
  } else {
    this.struct = read.csv(paste0(out_dir, "structure_296_", toupper(substr(s, 3, 3)), "1_vdc.csv"))
  }
  
  this.struct[,1] = paste0(s, "_", this.struct[,1])
  all_merge$structure[this.struct[,1]] = this.struct[,2]
}
all_merge$structure[which(all_merge$structure == "Olf. tract")] = "tract"
all_merge$structure[which(all_merge$structure == "Dc-1")] = "Dc-1/2"
all_merge$structure[which(all_merge$structure == "Dc-2")] = "Dc-1/2"
all_merge$structure[which(all_merge$structure == "Vs-m")] = "Vs-m"
all_merge$structure[which(all_merge$structure == "Vs-l")] = "Vs-l"
all_merge$structure[which(all_merge$structure == "Dm")] = "Dm-3"
all_merge$structure[which(all_merge$structure == "Dl-d")] = "Dl-v"
all_merge$structure[which(all_merge$structure == "tract")] = "ON"
all_merge$structure[which(all_merge$structure == "NT")] = "Dl-vv"
all_merge$structure[which(all_merge$structure == "SP-u")] = "Vx"
b_order = read.csv("~/Downloads/b2c2_struct_order.csv")[,1]
all_merge$structure_no_order = all_merge$structure
all_merge$structure = factor(all_merge$structure, levels = b_order)
all_merge$struct = all_merge$structure
# all_merge$struct_b2_vdc = all_merge$struct

pdf(paste0(out_dir, "b2_brain_structures_annotated.pdf"), width = 16, height = 2, onefile = F)
print(myB2SFP( all_merge, "structure", pal = scales::hue_pal()(length(unique(all_merge$structure))), showLegend = T ))
dev.off()

# Map the brain structures to spatial clusters
struct_clust_table = table(all_merge$struct, all_merge$seurat_clusters)
struct_clust = data.frame(matrix(struct_clust_table, ncol = ncol(struct_clust_table), dimnames = dimnames(struct_clust_table)))
struct_clust = struct_clust / rowSums(struct_clust)
colnames(struct_clust) = 0:(ncol(struct_clust)-1)

struct_clust_melt = reshape2::melt(as.matrix(struct_clust))
struct_clust_melt$Var2 = factor(struct_clust_melt$Var2)
pdf(paste0( "c2b2_brain_stuct_to_cluster.pdf"), width = 5, height = 4.5)
ggplot(struct_clust_melt, aes(x = Var2, y = Var1, fill = value)) + geom_raster() + coord_fixed() + scale_x_discrete(expand = c(0,0), name = "") + scale_y_discrete(expand = c(0,0), name = "") + theme_classic() + scale_fill_viridis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank())
dev.off()

# Map the brain structures to bb clusters
cells.struct.clust = t(rowsum(c2l_mean, group = all_merge$struct, na.rm = T))
cells.struct.clust = cells.struct.clust / colSums(cells.struct.clust)[col(cells.struct.clust)]
cells.struct.clust.melt = reshape2::melt(cells.struct.clust)
ggplot(cells.struct.clust.melt, aes(x = Var2, y = Var1, fill = value)) + geom_raster() + scale_fill_viridis() + xlab("") + ylab("") + coord_fixed() + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))

# Dotplot of the top 2 DEGs of brain structures
# Making brain structure DEGs not shown
top.deg = st.clust.deg %>% group_by(cluster) %>% slice(1:2)
Idents(all_merge) = all_merge$struct
pdf(paste0(out_dir, "b2_stuct_markers.pdf"), width = 12, height = 6, onefile = F)
DotPlot(all_merge, features = top.deg$gene) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_color_viridis() + scale_x_discrete(labels = top.deg$label)
dev.off()

# Dotplot of the top 2 DEGs of brain structures (Brianna)
# top.deg = read.csv('~/Downloads/ST_fig2_B2struct_topgenes_121422.csv')
top.deg = read.csv('~/Downloads/b2c2_fig2_struct_markers_013023.csv')
colnames(top.deg)[1] = "Structure"
top.deg$Structure = trimws(top.deg$Structure)
# top.deg[which(top.deg$seurat_name == "LOC101477225"), c("seurat_name", "label")] = c("st18", "st18")
# top.deg = top.deg[which(top.deg$Structure != "Dc-5"),]
# top.deg$Structure[which(top.deg$Structure == "Dc-4")] = "Dc-4/5"
# top.deg$label[which(top.deg$label == "LOC101464037")] = "sccl21"
# top.deg$label[which(top.deg$label == "LOC101473394")] = "nostrin"
# top.deg$label[which(top.deg$label == "LOC105941321")] = "sdc1"
top.deg$label[which(top.deg$label == "LOC101467413")] = "bcl11ab"
top.deg$label[which(top.deg$label == "LOC101486618")] = "her.2"
top.deg$label[which(top.deg$label == "LOC101464413")] = "kif5c"
top.deg$label[which(top.deg$label == "LOC101471517")] = "CACNA2D4"
# this.levels[which(this.levels == "Dc-4")] = "Dc-4/5"
top.deg$Structure = factor(top.deg$Structure, levels = levels(all_merge$struct))
top.deg = top.deg[order(top.deg$Structure),]
# this.ident = as.character(as.vector(all_merge$struct))
# this.ident[which(this.ident == "Dc-4" | this.ident == "Dc-5")] = "Dc-4/5"
# Idents(all_merge) = factor(this.ident, levels = this.levels)
pdf(paste0(out_dir, "c2b2_stuct_markers_bri.pdf"), width = 13, height = 5.25, onefile = F)
DotPlot(all_merge, features = top.deg$gene) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_color_viridis() + scale_x_discrete(labels = top.deg$label) + xlab("") + ylab("")
dev.off()

pdf(paste0( "c2b2_painted_structure_F.pdf"), width = 3, height = 3)
Idents(all_merge) = all_merge$struct
print(DimPlot(all_merge, cols = this.pal, order = F) + theme_void() + NoLegend() + ggtitle("") + coord_fixed())
dev.off()

# Plot brain structure specific genes on UMAP
bs.genes = read.csv("~/Downloads/struct_umap_list_firstpass.csv")
bs.genes$col = scales::hue_pal()(nrow(bs.genes))
all_merge$struct_gene = "gray90"
all_merge$struct_gene2 = "none"
this.threshold = 0.4
for (i in 1:nrow(bs.genes)) {
  this.score = range01(all_merge@assays$Spatial@counts[bs.genes$seurat_name[i],])
  all_merge$struct_gene[which(this.score > this.threshold)] = bs.genes$col[i]
  all_merge$struct_gene2[which(this.score > this.threshold)] = bs.genes$structure
}
df = data.frame(all_merge@reductions$umap@cell.embeddings)
df$struct_gene = all_merge$struct_gene
df$struct_gene = factor(df$struct_gene, levels = c(unique(df$struct_gene)[which(unique(df$struct_gene) != "gray90")], "gray90"))
df = df[order(df$struct_gene, decreasing = T),]
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = struct_gene)) + geom_point() + scale_color_identity(guide = "legend", labels = bs.genes$structure[match(levels(df$struct_gene), bs.genes$col)]) + theme_void() + ggtitle(paste0("Threshold = ", this.threshold))

#*******************************************************************************
# CellChat =====================================================================
#*******************************************************************************
cc = read.csv("~/scratch/st/results/cellchat/cellchat_st_ct_weights.csv")
convert53 = read.csv("~/scratch/st/data/convert53.csv")
cc$Sender   = unlist(lapply(1:nrow(cc), function(x) paste(strsplit(cc$X[x], "\\.", perl = T)[[1]][1:2], collapse = ".") ))
cc$Sender   = unlist(lapply(1:nrow(cc), function(x) { this.split = strsplit(cc$Sender[x], "_", perl = T)[[1]]; this.start = this.split[1]; this.end = strsplit(this.split[2], "\\.", perl = T)[[1]][1]; paste(c(this.start, this.end), collapse = "_") }  ))
cc$Receiver = unlist(lapply(1:nrow(cc), function(x) strsplit(cc$X[x], "\\d_[A-z0-9/-]+\\.", perl = T)[[1]][2] ))
cc$Sender   = factor(cc$Sender,   levels = convert53$new)
cc$Receiver = factor(cc$Receiver, levels = convert53$new)
pdf("~/scratch/st/results/cellchat_st_ct_weights.pdf", width = 10, height = 10)
ggplot(cc, aes(x = Sender, y = Receiver, fill = x)) + geom_raster() + scale_fill_viridis() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 10), axis.line=element_blank()) + force_panelsizes(rows = unit(length(unique(cc$Sender))/8, "in"), cols = unit(length(unique(cc$Receiver))/8, "in"))
dev.off()
system("rclone copy ~/scratch/st/results/cellchat_st_ct_weights.pdf dropbox:BioSci-Streelman/George/Brain/spatial/analysis/cellchat/")

cc = read.csv("~/scratch/st/results/cellchat/neuronchat_st_structure_weights.csv")
# cc = read.csv("C:/Users/miles/Downloads/neuronchat_st_structure_weights.csv")
cc[,c("Sender", "Receiver")] = reshape2::colsplit(cc$X, "\\.", c('1', '2'))
pdf("~/scratch/st/results/neuronchat_st_structure_weights.pdf", width = 10, height = 10)
ggplot(cc, aes(x = Sender, y = Receiver, fill = x)) + geom_raster() + scale_fill_viridis() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 10), axis.line=element_blank()) + force_panelsizes(rows = unit(length(unique(cc$Sender))/8, "in"), cols = unit(length(unique(cc$Receiver))/8, "in"))
dev.off()
system("rclone copy ~/scratch/st/results/neuronchat_st_structure_weights.pdf dropbox:BioSci-Streelman/George/Brain/spatial/analysis/cellchat/")

# Graph the network
library("igraph")
cc.backup = cc

cc = cc.backup
cc = cc[which(cc$x > quantile(cc$x, 0.9)),c("Sender", "Receiver", "x")]
colnames(cc)[3] = c("weight")
my.nodes = data.frame(label = unique(c(cc$Sender, cc$Receiver)))
my.nodes$sum = unlist(lapply(my.nodes$label, function(x) sum(cc.backup$x[which(cc$Sender == x | cc$Receiver == x)]) ))
my.nodes$color = scales::hue_pal()(nrow(my.nodes))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
my.nodes$sum = 15 * (range01(my.nodes$sum) + 0.5)
cc$weight = 2.5 * (range01(cc$weight) + 0.25)
g1 = graph_from_data_frame(cc, vertices = my.nodes)
V(g1)$color = my.nodes$color
V(g1)$size = my.nodes$sum
V(g1)$label.color = "black"
V(g1)$frame.color = NA
# E(g1)$color = E(g1)$col
E(g1)$width = E(g1)$weight
lfr = layout_with_fr(g1)
plot.igraph(g1, edge.arrow.size=.2)

tkid <- tkplot(g1, vertex.label=my.nodes$label, vertex.label.dist=1)
l <- tkplot.getcoords(tkid)
tk_close(tkid, window.close = T)

pdf("C:/Users/miles/Downloads/st_cc_structure_network.pdf", width = 6, height = 6)
plot.igraph(g1, layout = l, edge.arrow.size = 0.4)
dev.off()

#*******************************************************************************
# CCI ==========================================================================
#*******************************************************************************

#------------------------------------#
# Cell Type to Cell Type Connections #
#------------------------------------#
cci.p.list = list()
full.mat = expand.grid(convert53$new, convert53$new)
colnames(full.mat) = c("Sender", "Receiver")
full.mat$id = paste0(full.mat$Sender, ".", full.mat$Receiver)
mode = "sum_logp"
my.sh = c("c2dr", "c2dl", "c2cr", "c2cl", "c2br", "c2bl", "b2ar", "b2al", "b2br", "b2cr", "b2dr", "b2bl", "b2cl", "b2dl")
for (s in my.sh) {
  this.spa = readRDS(paste0("~/research/st/data/cci/", s, "_human.rds"))
  this.cci = this.spa@lrpair
  if (mode == "sum_logp") {
    this.cci$logp = -log10(this.cci$lr_co_ratio_pvalue)
    this.cci$logp[which(is.infinite(this.cci$logp))] = 4
    this.cci$id = paste0(this.cci$celltype_sender, ".", this.cci$celltype_receiver)
    this.cci.mat.melt = data.frame(id = unique(this.cci$id), value = unlist(lapply(unique(this.cci$id), function(x) sum(this.cci$logp[which(this.cci$id == x)]))))
    this.cci.mat.melt[, c("Sender", "Receiver")] = reshape2::colsplit(this.cci.mat.melt$id, "\\.", c('1', '2'))
    this.cci.mat.melt = this.cci.mat.melt[, c("Sender", "Receiver", "value", "id")]
  } else if (mode == "count") {
    this.cci.table = table(this.cci$celltype_sender, this.cci$celltype_receiver)
    this.cci.mat = matrix(this.cci.table, ncol = length(unique(this.cci$celltype_sender)), dimnames = dimnames(this.cci.table))
    this.cci.mat.melt = reshape2::melt(as.matrix(this.cci.mat))
  }
  this.cci.mat.melt[,1]   = convert53$new[match(this.cci.mat.melt[,1], convert53$old)]
  this.cci.mat.melt[,2] = convert53$new[match(this.cci.mat.melt[,2], convert53$old)]
  this.cci.mat.melt$id = paste0(this.cci.mat.melt[,1], ".", this.cci.mat.melt[,2])
  this.full.mat = full.mat[,1:3]
  this.full.mat$Sender   = factor(this.full.mat$Sender, levels = convert53$new)
  this.full.mat$Receiver = factor(this.full.mat$Receiver, levels = convert53$new)
  this.full.mat$value = this.cci.mat.melt[match(this.full.mat$id, this.cci.mat.melt$id),3]
  this.full.mat$value[which(is.na(this.full.mat$value))] = 0
  full.mat[,s] = this.full.mat$value
  # this.full.mat$value[which(this.full.mat$value > 60)] = 60
  p=ggplot(this.full.mat, aes(x = Receiver, y = Sender, fill = value)) + geom_raster() + scale_fill_viridis() + theme_classic() + scale_x_discrete(expand = c(0,0), drop = F) + scale_y_discrete(expand = c(0,0), drop = F) + coord_fixed() + ggtitle(paste0(s, ": 60+")) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  cci.p.list[[s]] = p + xlab("") + ylab("")
  print(p)
}

# Correlations using specific LR interactions
# all.lr.pairs = c()
# for (s in my.sh) {
#   this.spa = readRDS(paste0("~/research/st/data/cci/", s, "_human.rds"))
#   this.cci = this.spa@lrpair
#   this.cci$lr.pair = paste0(this.cci$ligand, "_", this.cci$receptor)
#   all.lr.pairs = c(all.lr.pairs, this.cci$lr.pair)
# }
# all.lr.pairs = sort(unique(all.lr.pairs))
# full.mat.lr = expand.grid(convert53$new, convert53$new, all.lr.pairs)
# colnames(full.mat.lr) = c("Sender", "Receiver", "lr.id")
# full.mat.lr[, c("Ligand", "Receptor")] = reshape2::colsplit(full.mat.lr$lr.id, "_", c('1', '2'))
# full.mat.lr$pop.id = paste0(full.mat.lr$Sender, ".", full.mat.lr$Receiver)
# full.mat.lr$full.id = paste0(full.mat.lr$pop.id, "_", full.mat.lr$lr.id)
# rownames(full.mat.lr) = full.mat.lr$full.id
# my.sh = c("c2dr", "c2dl", "c2cr", "c2cl", "c2br", "c2bl", "b2ar", "b2al", "b2br", "b2cr", "b2dr", "b2bl", "b2cl", "b2dl")
# full.mat.lr[, my.sh] = 0
# for (s in my.sh) {
#   this.spa = readRDS(paste0("~/research/st/data/cci/", s, "_human.rds"))
#   this.cci = this.spa@lrpair
#   this.cci$logp = -log10(this.cci$lr_co_ratio_pvalue)
#   this.cci$logp[which(is.infinite(this.cci$logp))] = 4
#   this.cci$celltype_sender_name   = convert53$new[match(this.cci$celltype_sender,    convert53$old)]
#   this.cci$celltype_receiver_name = convert53$new[match(this.cci$celltype_receiver , convert53$old)]
#   this.cci$pop.id = paste0(this.cci$celltype_sender_name, ".", this.cci$celltype_receiver_name)
#   this.cci$lr.id   = paste0(this.cci$ligand, "_", this.cci$receptor)
#   this.cci$full.id = paste0(this.cci$pop.id, "_", this.cci$lr.id)
#   full.mat.lr[this.cci$full.id, s] = this.cci$logp
# }

pdf(paste0(out_dir, "c2b2_bb53_cci.pdf"), width = 16, height = 16, onefile = F)
print(cowplot::plot_grid(plotlist = cci.p.list))
dev.off()

# correlation of samples with each other
cor.mat = cor(as.matrix(full.mat[,4:ncol(full.mat)]))
cor.hcl = hclust(as.dist(1 - abs(cor.mat)))
test = cutree(cor.hcl, k = 2)
full.mat.cor = reshape2::melt(cor.mat)
full.mat.cor$value[which(full.mat.cor$value < 0)] = 0
# full.mat.cor$value[which(full.mat.cor$value == 1)] = NA
pdf(paste0(out_dir, "c2b2_halves_bb53_cor.pdf"), width = 2.5, height = 2.5, onefile = F)
print(ggplot(full.mat.cor, aes(x = Var1, y = Var2, fill = value)) + geom_raster() + scale_fill_viridis(limits = c(0, 1)) + theme_classic() + scale_x_discrete(expand = c(0,0), drop = F) + scale_y_discrete(expand = c(0,0), drop = F) + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("") + NoLegend())
dev.off()

# Are the anterior-posterior correlations greater than permutations?
library(combinat)
perm.sh = data.frame(do.call('rbind', permn(my.sh)))
colnames(perm.sh) = my.sh
perm.sh$ant.mean.r = unlist(lapply(1:nrow(perm.sh), function(x) mean(cor.mat[unlist(perm.sh[x,1:5]),unlist(perm.sh[x,1:5])]) )) 
perm.sh$pos.mean.r = unlist(lapply(1:nrow(perm.sh), function(x) mean(cor.mat[unlist(perm.sh[x,6:8]),unlist(perm.sh[x,6:8])]) )) 
perm.sh$sum.mean.r = perm.sh$ant.mean.r + perm.sh$pos.mean.r
perm.sh$outside.mean.r = unlist(lapply(1:nrow(perm.sh), function(x) mean(cor.mat[unlist(perm.sh[x,6:8]),unlist(perm.sh[x,1:5])]) )) 

# correlation of sample halves by distance
half.mat.cor = cor.mat
half.mat.cor[lower.tri(cor.mat, diag = T)] = NA
half.mat.cor = reshape2::melt(half.mat.cor)
half.mat.cor = half.mat.cor[which(!is.na(half.mat.cor$value)),]
half.mat.cor$Var1.idx = match(half.mat.cor$Var1, my.sh)
half.mat.cor$Var2.idx = match(half.mat.cor$Var2, my.sh)
half.mat.cor$hop = factor(abs(half.mat.cor$Var1.idx - half.mat.cor$Var2.idx))
pdf(paste0(out_dir, "b2_halves_cci_hop.pdf"), width = 3, height = 2.5)
# ggplot(half.mat.cor, aes(x = hop, y = value, group = hop)) + geom_boxplot() + theme_bw() + ylab("Correlation") + xlab("Hop")
print(ggplot(half.mat.cor, aes(x = hop, y = value, group = 1)) + geom_point(alpha=0.2, position = position_jitter(width = 0.1)) + theme_bw() + ylab("Correlation") + xlab("Hop") + stat_summary(fun = "mean", geom = "point") + stat_summary(fun = "mean", geom = "line") + stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.4))
dev.off()
hop.lm = summary(lm(value ~ as.numeric(hop), data = half.mat.cor))
hop.lm$coefficients[2,4]

ap.changes     = full.mat %>% group_by(Sender)   %>% summarise(spop = sum(abs(lm.slope)))
ap.changes[,3] = full.mat %>% group_by(Receiver) %>% summarise(rpop = sum(abs(lm.slope))) %>% select("rpop")
ap.changes = reshape2::melt(ap.changes)
# full.mat2 = data.frame(values = abs(full.mat$lm.slope[which(full.mat$Sender == "8-9_Glut")]), pop = "8-9_Glut")
full.mat2 = data.frame()
for (i in rownames(ap.changes)[order(ap.changes[,3], decreasing = T)][1:5]) { this.col = ifelse(ap.changes[i,2] == "spop", "Sender", "Receiver"); this.df = data.frame(values = abs(full.mat$lm.slope[which(full.mat[,this.col] == ap.changes[i,1])]), pop = paste0(this.col, ".", ap.changes[i,1])); full.mat2 = rbind(full.mat2, this.df) }
full.mat2$pop = factor(full.mat2$pop, levels = unique(full.mat2$pop))
pdf(paste0(out_dir, "b2_halves_cci_top_changes.pdf"), width = 3, height = 4)
# print(ggplot(full.mat2, aes(x = pop, y = values)) + geom_point(position = position_jitter()) + stat_summary(fun = "mean", geom = "point", color = "red") + stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.4, color = "red"))
print(ggplot(full.mat2, aes(x = pop, y = values, color = pop, fill = pop)) + geom_bar(stat = "identity", alpha = 0.3) + scale_x_discrete(name="") + scale_y_continuous(expand=c(0,0), name="Change From Anterior-Posterior (Slope)") + theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + NoLegend())
dev.off()

# What are the specific connections that change the most in 1.2_RG
n.changes = 3
# rg.changes = reshape2::melt(full.mat[which(full.mat$Receiver == "1.2_RG"), c(3:11)], id.var = "id")
# rg.changes$slope = full.mat$lm.slope[match(rg.changes$id, full.mat$id)]
rg.changes = reshape2::melt(full.mat[which(full.mat$Receiver == "1.2_RG" & full.mat$b2cl != 0 & full.mat$b2dl != 0), c(3:11)], id.var = "id")
rg.changes$slope = full.mat$lm.slope[match(rg.changes$id, full.mat$id)]
# tmp.rg.changes = rg.changes[which(rg.changes$variable %in% c("b2ar", "b2al") & rg.changes$value != 0),]
rg.changes.p = rg.changes[order(rg.changes$slope, decreasing = T)[1:(8*n.changes)],]
rg.changes = reshape2::melt(full.mat[which(full.mat$Receiver == "1.2_RG" & full.mat$b2ar != 0 & full.mat$b2al != 0), c(3:11)], id.var = "id")
rg.changes$slope = full.mat$lm.slope[match(rg.changes$id, full.mat$id)]
rg.changes.p = rbind(rg.changes.p, rg.changes[order(rg.changes$slope, decreasing = F)[1:(8*n.changes)],])
rg.changes.p$dir = sign(rg.changes.p$slope)
pdf(paste0(out_dir, "b2_halves_cci_top_changes_in_rg.pdf"), width = 4, height = 4)
# print(ggplot(rg.changes.p, aes(x = variable, y = value, color = id, group = id)) + geom_point() + geom_line() + facet_wrap(~ dir, ncol = 1) + theme_bw())
print(ggplot(rg.changes.p, aes(x = variable, y = value, color = id, group = id)) + geom_point() + geom_line() + facet_wrap(~ dir, ncol = 1) + theme_bw() + NoLegend())
dev.off()

# Mean connection weights
full.mat$mean = rowMeans(full.mat[,4:ncol(full.mat)])
full.mat$mean[which(full.mat$mean > 150)] = 150
p=ggplot(full.mat, aes(x = Receiver, y = Sender, fill = mean)) + geom_raster() + scale_fill_viridis() + theme_classic() + scale_x_discrete(expand = c(0,0), drop = F) + scale_y_discrete(expand = c(0,0), drop = F) + coord_fixed() + ggtitle(": 60+") + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
pdf(paste0(out_dir, "b2_halves_bb53_cci_mean.pdf"), width = 8.5, height = 8.5, onefile = F)
print(p)
dev.off()

# How consistent are the weights across samples?
full.mat$sd = matrixStats::rowSds(as.matrix(full.mat[,4:11]))
full.mat$sd_plot = 1/full.mat$sd
full.mat$sd_plot[which(is.infinite(full.mat$sd_plot))] = max(full.mat$sd_plot[which(is.finite(full.mat$sd_plot))])
full.mat$sd_plot[which(full.mat$mean  <= 5)] = 0
pdf(paste0(out_dir, "b2_bb53_cci_mean.pdf"), width = 8, height = 8, onefile = F)
print(ggplot(full.mat, aes(x = Receiver, y = Sender, fill = mean)) + geom_raster() + scale_fill_viridis() + theme_classic() + scale_x_discrete(expand = c(0,0), drop = F) + scale_y_discrete(expand = c(0,0), drop = F) + coord_fixed() + ggtitle(": 60+") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()
pdf(paste0(out_dir, "b2_bb53_cci_inverse_std.pdf"), width = 8, height = 8, onefile = F)
print(ggplot(full.mat, aes(x = Receiver, y = Sender, fill = sd_plot)) + geom_raster() + scale_fill_viridis() + theme_classic() + scale_x_discrete(expand = c(0,0), drop = F) + scale_y_discrete(expand = c(0,0), drop = F) + coord_fixed() + ggtitle(": 60+") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
dev.off()

full.mat$b2a_dif = full.mat$b2a - full.mat$mean
p=ggplot(full.mat, aes(x = Receiver, y = Sender, fill = b2a_dif)) + geom_raster() + scale_fill_viridis() + theme_classic() + scale_x_discrete(expand = c(0,0), drop = F) + scale_y_discrete(expand = c(0,0), drop = F) + coord_fixed() + ggtitle(": 60+") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Are there significant differences across halves?
full.mat$num.zero = unlist(lapply(1:nrow(full.mat), function(x) length(which(full.mat[x,4:11] == 0))))
full.mat$lm.slope = unlist(lapply(1:nrow(full.mat), function(x) lm(value ~ order, data.frame(value = unlist(full.mat[x,4:11]), order = 1:8))[1]$coefficients[2] ))

# Plot some of the strongest and most consistent connections
full.mat$tmp = full.mat$sd_plot * full.mat$mean
full.mat.melt = reshape2::melt(full.mat[,c("id", "b2a", "b2b", "b2c", "b2d")], id.var = "id")
full.mat.melt$id = factor(full.mat.melt$id, levels = full.mat$id[order(full.mat$tmp, decreasing = T)])
ggplot(full.mat.melt[which(full.mat.melt$id %in% levels(full.mat.melt$id)[1:5]),], aes(x = id, y = value)) + geom_boxplot() + geom_point(position = position_jitter()) + ylab("Mean") + theme_bw()

# Graph the cell type -> cell type network
cc.res.mat.df = this.cci.mat.melt
cc.res.mat.df = cc.res.mat.df[which(cc.res.mat.df$value >= 30),]
cc.res.mat.df$col = convert53$col[match(cc.res.mat.df$Sender, convert53$new)]
my.nodes = data.frame(gene = unique(c(cc.res.mat.df$Sender, cc.res.mat.df$Receiver)), label = unique(c(cc.res.mat.df$Sender, cc.res.mat.df$Receiver)))
my.nodes$col = convert53[match(my.nodes$label, convert53$new), c("col")]
my.nodes$sum = unlist(lapply(my.nodes$label, function(x) sum(cc.res.mat.df$value[which(cc.res.mat.df$Sender == x | cc.res.mat.df$Receiver == x)]) ))
g1 = graph_from_data_frame(cc.res.mat.df, vertices = my.nodes)
V(g1)$color = my.nodes$col
V(g1)$size = log2(my.nodes$sum)
V(g1)$label.color = "black"
V(g1)$frame.color = NA
E(g1)$color = E(g1)$col
lfr = layout_with_fr(g1)
plot.igraph(g1, edge.arrow.size=.2)

tkid <- tkplot(g1, vertex.label=my.nodes$label, vertex.label.dist=1)
l <- tkplot.getcoords(tkid)
tk_close(tkid, window.close = T)

pdf("C:/Users/miles/Downloads/b2a_cci_thresh30.pdf", width = 7, height = 7)
plot.igraph(g1, layout = l, edge.arrow.size = 0.5)
dev.off()

#---------------------------------------#
# Specific Ligand-Receptor Interactions #
#---------------------------------------#
l.df  = data.frame()
r.df  = data.frame()
lr.df = data.frame()
mode = "count"
for (s in c("b2a", "b2b", "b2c", "b2d")) {
  this.spa = readRDS(paste0("~/research/st/data/cci/", s, "_human.rds"))
  this.cci = this.spa@lrpair
  if (mode == "sum_logp") {
    print("Not Yet Implemented")
  } else if (mode == "count") {
    l.df  = rbind(l.df, cbind(as.data.frame(table(this.cci$ligand)), s)  )
    r.df  = rbind(r.df,  cbind(as.data.frame(table(this.cci$receptor)), s) )
    lr.df = rbind(lr.df, cbind(as.data.frame(table(paste0( this.cci$ligand, ".", this.cci$receptor ))), s) )
  }
}

l.df.mat = reshape2::dcast(l.df, Var1 ~ s, value.var = "Freq")
l.df.mat[is.na(l.df.mat)] = 0
l.df.mat$mean = rowMeans(l.df.mat[,2:5])

r.df.mat = reshape2::dcast(r.df, Var1 ~ s, value.var = "Freq")
r.df.mat[is.na(r.df.mat)] = 0
r.df.mat$mean = rowMeans(r.df.mat[,2:5])

lr.df.mat = reshape2::dcast(lr.df, Var1 ~ s, value.var = "Freq")
lr.df.mat[is.na(lr.df.mat)] = 0
lr.df.mat$mean = rowMeans(lr.df.mat[,2:5])


#---------------------#
# Single cell spatial #
#---------------------#
stsc.meta = data.frame()
for (s in c("c2dr", "c2dl", "c2cr", "c2cl", "c2br", "c2bl", "b2ar", "b2al", "b2br", "b2cr", "b2dr", "b2bl", "b2cl", "b2dl")) {
  this.obj = readRDS(paste0("~/research/st/data/infer_cell/", s, ".rds"))
  if (s == "c2dr") { stsc.mat = this.obj@data$newdata } else { stsc.mat = cbind(stsc.mat, this.obj@data$newdata) }
  stsc.meta = rbind(stsc.meta, this.obj@meta$newmeta)
}
stsc.meta$cell = paste0("C", 1:nrow(stsc.meta))
colnames(stsc.mat) = paste0("C", 1:nrow(stsc.meta))

pdf("~/research/st/results/testing_stsc_mini.pdf", width = 16, height = 2)
print(myB2SFP(all_merge, "egr1", stsc.list = list(stsc.mat, stsc.meta), pal = colorRampPalette(viridis(100))))
dev.off()
Cairo::CairoPNG("~/research/st/results/testing_stsc_mini.png", width = 4000, height = 500, res = 240)
print(myB2SFP(all_merge, "egr1", stsc.list = list(stsc.mat, stsc.meta), pal = colorRampPalette(viridis(100))))
dev.off()

#*******************************************************************************
# GWAS =========================================================================
#*******************************************************************************
# 5x10-8, the traditional level of genome-wide significance
gwas = read.csv("~/Downloads/efotraits_MONDO_0005180-associations-2022-11-29.csv")
gwas = read.csv("C:/Users/miles/Downloads/efotraits_MONDO_0004975-associations-2022-11-28.csv") # gwas = read.csv("~/Downloads/efotraits_MONDO_0004975-associations-2022-11-28.csv")
deg15 = read.csv("C:/Users/miles/Downloads/bdeg_gdeg_qdeg_15cluster_summary.csv") # deg15 = read.csv("~/Downloads/bdeg_gdeg_qdeg_15cluster_summary.csv")
deg53 = read.csv("C:/Users/miles/Downloads/bdeg_gdeg_qdeg_53cluster_summary.csv") # deg53 = read.csv("~/Downloads/bdeg_gdeg_qdeg_53cluster_summary.csv")
gwas = as.data.frame(tidyr::separate_rows(gwas, Mapped.gene, sep=",\\s+"))
gwas$P.value.e = as.numeric(reshape2::colsplit(gwas$P.value, "-", c('1', '2'))[,2])
this.split = strsplit(gwas$P.value, ' ')
split.df = as.data.frame(do.call(rbind, this.split))
split.df$e = reshape2::colsplit(split.df[,3], '-', c('1', '2'))[,2]
split.df$num = as.numeric(split.df[,1]) * 10^-as.numeric(split.df$e)
gwas$p = split.df$num
gwas = gwas[which(gwas$p < 5e-8),]
gwas.genes = sort(unique(unlist(as.list(strsplit(gwas$Mapped.gene, ', ')))))
gwas.genes = gwas.genes[2:length(gwas.genes)]
gwas.genes.df = data.frame(human = gwas.genes, mzebra = gene_info$seurat_name[match(gwas.genes, gene_info$human)])
gwas.genes.df = gwas.genes.df[which(!is.na(gwas.genes.df$mzebra)),]
gwas.genes.df$mzebra.num.spots = rowSums(all_merge@assays$Spatial@counts[gwas.genes.df$mzebra,] > 0)
gwas.genes.df = gwas.genes.df[which(gwas.genes.df$mzebra.num.spots >= 3),]
gwas.deg15 = deg15[which(deg15$human %in% gwas.genes),]
gwas.deg53 = deg53[which(deg53$human %in% gwas.genes),]
deg15 = deg15[order(deg15$category, deg15$hmp),]
deg15$rank = 1:nrow(deg15)
tmp = reshape2::melt(deg15[which(deg15$category == "building"), c("rank", "p_min", "p_max", "hmp")], id.var = "rank")
ggplot(tmp, aes(x = rank, y = value, color = variable)) + geom_point()

my.cor = cor(t(as.matrix(all_merge@assays$Spatial[gwas.genes.df$mzebra,])))
diag(my.cor) = NA
my_callback = function(hcl, mat) { print(hcl); hcl <<- hcl; return(hcl) }
pheatmap::pheatmap(my.cor, clustering_callback = my_callback, show_rownames = F, show_colnames = F)
hcl.genes = colnames(my.cor)[hcl$order]
test = stats::cutree(hcl, k = 2)
test[hcl$labels[hcl$order]]

#*******************************************************************************
# Allen Brain Atlas (ABA) ======================================================
#*******************************************************************************
library("cocoframer")
library(purrr)
library(viridisLite) # optional - nice color palettes
library(dplyr)
library(reshape2)

# Get structure ontology annotations
ga <- get_ccf_grid_annotation()
ontology <- get_mba_ontology()
ontology_df <- flatten_mba_ontology(ontology)
# filtered_ontology_df = filter_mba_ontology_children(ontology_df, "CTX")

# build a 3d array of ontology structure acronyms - easier to deal with than IDs
oa <- array(ontology_df$acronym[match(ga, ontology_df$id)], dim = dim(ga))
# oa <- array(filtered_ontology_df$acronym[match(ga, filtered_ontology_df$id)], dim = dim(ga))

deg = read.csv("~/research/st/results/b2_cluster_markers_standard_111622.csv")
all.gene.ish = readRDS("~/research/st/data/aba_all_gene_list2.rds")
all.annot = sort(unique(as.vector(oa)))
all.annot = all.annot[which(!is.na(all.annot))]

# Method 2: Correlation of Gene Expression in Spatial w/ ABA
# myCalcPerSlice = function(this.slice) {
#   message(this.slice)
#   this.mat <- data.frame(matrix(0L, nrow = length(all.gene.ish), ncol = length(all.annot), dimnames = list(names(all.gene.ish), all.annot)))
#   this.annot = as.vector(slice_ccf_arr(oa, this.slice, "coronal"))
#   this.annot = factor(this.annot, levels = all.annot)
#   for (this.gene in names(all.gene.ish)) {
#     this.gene.exp = all.gene.ish[[this.gene]][this.slice,,]
#     this.gene.exp = as.vector(all.gene.ish[[this.gene]][this.slice,,])
#     this.gene.exp[which(this.gene.exp < 0)] = 0
#     this.sum = tapply(this.gene.exp, this.annot, sum)
#     this.sum[which(is.na(this.sum))] = 0
#     # print(length(which(this.sum)))
#     this.mat[this.gene,] = this.mat[this.gene,] + this.sum
#   }
#   print(length(which(is.na(this.mat))))
#   return(this.mat)
# }
# myCalcAnnotPerSlice = function(this.slice) {
#   this.annot = as.vector(slice_ccf_arr(oa, this.slice, "coronal"))
#   this.annot = factor(this.annot, levels = all.annot)
#   return(as.vector(table(this.annot))*length(all.gene.ish))
# }
# slice.mats = mclapply(1:58, function(x) myCalcPerSlice(x), mc.cores = 5)
# aba.gene.annot = as.data.frame(purrr::map(transpose(slice.mats), reduce, `+`))
# rownames(aba.gene.annot) = names(all.gene.ish)
# colnames(aba.gene.annot) = all.annot
# # aba.gene.annot[which(is.na(aba.gene.annot))] = 0
# annot.to.add = mclapply(1:58, function(x) myCalcAnnotPerSlice(x), mc.cores = 5)
# annot.sum.vect = do.call('rbind', annot.to.add)
# annot.sum.vect = colSums(annot.sum.vect)
# names(annot.sum.vect) = all.annot
# 
# aba.gene.annot.pct = aba.gene.annot/annot.sum.vect
# aba.gene.annot.pct = as.matrix(aba.gene.annot / annot.sum.vect[col(aba.gene.annot)])
# aba.gene.annot.pct[which(is.na(aba.gene.annot.pct))] = 0
# aba.gene.annot.pct = as.data.frame(aba.gene.annot.pct)
# saveRDS(aba.gene.annot, "~/research/st/data/aba_gene_by_annot.rds")
# saveRDS(aba.gene.annot.pct, "~/research/st/data/aba_gene_by_annot_pct.rds")

# Method 2: Start Here to Save Time
ga <- get_ccf_grid_annotation()
ontology <- get_mba_ontology()
ontology_df <- flatten_mba_ontology(ontology)
ontology_df$acronym[which(grepl("-", ontology_df$acronym))] = str_replace(ontology_df$acronym[which(grepl("-", ontology_df$acronym))], "-", ".")
ontology_df$acronym[which(grepl("/", ontology_df$acronym))] = str_replace(ontology_df$acronym[which(grepl("/", ontology_df$acronym))], "/", ".")
oa <- array(ontology_df$acronym[match(ga, ontology_df$id)], dim = dim(ga))
aba.gene.annot = readRDS("~/research/st/data/aba_gene_by_annot.rds")
all.gene.ish = readRDS("~/research/st/data/aba_all_gene_list2.rds")
all.annot = sort(unique(as.vector(oa)))
all.annot = all.annot[which(!is.na(all.annot))]
annot.to.add = mclapply(1:58, function(x) myCalcAnnotPerSlice(x), mc.cores = 5)
annot.sum.vect = do.call('rbind', annot.to.add)
annot.sum.vect = colSums(annot.sum.vect)
names(annot.sum.vect) = all.annot

diveOntology = function(id) {
  child.id = ontology_df$id[which(ontology_df$parent_structure_id == id)]
  deep.child.id = c(id)
  for (id.i in child.id) {
    deep.child.id = c(deep.child.id, diveOntology(id.i))
  }
  return (deep.child.id)
}
botUpOntology = function(id) {
  this.level   = ontology_df$st_level[which(ontology_df$id == id)]
  parent.id    = ontology_df$parent_structure_id[which(ontology_df$id == id)]
  parent.level = ontology_df$st_level[which(ontology_df$id == parent.id)]
  if (parent.level < 6) { 
    return(id) 
  } else if (parent.level == 6) {
    return(parent.id)
  } else {
    return(botUpOntology(parent.id))
  }
}
# diveOntology6 = function(id) {
#   child.id = ontology_df$id[which(ontology_df$parent_structure_id == id)]
#   deep.child.id = c(id)
#   for (id.i in child.id) {
#     deep.child.id = c(deep.child.id, diveOntology(id.i))
#   }
#   return (deep.child.id)
# }
child.ids = c(567, diveOntology(567))
my_ontology_df = ontology_df[which(ontology_df$id %in% child.ids),]
my_ontology_df6 = my_ontology_df[which(my_ontology_df$st_level >= 6),]
# for (id6 in my_ontology_df$id[which(my_ontology_df$st_level == 6)]) {
#   child6 = diveOntology(id6)
#   my_ontology_df6$parent_structure_id[which(my_ontology_df6$id %in% child6)] == id6
# }
# my_ontology_df6$parent_structure_id[which(my_ontology_df6$st_level == 6)] = my_ontology_df6$id[which(my_ontology_df6$st_level == 6)]
my_ontology_df6$parent_structure_id = unlist(mclapply(my_ontology_df6$id, function(x) botUpOntology(x), mc.cores = 5))
my_ontology_df6$parent_st_level = ontology_df$st_level[match(my_ontology_df6$parent_structure_id, ontology_df$id)]
# my_ontology_df6$this_sum = annot.sum.vect[my_ontology_df6$acronym]
my_ontology_df6$this_sum = cnu.aba$avg[match(my_ontology_df6$atlas_id, cnu.aba$atlas_id)]
my_ontology_df6 = data.frame(my_ontology_df6 %>% group_by(parent_structure_id) %>% mutate(parent_sum = sum(this_sum, na.rm = T)))

# *** Getting Brianna's Children ***#
b_struct = read.csv("~/Downloads/aba_region_list1_gmod.csv")
b_struct$isSmall = T
b_struct = rbind(b_struct, data.frame(LIST.Old = c("OLF", "sAMY", "CTX"), LIST.New = c("OLF", "sAMY", "U_CTX"), isSmall = F))
b_ont_df = data.frame()
for (i in 1:nrow(b_struct)) {
  struct = b_struct$LIST.Old[i]
  this.id = ontology_df$id[which(ontology_df$acronym == struct)]
  if ( b_struct$isSmall[i] ) { child.ids = diveOntology(this.id) } else { child.ids = this.id }
  child.df = ontology_df[match(child.ids, ontology_df$id),]
  child.df$parent_id = this.id
  child.df$parent_acr_old = struct
  child.df$parent_acr_new = b_struct$LIST.New[i]
  if ( struct == "CTX" ) { print(i); print(child.df) }
  b_ont_df = rbind(b_ont_df, child.df)
}
write.csv(b_ont_df, "~/Downloads/b_ont_df.csv")
length(which(oritz$ABA_acronym %in% b_ont_df$acronym))
oritz$b_parent = b_ont_df$parent_acr_new[match(oritz$ABA_acronym, b_ont_df$acronym)]
oritz2 = subset(oritz, cells = colnames(oritz)[which(!is.na(oritz$b_parent))])
# ***          End               ***#

# annot.sum.vect.small = as.vector(tapply(annot.sum.vect, my_ontology_df6$parent_structure_id[match(names(annot.sum.vect), my_ontology_df6$acronym)], sum))
aba.gene.annot.small = t(rowsum(t(aba.gene.annot), group = my_ontology_df6$parent_structure_id[match(colnames(aba.gene.annot), my_ontology_df6$acronym)], na.rm = T))
aba.gene.annot.small = aba.gene.annot.small[,which(!is.na(colnames(aba.gene.annot.small)))]

annot.sum.vect.small = my_ontology_df6$parent_sum[match(colnames(aba.gene.annot.small), as.character(my_ontology_df6$parent_structure_id))]
aba.gene.annot.small.pct = aba.gene.annot.small/annot.sum.vect.small[col(aba.gene.annot.small)]
aba.gene.annot.small.pct[which(is.na(aba.gene.annot.small.pct))] = 0
aba.gene.annot.small.pct = as.data.frame(aba.gene.annot.small.pct)
colnames(aba.gene.annot.small.pct) = ontology_df$acronym[match(colnames(aba.gene.annot.small.pct), ontology_df$id)]

# # No Log Transformation, Using Raw Counts Values
mz.homology.df = data.frame(mouse = rownames(aba.gene.annot), human = toupper(rownames(aba.gene.annot)), mz = gene_info$seurat_name[match(toupper(rownames(aba.gene.annot)), gene_info$one_to_one_human)])
Idents(all_merge) = all_merge$struct
# mz.exp = AverageExpression(all_merge, assay = "Spatial", features = mz.ch.names)[[1]]
# mz.aba.cor = cor(mz.exp, aba.gene.annot.pct)
# mz.aba.cor[which(is.na(mz.aba.cor))] = 0
# mz.aba.cor = mz.aba.cor[,which(colSums(mz.aba.cor) > 0)]
# pheatmap::pheatmap(mz.aba.cor, border_color = NA, scale = "row", cellheight = 5, cellwidth = 5, filename = "~/Downloads/spatial_aba_first_pass.pdf")

# aba.gene.annot.small.pct = aba.gene.annot / annot.sum.vect[col(aba.gene.annot)]
# aba.gene.annot.small.pct = aba.gene.annot.small.pct[,which(colnames(aba.gene.annot.small.pct) %in% my_ontology_df$acronym)]

# aba.gene.annot.small.pct.means = rowMeans(aba.gene.annot.small.pct)
mz.homology.df$isMzHVG = mz.homology.df$mz %in% all_merge@assays$SCT@var.features
mz.homology.df = mz.homology.df[which(mz.homology.df$isMzHVG),]
mz.exp = AverageExpression(all_merge, assay = "SCT", features = mz.homology.df$mz)[[1]]
mz.input = log(mz.exp+1)+0.1
mz.input = mz.input /  rowMeans(mz.input)
# mz.input = t(scale(t(mz.exp)))
mz.input[which(is.infinite(mz.input))] = 0
# aba.input = as.matrix(log(aba.gene.annot.small.pct[which(mz.ch.names.hvg),] / aba.gene.annot.small.pct.means[which(mz.ch.names.hvg)]))
aba.input = as.matrix(aba.gene.annot.small.pct[mz.homology.df$mouse,])
# aba.gene.annot.small.tmp = aba.gene.annot.small
# colnames(aba.gene.annot.small.tmp) = ontology_df$acronym[match(colnames(aba.gene.annot.small.tmp), ontology_df$id)]
# aba.input = as.matrix(aba.gene.annot.small.tmp[mz.homology.df$mouse,])
# aba.input = log(aba.input+1)+0.1
aba.input = aba.input /  rowMeans(aba.input)
# aba.input = t(scale(t(aba.input)))
aba.input[which(is.infinite(aba.input))] = 0
mz.aba.cor = cor(mz.input, aba.input, method = "spearman")
# mz.aba.cor[which(mz.aba.cor < 0)] = 0 # Colquitt Fig 3F: Negative correlations set to zero
mz.aba.cor[which(is.na(mz.aba.cor))] = 0
pheatmap::pheatmap(mz.aba.cor, color = viridis(100), border_color = NA, scale = "none", cellheight = 10, cellwidth = 10)
which(mz.aba.cor==max(mz.aba.cor[which(rownames(mz.aba.cor) != "Dc-4")]), arr.ind = T)

annot.sum.vect.big = annot.sum.vect[colnames(aba.gene.annot)]
aba.gene.annot.big.pct = aba.gene.annot / annot.sum.vect.big[col(aba.gene.annot)]
filtered_ontology_df = filter_mba_ontology_children(ontology_df, "CH")
filtered_ontology_df$acronym[which(grepl("-", filtered_ontology_df$acronym))] = str_replace(filtered_ontology_df$acronym[which(grepl("-", filtered_ontology_df$acronym))], "-", ".")
filtered_ontology_df$acronym[which(grepl("/", filtered_ontology_df$acronym))] = str_replace(filtered_ontology_df$acronym[which(grepl("/", filtered_ontology_df$acronym))], "/", ".")
aba.gene.annot.big.pct = aba.gene.annot.big.pct[,which(colnames(aba.gene.annot.big.pct) %in% filtered_ontology_df$acronym)]
aba.input.big = as.matrix(aba.gene.annot.big.pct[mz.homology.df$mouse,])
# aba.input.big = log(aba.input.big+1)+0.1
aba.input.big = aba.input.big /  rowMeans(aba.input.big)
aba.input.big[which(is.infinite(aba.input.big))] = 0
mz.aba.cor.big = cor(mz.input, aba.input.big, method = "spearman")

this.cluster = "Dc-3"
# high.res.annot = get_ccf_annotation()
annot.mat = reshape2::melt(high.res.annot[270,,]) # coronal (could try 285)
annot.mat$st_level = my_ontology_df$st_level[match(annot.mat$value, my_ontology_df$id)]
annot.mat$my.annot.id = NA
# annot.mat$my.annot.id[which(annot.mat$value %in% my_ontology_df6$id)] = my_ontology_df6$parent_structure_id[match(annot.mat$value[which(annot.mat$value %in% my_ontology_df6$id)], my_ontology_df6$id)]
annot.mat$my.annot.id[which(annot.mat$value %in% ontology_df$id[match(colnames(mz.aba.cor.big), ontology_df$acronym)] )] = annot.mat$value[which(annot.mat$value %in% ontology_df$id[match(colnames(mz.aba.cor.big), ontology_df$acronym)] )]
annot.mat$my.annot = ontology_df$acronym[match(annot.mat$my.annot.id, ontology_df$id)]
annot.mat$cor = NA
annot.mat$cor[which(annot.mat$st_level>0)] = 0
# annot.mat$cor[which(annot.mat$value > 10)] = 0
annot.mat$cor[which(!is.na(annot.mat$my.annot.id ))] = mz.aba.cor.big[this.cluster, annot.mat$my.annot[which(!is.na(annot.mat$my.annot.id ))]]
annot.mat$cor[which(annot.mat$cor < 0)] = 0
pdf("mz_aba_dc3_013023.pdf", width = 3, height = 2.5)
# ggplot(annot.mat, aes(x = Var2, y = -Var1, fill = cor)) + geom_raster() + coord_fixed() + scale_fill_gradientn(colors = c("lightgrey", "#0f326cff"), na.value = "white") + theme_void() # coronal
ggplot(annot.mat, aes(x = Var2, y = -Var1, fill = cor)) + geom_raster() + coord_fixed() + scale_fill_viridis(na.value = "white") + theme_void() # coronal
dev.off()

# Perms
# shuffleMZExp = function(x) {
#   Idents(all_merge) = sample(all_merge$struct)
#   this.mz.exp = AverageExpression(all_merge, assay = "SCT", features = mz.ch.names[which(mz.ch.names.hvg)])[[1]]
#   this.mz.exp.means = rowMeans(this.mz.exp)
#   this.mz.input = log(this.mz.exp / this.mz.exp.means)
#   this.mz.input[which(is.infinite(this.mz.input))] = 0
#   return(this.mz.input)
# }
permAbaCorTosches = function(old.mat) {
  new.mat.list = lapply(1:nrow(old.mat), function(x) sample(old.mat[x,]))
  new.mat = do.call('rbind', new.mat.list)
  perm.cor = cor(new.mat, aba.input, method = "spearman")
  perm.cor.melt = reshape2::melt(perm.cor)
  return(perm.cor.melt[,3])
}
spearmanCorAba = function(x) { return(cor.test(mz.input[,real.cor.melt$MZ[x]], aba.input[,real.cor.melt$ABA[x]], method = "spearman", alternative = "greater")$p.value) }

n.perms = 1000
perm.raw =  do.call('cbind', mclapply(1:n.perms, function(x) permAbaCorTosches(mz.input), mc.cores = 5))
real.cor = mz.aba.cor
mz.aba.cor.0 = mz.aba.cor
mz.aba.cor.0[which(mz.aba.cor.0 < 0)] = 0
real.cor.melt = reshape2::melt(real.cor)
colnames(real.cor.melt) = c("MZ", "ABA", "cor")
real.cor.melt$MZ = as.character(as.vector(real.cor.melt$MZ))
real.cor.melt$ABA = as.character(as.vector(real.cor.melt$ABA))
real.cor.melt$num.perm.greater = rowSums(perm.raw > real.cor.melt$cor)
real.cor.melt$p.perm = real.cor.melt$num.perm.greater / n.perms
real.cor.melt$bon.perm = p.adjust(real.cor.melt$p.perm, method = "bonferroni")
real.cor.melt$p.cor = unlist(mclapply(1:nrow(real.cor.melt), function(x) spearmanCorAba(x), mc.cores = 5))
real.cor.melt$bon.cor = p.adjust(real.cor.melt$p.cor, method = "bonferroni")
real.cor.melt$isSig = real.cor.melt$bon.perm < 0.05 & real.cor.melt$bon.cor < 0.05
# real.cor.melt$maxed = real.cor.melt$cor
# real.cor.melt$maxed[which(real.cor.melt$maxed >  0.2)] =  0.2
# real.cor.melt$maxed[which(real.cor.melt$maxed < -0.2)] = -0.2
real.cor.melt$maxed = reshape2::melt(mz.aba.cor.0)[,3]
real.cor.melt$maxed[which(real.cor.melt$maxed > 0.2)] = 0.2

mz.order.hclust  = hclust(dist(real.cor), method = "complete")
mz.order = mz.order.hclust$labels[mz.order.hclust$order]
aba.order.hclust  = hclust(dist(t(real.cor)), method = "complete")
aba.order = aba.order.hclust$labels[aba.order.hclust$order]
real.cor.melt$MZ = factor(real.cor.melt$MZ, levels=rev(mz.order))
real.cor.melt$ABA = factor(real.cor.melt$ABA, levels=aba.order)
pdf("mz_aba_w_sig_013023.pdf", width = 6.5, height = 4)
print(ggplot(real.cor.melt, aes(x = ABA, y = MZ, fill = maxed)) + geom_raster() + geom_point(data = real.cor.melt[which(real.cor.melt$isSig),], size = 0.5) + scale_fill_viridis() + coord_fixed() + scale_x_discrete(expand = c(0,0), name = NULL) + scale_y_discrete(expand = c(0,0), name = NULL) + theme(axis.text.x = element_text(angle=270, vjust=0.5, hjust=0)))
# print(ggplot(real.cor.melt, aes(x = ABA, y = MZ, fill = maxed)) + geom_raster() + geom_point(data = real.cor.melt[which(real.cor.melt$isSig),], size = 0.5) + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), limits = c(-0.2, 0.2)) + coord_fixed() + scale_x_discrete(expand = c(0,0), name = NULL) + scale_y_discrete(expand = c(0,0), name = NULL) + theme(axis.text.x = element_text(angle=270, vjust=0.5, hjust=0)))
dev.off()

# ggplot(reshape2::melt(as.matrix(aba.gene.annot[1:100])),        aes(x = Var1, y = Var2, fill = value)) + geom_raster() + scale_fill_viridis()
# ggplot(reshape2::melt(as.matrix(scale(aba.gene.annot[1:100]))), aes(x = Var1, y = Var2, fill = value)) + geom_raster() + scale_fill_viridis()

# ctx.aba = get_aba_ish_structure_data(mba_structure_id("CTX"))
ctx.aba = get_aba_ish_data(mba_structure_id("CTX"))
cnu.aba = get_aba_ish_structure_data(mba_structure_id("CNU"))

# ABA v2 =======================================================================
mz.homology.df = read.csv("~/research/st/results/aba_ish_mz_homology.csv")
# mz.homology.df = read.csv("~/research/st/results/aba_ish_bb_homology.csv")
# fetch_start_time <- proc.time()[[3]]
# fetchGeneFunc = function(this.gene) {
#   cat(".")
#   this.gene.ids <- get_gene_aba_ish_ids(this.gene, plane = "coronal")
#   return(get_aba_ish_structure_data(this.gene.ids[1]))
# }
# all.gene.ish = mclapply(mz.homology.df$mouse, function(x) fetchGeneFunc(x), mc.cores = 10)
# names(all.gene.ish) = mz.homology.df$mouse
# fetch_stop_time <- proc.time()[[3]]
# message(paste0("Fetching Genes Took: ", format(round(fetch_stop_time-fetch_start_time, 2), nsmall=2), " seconds"))
# saveRDS(all.gene.ish, "~/research/st/data/aba_ish_data3.rds")

all.gene.ish = readRDS("~/research/st/data/aba_ish_data3.rds")
# all.gene.ish = readRDS("~/research/st/data/aba_ish_data_bb.rds")
names(all.gene.ish) = mz.homology.df$mouse
ga <- get_ccf_grid_annotation()
ontology <- get_mba_ontology()
ontology_df <- flatten_mba_ontology(ontology)
filtered_ontology_df = filter_mba_ontology_children(ontology_df, "CH")
filtered_ontology_df$n_taxon = stringr::str_count(filtered_ontology_df$taxons, ";")
filtered_ontology_df$my_parent_id = filtered_ontology_df$id
filtered_ontology_df$my_parent_id[which(filtered_ontology_df$n_taxon > 5)] = unlist(lapply(which(filtered_ontology_df$n_taxon > 5), function(x) str_split(filtered_ontology_df$taxons[x], pattern = ';')[[1]][5] ))
filtered_ontology_df$my_parent_atlas_id = ontology_df$atlas_id[match(filtered_ontology_df$my_parent_id, as.character(ontology_df$id))]
# filtered_ontology_df$my_parent_atlas_id[which(filtered_ontology_df$n_children != 0)] = NA

oa.mat = matrix(0L, nrow = length(all.gene.ish), ncol = length(sort(unique(ontology_df$atlas_id))), dimnames = list(names(all.gene.ish), sort(unique(ontology_df$atlas_id))))
# for (i in 1:length(all.gene.ish)) { oa.mat[names(all.gene.ish)[i], all.gene.ish[[i]]$atlas_id] = all.gene.ish[[i]]$sum_pixel_intensity / all.gene.ish[[i]]$sum_pixels }
for (i in 1:length(all.gene.ish)) { oa.mat[names(all.gene.ish)[i], all.gene.ish[[i]]$atlas_id] = all.gene.ish[[i]]$voxel_energy_mean }
# for (i in 1:length(all.gene.ish)) { oa.mat[names(all.gene.ish)[i], all.gene.ish[[i]]$atlas_id] = all.gene.ish[[i]]$sum_pixel_intensity }
oa.mat.acr = oa.mat
colnames(oa.mat.acr) = ontology_df$acronym[match(colnames(oa.mat) , ontology_df$atlas_id)]
oa.mat.acr = oa.mat.acr[,which(colSums(oa.mat.acr) != 0)]

oa.mat.small = t(rowsum(t(oa.mat), group = filtered_ontology_df$my_parent_atlas_id[match(colnames(oa.mat), as.character(filtered_ontology_df$atlas_id))], na.rm = T))
oa.mat.small = oa.mat.small[,which(colSums(oa.mat.small) != 0 & !is.na(colnames(oa.mat.small)) )]
oa.mat.small.acr = oa.mat.small
colnames(oa.mat.small.acr) = ontology_df$acronym[match(colnames(oa.mat.small.acr), as.character(ontology_df$atlas_id))]

mz.exp = AverageExpression(all_merge, assay = "SCT", features = mz.homology.df$mz)[[1]]
# mz.exp = AverageExpression(bb, assay = "RNA", features = mz.homology.df$mz)[[1]]
mz.input = log(mz.exp+1)+0.1
mz.input = mz.input /  rowMeans(mz.input)
mz.input[which(is.infinite(mz.input) | is.na(mz.input))] = 0
# aba.input = oa.mat.acr
aba.input = oa.mat.small.acr
aba.input = log(aba.input+1)+0.1
aba.input = aba.input /  rowMeans(aba.input)
aba.input[which(is.infinite(aba.input) | is.na(aba.input))] = 0
mz.aba.cor = cor(mz.input, aba.input, method = "spearman")
# mz.aba.cor[which(mz.aba.cor < 0)] = 0 # Colquitt Fig 3F: Negative correlations set to zero
mz.aba.cor[which(is.na(mz.aba.cor))] = 0
pheatmap::pheatmap(mz.aba.cor, color = viridis(100), border_color = NA, scale = "none", cellheight = 10, cellwidth = 10)
mz.aba.cor.vis = reshape2::melt(mz.aba.cor)
mz.aba.cor.vis$scale.mz = reshape2::melt(t(scale(t(mz.aba.cor))))[,3]
ggplot(mz.aba.cor.vis, aes(x = value, y = scale.mz, color = value)) + geom_point(alpha = 0.5) + scale_color_viridis()

permAbaCorTosches = function(old.mat) {
  new.mat.list = lapply(1:nrow(old.mat), function(x) sample(old.mat[x,]))
  new.mat = do.call('rbind', new.mat.list)
  perm.cor = cor(new.mat, aba.input, method = "spearman")
  perm.cor.melt = reshape2::melt(perm.cor)
  return(perm.cor.melt[,3])
}
spearmanCorAba = function(x) { return(cor.test(mz.input[,real.cor.melt$MZ[x]], aba.input[,real.cor.melt$ABA[x]], method = "spearman", alternative = "greater")$p.value) }

n.perms = 1000
perm.raw =  do.call('cbind', mclapply(1:n.perms, function(x) permAbaCorTosches(mz.input), mc.cores = 5))
real.cor = mz.aba.cor
mz.aba.cor.0 = mz.aba.cor
mz.aba.cor.0[which(mz.aba.cor.0 < 0)] = 0
real.cor.melt = reshape2::melt(real.cor)
colnames(real.cor.melt) = c("MZ", "ABA", "cor")
real.cor.melt$MZ = as.character(as.vector(real.cor.melt$MZ))
real.cor.melt$ABA = as.character(as.vector(real.cor.melt$ABA))
real.cor.melt$num.perm.greater = rowSums(perm.raw > real.cor.melt$cor)
real.cor.melt$p.perm = real.cor.melt$num.perm.greater / n.perms
real.cor.melt$bon.perm = p.adjust(real.cor.melt$p.perm, method = "bonferroni")
real.cor.melt$p.cor = unlist(mclapply(1:nrow(real.cor.melt), function(x) spearmanCorAba(x), mc.cores = 5))
real.cor.melt$bon.cor = p.adjust(real.cor.melt$p.cor, method = "bonferroni")
real.cor.melt$isSig = real.cor.melt$bon.perm < 0.05 & real.cor.melt$bon.cor < 0.05
# real.cor.melt$maxed = real.cor.melt$cor
# real.cor.melt$maxed[which(real.cor.melt$maxed >  0.2)] =  0.2
# real.cor.melt$maxed[which(real.cor.melt$maxed < -0.2)] = -0.2
real.cor.melt$maxed = reshape2::melt(mz.aba.cor.0)[,3]
real.cor.melt$maxed[which(real.cor.melt$maxed > 0.2)] = 0.2

mz.order.hclust  = hclust(dist(real.cor), method = "complete")
mz.order = mz.order.hclust$labels[mz.order.hclust$order]
aba.order.hclust  = hclust(dist(t(real.cor)), method = "complete")
aba.order = aba.order.hclust$labels[aba.order.hclust$order]
real.cor.melt$MZ = factor(real.cor.melt$MZ, levels=rev(mz.order))
real.cor.melt$ABA = factor(real.cor.melt$ABA, levels=aba.order)
pdf("mz_aba_w_sig_020623.pdf", width = 6.5, height = 4)
print(ggplot(real.cor.melt, aes(x = ABA, y = MZ, fill = maxed)) + geom_raster() + geom_point(data = real.cor.melt[which(real.cor.melt$isSig),], size = 0.5) + scale_fill_viridis() + coord_fixed() + scale_x_discrete(expand = c(0,0), name = NULL) + scale_y_discrete(expand = c(0,0), name = NULL) + theme(axis.text.x = element_text(angle=270, vjust=0.5, hjust=0)))
# print(ggplot(real.cor.melt, aes(x = ABA, y = MZ, fill = maxed)) + geom_raster() + geom_point(data = real.cor.melt[which(real.cor.melt$isSig),], size = 0.5) + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), limits = c(-0.2, 0.2)) + coord_fixed() + scale_x_discrete(expand = c(0,0), name = NULL) + scale_y_discrete(expand = c(0,0), name = NULL) + theme(axis.text.x = element_text(angle=270, vjust=0.5, hjust=0)))
dev.off()

#*******************************************************************************
# Noise Zeisel =================================================================
#*******************************************************************************
zei = readRDS("~/scratch/bcs/data/l5_tel_norm.rds")
zei = subset(zei, cells = colnames(zei)[which(!zei$ClusterName %in% c("TEGLU8", "TEINH16", "RGDG"))])
zei.counts = zei@assays$RNA@counts
zei.num.to.remove = round(zei$nCount_RNA * 0.10)
# mclapply(1:ncol(zei.counts), function(x) { this.genes = rep(rownames(zei@assays$RNA@counts),table(zei@assays$RNA@counts[,1])); zei.num.to.remove[x] } )
# rm.idx.mat = matrix(0L, nrow = nrow(zei.counts), ncol = ncol(zei.counts), dimnames = list(rownames(zei.counts), colnames(zei.counts)))
# for (x in 1:ncol(zei.counts)) {
#   pos.idx = which(zei.counts[,x] > 0)
#   rm.idx = sample(pos.idx, zei.num.to.remove[x])
#   rm.idx.mat[rm.idx,x] = zei.counts[rm.idx,x]-1
# }
# mclapply(1:ncol(zei.counts), function(x) { pos.idx = which(zei.counts[,x] > 0); rm.idx = sample(pos.idx, zei.num.to.remove[x]); zei.counts[rm.idx,x] = zei.counts[rm.idx,x]-1 } )
idx.to.rm = mclapply(1:ncol(zei.counts), function(x) { pos.idx = which(zei.counts[,x] > 0); rm.idx = sample(pos.idx, zei.num.to.remove[x]); return(1:nrow(zei.counts) %in% rm.idx) }, mc.cores = 24 )
idx.rm.mat = do.call('cbind', idx.to.rm)
class(idx.rm.mat) = "numeric"
zei.counts.down = zei.counts - idx.rm.mat
down.zei = CreateSeuratObject(counts = zei.counts.down, meta.data = zei@meta.data)
saveRDS(down.zei, "~/scratch/bcs/data/zei_down_022823.rds")

mz.mouse.cor2 = mz.mouse.cor
for (i in rownames(mz.mouse.cor2)) {
  for (j in colnames(mz.mouse.cor2)) {
    mz.mouse.cor2[i, j] = i == j
  }
}
cor(as.vector(mz.mouse.cor), as.vector(mz.mouse.cor2))**2
mz.mouse.cor3 = mz.mouse.cor
mz.mouse.cor3[which(mz.mouse.cor3 < 0)] = 0
cor(as.vector(mz.mouse.cor3), as.vector(mz.mouse.cor2))**2

#*******************************************************************************
# Gene-Cell Heatmap ============================================================
#*******************************************************************************
# Load Libraries
message("Loading Libraries")
suppressMessages(source("~/scratch/st/st_scripts/st_f.R"))
suppressMessages(source("~/scratch/bcs/bcs_scripts/bcs_f.R"))
library("ggh4x")

# Mouse
mouse.path ="~/scratch/bcs/data/oritz_b_raw.rds"
mouse = readRDS(mouse.path)
Idents(mouse) = mouse$b_parent

# Cichlid
gene_info = read.csv("~/scratch/m_zebra_ref/gene_info_3.csv")
mz = readRDS("~/scratch/st/data/st_c2b2_hi_022023.rds") 
Idents(mz) = "structure"

gp = read.csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/st_oritzb_genes.csv")
gp$X = NULL
gp_long = data.frame()
for (i in seq(1, ncol(gp), by = 3)) { 
  this.gp = gp[,i:(i+2)]
  combo = colnames(this.gp)[1]
  # combo = stringr::str_replace_all(combo, "U_CTX","U_CTX")
  combo = stringr::str_replace_all(combo, "sAMY.other","sAMY-other")
  combo = stringr::str_replace_all(combo, "STRv.other","STRv-other")
  combo = stringr::str_replace_all(combo, "HIP.other","HIP-other")
  combo = stringr::str_replace_all(combo, "fiber.tracts","fiber tracts")
  mm_name = stringr::str_split(stringr::str_split(combo, "\\.")[[1]][1], "_")[[1]][2]
  mz_name = stringr::str_split(stringr::str_split(combo, "\\.")[[1]][2], "_")[[1]][2]
  if (TRUE) { mz_name = levels(mz)[as.numeric(mz_name)+1] }
  this.gp$mm_name = mm_name
  this.gp$mz_name = mz_name
  this.gp$combo = paste0(mm_name, "_", mz_name)
  this.gp = this.gp[which( !is.na(this.gp[,1]) & !is.na(this.gp[,2]) & !is.na(this.gp[,3]) ),]
  colnames(this.gp) = c("genes", "mm_pval", "mz_pval", "mm_name", "mz_name", "combo")
  this.gp[,c("mm_gene", "mz_gene")] = reshape2::colsplit(this.gp[,1], ";", c('1', '2')) 
  this.gp$mm_gene =  reshape2::colsplit(this.gp$mm_gene, "_", c('1', '2'))[,2]
  this.gp$mz_gene =  reshape2::colsplit(this.gp$mz_gene, "_", c('1', '2'))[,2]
  gp_long = rbind(gp_long, this.gp)
}
gp_long = gp_long[which(gp_long$combo %in% c("MOB_OB gc", "ACB_Vd-r", "LSX_Vv", 'DG_Dl-v', "Isocortex_Dl-g")),]
table(gp_long$combo)

# DEGs in only that cluster
gp_long_mz_table = as.data.frame(table(gp_long$mz_gene))
gp_long_mz_table = gp_long_mz_table[which(gp_long_mz_table$Freq == 1),]
gp_long_mm_table = as.data.frame(table(gp_long$mm_gene))
gp_long_mm_table = gp_long_mm_table[which(gp_long_mm_table$Freq == 1),]
gp_long_unique = gp_long[which(gp_long$mz_gene %in% gp_long_mz_table$Var1 & gp_long$mm_gene %in% gp_long_mm_table$Var1),]
table(gp_long_unique$combo)
gp_long_unique = gp_long_unique %>% group_by(combo) %>% slice(1:30)

# P-val threshold
gp_long$mm_neg_log = -log10(gp_long$mm_pval)
gp_long$mz_neg_log = -log10(gp_long$mz_pval)
gp_long_high = gp_long[which(gp_long$mm_neg_log > 10 & gp_long$mz_neg_log > 10),]
table(gp_long_high$combo)
gp_long_high = gp_long_high %>% group_by(combo) %>% slice(1:35)

mz.gp = subset(mz, idents = unique(gp_long_high$mz_name), features = unique(gp_long_high$mz_gene))
mz.gp = ScaleData(mz.gp, features = unique(gp_long_high$mz_gene))
# mz.gp.mat = mz.gp@assays$SCT@var.features

mouse.gp = subset(mouse, idents = unique(gp_long_high$mm_name), features = unique(gp_long_high$mm_gene))
mouse.gp = ScaleData(mouse.gp, features = unique(gp_long_high$mm_gene))

# mz.gp.fc = data.frame()
# for (i in unique(gp_long_high$mz_name)) { this.fc = FoldChange(mz, features = gp_long_high$mz_gene, ident.1 = i); this.fc$gene = gp_long_high$mz_gene; this.fc$cluster = i; mz.gp.fc = rbind(mz.gp.fc, this.fc); }
# mouse.gp.fc = data.frame()
# for (i in unique(gp_long_high$mm_name)) { this.fc = FoldChange(mouse, features = gp_long_high$mm_gene, ident.1 = i); this.fc$gene = gp_long_high$mm_gene; this.fc$cluster = i; mouse.gp.fc = rbind(mouse.gp.fc, this.fc); }
gp_long_high_tmp = gp_long[which(gp_long$mm_neg_log > 10 & gp_long$mz_neg_log > 10),]
mz.gp.fc.mat.small.tmp    = AverageExpression(mz,    features = gp_long_high_tmp$mz_gene, ident.1 = i)$SCT
mouse.gp.fc.mat.small.tmp = AverageExpression(mouse, features = gp_long_high_tmp$mm_gene, ident.1 = i)$SCT
mz.gp.fc.mat.tmp   = mz.gp.fc.mat.small.tmp[gp_long_high_tmp$mz_gene,   ]
mouse.gp.fc.mat.tmp = mouse.gp.fc.mat.small.tmp[gp_long_high_tmp$mm_gene,]
mz.gp.fc.mat.tmp    = as.matrix(t(scale(t(mz.gp.fc.mat.tmp))))
mouse.gp.fc.mat.tmp = as.matrix(t(scale(t(mouse.gp.fc.mat.tmp))))
gp_long_high_tmp$mz_top = unlist(lapply(1:nrow(gp_long_high_tmp), function(x) mz.gp.fc.mat.tmp[gp_long_high_tmp$mz_gene[x],    gp_long_high_tmp$mz_name[x]]))
gp_long_high_tmp$mm_top = unlist(lapply(1:nrow(gp_long_high_tmp), function(x) mouse.gp.fc.mat.tmp[gp_long_high_tmp$mm_gene[x], gp_long_high_tmp$mm_name[x]]))
gp_long_high_tmp$mean_top = rowMeans(gp_long_high_tmp[, c("mz_top", "mm_top")])
gp_long_high = data.frame(gp_long_high_tmp %>% group_by(combo) %>% arrange(-mean_top) %>% slice(1:35))

mz.gp.fc.mat.small    = AverageExpression(mz,    features = gp_long_high$mz_gene, ident.1 = i)$SCT
mouse.gp.fc.mat.small = AverageExpression(mouse, features = gp_long_high$mm_gene, ident.1 = i)$SCT
mz.gp.fc.mat    = mz.gp.fc.mat.small[gp_long_high$mz_gene,   ]
mouse.gp.fc.mat = mouse.gp.fc.mat.small[gp_long_high$mm_gene,]
mz.gp.fc.mat    = as.matrix(t(scale(t(mz.gp.fc.mat))))
mouse.gp.fc.mat = as.matrix(t(scale(t(mouse.gp.fc.mat))))
mz.gp.fc.mat    = mz.gp.fc.mat[,unique(gp_long_high$mz_name)]
mouse.gp.fc.mat = mouse.gp.fc.mat[,unique(gp_long_high$mm_name)]
mz.gp.fc = reshape2::melt(mz.gp.fc.mat)
mouse.gp.fc = reshape2::melt(mouse.gp.fc.mat)
colnames(mz.gp.fc) = colnames(mouse.gp.fc) = c("gene", "cluster", "avg_log2FC")

gp_long_high$gene_id = 1:nrow(gp_long_high)
mz.gp.fc$gene_id    = rep(gp_long_high$gene_id, length(unique(gp_long_high$mz_name)))
mouse.gp.fc$gene_id = rep(gp_long_high$gene_id, length(unique(gp_long_high$mm_name)))
mz.gp.fc$cluster    = factor(mz.gp.fc$cluster, levels = unique(gp_long_high$mz_name))
mouse.gp.fc$cluster = factor(mouse.gp.fc$cluster, levels = unique(gp_long_high$mm_name))
mz.maxed = 4
mz.gp.fc$maxed = mz.gp.fc$avg_log2FC
mz.gp.fc$maxed[which(mz.gp.fc$maxed >  mz.maxed)] =  mz.maxed
mz.gp.fc$maxed[which(mz.gp.fc$maxed < 0)] = 0
mouse.maxed = 4
mouse.gp.fc$maxed = mouse.gp.fc$avg_log2FC
mouse.gp.fc$maxed[which(mouse.gp.fc$maxed >  mouse.maxed)] =  mouse.maxed
mouse.gp.fc$maxed[which(mouse.gp.fc$maxed < 0)] = 0

# mz.gp.fc.mat = reshape2::dcast(mz.gp.fc, )
pdf('~/scratch/bcs/results/st_oritzb_gene_heatmap_avgexp.pdf', width = 12, height = 6)
p1 = ggplot(mz.gp.fc,    aes(x = gene_id, y = cluster, fill = maxed)) + geom_raster() + scale_fill_viridis() + theme_classic() + scale_y_discrete(name = "", expand=c(0,0)) + scale_x_discrete(name = "", expand=c(0,0)) + force_panelsizes(rows = unit(length(unique(gp_long_high$mz_name))/8, "in"), cols = unit(length(gp_long_high$mz_gene)/20, "in")) + theme(axis.line=element_blank())
p2 = ggplot(mouse.gp.fc, aes(x = gene_id, y = cluster, fill = maxed)) + geom_raster() + scale_fill_viridis() + theme_classic() + scale_y_discrete(name = "", expand=c(0,0)) + scale_x_discrete(name = "", expand=c(0,0)) + force_panelsizes(rows = unit(length(unique(gp_long_high$mm_name))/8, "in"), cols = unit(length(gp_long_high$mm_gene)/20, "in")) + theme(axis.line=element_blank())
cowplot::plot_grid(plotlist = list(p1, p2), ncol = 1)
dev.off()

Idents(mz.gp)    = factor(as.vector(mz.gp$structure),   levels = unique(gp_long_high$mz_name))
Idents(mouse.gp) = factor(as.vector(mouse.gp$b_parent), levels = unique(gp_long_high$mm_name))
pdf('~/scratch/bcs/results/st_oritzb_gene_heatmap.pdf', width = 16, height = 8)
p1 = DoHeatmap(mz.gp,    features = gp_long_high$mz_gene)
p2 = DoHeatmap(mouse.gp, features = gp_long_high$mm_gene)
cowplot::plot_grid(plotlist = list(p1, p2))
dev.off()

#*******************************************************************************
# SAMap Half Clusters ==========================================================
#*******************************************************************************
lapply(c("c2dr", "c2dl", "c2cr", "c2cl", "c2br", "c2bl", "b2ar", "b2al", "b2br", "b2cr", "b2dr", "b2bl", "b2cl", "b2dl"), sh_map)
sh_map = function(sh) {
  mzmm = as.matrix(read.csv(paste0("~/research/st/results/sh_samap/clusters/", sh, "_mapping_mine3.csv"), row.names = 1))
  colnames(mzmm) = str_sub(colnames(mzmm), 2, 50)
  mzmm[which(mzmm > quantile(mzmm, 0.99))] = quantile(mzmm, 0.99)
  mzmm.melt = reshape2::melt(mzmm)
  mzmm.melt = mzmm.melt[which(!is.na(mzmm.melt$value) & mzmm.melt$Var2 != ""),]
  colnames(mzmm.melt) = c("mm_name", "mz_name", "Score")
  mzmm.melt$id = paste0(mzmm.melt$mm_name, "_", mzmm.melt$mz_name)
  
  mzmm.p = as.matrix(read.csv(paste0("~/research/st/results/sh_samap/clusters/", sh, "_mapping_mine_p3.csv"), row.names = 1))
  colnames(mzmm.p) = str_sub(colnames(mzmm.p), 2, 50)
  mzmm.p.melt = reshape2::melt(mzmm.p)
  colnames(mzmm.p.melt) = c("mm_name", "mz_name", "p")
  mzmm.p.melt$id = paste0(mzmm.p.melt$mm_name, "_", mzmm.p.melt$mz_name)
  mzmm.melt$p_perm = mzmm.p.melt$p[match(mzmm.melt$id, mzmm.p.melt$id)]
  mzmm.melt$bh_perm = p.adjust(mzmm.melt$p_perm, method = "BH")
  mzmm.melt$bh_sig = mzmm.melt$bh_perm < 0.05
  mzmm.melt$p_sig  = mzmm.melt$p_perm < 0.05
  mzmm.melt$p0     = mzmm.melt$p_perm == 0
  
  mz.order  = hclust(dist(t(mzmm)), method = "complete")
  mz_order = reshape2::colsplit(mz.order$labels[mz.order$order], "\\.", c('1', '2'))[,1]
  mzmm.melt$mz_name = factor(mzmm.melt$mz_name, levels = mz_order)
  mouse.order = hclust(dist(mzmm), method = "complete")
  mzmm.melt$mm_name = factor(mzmm.melt$mm_name, levels = mouse.order$labels[mouse.order$order])
  
  pdf(paste0("~/research/st/results/sh_samap/", sh, ".pdf"), width = (nrow(mzmm)/5) + 2, height = (ncol(mzmm)/5) + 2)
  # ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = brewer.pal(9, "Greens"), breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
  # ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = rev(brewer.pal(11, "PiYG")[1:6]), breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
  # ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = rev(brewer.pal(11, "PiYG")[1:6]), breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="", labels=mzmm.melt$mm_num[match(levels(mzmm.melt$mm_name), mzmm.melt$mm_name)]) + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
  # ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = brewer.pal(9, "Blues"), breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
  # ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_viridis(breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$p_sig),], size = 0.6, color = "black") + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
  print(ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_viridis(breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$p_sig),], size = 0.6, color = "black") + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white"))
  dev.off()
}


#*******************************************************************************
# SAMap ========================================================================
#*******************************************************************************
library("jsonlite")
proteinToGene = function(x) {
  if ("transcripts" %in% names(x)) {
    if ("protein" %in% names(x$transcripts)) {
      if ("accessionVersion" %in% names(x$transcripts$protein)) {
        return(data.frame(protein = x$transcripts$protein$accessionVersion, gene_id = x$geneId, symbol = x$symbol))
      }
    }
  }
}
raw = readLines("~/scratch/m_zebra_ref/mouse_protein.jsonl")
json_raw = lapply(raw, fromJSON)
json_me = lapply(json_raw, proteinToGene)
json_df = do.call('rbind', json_me)
# gene_df = data.table::fread("~/scratch/m_zebra_ref/bird_genes.tsv", data.table = F)
# org = qs::qread(paste0(data_dir, "bird.qs"))
# gene_df$syb_in_sc = gene_df$Symbol %in% rownames(org@assays$RNA@counts)
# gene_df$syn_in_sc = gene_df$Synonyms %in% rownames(org@assays$RNA@counts)
# gene_df = gene_df[which( !(gene_df$syb_in_sc & gene_df$syn_in_sc) & (gene_df$syb_in_sc | gene_df$syn_in_sc) ),]
# gene_df$gene = gene_df$Symbol
# gene_df$gene[which(gene_df$syn_in_sc)] = gene_df$Synonyms[which(gene_df$syn_in_sc)]
# json_df$gene = gene_df$gene[match(json_df$gene_id, gene_df[,1])]
# write.csv(json_df, "~/scratch/m_zebra_ref/bird_gene_prot_table2.csv")
# length(which(unique(json_df$gene) %in% rownames(org@assays$RNA@counts)))

gtf = data.table::fread("~/scratch/m_zebra_ref/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gtf", data.table = F)
# gtf = data.table::fread("~/scratch/m_zebra_ref/GCF_000001635.27.gtf", data.table = F)
gtf$loc_id = paste0("LOC", reshape2::colsplit( reshape2::colsplit(gtf$V9, 'db_xref "GeneID:', c('1', '2'))[,2], '"', c('1', '2'))[,1])
gtf$gene_id = reshape2::colsplit( reshape2::colsplit(gtf$V9, '";', c('1', '2'))[,1], 'gene_id "', c('1', '2'))[,2]
gtf$protein_id = reshape2::colsplit( reshape2::colsplit(gtf$V9, 'protein_id "', c('1', '2'))[,2], '";', c('1', '2') )[,1]
gtf_small = unique(gtf[,c("loc_id", "gene_id", "protein_id")])
gtf_small = gtf_small[which(gtf_small$gene_id != "" & gtf_small$protein_id != ""),]
write.csv(gtf_small, "~/scratch/m_zebra_ref/turtle_prot_table3.csv")

gtf_small$gene2_id = stringr::str_replace_all(gtf_small$gene_id, "_", "-")
gtf_small$gene2_id_in_obj = gtf_small$gene2_id %in% rownames(org@assays$RNA@counts)
gtf_small$gene_id_in_obj = gtf_small$gene_id %in% rownames(org@assays$RNA@counts)
gtf_small$loc_id_in_obj  = gtf_small$loc_id %in% rownames(org@assays$RNA@counts)
gtf_small$gene = gtf_small$gene_id
gtf_small$gene[which(gtf_small$loc_id_in_obj)] = gtf_small$loc_id[which(gtf_small$loc_id_in_obj)]
gtf_small$gene[which(gtf_small$gene2_id_in_obj)] = gtf_small$gene2_id[which(gtf_small$gene2_id_in_obj)]

length(which(! unique(gtf_small$gene) %in% rownames(org@assays$RNA@counts) ))
length(which(! rownames(org@assays$RNA@counts) %in% unique(gtf_small$gene) ))
length(which(! unique(gtf_small$gene_id) %in% rownames(org@assays$RNA@counts) ))
length(which(! rownames(org@assays$RNA@counts) %in% unique(gtf_small$gene_id) ))

gene_df = data.table::fread("~/scratch/m_zebra_ref/mouse_genes.tsv", data.table = F)
gene_df$loc_id = paste0("LOC", gene_df[,1])
gtf_small[,c("Symbol", "Synonyms")] = gene_df[match(gtf_small$loc_id, gene_df$loc_id), c("Symbol", "Synonyms")]
gtf_small$gene[which(gtf_small$Symbol %in% rownames(org@assays$RNA@counts))] = gtf_small$Symbol[which(gtf_small$Symbol %in% rownames(org@assays$RNA@counts))]
for (i in which(!gtf_small$gene %in% rownames(org@assays$RNA@counts) & !is.na(gtf_small$Synonyms))) {
  syns = gtf_small$Synonyms[i]
  syns = stringr::str_split(syns, ",")[[1]]
  syns_logical = syns %in% rownames(org@assays$RNA@counts)
  if (any(syns_logical)) {
    gtf_small$gene[i] = syns[which(syns_logical)[1]]
  }
}

gene_info = read.csv("~/scratch/m_zebra_ref/gene_info_3.csv")
json_df$loc = paste0("LOC", json_df$gene_id)
json_df$seurat_name = gene_info$seurat_name[match(json_df$loc, gene_info$loc)]
write.csv(json_df, "~/scratch/m_zebra_ref/mz_gene_prot_table2.csv")

# conda activate SeuratDisk
.libPaths("/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/conda_envs/SeuratDisk/lib/R/library")
library("Seurat")
library("SeuratObject")
library("SeuratDisk")
library("ggplot2")
Convert("~/scratch/bcs/data/mz_mm_zei.h5ad", dest = "h5seurat", overwrite = TRUE)
merged = LoadH5Seurat("~/scratch/bcs/data/mz_mm_zei.h5seurat", meta.data = FALSE, misc = FALSE)
# library(rhdf5)
# merged[["mised_meta_value"]] <- h5read("~/scratch/bcs/data/mz_mm_zei_keys.h5ad", "/obs/mised_meta_value")
merged = readRDS("~/scratch/brain/data/mz_mm_zei.rds")
mz = readRDS("~/scratch/brain/data/bb_demux_102021.rds")
mm = readRDS("~/scratch/bcs/data/l5_tel_norm.rds")
# mm = readRDS("~/scratch/bcs/data/mouse_w_pc_down_norm.rds")
# mm = readRDS("~/scratch/bcs/data/oritz_b_raw.rds")
convert53 = read.csv("~/scratch/st/data/convert53.csv")
merged$mz_seuratclusters53 = mz$seuratclusters53[match(colnames(merged), colnames(mz))]
merged$mz_cluster = factor(convert53$new[match(merged$mz_seuratclusters53, convert53$old)], levels = convert53$new)
merged$mm_cluster = mm$ClusterName[match(colnames(merged), colnames(mm))]
# merged$mm_cluster = mm$tissue_cluster[match(colnames(merged), colnames(mm))]
# merged$mm_cluster = mm$b_parent[match(colnames(merged), colnames(mm))]
merged$cluster = as.vector(merged$mz_cluster)
merged$cluster[which(is.na(merged$cluster))] = merged$mm_cluster[which(is.na(merged$cluster))]
merged$species = plyr::revalue(as.character(colnames(merged) %in% colnames(mz)), replace = c("TRUE" = "mz", "FALSE" = "mm"))

# Plot the merged data
Idents(merged) = merged$species
pdf("~/scratch/bcs/results/mz_mm_zei_umap.pdf", width = 5, height = 5)
print(DimPlot(merged) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
dev.off()
Idents(merged) = merged$cluster
pdf("~/scratch/bcs/results/mz_mm_zei_split.pdf", width = 10, height = 5)
print(DimPlot(merged, split.by = "species", label = F) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend())
dev.off()
pdf("~/scratch/bcs/results/mz_mm_zei_split_label.pdf", width = 20, height = 5)
print(DimPlot(merged, split.by = "species", label = T) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
dev.off()

# pdf("~/scratch/bcs/results/mz_mm_zei_umap_black.pdf", width = 4, height = 4)
# print(DimPlot(merged, cols = c("black", RColorBrewer::brewer.pal(9, "Greens")[8])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
# dev.off()
# pdf("~/scratch/bcs/results/mz_mm_zei_split_black.pdf", width = 8, height = 4)
# print(DimPlot(merged, split.by = "species", label = F, cols = c("black", RColorBrewer::brewer.pal(9, "Greens")[8])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend())
# dev.off()

Cairo::Cairo("~/scratch/bcs/results/mz_mm_zei_umap_black.png", width = 1600, height = 1600, res = 300)
print(DimPlot(merged, cols = c("black", RColorBrewer::brewer.pal(9, "Greens")[8])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
dev.off()
Cairo::Cairo("~/scratch/bcs/results/mz_mm_zei_split_black.png", width = 3200, height = 1600, res = 300)
print(DimPlot(merged, split.by = "species", label = F, cols = c("black", RColorBrewer::brewer.pal(9, "Greens")[8])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend())
dev.off()

Cairo::Cairo("~/scratch/bcs/results/mz_mm_saunders_umap_black.png", width = 1600, height = 1600, res = 300)
print(DimPlot(merged, raster = F, cols = c("black", RColorBrewer::brewer.pal(11, "PiYG")[2])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
dev.off()
Cairo::Cairo("~/scratch/bcs/results/mz_mm_saunders_split_black.png", width = 3200, height = 1600, res = 300)
print(DimPlot(merged, split.by = "species", raster = F, label = F, cols = c("black", RColorBrewer::brewer.pal(11, "PiYG")[2])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend())
dev.off()

Cairo::Cairo("~/scratch/bcs/results/mz_mm_oritzb_umap_black.png", width = 1600, height = 1600, res = 300)
print(DimPlot(merged, raster = F, cols = c("black", RColorBrewer::brewer.pal(9, "Blues")[8])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
dev.off()
Cairo::Cairo("~/scratch/bcs/results/mz_mm_oritzb_split_black.png", width = 3200, height = 1600, res = 300)
print(DimPlot(merged, split.by = "species", label = F, raster = F, cols = c("black", RColorBrewer::brewer.pal(9, "Blues")[8])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend())
dev.off()

# pdf("~/scratch/bcs/results/mz_mm_saunders_umap_cols.pdf", width = 4, height = 4)
# print(DimPlot(merged, cols = c("#FBCBC7", RColorBrewer::brewer.pal(11, "PiYG")[2])) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()))
# dev.off()
# pdf("~/scratch/bcs/results/mz_mm_saunders_split_cols.pdf", width = 4, height = 8)
# print(DimPlot(merged, split.by = "species", label = F, cols = c("#FBCBC7", RColorBrewer::brewer.pal(11, "PiYG")[2]), ncol = 1) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + NoLegend())
# dev.off()

# Celltype-Celltype mapping
mzmm = read.csv("~/Downloads/mz_mm_oritz_mapping_table.csv")
mzmm = mzmm[54:nrow(mzmm),1:54]
# mzmm = mzmm[27:nrow(mzmm),1:27]
# colnames(mzmm) = levels(all_merge$struct)[match(as.numeric(reshape2::colsplit(colnames(mzmm), "_", c('1','2'))[,2])+1, 1:26)]
colnames(mzmm) = convert53$new[match(as.numeric(reshape2::colsplit(colnames(mzmm), "_", c('1','2'))[,2]), convert53$old)]
mzmm[,1] = reshape2::colsplit(mzmm[,1], "_", c('1', '2'))[,2]
mzmm.melt = reshape2::melt(mzmm)
mzmm.melt = mzmm.melt[which(!is.na(mzmm.melt$value)),]
colnames(mzmm.melt) = c("mm_name", "mz_name", "Score")
# mzmm.melt$mz_name = convert53$new[match(as.numeric(reshape2::colsplit(mzmm.melt$variable, "_", c('1','2'))[,2]), convert53$old)] 
# mzmm.melt$mz_name = factor(mzmm.melt$mz_name, levels = convert53$new)
# mzmm.melt$mm_name = reshape2::colsplit(mzmm.melt$X, "_", c('1', '2'))[,2]
mzmm.mat = mzmm; rownames(mzmm.mat) = mzmm.mat[,1]; mzmm.mat[,1] = NULL;
mz.order  = hclust(dist(t(mzmm.mat)), method = "complete")
mzmm.melt$mz_name = factor(mzmm.melt$mz_name, levels = mz.order$labels[mz.order$order])
mouse.order = hclust(dist(mzmm.mat), method = "complete")
mzmm.melt$mm_name = factor(mzmm.melt$mm_name, levels = mouse.order$labels[mouse.order$order])

pdf("~/research/st/results/mz_mm_oritz_mapping_keys.pdf", width = 4, height = 4.5)
# pdf("~/research/st/results/mz_mm_oritz_mapping.pdf", width = 4, height = 8)
# pdf("~/research/st/results/mz_mm_zei_broad_mapping.pdf", width = 4.5, height = 12)
# pdf("~/research/st/results/mz_mm_zei_mapping.pdf", width = 10, height = 10)
# ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = value)) + geom_raster() + scale_fill_gradientn(colors = c("#ffffff", rev(brewer.pal(11, "RdBu"))[7:11])) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("") + coord_fixed()
ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_viridis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("") + coord_fixed()
dev.off()

# My mapping
mzmm = as.matrix(read.csv("~/Downloads/bb_bird_mapping_mine3.csv", row.names = 1))
# colnames(mzmm) = str_replace_all(colnames(mzmm), "\\.", "-")
# colnames(mzmm) = plyr::revalue(colnames(mzmm), c("Dc-1-2" = "Dc-1/2"))
colnames(mzmm) = str_sub(colnames(mzmm), 2, 50)
colnames(mzmm) = str_replace(colnames(mzmm), "Astro", "RG")
colnames(mzmm) = plyr::revalue(colnames(mzmm), c("8.9_Glut"="8-9_Glut", "8.9_Glut.1"="8.9_Glut", "8.9_Glut.1"="8.9_Glut", "15.1_GABA.Glut"="15.1_GABA/Glut", "15.5_GABA.Glut"="15.5_GABA/Glut"))
# mzmm = scale(mzmm) # scale by cichlid cluster
mzmm[which(mzmm > quantile(mzmm, 0.99))] = quantile(mzmm, 0.99)
# mzmm[which(mzmm > 0.5)] = 0.5
mzmm.melt = reshape2::melt(mzmm)
mzmm.melt = mzmm.melt[which(!is.na(mzmm.melt$value) & mzmm.melt$Var2 != ""),]
colnames(mzmm.melt) = c("mm_name", "mz_name", "Score")
mzmm.melt$id = paste0(mzmm.melt$mm_name, "_", mzmm.melt$mz_name)

mzmm.p = as.matrix(read.csv("~/Downloads/bb_bird_mapping_mine_p3.csv", row.names = 1))
# colnames(mzmm.p) = str_replace_all(colnames(mzmm.p), "\\.", "-")
# colnames(mzmm.p) = plyr::revalue(colnames(mzmm.p), c("Dc-1-2" = "Dc-1/2"))
colnames(mzmm.p) = str_sub(colnames(mzmm.p), 2, 50)
colnames(mzmm.p) = str_replace(colnames(mzmm.p), "Astro", "RG")
colnames(mzmm.p) = plyr::revalue(colnames(mzmm.p), c("8.9_Glut"="8-9_Glut", "8.9_Glut.1"="8.9_Glut", "8.9_Glut.1"="8.9_Glut", "15.1_GABA.Glut"="15.1_GABA/Glut", "15.5_GABA.Glut"="15.5_GABA/Glut"))
mzmm.p.melt = reshape2::melt(mzmm.p)
colnames(mzmm.p.melt) = c("mm_name", "mz_name", "p")
mzmm.p.melt$id = paste0(mzmm.p.melt$mm_name, "_", mzmm.p.melt$mz_name)
mzmm.melt$p_perm = mzmm.p.melt$p[match(mzmm.melt$id, mzmm.p.melt$id)]
# mzmm.melt$bh_perm = p.adjust(mzmm.melt$p_perm, method = "BH")
mzmm.melt$bh_perm = p.adjust(mzmm.melt$p_perm, method = "BH")
mzmm.melt$bh_sig = mzmm.melt$bh_perm < 0.05
mzmm.melt$p_sig  = mzmm.melt$p_perm < 0.05
mzmm.melt$p0     = mzmm.melt$p_perm == 0
# mzmm.melt$sig    = mzmm.melt$bh_perm < 0.05
# mzmm.melt$sig    = mzmm.melt$p_perm == 0

mz.order  = hclust(dist(t(mzmm)), method = "complete")
mzmm.melt$mz_name = factor(mzmm.melt$mz_name, levels = mz.order$labels[mz.order$order])
mouse.order = hclust(dist(mzmm), method = "complete")
mzmm.melt$mm_name = factor(mzmm.melt$mm_name, levels = mouse.order$labels[mouse.order$order])
# mzmm.melt$mm_num = reshape2::colsplit(mzmm.melt$mm_name, "_", c('1', '2'))[,2]

pdf("~/research/st/results/bb_bird_mine3.pdf", width = (nrow(mzmm)/5) + 2, height = (ncol(mzmm)/5) + 2)
# ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = brewer.pal(9, "Greens"), breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
# ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = rev(brewer.pal(11, "PiYG")[1:6]), breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
# ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = rev(brewer.pal(11, "PiYG")[1:6]), breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="", labels=mzmm.melt$mm_num[match(levels(mzmm.melt$mm_name), mzmm.melt$mm_name)]) + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
# ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = brewer.pal(9, "Blues"), breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
# ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_viridis(breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$p_sig),], size = 0.6, color = "black") + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_viridis(breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) + geom_point(data = mzmm.melt[which(mzmm.melt$p_sig),], size = 0.6, color = "black") + geom_point(data = mzmm.melt[which(mzmm.melt$bh_sig),], size = 1.2, color = "gray") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.2, color = "white")
dev.off()

maxed.num = 1e-4
mzmm.melt$maxed = mzmm.melt$Score
mzmm.melt$maxed[which(mzmm.melt$maxed > maxed.num)]  =  maxed.num
mzmm.melt$maxed[which(mzmm.melt$maxed < -maxed.num)] = -maxed.num
pdf("~/research/st/results/mz_mm_saunders_mapping_minemax.pdf", width = 9, height = 9)
ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = maxed)) + geom_raster() + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), n.breaks = 6, limits = c(-maxed.num, maxed.num)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("") + coord_fixed()
dev.off()

# SAMap Mine vs Cors
cor.res = read.csv("~/Downloads/mz_saunders_cor.csv")
cor.res.mat = reshape2::acast(cor.res, mouse.cluster ~ mz.cluster, value.var = "cor")
# cor.res.mat = scale(cor.res.mat)
cor(as.vector(cor.res.mat), as.vector(mzmm[rownames(cor.res.mat), colnames(cor.res.mat)]))
cor(as.vector(cor.res.mat), as.vector(mzmm[rownames(cor.res.mat), colnames(cor.res.mat)]), method = "spearman")
cor.res.mat.melt = reshape2::melt(cor.res.mat)
cor.res.mat.melt$Var1 = factor(cor.res.mat.melt$Var1, levels = mouse.order$labels[mouse.order$order])
cor.res.mat.melt$Var2 = factor(cor.res.mat.melt$Var2, levels = mz.order$labels[mz.order$order])
maxed.num = 0.35
cor.res.mat.melt$maxed = cor.res.mat.melt$value
cor.res.mat.melt$maxed[which(cor.res.mat.melt$maxed > maxed.num)]   =  maxed.num
cor.res.mat.melt$maxed[which(cor.res.mat.melt$maxed < -maxed.num)]  = -maxed.num
pdf("~/research/st/results/mz_mm_saunders_cor_w_samap_hclust.pdf", width = 9, height = 9)
# ggplot(cor.res.mat.melt, aes(x = Var1, y = Var2, fill = maxed)) + geom_raster() + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), n.breaks = 6, limits = c(-maxed.num, maxed.num)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("") + coord_fixed()
ggplot(cor.res.mat.melt, aes(x = Var1, y = Var2, fill = value)) + geom_raster() + scale_fill_viridis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("") + coord_fixed()
dev.off()
mzmm.melt$cor = cor.res.mat.melt$value[match(paste0(mzmm.melt$mm_name, "_", mzmm.melt$mz_name), paste0(cor.res.mat.melt$Var1, "_", cor.res.mat.melt$Var2))]
ggplot(mzmm.melt, aes(x = Score, y = cor)) + geom_point(alpha = 0.5) 

cor.res.glut = read.csv("~/Downloads/bb_saunders_subcluster_glut_cor2.csv")
cor.res.gaba = read.csv("~/Downloads/bb_saunders_subcluster_gaba_cor.csv")
cor.res.nn   = read.csv("~/Downloads/bb_saunders_subcluster_nn_cor.csv")
# cor.res.stich = (cor.res.mat == "hi") * 1 # trying to make a matrix of 0's
mouse_clusters = unique(c(cor.res.glut$mouse.cluster, cor.res.gaba$mouse.cluster, cor.res.nn$mouse.cluster))
mz_clusters = unique(c(cor.res.glut$mz.cluster, cor.res.gaba$mz.cluster, cor.res.nn$mz.cluster))
cor.res.stich = reshape2::melt(matrix(0L, nrow = length(mouse_clusters), ncol = length(mz_clusters), dimnames = list(mouse_clusters, mz_clusters)))
cor.res.stich$value = 0
cor.res.stich$id = paste0(cor.res.stich$Var1, "_", cor.res.stich$Var2)
cor.res.glut$id  = paste0(cor.res.glut$mouse.cluster,  "_", cor.res.glut$mz.cluster)
cor.res.gaba$id  = paste0(cor.res.gaba$mouse.cluster,  "_", cor.res.gaba$mz.cluster)
cor.res.nn$id    = paste0(cor.res.nn$mouse.cluster,    "_", cor.res.nn$mz.cluster)
cor.res.stich$value[which(cor.res.stich$id %in% cor.res.glut$id)] = cor.res.glut$cor[match(cor.res.stich$id[which(cor.res.stich$id %in% cor.res.glut$id)], cor.res.glut$id)]
cor.res.stich$value[which(cor.res.stich$id %in% cor.res.gaba$id)] = cor.res.gaba$cor[match(cor.res.stich$id[which(cor.res.stich$id %in% cor.res.gaba$id)], cor.res.gaba$id)]
cor.res.stich$value[which(cor.res.stich$id %in% cor.res.nn$id)]   = cor.res.nn$cor[match(cor.res.stich$id[which(cor.res.stich$id   %in% cor.res.nn$id)],   cor.res.nn$id)]
maxed.num = 0.2
cor.res.stich$maxed = cor.res.stich$value
cor.res.stich$maxed[which(cor.res.stich$maxedd >  maxed.num)] =  maxed.num
cor.res.stich$maxed[which(cor.res.stich$maxedd < -maxed.num)] = -maxed.num
ggplot(cor.res.stich, aes(x = Var1, y = Var2, fill = maxed)) + geom_raster() + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), n.breaks = 6, limits = c(-maxed.num, maxed.num)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("") + coord_fixed()
cor.res.stich$samap = mzmm.melt$Score[match(cor.res.stich$id, mzmm.melt$id)]
cor(cor.res.stich$value[which(!is.na(cor.res.stich$samap) )], cor.res.stich$samap[which(!is.na(cor.res.stich$samap) )])

# Celltype-Celltype mapping gene pairs
gp = read.csv("~/Downloads/mz_mm_zei_broad_gene_pair.csv")
gp$X = NULL
gp.long = data.frame()
for (this.col in colnames(gp)[c(T,F,F)]) {
  print(which(colnames(gp) == this.col))
  this.split = stringr::str_split(this.col, "_")[[1]]
  this.mm.cluster =  stringr::str_split(this.split[3], "\\.")[[1]]
  this.mm.cluster = paste(this.mm.cluster[1:(length(this.mm.cluster)-1)], collapse = " ")
  this.mm.cluster = paste( c(this.split[2], "_", this.mm.cluster), collapse = "" )
  this.mz.cluster = convert53$new[match(as.numeric(this.split[4]), convert53$old)]
  for (this.row in 1:nrow(gp)) {
    this.genes = gp[this.row,this.col]
    if (this.genes == "" || is.na(this.genes)) { break; }
    this.genes = stringr::str_split(this.genes, "_")[[1]]
    mm.gene = stringr::str_split(this.genes[2], ";")[[1]][1]
    mz.gene = this.genes[3]
    this.p1 = gp[this.row, which(colnames(gp) == this.col)+1]
    this.p2 = gp[this.row, which(colnames(gp) == this.col)+2]
    this.df = data.frame(mz_cluster = this.mz.cluster, mm_cluster = this.mm.cluster, mz_gene = mz.gene, mm_gene = mm.gene, p1 = this.p1, p2 = this.p2)
    gp.long = rbind(gp.long, this.df)
  }
}
gp.long$neg_log_p1 = -log10(gp.long$p1)
gp.long$neg_log_p2 = -log10(gp.long$p2)
gp.long$neg_log_p1[which(gp.long$neg_log_p1 > 500)] = 500
gp.long$neg_log_p2[which(gp.long$neg_log_p2 > 500)] = 500
num.genes = as.data.frame(table(gp.long$mz_cluster, gp.long$mm_cluster))
ggplot(num.genes, aes(x = Var2, y = Var1, fill = Freq)) + geom_raster() + xlab("") + ylab("") + scale_fill_viridis() + coord_fixed() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf("~/research/st/results/test.pdf", width = 20, height = 20)
ggplot(gp.long, aes(x = neg_log_p1, y = neg_log_p2)) + geom_point() + facet_grid(vars(mz_cluster), vars(mm_cluster))
dev.off()


# SAMap vs Cors
# cor.res = read.csv("~/Downloads/mz_zei_cor_cluster.csv")
# samap.res = read.csv("~/Downloads/mz_mm_zei_mapping_table.csv")
cor.res = read.csv("~/Downloads/mz_saunders_cor.csv")
samap.res = read.csv("~/Downloads/mz_mm_saunders_cluster_mapping_table.csv")
cor.res.mat = reshape2::acast(cor.res, mouse.cluster ~ mz.cluster, value.var = "cor")
# cor.res.mat[which(cor.res.mat < 0)] = 0
samap.res = samap.res[54:nrow(samap.res),1:54]
rownames(samap.res) = reshape2::colsplit(samap.res[,1], "_", c('1', '2'))[,2]
colnames(samap.res) = reshape2::colsplit(colnames(samap.res), "_", c('1', '2'))[,2]
samap.res[,1] = NULL
colnames(samap.res) = convert53$new[match(as.numeric(colnames(samap.res)), convert53$old)]
samap.res = as.matrix(samap.res[,colnames(cor.res.mat)])
cor(as.vector(cor.res.mat),  as.vector(samap.res) )
cor(as.vector(cor.res.mat),  as.vector(samap.res), method = "spearman")
# cor(as.vector(cor.res.mat^2),  as.vector(samap.res^(1/2)) )
# cor(as.vector(cor.res.mat^2),  as.vector(samap.res^(1/2)) )^2

top.x = 5
cor.combos.top   = unlist(lapply(colnames(cor.res.mat), function(x) paste(rownames(cor.res.mat)[order(cor.res.mat[,x], decreasing = T)[1:top.x]], x) ) )
samap.combos.top = unlist(lapply(colnames(samap.res),   function(x) paste(rownames(samap.res)[order(samap.res[,x],     decreasing = T)[1:top.x]],     x) ) )
length(which(cor.combos.top %in% samap.combos.top)) / length(unique(c(cor.combos.top, samap.combos.top)))

cor.combos.top   = unlist(lapply(rownames(cor.res.mat), function(x) paste(colnames(cor.res.mat)[order(cor.res.mat[x,], decreasing = T)[1:top.x]], x) ) )
samap.combos.top = unlist(lapply(rownames(samap.res),   function(x) paste(colnames(samap.res)[order(samap.res[x,],     decreasing = T)[1:top.x]],     x) ) )
length(which(cor.combos.top %in% samap.combos.top)) / length(unique(c(cor.combos.top, samap.combos.top)))

pheatmap::pheatmap(samap.res^(1/2), color=viridis(100), border_color = NA, cluster_rows = F, cluster_cols = F, cellwidth = 6, cellheight = 6, show_rownames = F, show_colnames = F)
# pheatmap::pheatmap(cor.res.mat^2, color=viridis(100), border_color = NA, cluster_rows = F, cluster_cols = F)

cor.res.mat2 = t(scale(t( cor.res.mat^5 )))
cor.res.mat2[which(cor.res.mat2 < 0)] = 0
pheatmap::pheatmap(cor.res.mat2, color=viridis(100), border_color = NA, cluster_rows = F, cluster_cols = F, cellwidth = 6, cellheight = 6, show_rownames = F, show_colnames = F)
cor(as.vector(cor.res.mat2),  as.vector(samap.res^(1/2)) )
# cowplot::plot_grid(plist = list(p1, p2)) 

mzmm = read.csv("~/scratch/bcs/results/mz_mm_zei_broad_mapping_table.csv")
mzmm = mzmm[54:nrow(mzmm),1:54]
colnames(mzmm) = convert53$new[match(as.numeric(reshape2::colsplit(colnames(mzmm), "_", c('1','2'))[,2]), convert53$old)]
mzmm[,1] = reshape2::colsplit(mzmm[,1], "_", c('1', '2'))[,2]
rownames(mzmm) = mzmm[,1]; mzmm[,1] = NULL;
mz_mouse_cluster = data.frame(mz = colnames(mzmm), mm = unlist(lapply(1:ncol(mzmm), function(x) rownames(mzmm)[which.max(mzmm[,x])])))
mz_mouse_cluster$value = unlist(lapply(1:ncol(mzmm), function(x) max(mzmm[,x])))
mz_mouse_cluster = mz_mouse_cluster[which(mz_mouse_cluster$value >= 0.2),]

mz.min.size = min(table(mz$seuratclusters53))
mz.cells = unlist(lapply(mz_mouse_cluster$mz, function(x) sample(colnames(mz)[which(mz$good_names53 == x)], mz.min.size) ))

mm$broad = paste0(mm$Region, "_", mm$TaxonomyRank4)
mm.min.size = min(table(mm$broad))
mm.cells = unlist(lapply(unique(mz_mouse_cluster$mm), function(x) sample(colnames(mm)[which(mm$broad == x)], mm.min.size) ))

dp.mz = subset(mz, cells = mz.cells)
dp.mm = subset(mm, cells = mm.cells)

# Testing making my own SAMap score
xsim = data.table::fread("~/scratch/bcs/results/test.csv", data.table = F, header = T, drop = 1)
meta = data.table::fread("~/scratch/bcs/results/test2.csv", data.table = F, header = T)
xsim.comp = xsim[which(meta$species == "mz"), which(meta$species == "mm")]
# xsim.comp.small   = t(rowsum(t(xsim.comp), group = meta$mm_ABA_parent[which(meta$species == "mm")], na.rm = T))
# xsim.comp.smaller = t(rowsum(xsim.comp.small, group = meta$mz_struct_b2_vdc[which(meta$species == "mz")], na.rm = T))
all_mm_clusters = sort(unique(meta$mm_ABA_parent[which(meta$species == "mm")]))
all_mz_clusters = sort(unique(meta$mz_struct_b2_vdc[which(meta$species == "mz")]))
xsim.comp.smaller = matrix(0L, nrow = length(all_mm_clusters), ncol = length(all_mz_clusters), dimnames = list(all_mm_clusters, all_mz_clusters))
for (mm.cluster in all_mm_clusters) {
  for (mz.cluster in all_mz_clusters) {
    xsim.comp.smaller[mm.cluster, mz.cluster] = mean(as.matrix( xsim[which(meta$mm_ABA_parent == mm.cluster),which(meta$mz_struct_b2_vdc == mz.cluster)] ))
  }
}
pheatmap::pheatmap(xsim.comp.smaller, color = viridis(100), cellwidth = 10, cellheight = 10, filename = "~/scratch/bcs/results/test.pdf")
xsim.comp.smaller2 = matrix(0L, nrow = length(all_mz_clusters), ncol = length(all_mm_clusters), dimnames = list(all_mz_clusters, all_mm_clusters))
for (mz.cluster in all_mz_clusters) {
  for (mm.cluster in all_mm_clusters) {
    xsim.comp.smaller2[mz.cluster, mm.cluster] = mean(as.matrix( xsim[which(meta$mz_struct_b2_vdc == mz.cluster),which(meta$mm_ABA_parent == mm.cluster)] ))
  }
}
pheatmap::pheatmap(xsim.comp.smaller, color = viridis(100), cellwidth = 10, cellheight = 10, filename = "~/scratch/bcs/results/test2.pdf")

mz.df = data.frame(table(meta$mz_struct_b2_vdc[which(meta$mz_struct_b2_vdc != "unassigned")]))
mm.df = data.frame(table(meta$mm_ABA_parent[which(meta$mm_ABA_parent != "unassigned")]))
xsim.comp.smaller.norm = xsim.comp.smaller / mm.df[,2]
xsim.comp.smaller.norm = xsim.comp.smaller.norm / mz.df[col(xsim.comp.smaller.norm),2]
pheatmap::pheatmap(xsim.comp.smaller.norm, color = viridis(100), cellwidth = 10, cellheight = 10, filename = "~/scratch/bcs/results/test.pdf")

# Single cell vs spatial results
library(ggh4x)
mzmm.st = as.matrix(read.csv("~/Downloads/st_oritz_b_mapping_mine.csv", row.names = 1))
colnames(mzmm.st) = str_replace(colnames(mzmm.st), "\\.", "-")
mzmm.sc = as.matrix(read.csv("~/Downloads/bb_oritz_b_mapping_mine.csv", row.names = 1))
colnames(mzmm.sc) = str_sub(colnames(mzmm.sc), 2, 50)
colnames(mzmm.sc) = str_replace(colnames(mzmm.sc), "Astro", "RG")
colnames(mzmm.sc) = plyr::revalue(colnames(mzmm.sc), c("8.9_Glut"="8-9_Glut", "8.9_Glut.1"="8.9_Glut", "8.9_Glut.1"="8.9_Glut", "15.1_GABA.Glut"="15.1_GABA/Glut", "15.5_GABA.Glut"="15.5_GABA/Glut"))
mzmm.st = scale(mzmm.st)
mzmm.sc = scale(mzmm.sc)
mzmm.sc = mzmm.sc[,colnames(c2l_mean)]
spot.mm = as.data.frame(as.matrix(c2l_mean) %*% t(mzmm.sc))
# st.mm = matrix(0L, nrow = ncol(mzmm.st), ncol = ncol(spot.mm), dimnames = list(colnames(mzmm.st), colnames(spot.mm)))
spot.mm$struct = all_merge$struct
spot.mm = spot.mm %>% group_by(struct) %>% dplyr::summarise(across(everything(), mean, na.rm = TRUE))
spot.mm.mat = as.matrix(spot.mm); rownames(spot.mm.mat) = spot.mm.mat[,1]; spot.mm.mat = spot.mm.mat[,-1];
# class(spot.mm.mat) = "numeric"; spot.mm.mat = scale(spot.mm.mat); spot.mm = spot.mm.mat;
spot.mm.melt = reshape2::melt(spot.mm, id.var = "struct")
# colnames(spot.mm.melt) = c("struct", "variable", "value")
mz.order  = hclust(dist(spot.mm.mat), method = "complete")
spot.mm.melt$struct = factor(spot.mm.melt$struct, levels = mz.order$labels[mz.order$order])
mouse.order = hclust(dist(t(spot.mm.mat)), method = "complete")
spot.mm.melt$variable = factor(spot.mm.melt$variable, levels = mouse.order$labels[mouse.order$order])
maxed.num = 30
spot.mm.melt$maxed = spot.mm.melt$value
spot.mm.melt$maxed[which(spot.mm.melt$maxed >  maxed.num)] =  maxed.num
spot.mm.melt$maxed[which(spot.mm.melt$maxed < -maxed.num)] = -maxed.num
pdf("~/research/st/results/bb_to_oritz_to_spatial_max.pdf", width = 6, height = 4.5)
ggplot(spot.mm.melt, aes(x = variable, y = struct, fill = maxed)) + geom_raster() + scale_fill_viridis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("") + coord_fixed() + force_panelsizes(rows = unit(nrow(spot.mm.mat)/8, "in"), cols = unit(ncol(spot.mm.mat)/8, "in"))
dev.off()


#*******************************************************************************
# Trash Can ====================================================================
#*******************************************************************************

# Figure 2 B
pdf(paste0( "~/research/st/results/c2b2_eomesa.pdf"), width = 2.5, height = 2.5)
print(FeaturePlot(all_merge, "LOC101480282", order = T, cols=c("lightgrey", "#0f9783ff")) + theme_void() + NoLegend() + coord_fixed() + ggtitle(""))
dev.off()
pdf(paste0( "~/research/st/results/c2b2_dlx5.pdf"), width = 2.5, height = 2.5)
print(FeaturePlot(all_merge, "dlx5", order = T, cols=c("lightgrey", "#ED256E")) + theme_void() + NoLegend() + coord_fixed() + ggtitle(""))
dev.off()
pdf(paste0( "~/research/st/results/c2b2_slc17a6.pdf"), width = 2.5, height = 2.5)
print(FeaturePlot(all_merge, "slc17a6", order = T, cols=c("lightgrey", "#FDE725FF")) + theme_void() + NoLegend() + coord_fixed() + ggtitle(""))
dev.off()
pdf(paste0( "~/research/st/results/c2b2_gad1.pdf"), width = 2.5, height = 2.5)
print(FeaturePlot(all_merge, "gad1", order = T, cols=c("lightgrey", "#440154FF")) + theme_void() + NoLegend() + coord_fixed() + ggtitle(""))
dev.off()

# Figure Cross Species Drivers
pdf("~/research/st/results/paintings-cross-species/c2b2_zic1.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "zic1", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_prdm16.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "prdm16", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_penk_LOC101474395.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "LOC101474395", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_six3.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "six3", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_sp8_LOC101468530.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "LOC101468530", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_sall1.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "sall1", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_gpsm1.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "gpsm1", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_lamp5.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "lamp5", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_cacnb3.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "cacnb3", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_bhlhe22.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "bhlhe22", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()
pdf("~/research/st/results/paintings-cross-species/c2b2_neurod2.pdf", width = 14, height = 2, onefile = F)
print(myC2B2SFPFew(all_merge, "neurod2", pt.size.multiplier = 0.9, pal = colorRampPalette(c("#DCDCDC", brewer.pal(9, "Reds")[2:8])), rm.zero = F, scale.alpha = F, same.col.scale = F, showLegend = F))
dev.off()

# Figure 2 C
pdf(paste0(out_dir, "c2b2_npy_dark.pdf"), width = 14, height = 2, onefile = F)
this.col = "#E80975"
print(myC2B2SFPFew(all_merge, "npy", pt.size.multiplier = 1.1, pal = colorRampPalette(c(darken(this.col, amount = 0.2), this.col)), rm.zero = T, scale.alpha = T, same.col.scale = F, showLegend = F))
dev.off()
pdf(paste0(out_dir, "c2b2_npy_light.pdf"), width = 14, height = 2, onefile = F)
this.col = "#E80975"
print(myC2B2SFPFew(all_merge, "npy", pt.size.multiplier = 1.1, pal = colorRampPalette(c(lighten(this.col, amount = 0.4), this.col)), rm.zero = T, scale.alpha = T, same.col.scale = F, showLegend = F))
dev.off()
pdf(paste0(out_dir, "c2b2_sst_light.pdf"), width = 14, height = 2, onefile = F)
this.col = "#DD19E0"
print(myC2B2SFPFew(all_merge, "LOC101485677", pt.size.multiplier = 1.1, pal = colorRampPalette(c(lighten(this.col, amount = 0.4), this.col)), rm.zero = T, scale.alpha = T, same.col.scale = F, showLegend = F))
dev.off()
pdf(paste0(out_dir, "c2b2_th_light.pdf"), width = 14, height = 2, onefile = F)
this.col = "#BD4AFF"
print(myC2B2SFPFew(all_merge, "th", pt.size.multiplier = 1.1, pal = colorRampPalette(c(lighten(this.col, amount = 0.05), this.col)), rm.zero = T, scale.alpha = T, same.col.scale = F, showLegend = F))
dev.off()
pdf(paste0(out_dir, "c2b2_penk_light.pdf"), width = 14, height = 2, onefile = F)
this.col = "#653FE0"
print(myC2B2SFPFew(all_merge, "LOC101474395", pt.size.multiplier = 1.1, pal = colorRampPalette(c(lighten(this.col, amount = 0.4), this.col)), rm.zero = T, scale.alpha = T, same.col.scale = F, showLegend = F))
dev.off()
pdf(paste0(out_dir, "c2b2_pvalb_light.pdf"), width = 14, height = 2, onefile = F)
this.col = "#29A2AB"
this.col = "#43d0e6"
print(myC2B2SFPFew(all_merge, "LOC101487165", pt.size.multiplier = 1.1, pal = colorRampPalette(c(lighten(this.col, amount = 0.4), this.col)), rm.zero = T, scale.alpha = T, same.col.scale = F, showLegend = F))
dev.off()
pdf(paste0(out_dir, "c2b2_unc_light.pdf"), width = 14, height = 2, onefile = F)
this.col = "#22E38F"
print(myC2B2SFPFew(all_merge, "ucn", pt.size.multiplier = 1.1, pal = colorRampPalette(c(lighten(this.col, amount = 0.4), this.col)), rm.zero = T, scale.alpha = T, same.col.scale = F, showLegend = F))
dev.off()

fake.stsc.meta = data.frame()
for (i in 1:nrow(c2l_mean_melt)) {
  this.num = c2l_mean_melt$value[i]
  if (this.num > 0) {
    fake.stsc.meta=rbind(fake.stsc.meta, data.frame(cell="cell1", spot = rep(c2l_mean_melt$spot[i],this.num), ct = c2l_mean_melt$variable[i] ))
  }
}
fake.stsc.meta$sh = my.b2$sh[match(fake.stsc.meta$spot, colnames(my.b2))]
pdf(paste0(out_dir, "c2b2_mini_cells.pdf"), width = 25, height = 3, onefile = F)
print(myC2B2SFPFew( my.b2, "ct", stsc.list = list(fake.stsc.mat, fake.stsc.meta), pal = scales::hue_pal()(length(levels(all_merge$ct))), showLegend = T, pt.size.multiplier = 1.1 ))
dev.off()

neg_log_bon_thresh = 200
avg_logFC_thresh = 3
bri.deg = read.csv("~/research/st/results/st_c2b2_standard_cluster_markers_new_nums_011223.csv")
bri.deg.genes = read.csv("~/Downloads/b2c2_fig1_markers_011323.csv")
bri.deg$label[which(bri.deg$gene == "LOC101479570")] = "sox6"
bri.deg$label[which(bri.deg$gene == "LOC101486618")] = "her4.2"
bri.deg$label[which(bri.deg$gene == "LOC101464413")] = "kif5c"
bri.deg$label[which(bri.deg$gene == "LOC101486743")] = "homer3"
bri.deg$label[which(bri.deg$gene == "LOC101465831")] = "sema3f"
bri.deg.genes$seurat_name[which(bri.deg.genes$seurat_name == "crhbo")] = "crhbp"
bri.deg$cluster_gene = paste0(bri.deg$cluster, "_", bri.deg$gene)
bri.deg.genes$cluster_gene = paste0(bri.deg.genes$cluster, "_", bri.deg.genes$seurat_name)
bri.deg$neg_log_bon = -log10(bri.deg$p_val_adj)
bri.deg$neg_log_bon[which(bri.deg$neg_log_bon > neg_log_bon_thresh)] = neg_log_bon_thresh
bri.deg$avg_logFC[which(bri.deg$avg_logFC > avg_logFC_thresh)] = avg_logFC_thresh
bri.deg = bri.deg[which(bri.deg$cluster_gene %in% bri.deg.genes$cluster_gene),]
pdf("~/research/st/results/st_c2b2_dotplotlike.pdf", width = 15, height = 3)
# ggplot(top.deg, aes(x = gene, y = neg_log_bon, color = avg_logFC, size = pct.1)) + geom_point() + scale_color_viridis(limits = c(0, 3)) + scale_size_continuous(limits = c(0, 1), range = c(0.5, 4)) + facet_wrap(~ cluster, scales = "free_x",  ncol = 33) + theme_light() + theme(panel.spacing = unit(0, "lines"), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("")
ggplot(bri.deg, aes(x = label, y = neg_log_bon, color = avg_logFC, size = pct.1)) + geom_point() + scale_color_viridis(limits = c(0, 3)) + scale_size_continuous(limits = c(0, 1), range = c(0.5, 4)) + facet_wrap(~ cluster, scales = "free_x",  ncol = 33) + theme_light() + theme(strip.text = element_text(size=15), panel.spacing = unit(0, "lines"), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("")
dev.off()

feature.df = my.b2@meta.data[, c("sample", "nFeature_Spatial"),]
pdf("~/research/st/results/c2b2_nFeature_vln.pdf", width = 5.75, height = 2)
ggplot(feature.df, aes(x = sample, y = nFeature_Spatial, color = sample, fill = sample)) + geom_violin(alpha = 0.1) + stat_summary(fun = "mean", fun.min = min, fun.max = max, geom="pointrange", stroke = 0.1, size = 0.7) + theme_classic() + scale_y_continuous(expand = c(0,0)) + NoLegend() + ylab("")
dev.off()
count.df = my.b2@meta.data[, c("sample", "nCount_Spatial"),]
pdf("~/research/st/results/c2b2_nCount_vln.pdf", width = 5.75, height = 2)
ggplot(count.df, aes(x = sample, y = nCount_Spatial, color = sample, fill = sample)) + geom_violin(alpha = 0.1) + stat_summary(fun = "mean", fun.min = min, fun.max = max, geom="pointrange", stroke = 0.1, size = 0.7) + theme_classic() + scale_y_continuous(expand = c(0,0)) + NoLegend() + ylab("")
dev.off()

bri_tmp = read.csv("~/Downloads/295_A1_E2.csv")
all_merge_hi$bri_annot = "other"
all_merge_hi$bri_annot[paste0("c2a_", bri_tmp$Barcode)] = bri_tmp$Structure
Idents(all_merge_hi) = all_merge_hi$bri_annot
bri_deg = FindAllMarkers(all_merge_hi)
bri_deg_sig = bri_deg[which(bri_deg$p_val_adj < 0.05 & bri_deg$cluster != "other"),]
bri_deg_sig$pct.dif = bri_deg_sig$pct.1 - bri_deg_sig$pct.2

this.ident = "E"
pdf(paste0(out_dir, "bri_", this.ident, ".pdf"), width = 12, height = 12, onefile = F)
this.deg = bri_deg_sig[which(bri_deg_sig$pct.dif > 0.5 & bri_deg_sig$avg_log2FC > 0.25 & bri_deg_sig$cluster == this.ident), ]
print(nrow(this.deg))
all_merge_hi@meta.data[,paste0(this.ident, ".score.", nrow(this.deg))] = colSums(all_merge_hi@assays$Spatial@counts[this.deg$gene,])
myMultiSFP(all_merge_hi, feature = paste0(this.ident, ".score.", nrow(this.deg)), pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(viridis(100)), high.res = T, rm.zero = F, scale.alpha = F)
dev.off()

# Old Attempt to Label Brain Structures
label_spatial = function() {
  plot(this.df2)
  selectedPoints <- gatepoints::fhs(this.df2)
  return(as.numeric(as.vector(selectedPoints)))
}
# b2a_labels = data.frame(spot = colnames())
this.sample = "b2a"
this.df = GetTissueCoordinates(subset(all_merge, cells = colnames(all_merge)[which(all_merge$sample == this.sample)]))
this.df$label = ""
this.df2 = data.frame(x = this.df$imagecol, y = -this.df$imagerow)
this.idx = label_spatial()
this.df$label[which(this.df$label == "" & 1:nrow(this.df) %in% this.idx)] = "Dd"
# this.df$label[this.idx] = "Dl-v"
ggplot(this.df, aes(x = imagecol, y = -imagerow, color = label)) + geom_point(size = 3) + theme_bw()

# Plot Cell2location results
all_merge_hi = qs::qread(paste0(data_dir, "all_merge_hi.qs"))
means = read.csv(paste0(out_dir, "cell2location/bb_secondary/cell2location_spatial_output_means.csv"))
rownames(means) = means$X
means$X = NULL
colnames(means) = str_replace(colnames(means), "meanscell_abundance_w_sf_", "")
means_round = round(means)
all_merge_hi$sum = rowSums(means)

for (i in as.character(0:52)) {
  print(paste0("---- ", i, " ----"))
  this.name = bb_convert53$new[match(i, bb_convert53$old)]
  this.name = str_replace(this.name, "/", "_")
  all_merge_hi$tmp = means_round[,i]
  all_merge$tmp = all_merge_hi$tmp[colnames(all_merge)]
  # grDevices::svg(paste0(out_dir, "cell2location/bb_secondary/rounded_cell/", this.name, ".svg"), width = 7, height = 7, onefile = F)
  # myMultiSFP(all_merge_hi, feature = "tmp", pt.size.multiplier = 1, pal = colorRampPalette(viridis(100)), high.res = T)
  # dev.off()
  
  grDevices::svg(paste0(out_dir, "cell2location/bb_secondary/rounded_cell/", this.name, ".svg"), width = 16, height = 2, onefile = F)
  myB2SFP(all_merge, feature = "tmp", pt.size.multiplier = 1.3, pal = colorRampPalette(viridis(100)))
  dev.off()
  
  # pct_sample = unlist(lapply(levels(all_merge_hi$sample), function(x) range01(all_merge_hi$tmp[which(all_merge_hi$sample ==x)])*100))
  # names(pct_sample) = unlist(lapply(levels(all_merge_hi$sample), function(x) colnames(all_merge_hi)[which(all_merge_hi$sample ==x)]))
  # all_merge_hi$tmp2 = 0
  # all_merge_hi$tmp2[names(pct_sample)] = pct_sample
  # grDevices::svg(paste0(out_dir, "cell2location/bb_secondary/rounded_cell_max/", this.name, ".svg"), width = 7, height = 7, onefile = F)
  # print(myMultiSFP(all_merge_hi, feature = "tmp2", pt.size.multiplier = 1, pal = colorRampPalette(viridis(100)), high.res = T))
  # dev.off()
  # 
  # all_merge_hi$tmp = (means_round[,i] / all_merge_hi$sum) * 100
  # grDevices::svg(paste0(out_dir, "cell2location/bb_secondary/rounded_pct/", this.name, ".svg"), width = 7, height = 7, onefile = F)
  # myMultiSFP(all_merge_hi, feature = "tmp", pt.size.multiplier = 1, pal = colorRampPalette(viridis(100)), high.res = T)
  # dev.off()
  # 
  # pct_sample = unlist(lapply(levels(all_merge_hi$sample), function(x) range01(all_merge_hi$tmp[which(all_merge_hi$sample ==x)])*100))
  # names(pct_sample) = unlist(lapply(levels(all_merge_hi$sample), function(x) colnames(all_merge_hi)[which(all_merge_hi$sample ==x)]))
  # all_merge_hi$tmp2 = 0
  # all_merge_hi$tmp2[names(pct_sample)] = pct_sample
  # grDevices::svg(paste0(out_dir, "cell2location/bb_secondary/rounded_pct_max/", this.name, ".svg"), width = 7, height = 7, onefile = F)
  # myMultiSFP(all_merge_hi, feature = "tmp2", pt.size.multiplier = 1, pal = colorRampPalette(viridis(100)), high.res = T)
  # dev.off()
}

for (i in 0:52) {
  print(paste0("---- ", i, " ----"))
  this.name = convert53$new[match(i, convert53$old)]
  this.name.clean = str_replace(this.name, "/", "_")
  all_merge$tmp = c2l_mean[colnames(all_merge),this.name]
  
  # grDevices::svg(paste0(out_dir, "cell2location/bb_secondary/rounded_cell/", this.name.clean, ".svg"), width = 14, height = 4, onefile = F)
  # myC2B2SFPAll(all_merge, feature = "tmp", pt.size.multiplier = 1.3, pal = colorRampPalette(viridis(100)))
  # dev.off()
  
  Cairo::Cairo(paste0(out_dir, "cell2location/bb_secondary/rounded_cell/", this.name.clean, ".png"), width = 1400, height = 400, onefile = F)
  myC2B2SFPAll(all_merge, feature = "tmp", pt.size.multiplier = 1.3, pal = colorRampPalette(viridis(100)))
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

pdf(paste0(out_dir, "npy_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#5cc456"
myMultiSFP(all_merge_hi, feature = "npy", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "LOC101480282_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#a03fba"
myMultiSFP(all_merge_hi, feature = "LOC101480282", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "ucn_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#a30e07"
myMultiSFP(all_merge_hi, feature = "ucn", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "LOC101482567_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#7656c4"
myMultiSFP(all_merge_hi, feature = "LOC101482567", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "nkx2-1_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#0423bf"
myMultiSFP(all_merge_hi, feature = "nkx2-1", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "dlx5_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#56c49c"
myMultiSFP(all_merge_hi, feature = "dlx5", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "gad1_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#a661c2"
myMultiSFP(all_merge_hi, feature = "gad1", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "slc17a6_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#5cc456"
myMultiSFP(all_merge_hi, feature = "slc17a6", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "LOC101475168_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#f2de41"
myMultiSFP(all_merge_hi, feature = "LOC101475168", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "LOC101468574_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#23b5de"
myMultiSFP(all_merge_hi, feature = "LOC101468574", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "LOC101487165_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#a03fba"
myMultiSFP(all_merge_hi, feature = "LOC101487165", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "LOC106675461_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#db7a39"
myMultiSFP(all_merge_hi, feature = "LOC106675461", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "LOC101485677_dif_remove_zero_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#288f62"
myMultiSFP(all_merge_hi, feature = "LOC101485677", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
pdf(paste0(out_dir, "th_dif_remove_zero_alpha_alpha.pdf"), width = 12, height = 12, onefile = F)
this.col = "#dbb43d"
myMultiSFP(all_merge_hi, feature = "th", pt.size.multiplier = 1.6, same.col.scale = F, pal = colorRampPalette(c(darken(this.col, amount = 0.1), this.col)), high.res = T, rm.zero = F, scale.alpha = T)
dev.off()
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

# ABA **************************************************************************
ish_slice_heatmap <- function(mat,
                              anno = NULL,
                              slice_num,
                              plane = "coronal",
                              normalize = "slice",
                              colorset = c("darkblue","gray90","red")) {
  
  library(rbokeh)
  library(dplyr)
  library(reshape2)
  library(scrattch.vis)
  
  slice_mat <- slice_ccf_arr(mat, slice_num, plane)
  
  if(normalize == "slice") {
    max_val <- max(slice_mat)
  } else if(normalize == "all") {
    max_val <- max(mat)
  }
  
  slice_flat <- reshape2::melt(slice_mat)
  names(slice_flat) <- c("y","x","value")
  
  slice_flat$value[slice_flat$value < 0] <- 0
  slice_flat$color <- scrattch.vis::values_to_colors(slice_flat$value,
                                                     colorset = colorset,
                                                     min_val = 0,
                                                     max_val = max_val)
  
  if(is.null(anno)) {
    hover_list <- list("Value" = "value")
  } else {
    anno_flat <- reshape2::melt(slice_ccf_arr(anno, slice_num, plane))
    names(anno_flat) <- c("y","x","annotation")
    slice_flat <- dplyr::left_join(slice_flat, anno_flat, by = c("x","y"))
    hover_list <- list("Value" = "value",
                       "Annotation" = "annotation")
  }
  
  if(plane == "coronal") {
    f <- rbokeh::figure(width = dim(slice_mat)[2]*10,
                        height = dim(slice_mat)[1]*10) %>%
      rbokeh::ly_crect(data = slice_flat,
                       x = x,
                       y = -y,
                       fill_color = color,
                       line_color = NA,
                       fill_alpha = 1,
                       hover = hover_list)
  } else if(plane == "horizontal") {
    f <- rbokeh::figure(width = dim(slice_mat)[1]*10,
                        height = dim(slice_mat)[2]*10) %>%
      rbokeh::ly_crect(data = slice_flat,
                       x = y,
                       y = x,
                       fill_color = color,
                       line_color = NA,
                       fill_alpha = 1,
                       hover = hover_list)
  } else if(plane == "saggital") {
    f <- rbokeh::figure(width = dim(slice_mat)[1]*10,
                        height = dim(slice_mat)[2]*10) %>%
      rbokeh::ly_crect(data = slice_flat,
                       x = y,
                       y = -x,
                       fill_color = color,
                       line_color = NA,
                       fill_alpha = 1,
                       hover = hover_list)
  }
  
  return(f)
}
Slc17a7_ids <- get_gene_aba_ish_ids("Slc17a7", plane = "coronal")
Pvalb_ids <- get_gene_aba_ish_ids("Pvalb", plane = "coronal")
Sst_ids <- get_gene_aba_ish_ids("Sst", plane = "coronal")

# all_ids <- c(Slc17a7_ids[1], Pvalb_ids[1], Sst_ids[1])
all_ids <- c(Slc17a7_ids[2], Pvalb_ids[2])
all_data <- map(all_ids, get_aba_ish_data)

ish_slice_heatmap(all_data[[1]],
                  anno = oa,
                  plane = "coronal",
                  colorset = viridis(10),
                  slice = 42)

this.mat = reshape2::melt(all_data[[2]][,,42]) # saggital
ggplot(this.mat, aes(x = Var1, y = -Var2, fill = value)) + geom_raster() + scale_fill_viridis() + coord_fixed() # saggital
this.mat = reshape2::melt(all_data[[2]][42,,]) # coronal
ggplot(this.mat, aes(x = Var2, y = -Var1, fill = value)) + geom_raster() + scale_fill_viridis() + coord_fixed() # coronal

st.gene.num.spots = rowSums(all_merge@assays$Spatial@counts > 0)
st.gene.num.spots = names(st.gene.num.spots)[which(st.gene.num.spots > 0.01*ncol(all_merge))]
ch.names = stringr::str_to_title( unique(gene_info$human[match(st.gene.num.spots, gene_info$seurat_name)]) )
ch.names = ch.names[which(!is.na(ch.names) & ch.names != "")]
ch.names = sort(ch.names)
all.gene.ish = list()
fetch_start_time <- proc.time()[[3]]
for (this.gene in ch.names) {
  if (which(ch.names == this.gene) %% 100 == 0) { print(which(ch.names == this.gene) %% 100) }
  this.gene.ids <- get_gene_aba_ish_ids(this.gene, plane = "coronal")
  if (length(this.gene.ids) > 0) { all.gene.ish[[this.gene]] = get_aba_ish_data(this.gene.ids[1], values = "energy")  }
}
fetch_stop_time <- proc.time()[[3]]
message(paste0("Fetching Genes Took: ", format(round(fetch_stop_time-fetch_start_time, 2), nsmall=2), " seconds"))

# Method 1: Find Expression of DEG Cluster Markers in ABA
cluster.annot = data.frame(matrix(0L, nrow = length(all.annot), ncol = length(unique(deg$cluster))))
rownames(cluster.annot) = all.annot
for (this.clust in unique(deg$cluster)) {
  # this.clust = 2
  deg.clust.hgnc = deg$hgnc[which(deg$cluster == this.clust & -log10(deg$p_val_adj) > 10)[1:5]]
  deg.clust.hgnc = sort(unique( stringr::str_to_title(deg.clust.hgnc[which(deg.clust.hgnc != "" & !is.na(deg.clust.hgnc))]) ))
  deg.clust.hgnc = deg.clust.hgnc[which(deg.clust.hgnc %in% names(all.gene.ish))]
  mean.s = list()
  sum.df = annot.df = data.frame(matrix(0L, nrow = length(all.annot), ncol = 58))
  rownames(sum.df) = rownames(annot.df) = all.annot
  colnames(sum.df) = colnames(annot.df) = paste0("slice", 1:58)
  for (this.slice in 1:58) {
    this.annot = as.vector(slice_ccf_arr(oa, this.slice, "coronal"))
    this.annot = factor(this.annot, levels = all.annot)
    annot.df[,this.slice] = as.vector(table(this.annot))*length(deg.clust.hgnc)
    # mean.s[[this.slice]] = matrix(0, nrow = 41, ncol = 58)
    # for (this.gene in deg.clust.hgnc) { this.gene.exp = all.gene.ish[[this.gene]][this.slice,,]; this.gene.exp[which(this.gene.exp < 0)] = 0; mean.s[[this.slice]] = mean.s[[this.slice]] + this.gene.exp }
    # mean.s[[this.slice]] = mean.s[[this.slice]] / length(deg.clust.hgnc)
    # this.mat = reshape2::melt(mean.s[[this.slice]]) # coronal
    # print(ggplot(this.mat, aes(x = Var2, y = -Var1, fill = value)) + geom_raster() + scale_fill_viridis() + coord_fixed() + ggtitle(paste0("Slice = ", this.slice, ", Markers From Cluster ", this.clust))) # coronal
    # Sys.sleep(1)
    for (this.gene in deg.clust.hgnc) { 
      this.gene.exp = all.gene.ish[[this.gene]][this.slice,,]
      this.gene.exp = as.vector(all.gene.ish[[this.gene]][this.slice,,])
      this.gene.exp[which(this.gene.exp < 0)] = 0
      this.sum = tapply(this.gene.exp, this.annot, sum)
      this.sum[which(is.na(this.sum))] = 0
      sum.df[,this.slice]   = sum.df[,this.slice]   + this.sum
      # annot.df[,this.slice] = annot.df[,this.slice] + as.vector(table(this.annot))
    }
  }
  annot.mean = data.frame(value = rowSums(sum.df) / rowSums(annot.df), annot = rownames(annot.df))
  cluster.annot[,(this.clust+1)] = annot.mean$value
  annot.mean = annot.mean[order(annot.mean$value, decreasing = T),]
  annot.mean = annot.mean[1:20,]
  annot.mean$annot = factor(annot.mean$annot, levels = unique(tmp$annot))
  print(ggplot(annot.mean, aes(x = annot, y = value)) + geom_bar(stat = "identity"))
}
ggplot(reshape2::melt(as.matrix(cluster.annot)),        aes(x = Var1, y = Var2, fill = value)) + geom_raster() + scale_fill_viridis()
ggplot(reshape2::melt(as.matrix(scale(cluster.annot))), aes(x = Var1, y = Var2, fill = value)) + geom_raster() + scale_fill_viridis()

# My method
# mz.deg$hgnc = gene_info$human[match(mz.deg$gene, gene_info$seurat_name)]
# mouse.deg$hgnc = toupper(mouse.deg$gene)
# res = bigHeatmap(list('mz' = mz.deg, 'mm' = mouse.deg), pdf.name = "~/scratch/bcs/results/mz_zei_cluster_pct.pdf", rm.self = T)
# mz.deg.unique.genes = data.frame(table(mz.deg$gene))
# mz.deg.unique.genes = mz.deg.unique.genes[which(mz.deg.unique.genes[,2] == 1),1]
# mz.deg.unique = mz.deg[which(mz.deg$gene %in% mz.deg.unique.genes),]
# mouse.deg.unique.genes = data.frame(table(mouse.deg$gene))
# mouse.deg.unique.genes = mouse.deg.unique.genes[which(mouse.deg.unique.genes[,2] == 1),1]
# mouse.deg.unique = mouse.deg[which(mouse.deg$gene %in% mouse.deg.unique.genes),]
# res = bigHeatmap(list('mz' = mz.deg.unique, 'mm' = mouse.deg.unique), pdf.name = "~/scratch/bcs/results/mz_zei_cluster_small_pct_unique.pdf", rm.self = T)
# Testing
# Cichlid
findMZpct = function(this.ident) { this.pct = FoldChange(mz,    ident.1 = this.ident, features = mz.common.gene.set);    this.pct$gene = rownames(this.pct); this.pct$cluster = this.ident; return(data.frame(this.pct)) }
findMMpct = function(this.ident) { this.pct = FoldChange(mouse, ident.1 = this.ident, features = mouse.common.gene.set); this.pct$gene = rownames(this.pct); this.pct$cluster = this.ident; return(data.frame(this.pct)) }
mz.assay.to.use = ifelse("SCT" %in% names(mz@assays), "SCT", "RNA")
mz.avg.exp = AverageExpression(mz, features = mz.common.gene.set, assays = mz.assay.to.use, slot = "data")[[1]]
mz.avg.exp.norm = log(mz.avg.exp+1) + 0.1
mz.avg.exp.norm = mz.avg.exp.norm / rowMeans(mz.avg.exp.norm)
mz.pct.df = data.frame(data.table::rbindlist( parallel::mclapply(unique(Idents(mz)), function(x) findMZpct(x), mc.cores = 4) ))
mz.pct.df$pct.dif = mz.pct.df$pct.1 - mz.pct.df$pct.2
mz.pct1.mat = reshape2::acast(mz.pct.df, gene ~ cluster, value.var = "pct.1")
mz.pct2.mat = reshape2::acast(mz.pct.df, gene ~ cluster, value.var = "pct.2")
mz.pct.dif.mat = reshape2::acast(mz.pct.df, gene ~ cluster, value.var = "pct.dif")
mz.pct.dif.mat = t(scale(t(mz.pct.dif.mat[rownames(mz.avg.exp.norm), colnames(mz.avg.exp.norm)])))
# Mouse
mouse.avg.exp = AverageExpression(mouse, features = mouse.common.gene.set, assays = "SCT", slot = "data")[[1]]
mouse.avg.exp.norm = log(mouse.avg.exp+1) + 0.1
mouse.avg.exp.norm = mouse.avg.exp.norm / rowMeans(mouse.avg.exp.norm)
mouse.pct.df = data.frame(data.table::rbindlist( parallel::mclapply(unique(Idents(mouse)), function(x) findMMpct(x), mc.cores = 4) ))
mouse.pct.df$pct.dif = mouse.pct.df$pct.1 - mouse.pct.df$pct.2
mouse.pct1.mat = reshape2::acast(mouse.pct.df, gene ~ cluster, value.var = "pct.1")
mouse.pct2.mat = reshape2::acast(mouse.pct.df, gene ~ cluster, value.var = "pct.2")
mouse.pct.dif.mat = reshape2::acast(mouse.pct.df, gene ~ cluster, value.var = "pct.dif")
mouse.pct.dif.mat = t(scale(t(mouse.pct.dif.mat[rownames(mouse.avg.exp.norm), colnames(mouse.avg.exp.norm)])))
mz.mouse.cor = cor(mz.avg.exp.norm*mz.pct.dif.mat, mouse.avg.exp.norm*mouse.pct.dif.mat, method = "spearman")
pheatmap::pheatmap(mz.mouse.cor, cellwidth = 8, cellheight = 8, filename="~/scratch/bcs/results/spec_ind_test.pdf")