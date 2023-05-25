# Read Input ===================================================================
# this.run = 1; do.down = T; is.real = F; num.perms = 100;
# this.run = 1; do.down = F; is.real = T; num.perms = 1; ind = 0;
args = commandArgs(trailingOnly=TRUE)
my.dataset  = as.character(args[1])
meta.col = as.character(args[2])
set.seed(1)
message(paste0("Initializng run with: ", my.dataset))

# Load Libraries ===============================================================
suppressMessages(library('CellChat',  quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('patchwork', quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('stringr',   quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('parallel',  quietly = T, warn.conflicts = F, verbose = F))
options(stringsAsFactors = FALSE)
source("~/scratch/bcs/bcs_scripts/bcs_f.R")

# Load Data ====================================================================
message("Loading Object...")
gene_info = read.table("~/scratch/brain/cellchat/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
combined = readRDS("~/scratch/st/data/st_c2b2_hi_022023.rds")
if (my.dataset == "st.sc") { stsc = readRDS() }
message("Done.")

# Human Object =================================================================
# message("Creating a Human Object...")
# mz.df = data.frame(mz = rownames(combined@assays$RNA@counts), human = gene_info$human[match(rownames(combined@assays$RNA@counts), gene_info$mzebra)])
# mz.df$rowsums = rowSums(combined@assays$RNA@data)
# mz.df = mz.df[order(-mz.df$rowsums),]
# mz.df = mz.df[which(mz.df$rowsums != 0 & mz.df$human != "" & !is.na(mz.df$human)),]
# mz.df$human[which(mz.df$human == "NRG2")] = "NRG3"
# mz.df = rbind(data.frame(mz = "LOC101470250", human = "NRG2", rowsums = 5e6), mz.df)
# 
# mz.df = mz.df[!duplicated(mz.df$human),]
mz.df = read.csv("~/scratch/brain/cellchat/bb_cc_gene_name_converter.csv")
mz.df$X = NULL
# data.input = as.matrix(combined@assays$RNA@data[mz.df$mz,])
if (my.dataset == "st.sc") {
  
} else {
  combined = NormalizeData(combined)
  data.input = as.matrix(combined@assays$Spatial@data[mz.df$mz,])  
}
rownames(data.input) = mz.df$human
message("Done.")

# Cell Chat ====================================================================
CellChatWeights = function(x) {
  this.cells = colnames(data.input)
  this.meta = data.frame(label = meta.label, row.names = this.cells)
  
  cellchat = createCellChat(object = data.input[,this.cells], meta = this.meta, group.by = "label")
  cellchat = addMeta(cellchat, meta = this.meta)
  cellchat = setIdent(cellchat, ident.use = "label")
  
  cellchat@DB = CellChatDB.human
  cellchat = subsetData(cellchat)
  cellchat = identifyOverExpressedGenes(cellchat)
  cellchat = identifyOverExpressedInteractions(cellchat)
  
  cellchat = computeCommunProb(cellchat, type =  "truncatedMean", trim = 0.1, population.size = F)
  df.net_lig_recept = subsetCommunication(cellchat) 
  cellchat = aggregateNet(cellchat)
  cellchat = filterCommunication(cellchat, min.cells = 10)
  net_weight = data.frame(cellchat@net$weight)
  net_weight_vect = unlist(net_weight)
  name_rep = rep(rownames(net_weight), ncol(net_weight))
  names(net_weight_vect) = paste0(name_rep, ".", sort(name_rep))
  # test = reshape2::melt(as.matrix(net_weight))
  # colnames(test) = c("Sender", "Receiver", "value")
  # return(list(df.net_lig_recept, test))
  return(net_weight_vect)
}

# Getting correct labels
if (meta.col == "ct"        && my.dataset != "stsc") { meta.label = as.vector(combined@meta.data[, meta.col]) }
if (meta.col == "structure" && my.dataset != "stsc") { meta.label = combined@meta.data[, meta.col] }

rm(combined) # delete original Seurat object to save memory

message("Running cellchat (this while take awhile)...")
num.parallel.jobs = 1
message(paste0("Using ", num.parallel.jobs, " cores."))
sink(file="~/scratch/brain/cellchat_sink.txt")
run_outs = CellChatWeights(x)
sink()

n.success = length(run_outs)
# if (n.success != num.perms) { message(paste0("Not all runs were successful (", (num.perms - n.success), "/", num.perms, ")")) }
# out = as.data.frame(do.call('cbind', run_outs))
# colnames(out) = paste0("run", 1:n.success)
# out[, c("clust1", "clust2")] = reshape2::colsplit(names(run_outs[[1]]), "\\.", c("1", "2"))
# out = out[, c(n.success+1, n.success+2, 1:n.success)]
out = run_outs
# test = reshape2::melt(as.matrix(net_weight))
# colnames(test) = c("Sender", "Receiver", "value")
message("Done.")

# Save Output ==================================================================
message("Writing Output...")
todays.date = stringr::str_split(Sys.Date(), pattern = "-")[[1]]
todays.date = paste0(todays.date[2], todays.date[3], substr(todays.date[1], 3, 4))
out.str = paste0("~/scratch/st/results/cellchat/cellchat_", my.dataset, "_", meta.col, "_weights.csv")
write.csv(out, out.str)
message("Done.")
message("All Done.")
