# Read Input ===================================================================
args = commandArgs(trailingOnly=TRUE)
s = args[1]
my_n_cores = 24
if (length(args) > 1) { my_n_cores = args[2] }
message(paste0("Using paramters: s = ", s, ", my_n_cores = ", my_n_cores))

# Load Libraries and Data ======================================================
library("SpaTalk")
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

gene_info = read.table(paste0(main_path, "/gene_info_2.txt"), header = T, stringsAsFactors = F)
all_merge = qs::qread(paste0(data_dir, "st_070822.qs"))
spo = qs::qread(paste0(data_dir, "st_obj_list_070822.qs"))
bb = readRDS("~/scratch/brain/data/bb_demux_102021.rds")
# bb = readRDS("~/research/brain/data/bb_demux_102021.rds")

# SpaTalk Helpers ==============================================================
my_dec_celltype = function(object, sc_data, sc_celltype, st_celltype, n_cores = 4, if_doParallel = T, iter_num = 1000) {
  n.threads = n_cores
  st_coef = st_celltype
  st_meta = object@meta[['rawmeta']]
  st_data = object@data[["rawdata"]]
  st_dist = my_st_dist(st_meta)
  st_ndata = my_normalize_data(st_data)
  sc_ndata = sc_data@assays$RNA@data
  sc_celltype = data.frame(cell = colnames(sc_data), celltype = sc_celltype, stringsAsFactors = F)
  sc_celltype$celltype <- stringr::str_replace_all(sc_celltype$celltype, pattern = "-", replacement = "_")
  
  coef_name <- colnames(st_coef)
  coef_name <- coef_name[order(coef_name)]
  st_coef <- st_coef[ ,coef_name]
  st_coef = as.matrix(st_coef)
  object@coef <- st_coef
  # st_meta <- cbind(st_meta, my_coef_nor(st_coef))
  st_meta = cbind(st_meta, as.data.frame(st_coef))
  
  # st_dist <- object@dist
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  newmeta <- my_generate_newmeta_doParallel(st_meta, st_dist, 0.00000001, n_cores)
  doParallel::stopImplicitCluster()
  parallel::stopCluster(cl)
  
  print(head(newmeta))
  newmeta_cell <- my_generate_newmeta_cell(newmeta, st_ndata, sc_ndata, sc_celltype, iter_num, if_doParallel, n_cores)
  
  newdata <- sc_ndata[, newmeta_cell$cell_id]
  colnames(newdata) <- newmeta_cell$cell
  object@data$newdata <- methods::as(newdata, Class = "dgCMatrix")
  object@meta$newmeta <- newmeta_cell
  st_meta[st_meta$spot %in% newmeta_cell$spot, ]$celltype <- "sure"
  message("Calcualte distance matrix (George)")
  object@dist <- my_st_dist(newmeta_cell)
  message("Done.")
  object@meta$rawmeta <- st_meta
  return(object)
}

my_coef_nor <- function(st_coef) {
  for (i in 1:nrow(st_coef)) {
    st_coef1 <- as.numeric(st_coef[i, ])
    st_coef1 <- st_coef1/sum(st_coef1)
    st_coef[i, ] <- st_coef1
  }
  return(as.data.frame(st_coef))
}

my_get_weight1 <- function(st_meta_neighbor, sc_name){
  sc_name_ratio <- st_meta_neighbor[,sc_name]
  names(sc_name_ratio) <- st_meta_neighbor$w1
  if (!"I" %in% names(sc_name_ratio)) {
    sc_name_ratio <- c(sc_name_ratio, 0)
    names(sc_name_ratio)[length(sc_name_ratio)] <- "I"
  }
  if (!"II" %in% names(sc_name_ratio)) {
    sc_name_ratio <- c(sc_name_ratio, 0)
    names(sc_name_ratio)[length(sc_name_ratio)] <- "II"
  }
  if (!"III" %in% names(sc_name_ratio)) {
    sc_name_ratio <- c(sc_name_ratio, 0)
    names(sc_name_ratio)[length(sc_name_ratio)] <- "III"
  }
  if (!"IV" %in% names(sc_name_ratio)) {
    sc_name_ratio <- c(sc_name_ratio, 0)
    names(sc_name_ratio)[length(sc_name_ratio)] <- "IV"
  }
  sc_name_ratio <- sc_name_ratio[c("I", "II", "III", "IV")]
  sc_name_ratio <- sc_name_ratio + 1
  sc_name_ratio <- sc_name_ratio/sum(sc_name_ratio)
  sc_w1 <- rep(sc_name_ratio, each = 90)
  return(sc_w1)
}

my_get_weight2 <- function(st_meta_neighbor, sc_name, st_dist1, st_angle_new, spot_ratio){
  neighbor_weight <- 0
  if (st_angle_new > 0 & st_angle_new <= 90) {
    if ("I" %in% st_meta_neighbor$w1) {
      st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$w1 == "I",]
      neighbor_weight <- st_meta_neighbor1[ ,sc_name]
    }
  }
  if (st_angle_new > 90 & st_angle_new <= 180) {
    if ("II" %in% st_meta_neighbor$w1) {
      st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$w1 == "II",]
      neighbor_weight <- st_meta_neighbor1[ ,sc_name]
    }
  }
  if (st_angle_new > 180 & st_angle_new <= 270) {
    if ("III" %in% st_meta_neighbor$w1) {
      st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$w1 == "III",]
      neighbor_weight <- st_meta_neighbor1[ ,sc_name]
    }
  }
  if (st_angle_new > 270 & st_angle_new <= 360) {
    if ("IV" %in% st_meta_neighbor$w1) {
      st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$w1 == "IV",]
      neighbor_weight <- st_meta_neighbor1[ ,sc_name]
    }
  }
  spot_ratio <- c(spot_ratio, neighbor_weight)
  spot_ratio <- spot_ratio + 1
  spot_ratio <- spot_ratio/sum(spot_ratio)
  sc_w2 <- rep(spot_ratio, each = 5)
  sc_w2 <- sample(c(1:10), size = 1, prob = sc_w2)/10
  st_dist_new <- st_dist1[1]*sc_w2/2
  return(st_dist_new)
}

my_det_neighbor <- function(st_meta_neighbor, spot_x, spot_y, st_dist1){
  st_meta_neighbor$w1 <- "NA"
  # right-down
  st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x >= spot_x & st_meta_neighbor$y > spot_y,]
  if (nrow(st_meta_neighbor1)> 1) {
    neighbor_name <- st_meta_neighbor1$spot
    st_dist2 <- st_dist1[neighbor_name]
    st_dist3 <- names(which(st_dist2 == min(st_dist2)))
    if (length(st_dist3) > 1) {
      st_dist3 <- st_dist3[1]
    }
    st_dist2 <- names(st_dist2)
    st_dist2 <- st_dist2[!st_dist2 %in% st_dist3]
    st_meta_neighbor <- st_meta_neighbor[!st_meta_neighbor$spot %in% st_dist2,]
  }
  st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x >= spot_x & st_meta_neighbor$y > spot_y,]
  if (nrow(st_meta_neighbor1) > 0) {
    st_meta_neighbor[st_meta_neighbor$spot == st_meta_neighbor1$spot, ]$w1 <- "I"
  }
  # right-down
  st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x > spot_x & st_meta_neighbor$y <= spot_y,]
  if (nrow(st_meta_neighbor1)> 1) {
    neighbor_name <- st_meta_neighbor1$spot
    st_dist2 <- st_dist1[neighbor_name]
    st_dist3 <- names(which(st_dist2 == min(st_dist2)))
    if (length(st_dist3) > 1) {
      st_dist3 <- st_dist3[1]
    }
    st_dist2 <- names(st_dist2)
    st_dist2 <- st_dist2[!st_dist2 %in% st_dist3]
    st_meta_neighbor <- st_meta_neighbor[!st_meta_neighbor$spot %in% st_dist2,]
  }
  st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x > spot_x & st_meta_neighbor$y <= spot_y,]
  if (nrow(st_meta_neighbor1) > 0) {
    st_meta_neighbor[st_meta_neighbor$spot == st_meta_neighbor1$spot, ]$w1 <- "IV"
  }
  # left-up
  st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x < spot_x & st_meta_neighbor$y >= spot_y,]
  if (nrow(st_meta_neighbor1)> 1) {
    neighbor_name <- st_meta_neighbor1$spot
    st_dist2 <- st_dist1[neighbor_name]
    st_dist3 <- names(which(st_dist2 == min(st_dist2)))
    if (length(st_dist3) > 1) {
      st_dist3 <- st_dist3[1]
    }
    st_dist2 <- names(st_dist2)
    st_dist2 <- st_dist2[!st_dist2 %in% st_dist3]
    st_meta_neighbor <- st_meta_neighbor[!st_meta_neighbor$spot %in% st_dist2,]
  }
  st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x < spot_x & st_meta_neighbor$y >= spot_y,]
  if (nrow(st_meta_neighbor1) > 0) {
    st_meta_neighbor[st_meta_neighbor$spot == st_meta_neighbor1$spot, ]$w1 <- "II"
  }
  # left-down
  st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x <= spot_x & st_meta_neighbor$y < spot_y,]
  if (nrow(st_meta_neighbor1)> 1) {
    neighbor_name <- st_meta_neighbor1$spot
    st_dist2 <- st_dist1[neighbor_name]
    st_dist3 <- names(which(st_dist2 == min(st_dist2)))
    if (length(st_dist3) > 1) {
      st_dist3 <- st_dist3[1]
    }
    st_dist2 <- names(st_dist2)
    st_dist2 <- st_dist2[!st_dist2 %in% st_dist3]
    st_meta_neighbor <- st_meta_neighbor[!st_meta_neighbor$spot %in% st_dist2,]
  }
  st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x <= spot_x & st_meta_neighbor$y < spot_y,]
  if (nrow(st_meta_neighbor1) > 0) {
    st_meta_neighbor[st_meta_neighbor$spot == st_meta_neighbor1$spot, ]$w1 <- "III"
  }
  return(st_meta_neighbor)
}

my_generate_newmeta_doParallel <- function(st_meta, st_dist, min_percent, n_cores) {
  # generate new data
  st_meta <- st_meta[st_meta$label != "less nFeatures", ]
  cellname <- colnames(st_meta)[-c(1:7)]
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  newmeta <- foreach::foreach (i = 1:nrow(st_meta), .combine = rbind, .packages = "Matrix", .export = c("my_det_neighbor", "my_get_weight1", "my_get_weight2")) %dopar% {
    # for (i in 1:nrow(st_meta)) {
    newmeta_spot <- NULL
    newmeta_ratio <- NULL
    newmeta_cell <- NULL
    newmeta_x <- NULL
    newmeta_y <- NULL
    spot_name <- st_meta$spot[i]
    spot_x <- st_meta$x[i]
    spot_y <- st_meta$y[i]
    spot_percent <- as.numeric(st_meta[i, -c(1:7)])
    spot_percent[which(is.na(spot_percent))] = 0
    spot_cellnum <- st_meta$cell_num[i]
    spot_celltype <- which(spot_percent > 0)
    if (length(spot_celltype) > 0) {
      spot_percent <- spot_percent[spot_celltype]
      spot_celltype <- cellname[spot_celltype]
      spot_cell <- NULL
      for (j in 1:length(spot_celltype)) {
        newmeta_ratio <- c(newmeta_ratio, rep(spot_percent[j], spot_percent[j]))
        # newmeta_cell <- c(newmeta_cell, spot_celltype[j])
        spot_cell <- c(spot_cell, rep(spot_celltype[j], spot_percent[j]))
      }
    }
    newmeta_cell = spot_cell
    k = spot_cellnum
    if (k > 0) {
      newmeta_spot <- c(newmeta_spot, rep(spot_name, k))
      if (k == 1) {
        newmeta_x <- c(newmeta_x, spot_x)
        newmeta_y <- c(newmeta_y, spot_y)
      } else {
        n_neighbor <- 4
        st_dist1 <- st_dist[, spot_name]
        st_dist1 <- st_dist1[st_dist1 > 0]
        st_dist1 <- st_dist1[order(st_dist1)]
        st_dist1 <- st_dist1[1:n_neighbor]
        st_meta_neighbor <- st_meta[st_meta$spot %in% names(st_dist1), ]
        st_meta_neighbor <- my_det_neighbor(st_meta_neighbor, spot_x, spot_y, st_dist1)
        if (nrow(st_meta_neighbor) == 0) {
          for (j in 1:k) {
            st_angle_new <- sample(x = c(1:360), size = 1)
            st_dist_new <- sample(x = c(0:st_dist1), size = 1)
            newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * pi/180)
            newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * pi/180)
            newmeta_x <- c(newmeta_x, newmeta_x1)
            newmeta_y <- c(newmeta_y, newmeta_y1)
          }
        } else {
          for (j in 1:k) {
            sc_name <- newmeta_cell[j]
            sc_w1 <- my_get_weight1(st_meta_neighbor, sc_name)
            set.seed(j)
            st_angle_new <- sample(x = c(1:360), size = 1, prob = as.numeric(sc_w1))
            spot_ratio <- st_meta[st_meta$spot == spot_name, sc_name]
            st_dist_new <- my_get_weight2(st_meta_neighbor, sc_name, st_dist1, st_angle_new, spot_ratio)
            newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * pi/180)
            newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * pi/180)
            newmeta_x <- c(newmeta_x, newmeta_x1)
            newmeta_y <- c(newmeta_y, newmeta_y1)
          }
        }
      }
      data.frame(spot = newmeta_spot, cell_ratio = newmeta_ratio, celltype = newmeta_cell, x = as.numeric(newmeta_x), y = as.numeric(newmeta_y), stringsAsFactors = F)
    } else {
      data.frame(spot = "NA", cell_ratio = "NA", celltype = "NA", x = "NA", y = "NA", stringsAsFactors = F)
    }
  } # end foreach
  doParallel::stopImplicitCluster()
  parallel::stopCluster(cl)
  return(newmeta)
}

my_generate_newmeta_cell <- function(newmeta, st_ndata, sc_ndata, sc_celltype, iter_num, if_doParallel, n_cores) {
  newmeta_spotname <- unique(newmeta$spot)
  newmeta_cell <- NULL
  cat(crayon::cyan("Generating single-cell data for each spot", "\n"))
  if (if_doParallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    newmeta_cell <- foreach::foreach (i = 1:length(newmeta_spotname), .combine = "rbind", .packages = "Matrix", .export = "my_generate_newmeta_spot") %dopar% {
      spot_name <- newmeta_spotname[i]
      my_generate_newmeta_spot(spot_name, newmeta, st_ndata, sc_ndata, sc_celltype, iter_num)
    }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
  } else {
    for (i in 1:length(newmeta_spotname)) {
      spot_name <- newmeta_spotname[i]
      newmeta_spot <- my_generate_newmeta_spot(spot_name, newmeta, st_ndata, sc_ndata, sc_celltype, iter_num)
      newmeta_cell <- rbind(newmeta_cell, newmeta_spot)
    }
  }
  newmeta_cell$cell <- paste0("C", 1:nrow(newmeta))
  newmeta_cell <- newmeta_cell[, c(8, 4, 5, 3:1, 7, 6)]
  return(newmeta_cell)
}

my_generate_newmeta_spot <- function(spot_name, newmeta, st_ndata, sc_ndata, sc_celltype, iter_num) {
  newmeta_spot <- newmeta[newmeta$spot == spot_name, ]
  spot_ndata <- as.numeric(st_ndata[, spot_name])
  # random sampling
  score_cor <- NULL
  spot_cell_id <- list()
  for (k in 1:iter_num) {
    cat(paste0(k, "."))
    spot_cell_id_k <- NULL
    for (j in 1:nrow(newmeta_spot)) {
      spot_celltype <- newmeta_spot$celltype[j]
      sc_celltype1 <- sc_celltype[sc_celltype$celltype == spot_celltype, "cell"]
      sc_celltype1 <- sc_celltype1[sample(x = 1:length(sc_celltype1), size = 1)]
      spot_cell_id_k <- c(spot_cell_id_k, sc_celltype1)
    }
    if (length(spot_cell_id_k) == 1) {
      spot_ndata_pred <- as.numeric(sc_ndata[, spot_cell_id_k])
    } else {
      spot_ndata_pred <- as.numeric(rowSums(sc_ndata[, spot_cell_id_k]))
    }
    spot_ndata_cor <- cor(spot_ndata, spot_ndata_pred)
    score_cor <- c(score_cor, spot_ndata_cor)
    spot_cell_id[[k]] <- spot_cell_id_k
  }
  spot_cell_id <- spot_cell_id[[which.max(score_cor)]]
  newmeta_spot$cell_id <- spot_cell_id
  newmeta_spot$cor <- max(score_cor)
  return(newmeta_spot)
}

my_normalize_data <- function(rawdata) {
  rawdata <- Seurat::CreateSeuratObject(rawdata)
  rawdata <- Seurat::NormalizeData(rawdata, verbose = F)
  rawdata <- rawdata[["RNA"]]@data
  return(rawdata)
}

my_st_dist <- function(st_meta) {
  st_dist <- as.matrix(stats::dist(x = cbind(st_meta$x, st_meta$y)))
  rownames(st_dist) <- st_meta[, 1]
  colnames(st_dist) <- st_meta[, 1]
  return(st_dist)
}

# Main Body ====================================================================
st.counts =  spo[[s]]@assays$Spatial@counts
st.meta = data.frame(spot = colnames(st.counts), x = spo[[s]]@images$slice1@coordinates$imagecol, y = -spo[[s]]@images$slice1@coordinates$imagerow)

# Cell2location integration results
st.celltype = read.csv(paste0(out_dir, "cell2location_spatial_output_means.csv"))
st.celltype = st.celltype[match(colnames(st.counts), st.celltype$X),]
rownames(st.celltype) = st.celltype$X
st.celltype = st.celltype[,2:ncol(st.celltype)]
colnames(st.celltype) = str_replace(colnames(st.celltype), "meanscell_abundance_w_sf_", "")

# Get the celltype with the most cells per spot
st.celltype.char = unlist(lapply(1:nrow(st.celltype), function(x) colnames(st.celltype)[which.max(st.celltype[x,])]))
st.celltype = round(st.celltype)
zero.cell.st = which(rowSums(st.celltype) == 0)
for (i in zero.cell.st) { this.st.celltype.char = st.celltype.char[i]; st.celltype[i, this.st.celltype.char] = 1; }

# SpaTalk
message("Creating Spatalk object (George)")
obj = createSpaTalk(st_data = st.counts[,1:100], st_meta = st.meta[1:100,], species = "Human", if_st_is_sc = F, spot_max_cell = 1000, celltype = st.celltype.char[1:100])
# When creating the SpaTalk object, there's a line that removes genes with 0 expression: st_data <- st_data[which(rowSums(st_data) > 0), ]
# But I want a matrix in the end with every gene, so I'm going to modify the object's data to include the genes w/ 0 expression
obj = new("SpaTalk", data = list(rawdata = st.counts[,1:100]), meta = list(rawmeta = obj@meta[['rawmeta']]),
              para = list(species = 'Human', st_type = 'spot', spot_max_cell = 1000, if_skip_dec_celltype = T))
obj@meta[['rawmeta']]$cell_num = rowSums(st.celltype[1:100,])
message("Creating Single Cell Matrix (George)")
obj = my_dec_celltype(object = obj, sc_data = bb, sc_celltype = as.vector(bb$seuratclusters53), st_celltype = st.celltype[1:100,], n_cores = my_n_cores)
message("Saving object (George)")
saveRDS(obj, paste0(out_dir, "/spatalk/", s, ".rds"))

# Convert to Human Names
# In cases of multiple mz -> hgnc, keep the mz with the greatest mean ranked expression in st and snRNA-seq
message("Converting cichlid to human names (George)")
sn.rank = rowSums(bb@assays$RNA@data)
sn.rank = names(sort(sn.rank, decreasing = T))
st.rank = rowSums(my_normalize_data(st.counts))
st.rank = names(sort(st.rank, decreasing = T))

all.rank = data.frame(mz = gene_info$seurat_name, hgnc = gene_info$human)
all.rank$sn.rank = match(all.rank$mz, sn.rank)
all.rank$st.rank = match(all.rank$mz, st.rank)
all.rank$mean.rank = rowMeans(all.rank[, c("sn.rank", "st.rank")])
all.rank = all.rank[order(all.rank$mean.rank, decreasing = F),]
all.rank = all.rank[which(!duplicated(all.rank$hgnc) & !is.na(all.rank$hgnc)),]

hraw = obj@data[['rawdata']][all.rank$mz,]
hnew = obj@data[['newdata']][all.rank$mz,]
rownames(hraw) = rownames(hnew) = all.rank$hgnc
meta_w_celltype = obj@meta[['rawmeta']]
meta_w_celltype$celltype = st.celltype.char[1:100]
hobj = new("SpaTalk", data = list(rawdata = hraw, newdata = hnew), meta = list(rawmeta = meta_w_celltype, newmeta = obj@meta[['newmeta']]),
           dist = obj@dist,
           para = list(species = 'Human', st_type = 'spot', spot_max_cell = 1000, if_skip_dec_celltype = F))
hobj@data$rawndata = hobj@data$rawdata

# Cell-cell interactions
# hobj <- find_lr_path(object = hobj, lrpairs = lrpairs, pathways = pathways)
# hobj <- dec_cci(object = hobj, celltype_sender = '0', celltype_receiver = '1') # if_skip_dec_celltype has ta be F for this to work
# obj_lr_path <- get_lr_path(object = hobj, celltype_sender = '0', celltype_receiver = '1', ligand = 'SEMA3F', receptor = 'PLXNA3') # if you don't put a valid ligand/receptor pair, you will get an error

message("Finding CCI for all Cell Types and ligand/receptor pairs (George)")
hobj <- dec_cci_all(object = hobj, n_cores = my_n_cores)
message("Saving Human Object (George)")
saveRDS(hobj, paste0(out_dir, "/spatalk/", s, "_human.rds"))
message("All Done.")

