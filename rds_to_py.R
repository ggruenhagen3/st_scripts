args = commandArgs(trailingOnly=TRUE)
rds_path = args[1]
isRDS     = as.logical(args[2])
isSpatial = as.logical(args[3])

if (isRDS)     { obj = readRDS(rds_path) } else { obj = qs::qread(rds_path) }
if (isSpatial) { obj_assay = "Spatial"   } else { obj_assay = "RNA"         }

message("Dieting object")
DefaultAssay(obj) = obj_assay
obj = Seurat::DietSeurat(obj, dimreducs = "umap", assays = obj_assay)
base_name = stringr::str_split(rds_path, "\\.")[[1]][1]
message("Saving ")
SaveH5Seurat(obj, filename = paste0(base_name, ".h5seurat"))
Convert(paste0(base_name, ".h5seurat"), dest = "h5ad")
message("All Done Converting Seurat to Python")
