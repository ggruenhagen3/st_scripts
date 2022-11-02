

obj_lr_path2 <- get_lr_path(object = obj,
                           celltype_sender = 'SST',
                           celltype_receiver = 'PVALB',
                           ligand = 'Sst',
                           receptor = 'Sstr2')

# deconvolution makes obj@data$newdata
# the number of cells per spot can be found in obj@meta[[1]]$cell_num

# I could add the # of cells based on cell2location and the percentages based on cell2location

# Mine =========================================================================
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

gene_info = read.table(paste0(main_path, "/all_research/gene_info_2.txt"), header = T, stringsAsFactors = F)
all_merge = qs::qread(paste0(data_dir, "st_070822.qs"))
spo = qs::qread(paste0(data_dir, "st_obj_list_070822.qs"))

st.counts =  spo[["b1c"]]@assays$Spatial@counts
st.meta = data.frame(spot = colnames(st.counts), x = spo[["b1c"]]@images$slice1@coordinates$imagecol, y = -spo[["b1c"]]@images$slice1@coordinates$imagerow)

obj <- createSpaTalk(st_data = st.counts, st_meta = st.meta, species = "Human", if_st_is_sc = F, spot_max_cell = 1000, celltype = st.celltype)
