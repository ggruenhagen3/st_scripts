# MISTy
library(mistyR)
library(future)

# Seurat
library(Seurat)

# data manipulation
library(Matrix)
library(tibble)
library(purrr)
library(dplyr)

# normalization
library(sctransform)

# resource
library(progeny)

# setup parallel execution
options(future.globals.maxSize = 1024^3)
plan(multisession)

# Helper Functions =============================================================
run_misty_seurat <- function(visium.slide,
                             # Seurat object with spatial transcriptomics data.
                             view.assays,
                             # Named list of assays for each view.
                             view.features = NULL,
                             # Named list of features/markers to use.
                             # Use all by default.
                             view.types,
                             # Named list of the type of view to construct
                             # from the assay.
                             view.params,
                             # Named list with parameters (NULL or value)
                             # for each view.
                             spot.ids = NULL,
                             # spot IDs to use. Use all by default.
                             out.alias = "results"
                             # folder name for output
) {
  
  # Extracting geometry
  geometry <- GetTissueCoordinates(visium.slide,
                                   cols = c("row", "col"), scale = NULL
  )
  
  # Extracting data
  view.data <- map(view.assays,
                   extract_seurat_data,
                   geometry = geometry,
                   visium.slide = visium.slide
  )
  
  # Constructing and running a workflow
  build_misty_pipeline(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids,
    out.alias = out.alias
  )
}

extract_seurat_data <- function(visium.slide,
                                assay,
                                geometry) {
  data <- GetAssayData(visium.slide, assay = assay) %>%
    t() %>%
    as_tibble(rownames = NA)
  
  return(data %>% slice(match(rownames(.), rownames(geometry))))
}

# Filters data to contain only features of interest
filter_data_features <- function(data,
                                 features) {
  if (is.null(features)) features <- colnames(data)
  
  return(data %>% rownames_to_column() %>%
           select(rowname, all_of(features)) %>% rename_with(make.names) %>%
           column_to_rownames())
}

# Builds views depending on the paramaters defined
create_default_views <- function(data,
                                 view.type,
                                 view.param,
                                 view.name,
                                 spot.ids,
                                 geometry) {
  view.data.init <- create_initial_view(data)
  
  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }
  
  if (view.type == "intra") {
    data.red <- view.data.tmp$data %>%
      rownames_to_column() %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>%
      add_paraview(geometry, l = view.param)
    
    data.ix <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>%
      add_juxtaview(
        positions = geometry,
        neighbor.thr = view.param
      )
    
    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  }
  
  if (is.null(view.param) == TRUE) {
    misty.view <- create_view(
      paste0(view.name),
      data.red
    )
  } else {
    misty.view <- create_view(
      paste0(view.name, "_", view.param),
      data.red
    )
  }
  
  return(misty.view)
}

# Builds automatic MISTy workflow and runs it
build_misty_pipeline <- function(view.data,
                                 view.features,
                                 view.types,
                                 view.params,
                                 geometry,
                                 spot.ids = NULL,
                                 out.alias = "default") {
  
  # Adding all spots ids in case they are not defined
  if (is.null(spot.ids)) {
    spot.ids <- rownames(view.data[[1]])
  }
  
  # First filter the features from the data
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # Create initial view
  views.main <- create_initial_view(view.data.filt[[1]] %>%
                                      rownames_to_column() %>%
                                      filter(rowname %in% spot.ids) %>%
                                      select(-rowname))
  
  # Create other views
  view.names <- names(view.data.filt)
  
  all.views <- pmap(list(
    view.data.filt[-1],
    view.types[-1],
    view.params[-1],
    view.names[-1]
  ),
  create_default_views,
  spot.ids = spot.ids,
  geometry = geometry
  )
  
  pline.views <- add_views(
    views.main,
    unlist(all.views, recursive = FALSE)
  )
  
  
  # Run MISTy
  run_misty(pline.views, out.alias)
}

# My Helper Functions ==========================================================
my_collect_results <- function(folders) {
  samples <- R.utils::getAbsolutePath(folders)
  
  message("\nCollecting improvements")
  improvements <- samples %>%
    furrr::future_map_dfr(function(sample) {
      performance <- readr::read_table2(paste0(sample, .Platform$file.sep, "performance.txt"),
                                       na = c("", "NA", "NaN"), col_types = readr::cols()
      ) %>% dplyr::distinct()
      
      performance %>%
        dplyr::mutate(
          sample = sample,
          gain.RMSE = 100 * (.data$intra.RMSE - .data$multi.RMSE) / .data$intra.RMSE,
          gain.R2 = .data$multi.R2 - .data$intra.R2
        )
    }, .progress = TRUE) %>%
    tidyr::pivot_longer(-c(.data$sample, .data$target), names_to = "measure")
  
  
  message("\nCollecting contributions")
  contributions <- samples %>% furrr::future_map_dfr(function(sample) {
    coefficients <- readr::read_table2(paste0(sample, .Platform$file.sep, "coefficients.txt"),
                                      na = c("", "NA", "NaN"), col_types = readr::cols()
    ) %>% dplyr::distinct()
    
    coefficients %>%
      dplyr::mutate(sample = sample, .after = "target") %>%
      tidyr::pivot_longer(cols = -c(.data$sample, .data$target), names_to = "view")
  }, .progress = TRUE)
  
  message("\nCollecting importances")
  importances <- samples %>%
    furrr::future_map_dfr(function(sample) {
      targets <- contributions %>%
        dplyr::filter(.data$sample == !!sample) %>%
        dplyr::pull(.data$target) %>%
        unique() %>%
        sort()
      views <- contributions %>%
        dplyr::pull(.data$view) %>%
        unique() %>%
        stringr::str_subset("^p\\.", negate = TRUE) %>%
        stringr::str_subset("^intercept$", negate = TRUE)
      
      # one heatmap per view
      maps <- views %>%
        furrr::future_map_dfr(function(view) {
          all.importances <- targets %>% purrr::map(~ readr::read_csv(paste0(
            sample, .Platform$file.sep, "importances_", .x, "_", view, ".txt"
          ),
          col_types = readr::cols()
          ) %>%
            dplyr::distinct() %>%
            dplyr::rename(feature = target))
          
          features <- all.importances %>%
            purrr::map(~ .x$feature) %>%
            unlist() %>%
            unique() %>%
            sort()
          
          pvalues <- contributions %>%
            dplyr::filter(.data$sample == !!sample, view == paste0("p.", !!view)) %>%
            dplyr::mutate(value = 1 - .data$value)
          
          # importances are standardized for each target
          # and multiplied by 1-pval(view)
          all.importances %>%
            purrr::imap_dfc(~
                              tibble::tibble(feature = features, zero.imp = 0) %>%
                              dplyr::left_join(.x, by = "feature") %>%
                              dplyr::arrange(.data$feature) %>%
                              dplyr::mutate(
                                imp = scale(.data$imp)[, 1],
                                !!targets[.y] := .data$zero.imp + (.data$imp *
                                                                     (pvalues %>%
                                                                        dplyr::filter(target == targets[.y]) %>%
                                                                        dplyr::pull(.data$value)))
                              )
                            %>%
                              dplyr::select(targets[.y])) %>%
            dplyr::mutate(Predictor = features) %>%
            tidyr::pivot_longer(
              names_to = "Target",
              values_to = "Importance",
              -.data$Predictor
            ) %>%
            dplyr::mutate(Importance = replace(
              .data$Importance,
              is.nan(.data$Importance), 0
            )) %>%
            dplyr::mutate(view = view, .before = 1)
        }) %>%
        dplyr::mutate(sample = sample, .before = 1)
    }, .progress = TRUE)
  
  message("\nAggregating")
  
  misty.results <- c(
    list(
      improvements = improvements,
      contributions = contributions,
      importances = importances
    ),
    my_aggregate_results(improvements, contributions, importances)
  )
  
  return(misty.results)
}

my_aggregate_results_subset <- function(misty.results, folders) {
  assertthat::assert_that(("importances" %in% names(misty.results)),
                          msg = "The provided result list is malformed. Consider using collect_results()."
  )
  
  normalized.folders <- R.utils::getAbsolutePath(folders)
  # check if folders are in names of misty.results
  assertthat::assert_that(all(normalized.folders %in%
                                (misty.results$importances %>% dplyr::pull(.data$sample))),
                          msg = "The provided results list doesn't contain information about some of
    the requested result folders. Consider using collect_results()."
  )
  
  message("Aggregating subset")
  importances.aggregated.subset <- misty.results$importances %>%
    dplyr::filter(.data$sample %in% normalized.folders) %>%
    tidyr::unite(".PT", "Predictor", "Target", sep = "&") %>%
    dplyr::group_by(.data$view, .data$.PT) %>%
    dplyr::summarise(
      Importance = mean(.data$Importance),
      nsamples = dplyr::n(), .groups = "drop"
    ) %>%
    tidyr::separate(".PT", c("Predictor", "Target"), sep = "&")
  
  misty.results[["importances.aggregated.subset"]] <- importances.aggregated.subset
  
  return(misty.results)
}

my_aggregate_results <- function(improvements, contributions, importances) {
  improvements.stats <- improvements %>%
    dplyr::filter(!stringr::str_starts(.data$measure, "p\\.")) %>%
    dplyr::group_by(.data$target, .data$measure) %>%
    dplyr::summarise(
      mean = mean(.data$value), sd = stats::sd(.data$value),
      cv = .data$sd / .data$mean, .groups = "drop"
    )
  
  
  contributions.stats <- dplyr::inner_join(
    # mean coefficients
    (contributions %>%
       dplyr::filter(!stringr::str_starts(.data$view, "p\\.") &
                       .data$view != "intercept") %>%
       dplyr::group_by(.data$target, .data$view) %>%
       dplyr::summarise(mean = mean(.data$value), .groups = "drop_last") %>%
       dplyr::mutate(fraction = abs(.data$mean) / sum(abs(.data$mean))) %>%
       dplyr::ungroup()),
    # p values
    (contributions %>%
       dplyr::filter(stringr::str_starts(.data$view, "p\\.") &
                       !stringr::str_detect(.data$view, "intercept")) %>%
       dplyr::group_by(.data$target, .data$view) %>%
       dplyr::mutate(view = stringr::str_remove(.data$view, "^p\\.")) %>%
       dplyr::summarise(
         p.mean = mean(.data$value),
         p.sd = stats::sd(.data$value),
         .groups = "drop"
       )),
    by = c("target", "view")
  )
  
  importances.aggregated <- importances %>%
    tidyr::unite(".PT", "Predictor", "Target", sep = "&") %>%
    dplyr::group_by(.data$view, .data$.PT) %>%
    dplyr::summarise(
      Importance = mean(.data$Importance),
      nsamples = dplyr::n(), .groups = "drop"
    ) %>%
    tidyr::separate(".PT", c("Predictor", "Target"), sep = "&")
  
  return(list(
    improvements.stats = improvements.stats,
    contributions.stats = contributions.stats,
    importances.aggregated = importances.aggregated
  ))
}


# Main Body ====================================================================
folder <- "breast_A_1"

# Load the HDF5 object and normalize the expression
seurat.vs <-
  Load10X_Spatial(
    data.dir = folder,
    filename = "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
  )

sct.data <- vst(GetAssayData(
  object = seurat.vs,
  slot = "counts",
  assay = "Spatial"
),
verbosity = 0
)

seurat.vs[["SCT"]] <- CreateAssayObject(data = sct.data$y)

# Filtering genes that are expressed in at least 5% of spots
gene.expression <- GetAssayData(seurat.vs, assay = "SCT")
coverage <- rowSums(gene.expression > 0) / ncol(gene.expression)
slide.markers <- names(which(coverage >= 0.05))

# Defining Hypoxia and Estrogen responsive genes
estrogen.footprints <- getModel(top = 15) %>%
  tibble::rownames_to_column("gene") %>%
  filter(Estrogen != 0, gene %in% slide.markers) %>%
  pull(gene)

hypoxia.footprints <- getModel(top = 15) %>%
  tibble::rownames_to_column("gene") %>%
  filter(Hypoxia != 0, gene %in% slide.markers) %>%
  pull(gene)

# Define assay for each view
view.assays <- list(
  "main" = "SCT",
  "para.hypoxia" = "SCT",
  "para.estrogen" = "SCT"
)

# Define features for each view
view.features <- list(
  "main" = hypoxia.footprints,
  "para.hypoxia" = hypoxia.footprints,
  "para.estrogen" = estrogen.footprints
)

# Define spatial context for each view
view.types <- list(
  "main" = "intra",
  "para.hypoxia" = "para",
  "para.estrogen" = "para"
)

# Define additional parameters (l in the case of paraview)
view.params <- list(
  "main" = NULL,
  "para.hypoxia" = 10,
  "para.estrogen" = 10
)

misty.out <- "vignette_model_seurat"

# Run MISTy pipeline and collect results
misty.results <- run_misty_seurat(
  visium.slide = seurat.vs,
  view.assays = view.assays,
  view.features = view.features,
  view.types = view.types,
  view.params = view.params,
  spot.ids = NULL, # Using the whole slide
  out.alias = misty.out
) %>%
  my_collect_results()

misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  arrange(value)

misty.results %>% plot_view_contributions()

misty.results %>% plot_interaction_heatmap(view = "intra")

misty.results$importances.aggregated %>% 
  filter(view == "intra", Target == "PGK1") %>%
  arrange(-Importance)

misty.results %>% plot_interaction_heatmap(view = "para.hypoxia_10")

misty.results %>% plot_interaction_heatmap(view = "para.estrogen_10")

# My Run =======================================================================
spo = qs::qread(paste0(data_dir, "st_obj_list_070822.qs"))
obj = spo[["b1c"]]

coverage <- rowSums(obj@assays$Spatial@counts > 0) / ncol(obj)
slide.markers <- names(which(coverage >= 0.05))
slide.markers.hgnc = gene_info$human[match(slide.markers, gene_info$seurat_name)]

# estrogen.footprints <- getModel(top = 15) %>% tibble::rownames_to_column("gene") %>% filter(Estrogen != 0, gene %in% slide.markers) %>% pull(gene)
# hypoxia.footprints  <- getModel(top = 15) %>% tibble::rownames_to_column("gene") %>% filter(Hypoxia  != 0, gene %in% slide.markers) %>% pull(gene)
estrogen.footprints.hgnc <- getModel(top = 100) %>% tibble::rownames_to_column("gene") %>% filter(Estrogen != 0, gene %in% slide.markers.hgnc) %>% arrange(desc(Estrogen)) %>% slice(1:15) %>% pull(gene)
hypoxia.footprints.hgnc  <- getModel(top = 100) %>% tibble::rownames_to_column("gene") %>% filter(Hypoxia  != 0, gene %in% slide.markers.hgnc) %>% arrange(desc(Hypoxia)) %>% slice(1:15) %>% pull(gene)
hypoxia.footprints.mz  = slide.markers[which(slide.markers.hgnc %in% hypoxia.footprints.hgnc)]
estrogen.footprints.mz = slide.markers[which(slide.markers.hgnc %in% estrogen.footprints.hgnc)]

# Define assay, features, view type, and extra parameters for each view
view.assays <- list("main" = "SCT", "para.hypoxia" = "SCT", "para.estrogen" = "SCT")
view.features <- list("main" = hypoxia.footprints.mz, "para.hypoxia" = hypoxia.footprints.mz, "para.estrogen" = estrogen.footprints.mz)
view.types <- list("main" = "intra", "para.hypoxia" = "para", "para.estrogen" = "para")
view.params <- list("main" = NULL, "para.hypoxia" = 10, "para.estrogen" = 10)

misty.results <- run_misty_seurat(visium.slide = obj, view.assays = view.assays, 
                                  view.features = view.features, view.types = view.types,
                                  view.params = view.params, spot.ids = NULL, out.alias = misty.out) %>% my_collect_results()

misty.results %>% plot_improvement_stats("gain.R2") %>% plot_improvement_stats("gain.RMSE")
misty.results$improvements %>% filter(measure == "p.R2") %>% arrange(value)
misty.results %>% plot_view_contributions()
misty.results %>% plot_interaction_heatmap(view = "intra")
misty.results$importances.aggregated %>% filter(view == "intra", Target == "LOC101464960") %>% arrange(-Importance)
misty.results %>% plot_interaction_heatmap(view = "para.hypoxia_10")
misty.results %>% plot_interaction_heatmap(view = "para.estrogen_10")

# My Run Complex ===============================================================
get.mz.path.genes = function(path, top.x = 15) {
  this.hgnc = getModel(top = 100) %>% tibble::rownames_to_column("gene") %>% filter(!!as.name(path)  != 0, gene %in% slide.markers.hgnc) %>% arrange(desc(!!as.name(path))) %>% slice(1:top.x) %>% pull(gene)
  this.mz = slide.markers[which(slide.markers.hgnc %in% this.hgnc)]
  return(this.mz)
}
path.genes = lapply(colnames(getModel(top = 1)), function(x) get.mz.path.genes(x))
path.scores = lapply(1:length(path.genes), function(x) colSums(obj@assays$SCT[path.genes[[x]],]))
score.mat = do.call('rbind', path.scores)
rownames(score.mat) = colnames(getModel(top = 1))
new.assay = CreateAssayObject(counts = score.mat)
obj[["path.activity"]] = new.assay

ic <- import_omnipath_intercell()
ligand = ic[which(ic$category == "ligand"),]
ligand = sort(unique(ligand$genesymbol))
ligand = ligand[which(ligand != "")]
ligand.mz = slide.markers[which(slide.markers.hgnc %in% ligand)]

view.assays <- list("main" = "path.activity", "para.ligand" = "SCT", "para.path" = "path.activity")
view.features <- list("main" = NULL, "para.ligand" = ligand.mz, "para.path" = NULL) # TODO check if NULL = all
view.types <- list("main" = "intra", "para.ligand" = "para", "para.path" = "para")
view.params <- list("main" = NULL, "para.ligand" = 10, "para.path" = 10)

misty.results <- run_misty_seurat(visium.slide = obj, view.assays = view.assays, 
                                  view.features = view.features, view.types = view.types,
                                  view.params = view.params, spot.ids = NULL, out.alias = misty.out) %>% my_collect_results()

misty.results %>% plot_improvement_stats("gain.R2") %>% plot_improvement_stats("gain.RMSE")
misty.results$improvements %>% filter(measure == "p.R2") %>% arrange(value)
misty.results %>% plot_view_contributions()
misty.results %>% plot_interaction_heatmap(view = "intra")
misty.results$importances.aggregated %>% filter(view == "intra", Target == "NFkB") %>% arrange(-Importance)
misty.results %>% plot_interaction_heatmap(view = "para.ligand_10")
misty.results %>% plot_interaction_heatmap(view = "para.path_10")

test = as.data.frame(misty.results$importances[which(misty.results$importances$view == "para.ligand_10"),])
test = test[which(test$Importance > 2),]
test$Predictor_label = gene_info$label[match(test$Predictor, gene_info$seurat_name)]
ggplot(test, aes(x = Predictor_label, y = Target, fill = Importance)) + geom_tile() + scale_fill_viridis() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
