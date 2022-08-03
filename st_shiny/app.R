# =========== #
# Spatial App #
# =========== #
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2")
library("shinycssloaders")
library("shinyWidgets")
library("DT")
library("qs")
library("RColorBrewer")
library("patchwork")
# library("shinyjs")
library("grid")
# library("ggpubr")
library("cowplot")
options(warn=-1)

# Load object
# all_merge_diet = DietSeurat(all_merge, counts = F, data = T, assays = "Spatial", dimreducs = "umap")
obj = qs::qread("data/st_diet_080322.qs")
obj@active.assay = "Spatial"

# Set Some variables for the whole app
gene_info = read.table("data/gene_info_2.txt", sep="\t", stringsAsFactors = F, header = T)
gene_info = gene_info[which(! duplicated(gene_info$mzebra) ),]
gene_info$human_description[which(gene_info$human_description == "")] = NA
gene_info$mzebra_description = do.call('rbind', strsplit(as.character(gene_info$mzebra_description),' [Source',fixed=TRUE) )[,1]
gene_info$human_description = do.call('rbind', strsplit(as.character(gene_info$human_description),' [Source',fixed=TRUE) )[,1]

gene_info$human_orig = gene_info$human
gene_info$human = toupper(gene_info$label)
gene_info$human[which(startsWith(gene_info$human, "LOC"))] = gene_info$human_orig[which(startsWith(gene_info$human, "LOC"))]
human_genes = sort(unique(toupper(gene_info$human)))
human_genes = sort(unique(gene_info$human))
human_gene_names = human_genes
gene_names <- rownames(obj@assays$Spatial)

sample_pt_size = c(4, 3, 3, 4, 3, 3, 3, 4, 3, 4, 4, 4, 3)
names(sample_pt_size) = levels(obj$sample)

num_clusters = max(as.vector(obj$cluster))
clusters = 1:num_clusters
cluster_choices = c(clusters)

bhve_samples = c("b1", "b2")
ctrl_samples = c("c1", "c2")

discrete.colors = c("#f8e16c", "#00c49a", "#156064", "#ffecd100")

# Functions ====================================================================
# 07/25/22
mySingleSFP = function(obj = NULL, feature = NULL, assay = NULL, slot = NULL, coords = NULL, values = NULL, img.grob = NULL, points.as.text = F, rot.text = T, my.pt.size = 0.8, zoom.out.factor = 0.05, pal = colorRampPalette(colors = rev(brewer.pal(11, "Spectral"))), col.min = NULL, col.max = NULL, angle = NULL, discrete = F, rm.zero = F, col.ident = F) {
  #' My version of SpatialFeaturePlot for a single object.
  #' 
  #' Input either an object+feature+assay+slot or coords+values+image grob.
  #' 
  #' @param obj Seurat object
  #' @param feature a feature such as a gene
  #' @param assay assay of Seurat object to pull data from
  #' @param slot slot of Seurat object to pull data from
  #' 
  #' @param coords tissue coordinates
  #' @param values values to plot
  #' @param img.grob image grob of tissue (can be retrieved with GetImage)
  #' 
  #' @param my.pt.size point size of spatial dots
  #' @param zoom.out.factor zoom_out on tissue window (larger number -> more zoomed out)
  #' @param pal color pallette
  #' @param col.min min value for coloring (helpful for consistent color scale across multiple plots)
  #' @param col.max max value for coloring (helpful for consistent color scale across multiple plots)
  #' 
  #' @return ggplot object
  
  # Input Checking
  if (is.null(obj)) {
    if (is.null(coords) || is.null(values) || is.null(img.grob)) { message("Error. No object input, but coords, values, and/or img.grob were not specified."); return(NULL); }
    coords$value = values
  } else {
    if (is.null(feature) || is.null(assay) || is.null(slot)) { message("Error. Object was input, but no feature, assay, and/or slot were specified."); return(NULL) }
    coords = GetTissueCoordinates(object = obj)
    
    # See if feature is a gene or metadata
    if (feature %in% rownames(obj[[assay]])) {
      obj@active.assay = assay
      coords$value = FetchData(object = obj, vars = feature, slot = slot)[,1]
    } else if (feature %in% colnames(obj@meta.data)) {
      coords$value = obj@meta.data[,feature]
    } else {
      message(paste0("Error. Feature, ", feature, ", was not found in the assay nor in the metadata for the object."))
      return (NULL)
    }
    
    img.grob = GetImage(obj)
  }
  
  # Finding Coordinates of the Tissue
  my.imagerow.min = min(coords$imagerow)
  my.imagecol.min = min(coords$imagecol)
  my.imagerow.max = max(coords$imagerow)
  my.imagecol.max = max(coords$imagecol)
  myratio = (my.imagerow.max - my.imagerow.min) / (my.imagecol.max - my.imagecol.min)
  
  # Adjusting the tissue window: zoom out or zoom in
  my.x.width = my.imagerow.max - my.imagerow.min
  my.x.zoom = my.x.width * zoom.out.factor
  my.x.min = round(my.imagerow.min - my.x.zoom)
  my.x.max = round(my.imagerow.max + my.x.zoom)
  my.y.width = my.imagecol.max - my.imagecol.min
  my.y.zoom = my.y.width * zoom.out.factor
  my.y.min = round(my.imagecol.min - my.y.zoom)
  my.y.max = round(my.imagecol.max + my.y.zoom)
  img.grob.test = img.grob$raster[my.x.min:my.x.max, my.y.min:my.y.max]
  img.grob.test.grob = img.grob
  img.grob.test.grob$raster = img.grob.test
  
  # Set color range
  if (!discrete & is.null(col.min) && is.null(col.max)) { col.min = min(coords$value); col.max = max(coords$value); }
  
  # Remove spots with zero expression
  if (rm.zero) { coords = coords[which(coords$value != 0),] }
  
  # Plot
  if (points.as.text) {
    # This is for plotting the number of the cluster at each spot instead of plotting spots as points
    this.angle = 0
    if (rot.text) { this.angle = angle }
    p = ggplot(coords, aes(x=imagecol, y=-imagerow, color = value)) + annotation_custom(img.grob.test.grob) + geom_point(shape = 1, size = my.pt.size, stroke = 0) + geom_text(size = my.pt.size*0.8, angle = -this.angle, aes(label = value)) + scale_color_gradientn(colors=pal(100), limits = c(col.min, col.max)) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio)
  } else {
    if (discrete) {
      # Discrete Color Scales
      p = ggplot(coords, aes(x=imagecol, y=-imagerow, color = value)) + annotation_custom(img.grob.test.grob) + geom_point(size = my.pt.size, stroke = 0) + scale_color_manual(values=pal, drop = F)                             + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio)
    } else {
      # Continuous Color Scales
      p = ggplot(coords, aes(x=imagecol, y=-imagerow, color = value)) + annotation_custom(img.grob.test.grob) + geom_point(size = my.pt.size, stroke = 0) + scale_color_gradientn(colors=pal(100), limits = c(col.min, col.max)) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio)
    }
    if (col.ident) {
      # Identity Color Scale
      p = p + scale_color_identity() 
    }
  }
  
  return(p)
}

allSamplesSFP = function(obj, feature, assay = "SCT", slot = "data", points.as.text = F, rot.text = T, pt.size.multiplier = 1, zoom.out.factor = 0.05, pal = colorRampPalette(colors = rev(brewer.pal(11, "Spectral"))), rm.zero = F, col.ident = F) {
  #' My version of SpatialFeaturePlot for all samples. Plots are angled to the correct orientation.
  #' 
  #' @param obj Seurat object that contains all samples
  #' @param assay assay of Seurat object to pull data from
  #' @param slot slot of Seurat object to pull data from
  #' @param pt.size.multiplier multiplier to increase/decrease size of points for all samples proportionally
  #' @param zoom.out.factor zoom_out on tissue window (larger number -> more zoomed out)
  #' @param pal color pallette
  #' 
  #' @return nothing
  
  obj@active.assay = assay
  
  # Samples to Plot
  real.samples = c("c2a", "c2b", "c2c", "c2d", "b2a", "b2b", "b2c", "b2d", "c1a", "c1b", "c1c", "c1d", "b1c")
  all.samples = c(real.samples[1:(length(real.samples)-1)], "b1a", "b1b", "b1c", "b1d")
  
  # Get all values for feature
  value_list = list()
  isDiscrete = F
  if (length(feature) == 1) {
    for (s in real.samples) {
      # See if feature is a gene or metadata
      if (feature %in% rownames(obj[[assay]])) {
        value_list[[s]] = FetchData(object = obj, vars = feature, cells = colnames(obj)[which(obj$sample == s)], slot = slot)[,1]
      } else if (feature %in% colnames(obj@meta.data)) {
        value_list[[s]] = obj@meta.data[colnames(obj)[which(obj$sample == s)], feature]
      } else {
        message(paste0("Error. Feature, ", feature, ", was not found in the assay nor in the metadata for the object."))
        return (NULL)
      }
    }
    min.val = min(unlist(value_list))
    max.val = max(unlist(value_list))
  } else {
    isDiscrete = T
    message("Multiple features selected. Using binary expression values from Spatial counts.")
    
    for (f in feature) {
      if (! f %in% rownames(obj@assays$Spatial@data)) {
        message(paste0("Error. Feature, ", f, ", was not found in the object."))
        return (NULL)
      }
      all_value_sum = colSums(obj@assays$Spatial@data[c(feature[1], feature[2], feature[2]),] > 0)
      all_value_sum = plyr::revalue(as.character(all_value_sum), replace = c("0" = "none", "1" = feature[1], "2" = feature[2], "3" = "both"))
      all_value_sum = factor(all_value_sum, levels = c("both", feature, "none"))
      for (s in real.samples) {
        this.value = all_value_sum[which(obj$sample == s)]
        value_list[[s]] = this.value
      }
    }
    min.val = "dummy"
    max.val = "dummy"
  }
  
  
  # Point Sizes
  sample_pt_size = c(2.3, 1.75, 1.75, 2.1, 1.4, 1.5, 1.5, 1.8, 1.25, 2.3, 2, 2, 1.3)
  names(sample_pt_size) = real.samples
  
  # Angle to rotate plots
  angle.df = as.data.frame(c("c2a" = 155, "c2b" = 145, "c2c" = -115, "c2d" = 155,
                             "b2a" =  90, "b2b" =  95, "b2c" =  100, "b2d" =  85,
                             "c1a" =  90, "c1b" =  95, "c1c" =   98, "c1d" =  90,
                             "b1a" =   0, "b1b" =   0, "b1c" = -118, "b1d" =   0))
  colnames(angle.df) = "angle"
  
  # Create all the separate sample plots
  p_list = list()
  for (s in real.samples) {
    this.coords = GetTissueCoordinates(object = obj, image = s)
    this.values = value_list[[s]]
    this.img.grob = GetImage(obj, image = s)
    this.angle = angle.df[s, "angle"]
    p_list[[s]] = mySingleSFP(coords = this.coords, assay = assay, slot = slot, values = this.values, img.grob = this.img.grob, points.as.text = F, rot.text = T, my.pt.size = sample_pt_size[s]*pt.size.multiplier, zoom.out.factor = 0.05, pal = pal, col.min = min.val, col.max = max.val, angle = this.angle, discrete = isDiscrete, rm.zero = rm.zero, col.ident = col.ident)
  }
  
  # Get color legend
  if (isDiscrete) {
    leg_p = ggplot(data.frame(a = this.value, b = 1), aes(a, b, color = this.value)) + geom_point() + scale_color_manual(values=pal, name = NULL, drop = F) + theme(legend.position = 'bottom', legend.background = element_blank(), legend.key=element_blank()) + guides(color = guide_legend(ncol=2, by.row=T,  override.aes = list(size=3)))
  } else if (col.ident) {
    print(head(unique(unlist(value_list))))
    leg_p = ggplot(data.frame(a = unique(unlist(value_list)), b = 1), aes(a, b, color = a)) + geom_point() + scale_color_identity(name = NULL) + theme(legend.position = 'bottom', legend.background = element_blank())     
  } else {
    print(min.val)
    print(max.val)
    leg_p = ggplot(data.frame(a = 1, b = 1), aes(a, b, color = a)) + geom_point() + scale_color_gradientn(colors=pal(100), limits = c(min.val, max.val), name = NULL) + theme(legend.position = 'bottom', legend.background = element_blank()) 
  }
  leg_p_grob = get_legend(leg_p)
  
  # Create a grid of plots using viewports 
  all.samples.mtx = matrix(all.samples, ncol = 4, byrow = T)
  grid.newpage()
  pushViewport(viewport(width=1, height=1, xscale=c(0, nrow(all.samples.mtx)), yscale=c(0,ncol(all.samples.mtx))))
  for(i in 1:nrow(all.samples.mtx)){
    for(j in 1:ncol(all.samples.mtx)){
      s = as.character(all.samples.mtx[i,j])
      this.x = j-0.5 # coordinate for the position on the grid
      this.y = nrow(all.samples.mtx)-i+0.5 # coordinate for the position on the grid
      vp = viewport(x=unit(this.x,"native"), y=unit(this.y,"native"), width=unit(1,"native"), height=unit(1,"native"), clip=T) # current cell of the grid
      pushViewport(vp) # starts the viewport window
      
      # Determine background color and reorientation-angle of the sample
      this.back.color = switch(s, "c2a" = "#AEADAD", "c2b" = "#AEADAD", "c2c" = "#AEADAD", "c2d" = "#AEADAD",
                               "b2a" = "#AEADAD", "b2b" = "#AEADAD", "b2c" = "#AEADAD", "b2d" = "#AEADAD",
                               "c1a" = "#B6B6B5", "c1b" = "#B6B6B5", "c1c" = "#B6B6B5", "c1d" = "#B6B6B5",
                               "b1a" = "#B4B3B2", "b1b" = "#B4B3B2", "b1c" = "#B4B3B2", "b1d" = "#B4B3B2")
      this.angle = angle.df[s, "angle"]
      
      # Set the background color
      grid.rect(gp=gpar(fill=this.back.color, col=this.back.color)) # background color
      
      # Create a viewport within the current viewport (the one that represents the cell of the grid) 
      # This allows me to adjust the angle of the plots
      vp2 = viewport(width = unit(0.9, "npc"), height = unit(0.9, "npc"), angle = this.angle, clip = T)
      pushViewport(vp2)
      if (! s %in% c("b1a", "b1b", "b1d")) {
        grid.draw(ggplotGrob(p_list[[s]]))
      } else if (s == "b1a") {
        grid.text(as.character(feature), 0.5, 0.5, gp=gpar(cex=1.3))
      } else if (s == "b1b") {
        grid.draw(leg_p_grob)
      }
      popViewport() # pop back to the viewport that represents the cell of the grid
      
      # Sample Labels
      if (! s %in% c("b1a", "b1b", "b1d")) {
        grid.text(as.character(s), 0.1, 0.95, gp=gpar(cex=0.75))
      } 
      
      upViewport()
      
    }
  }
  popViewport(1)
}


# UI ===========================================================================
# Define UI for app
ui <- fluidPage(
  
  tags$head(
    tags$style(
      HTML(
        ".checkbox-inline { 
                    margin-left: 0px;
                    margin-right: 10px;
          }
         .checkbox-inline+.checkbox-inline {
                    margin-left: 0px;
                    margin-right: 10px;
          }
        "
      )
    )
  ),
  
  # App title
  # titlePanel("Paint Expression of Spots"),
  
  # Sidebar layout with input and output definitions
  # sidebarLayout(
  div(
    # Sidebar panel for inputs
    # sidebarPanel(
    wellPanel(
      tags$style(type="text/css", '#leftPanel { width:200px; float:left;}'),
      id = "leftPanel",
      
      # selectizeInput(inputId = "gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 1000)),
      
      selectizeInput(inputId = "all_gene_cichlid", label = "Cichlid Gene Name", choices = NULL, selected = "celsr1a", multiple = TRUE, options = list(maxOptions = 1000)),
      textOutput(outputId = "closest_human_gene"),
      selectizeInput(inputId = "all_gene_human", label = "Human Gene Name", choices = NULL, selected = "Celsr1", multiple = TRUE, options = list(maxOptions = 1000)),
      conditionalPanel(condition = "output.show_otho_cond_panel", id = "ortho.cond.panel",
                       # radioButtons(inputId = "cichlid.ortho", label = "MZ Orthologs",
                       checkboxGroupInput(inputId = "cichlid.ortho", label = "MZ Orthologs",
                                          choices = "junk",
                                          selected = "junk",
                                          inline = T),
                       tags$head(tags$style('#info{background-color: white;
                                    font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
                                    # border-style: solid;
                                    # border-wdith: thin;
                                    border-radius: 3px;
                                   }'))
      ),
      
      checkboxGroupInput(inputId = "samples", label = "Samples", 
                         choiceNames  = str_to_title(names(sample_pt_size)),
                         choiceValues = names(sample_pt_size),
                         selected = names(sample_pt_size), 
                         inline = T),
      checkboxGroupInput(inputId = "umap_add", label = "Disply w/ UMAP", 
                         choiceNames  = c("Comp. by Cluster", "Comp. by Sample"),
                         choiceValues = c("cluster_comp", "sample_comp"),
                         selected = c("cluster_comp", "sample_comp"), 
                         inline = F),
      sliderInput(inputId = "my.ncol", "Number of Columns", value = 4, min = 1, max = 13),
      sliderInput(inputId = "pt.size.multiplier", "Point Size", value = 0.5, min = 0, max = 1),
      checkboxInput(inputId = "toInfo", label = "Display Gene Info", value = T, width = NULL),
      
      conditionalPanel(condition = "input.toInfo && (input.all_gene_cichlid != null || input.all_gene_human != null)",
                       strong("Gene Info"),
                       textOutput(outputId = "info", container = pre),
                       tags$head(tags$style('#info{background-color: white;
                                    font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
                                    # border-style: solid;
                                    # border-wdith: thin;
                                    border-radius: 3px;
                                   }'))
      )
      
    ),
    
    # Main panel for displaying outputs
    # mainPanel(
    div(
      style = "flex-grow:1; resize:horizontal; overflow: hidden",
      tabsetPanel(id = "tabs", type = "tabs",
                  tabPanel("Spatial", value="sp_plot", plotOutput("sp_plot", width = "100%", height="90vh") %>% withSpinner(color="#FFA500")),
                  tabPanel("UMAP", value="umap_plot", plotOutput("umap_plot", width = "100%", height="90vh") %>% withSpinner(color="#FFA500"))
                  # tabPanel("UMAP w/ Info", value="umap_info_plot", plotOutput("umap_info_plot", width = "100%", height="100vh") %>% withSpinner(color="#FFA500"))
                  # tabPanel("By Fish", value="fishplot", plotOutput("fishplot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  # tabPanel("BHVE vs CTRL", value="bcplot", plotOutput("bcplot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  # tabPanel("BarPlot Summary", value="barplot", plotOutput("barplot", width = "100%", height="500px")  %>% withSpinner(color="#0dc5c1")),
                  # tabPanel("BarPlot %", value="barplotAvg", plotOutput("barplotAvg", width = "100%", height="500px")  %>% withSpinner(color="#0dc5c1")),
                  # tabPanel("Summary", value="summary", DT::dataTableOutput("summary") %>% withSpinner(color="#0dc5c1"))
      ),
      downloadButton(outputId = "down", label = "Download the plot")
      
      # conditionalPanel(condition = 'input.tabs == "summary || input.tabs == "barplot || input.tabs == "barplotAvg',
      #                  downloadButton(outputId = "sumDown", label = "Download Summary")
      # )
      
    )
  )
)

# Server =======================================================================
# Define server logic
server = function(input, output, session) {
  # updateSelectizeInput(session, "gene", choices = all_genes, server = TRUE)
  updateSelectizeInput(session, "all_gene_human",   choices = human_gene_names, server = TRUE)
  updateSelectizeInput(session, "all_gene_cichlid", choices = gene_names, server = TRUE)
  
  output$show_otho_cond_panel = reactive({ length(input$all_gene_human) > 0 })
  outputOptions(output, "show_otho_cond_panel", suspendWhenHidden = FALSE)
  
  # Update num MZ orthologs choices based on input Human gene - reactive
  observeEvent(input$all_gene_human, {
    if (length(input$all_gene_human) > 0) {
      myValues = gene_info$mzebra[which(gene_info$human == input$all_gene_human)]
      # myENS = gene_info$ens[which(gene_info$human == input$all_gene_human)]
      # myENS[which(startsWith(myENS, "ENSMZEG"))] = NA
      # myENS[which(!is.na(myENS))] = paste0(" (", myENS[which(!is.na(myENS))], ")")
      # myENS[which(is.na(myENS))] = ""
      this.symbol = gene_info$nd_symbol[which(gene_info$human == input$all_gene_human)]
      this.symbol[which(startsWith(this.symbol, "LOC") | startsWith(this.symbol, "zgc:") | startsWith(this.symbol, "si:"))] = NA
      this.symbol[which(!is.na(this.symbol))] = paste0(" (", this.symbol[which(!is.na(this.symbol))], ")")
      this.symbol[which(is.na(this.symbol))] = ""
      pat.mart.agree = gene_info$human_pat[which(gene_info$human == input$all_gene_human)] == gene_info$human_mart[which(gene_info$human == input$all_gene_human)]
      myNames = paste0(myValues, this.symbol)
      myNames[which(pat.mart.agree)] = paste0(myNames[which(pat.mart.agree)], " ✓")
      # updateRadioButtons(session = getDefaultReactiveDomain(), inputId = "cichlid.ortho", choiceNames = myNames, choiceValues = myValues, selected = findClosestCichlid(input$all_gene_human))
      updateCheckboxGroupInput(session = getDefaultReactiveDomain(), inputId = "cichlid.ortho", choiceNames = myNames, choiceValues = myValues, selected = findClosestCichlid(input$all_gene_human))
    }
  }, ignoreNULL = FALSE, ignoreInit = T)
  
  findClosestCichlid = function(human_gene) {
    this_rows = gene_info[which(gene_info$human == human_gene),]
    this_cichlid_genes = this_rows[, 1]
    if (length(this_cichlid_genes) > 1) {
      starts_with_hgnc = startsWith(tolower(this_cichlid_genes), tolower(human_gene))
      print(starts_with_hgnc)
      if (length(which(starts_with_hgnc)) > 0) {
        cichlid_gene <- this_cichlid_genes[which(starts_with_hgnc)[1]]
      } else {
        ens_starts_with_hgnc = startsWith(tolower(this_rows$ens), tolower(human_gene))
        print(ens_starts_with_hgnc)
        if (length(which(ens_starts_with_hgnc)) > 0) {
          cichlid_gene <- this_cichlid_genes[which(ens_starts_with_hgnc)[1]]
        } else {
          cichlid_gene <- this_cichlid_genes[1]
        }
      } # end bad multiple
    } else {
      cichlid_gene <- this_cichlid_genes
    }
    
    return(cichlid_gene)
  }
  
  findClosestHuman = function(cichlid_gene) {
    this_human_genes = gene_info[which(gene_info$mzebra == cichlid_gene),2]
    if (length(this_human_genes) > 1) {
      upper_human_gene <- this_human_genes[which(startsWith(tolower(this_human_genes), tolower(cichlid_gene)))]
      if (length(upper_human_gene) == 1) {
        human_gene <- upper_mouse_gene
      } else {
        human_gene <- this_human_genes[1]
      } # end bad multiple
    } else {
      human_gene <- this_human_genes
    }
    
    return(human_gene)
  }
  
  # oldCreateSpatialPlot = function(gene, split, samples) {
  #   # obj@active.assay <- "SCT"
  #   if (split == "cond") {
  #     print("Not yet implemented.")
  #   } else {
  #     this.obj <- obj[,which(obj$sample %in% samples)]
  #     plist = SpatialFeaturePlot(this.obj, images = samples, features = gene, ncol = 4, combine = F)
  #     names(plist) = samples
  #     plist2 = list()
  #     gene.values = FetchData(object = this.obj, vars = gene)
  #     max.val = max(gene.values)
  #     min.val = min(gene.values)
  #     for (s in names(plist)) { plist[[s]]$layers[[1]]$aes_params$point.size.factor = sample_pt_size[s]*(input$pt.size.multiplier*2) }
  #     for (s in names(plist)) { plist2[[s]] = plist[[s]] + ggplot2::scale_fill_gradientn(colors = rev(brewer.pal(11, "Spectral")), limits = c(min.val, max.val)) }
  #     wrap_plots(plist2, ncol = input$my.ncol, guides = "collect") & theme(legend.position = "bottom")
  #     # & plot_annotation(theme = theme(plot.background = element_rect(fill ="#b5b5b5")))
  #   }
  # }
  # 
  # createSpatialPlot = function(gene, split, samples) {
  #   # obj@active.assay <- "SCT"
  #   if (split == "cond") {
  #     print("Not yet implemented.")
  #   }
  #   if (all(levels(obj$sample)%in% samples)) {
  #     allSamplesSFP(obj, gene, "Spatial", "data", pt.size.multiplier = input$pt.size.multiplier*2)
  #   }
  # }
  
  createUmapPlot = function(gene, split, samples, cluster_comp = F, sample_comp = F) {
    
    this.obj <- obj[,which(obj$sample %in% samples)]
    plist = list()
    
    if (length(gene) > 1) {
      for (f in gene) {
        if (! f %in% rownames(this.obj@assays$Spatial@data)) {
          message(paste0("Error. Feature, ", f, ", was not found in the this.object."))
          return (NULL)
        }
      }
      this.obj$all_value_sum = colSums(this.obj@assays$Spatial@data[c(gene[1], gene[2], gene[2]),] > 0)
      this.obj$all_value_sum = plyr::revalue(as.character(this.obj$all_value_sum), replace = c("0" = "none", "1" = gene[1], "2" = gene[2], "3" = "both"))
      this.obj$all_value_sum = factor(this.obj$all_value_sum, levels = c("both", gene, "none"))
      
      this.df = data.frame(this.obj@reductions$umap@cell.embeddings[,1:2])
      this.df$all_value_sum = this.obj$all_value_sum
      this.df = this.df[order(this.df$all_value_sum, decreasing = T),]
      plist[[1]] = wrap_elements(plot = ggplot(this.df, aes(x = UMAP_1, y = UMAP_2, color = all_value_sum)) + geom_point(size = 0.9*(input$pt.size.multiplier*2)) + scale_color_manual(values = c(discrete.colors[1:(length(discrete.colors)-1)], "gray80"), name = NULL) + theme_void() + coord_fixed() + ggtitle(paste(gene, collapse = " + ")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
      if (cluster_comp) {
        cluster.df = data.frame(table(this.obj$cluster))
        cluster.df$pos = data.frame(table(this.obj$cluster[which( this.obj$all_value_sum == "both" )]))[,2]
      }
      if (sample_comp) {
        sample.df = data.frame(table(this.obj$sample))
        sample.df$pos = data.frame(table(this.obj$sample[which( this.obj$all_value_sum == "both" )]))[,2]
      }
    } else {
      plist[[1]] = FeaturePlot(this.obj, features = gene, order = T, pt.size = 0.9*(input$pt.size.multiplier*2)) + theme_void() + coord_fixed() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      if (cluster_comp) {
        cluster.df = data.frame(table(this.obj$cluster))
        cluster.df$pos = data.frame(table(this.obj$cluster[which( this.obj@assays$Spatial@data[gene,] > 0 )]))[,2]
      }
      if (sample_comp) {
        sample.df = data.frame(table(this.obj$sample))
        sample.df$pos = data.frame(table(this.obj$sample[which( this.obj@assays$Spatial@data[gene,] > 0 )]))[,2]
      }
    }
    
    if (cluster_comp) {
      cluster.df$pct = (cluster.df$pos/cluster.df$Freq) * 100
      cluster.df$round_pct = format(round(cluster.df$pct, 1), nsmall = 1)
      cluster.df$round_pct[which(cluster.df$pct == 0)] = 0
      plist[[2]] = ggplot(cluster.df, aes(x = Var1, y = pct, fill = Var1)) + geom_bar(stat  = "identity") + xlab("Cluster") + ylab("% of Spots in Cluster") + scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + ggtitle("Comp. by Cluster") + theme_classic() + NoLegend() + geom_text(aes(label=round_pct), vjust=1.25, hjust = 0.5, color="black", position = position_dodge(0.9), size=3)
    }
    if (sample_comp) {
      sample.df$pct = (sample.df$pos/sample.df$Freq) * 100
      sample.df$round_pct = format(round(sample.df$pct, 1), nsmall = 1)
      sample.df$round_pct[which(sample.df$pct == 0)] = 0
      plist[[length(plist)+1]] = ggplot(sample.df, aes(x = Var1, y = pct, fill = Var1)) + geom_bar(stat  = "identity") + xlab("Sample") + ylab("% of Spots in Sample") + scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + ggtitle("Comp. by Sample") + theme_classic() + NoLegend() + geom_text(aes(label=round_pct), vjust=1.25, hjust = 0.5, color="black", position = position_dodge(0.9), size=3)
    }

    if (! (cluster_comp && sample_comp)) {
      wrap_plots(plist, ncol = 1)
    } else {
      # (plist[[1]] | (plist[[2]] / plist[[3]]))
      wrap_plots(plist, ncol = 1)
    }
    
  }
  
  geneInfo = function(gene) {
    str = ""
    if (length(gene) > 0) {
      if (gene %in% gene_names) {
        human = gene_info$human[which(gene_info$mzebra == gene)]
        human_description = gene_info$human_description[which(gene_info$mzebra == gene)]
        mzebra_description = gene_info$mzebra_description[which(gene_info$mzebra == gene)] 
        str = paste0(str, "\nClosest Human Gene:\n", human, "\n")  
        str = paste0(str, "\nHuman Description:\n\"", human_description, '"\n')
        str = paste0(str, "\nMZ Description:\n\"", mzebra_description, '"\n') 
        
      } else {
        human_description = gene_info$human_description[which(gene_info$human == gene)]
        if (length(human_description) > 1) {
          str = paste0(str, "\n# of Close MZ Genes:\n", length(human_description), "\n")
        } else {
          mzebra = gene_info$mzebra[which(gene_info$human == gene)]
          mzebra_label = gene_info$label[which(gene_info$human == gene)]
          mzebra_ens = gene_info$ens[which(gene_info$human == gene)]
          mzebra_description = gene_info$mzebra_description[which(gene_info$human == gene)]
          str = paste0(str, "Closest MZ Gene:\n", mzebra, "\n")
          str = paste0(str, "\nClosest MZ Label:\n", mzebra_label, "\n")
          str = paste0(str, "\nClosest MZ ENS:\n", mzebra_ens, "\n")
          str = paste0(str, "\nMZ Description:\n\"", mzebra_description, '"\n')    
        }
        str = paste0(str, "\nHuman Description:\n\"", human_description[1], '"\n')
      }
    }
    return(str)
  }
  
  output$info = renderText({
    
    if (length(input$all_gene_cichlid) > 0) {
      geneInfo(input$all_gene_cichlid)
    } else if (length(input$all_gene_human) > 0) {
      gene <- input$all_gene_cichlid
      geneInfo(input$all_gene_human)
    }
    
  })
  
  output$closest_cichlid_gene <- renderText({
    if (length(input$all_gene_human) > 0) {
      str <- "Closest Cichlid Gene: "
      str <- paste(str, as.character(findClosestCichlid(input$all_gene_human)))
      return(str)
    }
  })
  
  output$closest_human_gene <- renderText({
    if (length(input$all_gene_cichlid) > 0) {
      str <- "Closest Human Gene: "
      str <- paste(str, as.character(findClosestHuman(input$all_gene_cichlid)))
      return(str)
    }
  })
  
  # Normal Plot (All Samples)
  output$sp_plot <- renderPlot({
    
    if (length(input$all_gene_cichlid) > 0 || length(input$all_gene_cichlid) == 0 && length(input$all_gene_human) > 0 && length(input$cichlid.ortho) > 0 && input$cichlid.ortho %in% gene_info$mzebra[which(gene_info$human == input$all_gene_human)]) {
      if (length(input$all_gene_cichlid) > 0) {
        gene = input$all_gene_cichlid
      } else {
        print(input$cichlid.ortho)
        gene = input$cichlid.ortho
      }

      if (length(gene) == 1) {
        suppressMessages(allSamplesSFP(obj, gene, "Spatial", "data", pt.size.multiplier = input$pt.size.multiplier*2))
      } else {
        suppressMessages(allSamplesSFP(obj, gene, "Spatial", "data", pt.size.multiplier = input$pt.size.multiplier*2, pal = discrete.colors))
      }
    }

  })
  
  output$umap_plot <- renderPlot({
    if (length(input$all_gene_cichlid) > 0) {
      gene <- input$all_gene_cichlid
      suppressMessages(createUmapPlot(gene, "split", input$samples, cluster_comp = "cluster_comp" %in% input$umap_add, sample_comp = "sample_comp" %in% input$umap_add))
    } else if (length(input$all_gene_cichlid) == 0 && length(input$all_gene_human) > 0 && length(input$cichlid.ortho) > 0 && input$cichlid.ortho %in% gene_info$mzebra[which(gene_info$human == input$all_gene_human)]) {
      print(input$cichlid.ortho)
      gene = input$cichlid.ortho
      suppressMessages(createUmapPlot(gene, "split", input$samples, cluster_comp = "cluster_comp" %in% input$umap_add, sample_comp = "sample_comp" %in% input$umap_add))
    }
  })
  
  # Plots cells that are overlapping if >1 genes provided (BHVE vs CTRL)
  output$bc_ovlp_plot <- renderPlot({
    gene = geneParser()
    if (length(gene) > 0) {
      createOvlpPlot(gene, "cond")
    }
  })
  
  output$barplot <- renderPlot({
    gene = geneParser()
    if (length(gene) > 0) {
      createBarPlot(gene, FALSE)
    }
  })
  
  output$barplotAvg <- renderPlot({
    gene = geneParser()
    if (length(gene) > 0) {
      createBarPlot(gene, TRUE)
    }
  })
  
  output$splClust <- renderPlot({
    if (length(input$cluster) > 0 || length(input$gene) > 0) {
      gene = geneParser()
      createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode, input$origClust)
    }
  })
  
  output$pntSplClust <- renderPlot({
    if (length(input$cluster) > 0 && length(input$gene) > 0) { 
      gene = geneParser()
      createSplFeaturePlot(gene, input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
    }
  })
  
  output$summary <- DT::renderDataTable({
    gene = geneParser()
    if (length(gene) > 0) {
      DT::datatable(createSummary(gene))
    }
  })
  
  # Download Functionality
  output$sumDown <- downloadHandler(
    filename = function() {
      gene = geneParser()
      paste0(paste0(c("summary", input$cluster, gene), collapse = "_"), ".tsv")
    }, 
    content = function(filename) {
      gene = geneParser()
      summary <- createSummary(gene)
      write.table(summary, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )
  
  output$degClustCond <- downloadHandler(
    filename = function() {
      gene = geneParser()
      paste0(paste0(c("degByClustCond", input$cluster, gene), collapse = "_"), ".tsv")
    }, 
    content = function(filename) {
      if (input$tabs == "splClust" || input$tabs == "pntSplClust") {
        obj <- createSpltObj(gene, input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
        deg <- data.frame()
        obj$clust.cond <- paste0(obj$seurat_clusters, obj$cond)
        obj_num_clusters <- as.numeric(tail(levels(obj@meta.data$seurat_clusters), n=1))
        for (i in 0:obj_num_clusters) {
          newRow <- FindMarkers(obj, ident.1 = paste0("CTRL_", i), ident.2 = past0("BHVE_", i))
          newRow$gene <- rownames(newRow)
          newRow$cluster <- i
          deg <- rbind(deg, newRow)
        }
        write.table(deg, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
      }
      
    }
  )
  output$degClust <- downloadHandler(
    filename = function() {
      gene = geneParser()
      paste0(paste0(c("degByClust", input$cluster, gene), collapse = "_"), ".tsv")
    }, 
    content = function(filename) {
      if (input$tabs == "splClust" || input$tabs == "pntSplClust") {
        gene = geneParser()
        obj <- createSpltObj(gene, input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
        deg <- FindAllMarkers(obj)
        deg$gene <- rownames(deg)
        write.table(deg, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
      }
      
    }
  )
  
  output$down <- downloadHandler(
    filename = function() {
      gene = geneParser()
      if (length(gene) > 1) {
        paste(paste(gene, collapse = '_'), ".png", sep="")
      } else {
        paste(gene, ".png", sep="")
      }
    }, 
    content = function(filename) {
      print(filename)
      print(input$gene)
      if (input$tabs == "bc_ovlp_plot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createOvlpPlot(gene, "cond")
      } else if (input$tabs == "b1c1_ovlp_plot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createOvlpPlot(gene, "sample")
      } else if (input$tabs == "allplot") {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot(gene, "sample", c("b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5"))
      } else if (input$tabs == "bcplot")  {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot(gene, "cond", c())
      } else if (input$tabs == "barplot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createBarPlot(gene, FALSE)
      } else if (input$tabs == "barplotAvg") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createBarPlot(gene, TRUE)
      } else if (intput$tabs == "splClust") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode, input$origClust)
      } else if (intput$tabs == "pntSplClust") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
      }
      print(p)
      dev.off()
    }
  )
  
}

shinyApp(ui = ui, server = server)