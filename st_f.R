# Load libraries
library("stringr")
library("ggplot2")
library("Seurat")
library("SeuratObject")
library("Matrix")
library("RColorBrewer")
library("viridis")
library("reshape2")
library("data.table")
library("colourvalues")
library("colorspace")
library("patchwork")
library("dplyr")
library("parallel")
library("grid")
library("ggpubr")

# Functions

mySingleSFP = function(obj = NULL, feature = NULL, assay = NULL, slot = NULL, coords = NULL, values = NULL, img.grob = NULL, points.as.text = F, rot.text = T, my.pt.size = 0.8, zoom.out.factor = 0.05, pal = colorRampPalette(colors = rev(brewer.pal(11, "Spectral"))), col.min = NULL, col.max = NULL, angle = NULL, discrete = F, rm.zero = F, col.ident = F, scale.alpha = F, doFlip = F, interactive = FALSE) {
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
  
  if (doFlip) { 
    coords$imagerow = nrow(img.grob$raster)-coords$imagerow
    img.grob.test.grob = img.grob
    img.grob.test.grob$raster = img.grob$raster[nrow(img.grob$raster):1,]
    img.grob = img.grob.test.grob
  }
  
  # Finding Coordinates of the Tissue
  my.imagerow.min = min(coords$imagerow)
  my.imagecol.min = min(coords$imagecol)
  my.imagerow.max = max(coords$imagerow)
  my.imagecol.max = max(coords$imagecol)
  myratio = (my.imagerow.max - my.imagerow.min) / (my.imagecol.max - my.imagecol.min) # locks the aspect ratio of the image to prevent distortion/elongation
  if (doFlip) { my.ratio = 20 }
    
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
  
  # Change factors to vector
  if (is.factor(coords$value)) { message("A factor was input. Now converting it to a numeric vector."); coords$value = as.numeric(as.vector(coords$value)) }
  
  # Set color range
  if (!discrete & is.null(col.min) && is.null(col.max)) { col.min = min(coords$value); col.max = max(coords$value); }
  
  # Remove spots with zero expression
  if (rm.zero) { coords = coords[which(coords$value != 0),] }
  
  if(scale.alpha) { scale.alpha = "value"; alpha.min = 0.1; } else { scale.alpha = "1"; alpha.min = 1; }
  
  # Plot
  if (points.as.text) {
    # This is for plotting the number of the cluster at each spot instead of plotting spots as points
    this.angle = 0
    if (rot.text) { this.angle = angle }
    p = ggplot(coords, aes_string(x="imagecol", y="-imagerow", color = "value", alpha = scale.alpha))   + annotation_custom(img.grob.test.grob) + geom_text(size = my.pt.size*0.8, fontface = "bold", angle = -this.angle, aes(label = value)) + scale_color_gradientn(colors=pal(100), limits = c(col.min, col.max)) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio) + scale_alpha_continuous(range = c(alpha.min, 1))
    # p = ggplot(coords, aes(x=imagecol, y=-imagerow, color = value))   + annotation_custom(img.grob.test.grob) + geom_point(shape = 1, size = my.pt.size, stroke = 0) + scale_color_gradientn(colors=pal(100), limits = c(col.min, col.max)) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio)
    # p = p + geom_text(size = my.pt.size*0.8, angle = -this.angle, aes(label = value))
  } else {
    if (discrete) {
      # Discrete Color Scales
      p = ggplot(coords, aes_string(x="imagecol", y="-imagerow", color = "value", alpha = scale.alpha)) + annotation_custom(img.grob.test.grob) + geom_point(size = my.pt.size, stroke = 0)            + scale_color_manual(values=pal, drop = F)                             + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio) + scale_alpha_continuous(range = c(alpha.min, 1))
    } else {
      # Continuous Color Scales
      p = ggplot(coords, aes_string(x="imagecol", y="-imagerow", color = "value", alpha = scale.alpha)) + annotation_custom(img.grob.test.grob) + geom_point(size = my.pt.size, stroke = 0)            + scale_color_gradientn(colors=pal(100), limits = c(col.min, col.max)) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio) + scale_alpha_continuous(range = c(alpha.min, 1))
      # coords$value = factor(coords$value, levels = sort(as.numeric(unique(coords$value))))
      # coords = coords[which(coords$value != 29),]
      # coords$spot = rownames(coords)
      # p = ggplot(coords, aes(x=imagecol, y=-imagerow, color = value)) + annotation_custom(img.grob.test.grob) + geom_point(size = my.pt.size, stroke = 0) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio) + ggforce::geom_mark_hull(aes(fill = value, group = spot), expand = unit(3, "mm"))
      # p = ggplot(coords, aes(x=imagecol, y=-imagerow, color = value)) + annotation_custom(img.grob.test.grob) + geom_point(size = my.pt.size, stroke = 0) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend() + theme(aspect.ratio = myratio) + ggforce::geom_voronoi_tile(aes(group = 1L), max.radius = 0.2)
    }
    if (col.ident) {
      # Identity Color Scale
      p = p + scale_color_identity() 
    }
  }
  
  # if (interactive) { p = ggplotly(ggplot(coords, aes_string(x="imagecol", y="-imagerow", color = "value")) + geom_point(size = my.pt.size, stroke = 0)            + scale_color_manual(values=pal, drop = F)                             + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + NoLegend()) %>% highlight("plotly_selected", dynamic = TRUE) }
  
  return(p)
}

myMultiSFP = function(obj, feature, samples = NULL, assay = "SCT", slot = "data", same.col.scale = T, points.as.text = F, rot.text = T, pt.size.multiplier = 1, zoom.out.factor = 0.05, pal = colorRampPalette(colors = rev(brewer.pal(11, "Spectral"))), rm.zero = F, col.ident = F, high.res = F, scale.alpha = F) {
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
  
  # Point Sizes
  sample_pt_size = c(2.3, 1.75, 2, 2.1, 1.4, 1.5, 1.5, 1.5, 1.25, 2.3, 2, 2, 1.3)
  names(sample_pt_size) = real.samples
  
  if (!is.null(samples)) { real.samples = real.samples[which(real.samples %in% samples)]; all.samples = all.samples[which(all.samples %in% samples)] }
  print(real.samples)
  
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
      if (! f %in% rownames(obj@assays$Spatial@counts)) {
        message(paste0("Error. Feature, ", f, ", was not found in the object."))
        return (NULL)
      }
    all_value_sum = colSums(obj@assays$Spatial@counts[c(feature[1], feature[2], feature[2]),] > 0)
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
  
  if (!same.col.scale) { min.val = NULL; max.val = NULL; }

  
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
    p_list[[s]] = mySingleSFP(coords = this.coords, assay = assay, slot = slot, values = this.values, img.grob = this.img.grob, points.as.text = points.as.text, rot.text = T, my.pt.size = sample_pt_size[s]*pt.size.multiplier, zoom.out.factor = 0.05, pal = pal, col.min = min.val, col.max = max.val, angle = this.angle, discrete = isDiscrete, rm.zero = rm.zero, col.ident = col.ident, scale.alpha = scale.alpha)
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
  this.ncol = 4
  # if (length(all.samples) < 4) { this.ncol = length(all.samples) }
  all.samples.mtx = matrix(all.samples, ncol = this.ncol, byrow = T)
  print(all.samples.mtx)
  grid.newpage()
  # pushViewport(viewport(width=1, height=1, xscale=c(0, nrow(all.samples.mtx)), yscale=c(0,ncol(all.samples.mtx))))
  pushViewport(viewport(width=1, height=1, xscale=c(0, ncol(all.samples.mtx)), yscale=c(0,nrow(all.samples.mtx))))
  for(i in 1:nrow(all.samples.mtx)){
    for(j in 1:ncol(all.samples.mtx)){
      s = as.character(all.samples.mtx[i,j])
      this.x = j-0.5 # coordinate for the position on the grid
      this.y = nrow(all.samples.mtx)-i+0.5 # coordinate for the position on the grid
      # print(paste0(s, ": x=", this.x, ", y=", this.y))
      vp = viewport(x=unit(this.x,"native"), y=unit(this.y,"native"), width=unit(1,"native"), height=unit(1,"native"), clip=T) # current cell of the grid
      pushViewport(vp) # starts the viewport window
      
      # Determine background color and reorientation-angle of the sample
      back.color.df = as.data.frame(c("c2a" = "#AEADAD", "c2b" = "#AEADAD", "c2c" = "#AEADAD", "c2d" = "#AEADAD",
                                      "b2a" = "#AEADAD", "b2b" = "#AEADAD", "b2c" = "#AEADAD", "b2d" = "#AEADAD",
                                      "c1a" = "#B6B6B5", "c1b" = "#B6B6B5", "c1c" = "#B6B6B5", "c1d" = "#B6B6B5",
                                      "b1a" = "#B4B3B2", "b1b" = "#B4B3B2", "b1c" = "#B4B3B2", "b1d" = "#B4B3B2"))
      zoom.df = as.data.frame(c("c2a" = 1, "c2b" = 0.9, "c2c" = 1, "c2d" = 0.95,
                                      "b2a" = 0.9, "b2b" = 1, "b2c" = 1, "b2d" = 1,
                                      "c1a" = 0.9, "c1b" = 1.05, "c1c" = 0.9, "c1d" = 1,
                                      "b1a" = 1, "b1b" = 1, "b1c" = 0.9, "b1d" = 1))
      colnames(back.color.df) = "back.color"
      colnames(zoom.df) = "zoom"
      
      this.zoom = zoom.df[s, "zoom"]
      this.back.color = back.color.df[s, "back.color"]
      this.angle = angle.df[s, "angle"]
      
      # Set the background color
      if (! high.res ) { grid.rect(gp=gpar(fill=this.back.color, col=this.back.color)) } # sets the background color
      # grid.rect(gp=gpar(fill=this.back.color, col=this.back.color)) 
      
      # Create a viewport within the current viewport (the one that represents the cell of the grid) 
      # This allows me to adjust the angle of the plots
      # vp2 = viewport(width = unit(0.9, "npc"), height = unit(0.9, "npc"), angle = this.angle, clip = T)
      vp2 = viewport(width = unit(this.zoom, "npc"), height = unit(this.zoom, "npc"), angle = this.angle, clip = T)
      pushViewport(vp2)
      if (s %in% real.samples) {
        grid.draw(ggplotGrob(p_list[[s]]))
      } else if (s == "b1a") {
        grid.text(as.character(feature), 0.5, 0.5, gp=gpar(cex=1.3))
      } else if (s == "b1b") {
        grid.draw(leg_p_grob)
      }
      popViewport() # pop back to the viewport that represents the cell of the grid
      
      # Sample Labels
      if (! s %in% c("b1a", "b1b", "b1d")) {
        grid.text(as.character(s), 0.1, 0.95, gp=gpar(cex=1))
      } 
      
      upViewport()
      
    }
  }
  popViewport(1)
}


myB2SFP = function(obj, feature, sample.halves = NULL, assay = "SCT", slot = "data", same.col.scale = T, points.as.text = F, rot.text = T, pt.size.multiplier = 1, zoom.out.factor = 0.05, pal = colorRampPalette(colors = rev(brewer.pal(11, "Spectral"))), rm.zero = F, col.ident = F, scale.alpha = F) {
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
  samples = sample.halves
  
  # Samples to Plot
  sample.df = data.frame(samples = c("b2al", "b2bl", "b2cl", "b2dl", "b2ar", "b2br", "b2cr", "b2dr"),
                         zoom1 =   c(0.05,    0.05,   0.1,    0.2,    0.05,   0.05,   0.1,    0.1),
                         angle =   c(90,      88,     82,     120,    90,     100,    115,    85),
                         zoom  =   c(1,       0.84,   0.8,    0.75,   0.98,   0.9,    0.9,    0.75),
                         pt.size = c(0.95,    1.15,   1.3,    1.2,    0.9,    0.85,   1.05,   1.10),
                         flip =    c(T,       T,      T,      T,      F,      F,      F,      F))
  sample.df$pt.size = sample.df$pt.size * 1.4
  rownames(sample.df) = sample.df$samples
  sample.df = sample.df[c("b2ar", "b2al", "b2br", "b2cr", "b2dr", "b2bl", "b2cl", "b2dl"),]
  
  if (!is.null(samples)) { sample.df = sample.df[which(sample.df$samples %in% samples)] }
  
  # Get all values for feature
  value_list = list()
  isDiscrete = F
  if (length(feature) == 1) {
    for (s in sample.df$samples) {
      # See if feature is a gene or metadata
      if (feature %in% rownames(obj[[assay]])) {
        value_list[[s]] = FetchData(object = obj, vars = feature, cells = colnames(obj)[which(obj$sh == s)], slot = slot)[,1]
      } else if (feature %in% colnames(obj@meta.data)) {
        value_list[[s]] = obj@meta.data[colnames(obj)[which(obj$sh == s)], feature]
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
      if (! f %in% rownames(obj@assays$Spatial@counts)) {
        message(paste0("Error. Feature, ", f, ", was not found in the object."))
        return (NULL)
      }
      all_value_sum = colSums(obj@assays$Spatial@counts[c(feature[1], feature[2], feature[2]),] > 0)
      all_value_sum = plyr::revalue(as.character(all_value_sum), replace = c("0" = "none", "1" = feature[1], "2" = feature[2], "3" = "both"))
      all_value_sum = factor(all_value_sum, levels = c("both", feature, "none"))
      for (s in sample.df$samples) {
        this.value = all_value_sum[which(obj$fish == "b2")]
        value_list[[s]] = this.value
      }
    }
    min.val = "dummy"
    max.val = "dummy"
  }
  
  if (!same.col.scale) { min.val = NULL; max.val = NULL; }
  
  # Create all the separate sample plots
  p_list = list()
  for (s in sample.df$samples) {
    this.coords = GetTissueCoordinates(object = obj, image = s)
    this.coords = this.coords[colnames(obj)[which(obj$sh == s)],]
    this.values = value_list[[s]]
    this.img.grob = GetImage(obj, image = s)
    this.angle = sample.df[s, "angle"]
    p_list[[s]] = mySingleSFP(coords = this.coords, assay = assay, slot = slot, 
                              values = this.values, img.grob = this.img.grob, 
                              points.as.text = points.as.text, rot.text = T, 
                              my.pt.size = sample.df[s, "pt.size"]*pt.size.multiplier, 
                              zoom.out.factor = sample.df[s, "zoom1"], pal = pal, 
                              col.min = min.val, col.max = max.val, angle = this.angle, 
                              discrete = isDiscrete, rm.zero = rm.zero, col.ident = col.ident, 
                              scale.alpha = scale.alpha, doFlip = sample.df[s, "flip"])
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
  this.ncol = 8
  # if (length(all.samples) < 4) { this.ncol = length(all.samples) }
  all.samples.mtx = matrix(sample.df$samples, ncol = this.ncol, byrow = T)
  print(all.samples.mtx)
  grid.newpage()
  # pushViewport(viewport(width=1, height=1, xscale=c(0, nrow(all.samples.mtx)), yscale=c(0,ncol(all.samples.mtx))))
  pushViewport(viewport(width=1, height=1, xscale=c(0, ncol(all.samples.mtx)), yscale=c(0,nrow(all.samples.mtx))))
  for(i in 1:nrow(all.samples.mtx)){
    for(j in 1:ncol(all.samples.mtx)){
      s = as.character(all.samples.mtx[i,j])
      this.x = j-0.5 # coordinate for the position on the grid
      this.y = nrow(all.samples.mtx)-i+0.5 # coordinate for the position on the grid
      # print(paste0(s, ": x=", this.x, ", y=", this.y))
      vp = viewport(x=unit(this.x,"native"), y=unit(this.y,"native"), width=unit(1,"native"), height=unit(1,"native"), clip=T) # current cell of the grid
      pushViewport(vp) # starts the viewport window
      
      # Determine background color and reorientation-angle of the sample
      this.zoom  = sample.df[s, "zoom"]
      this.angle = sample.df[s, "angle"]
      
      # Create a viewport within the current viewport (the one that represents the cell of the grid) 
      # This allows me to adjust the angle of the plots
      # vp2 = viewport(width = unit(0.9, "npc"), height = unit(0.9, "npc"), angle = this.angle, clip = T)
      vp2 = viewport(width = unit(this.zoom, "npc"), height = unit(this.zoom, "npc"), angle = this.angle, clip = T)
      pushViewport(vp2)
      if (s %in% sample.df$samples) {
        grid.draw(ggplotGrob(p_list[[s]]))
      } else if (s == "b1a") {
        grid.text(as.character(feature), 0.5, 0.5, gp=gpar(cex=1.3))
      } else if (s == "b1b") {
        grid.draw(leg_p_grob)
      }
      popViewport() # pop back to the viewport that represents the cell of the grid
      
      # Sample Labels
      if (! s %in% c("b1a", "b1b", "b1d")) {
        grid.text(as.character(s), 0.1, 0.95, gp=gpar(cex=1))
      } 
      
      upViewport()
      
    }
  }
  popViewport(1)
}
