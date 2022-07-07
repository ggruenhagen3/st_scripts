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
mySingleSFP = function(obj = NULL, feature = NULL, assay = NULL, slot = NULL, coords = NULL, values = NULL, img.grob = NULL, my.pt.size = 0.8, zoom.out.factor = 0.05, pal = colorRampPalette(colors = rev(brewer.pal(11, "Spectral"))), col.min = NULL, col.max = NULL) {
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
  #' 
  # TODO: assays
  # TODO: metadata
  
  # Input Checking
  if (is.null(obj)) {
    if (is.null(coords) || is.null(values) || is.null(img.grob)) { message("No object input, but coords, values, and/or img.grob were not specified."); return(NULL); }
    coords$value = values
  } else {
    if (is.null(feature) || is.null(assay) || is.null(slot)) { message("Object was input, but no feature, assay, and/or slot were specified.") }
    coords = GetTissueCoordinates(object = obj)
    coords$value = FetchData(object = obj, vars = feature, slot = slot)[,1]
    img.grob = GetImage(obj)
  }
  
  # Finding Cooridinates of the Tissue
  my.imagerow.min = min(coords$imagerow)
  my.imagecol.min = min(coords$imagecol)
  my.imagerow.max = max(coords$imagerow)
  my.imagecol.max = max(coords$imagecol)
  
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
  if (is.null(col.min) && is.null(col.max)) { col.min = min(coords$value); col.max = max(coords$value); }
  
  # Plot
  p = ggplot(coords, aes(x=imagecol, y=-imagerow, color = value)) + annotation_custom(img.grob.test.grob) + geom_point(size = my.pt.size) + scale_color_gradientn(colors=pal(100), limits = c(col.min, col.max)) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min)) + theme_void() + NoLegend()
  
  return(p)
}

allSamplesSFP = function(obj, feature, assay, slot, pt.size.multiplier = 1, zoom.out.factor = 0.05, pal = colorRampPalette(colors = rev(brewer.pal(11, "Spectral")))) {
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
  #' 
  #TODO assays
  #TODO metadata
  
  # Samples to Plot
  real.samples = c("c2a", "c2b", "c2c", "c2d", "b2a", "b2b", "b2c", "b2d", "c1a", "c1b", "c1c", "c1d", "b1c")
  all.samples = c(real.samples[1:(length(real.samples)-1)], "b1a", "b1b", "b1c", "b1d")
  
  # Get all values for feature
  value_list = list()
  for (s in real.samples) {
    value_list[[s]] = FetchData(object = obj, vars = feature, cells = colnames(obj)[which(obj$sample == s)], slot = slot)[,1]
  }
  min.val = min(unlist(value_list))
  max.val = max(unlist(value_list))
  
  # Point Sizes
  sample_pt_size = c(2.25, 1.5, 1.5, 2, 1.5, 1.5, 1.5, 2, 1.25, 2, 2, 2, 1.5)
  names(sample_pt_size) = real.samples
  
  # Create all the separate sample plots
  p_list = list()
  for (s in real.samples) {
    this.coords = GetTissueCoordinates(object = obj, image = s)
    this.values = value_list[[s]]
    this.img.grob = GetImage(obj, image = s)
    p_list[[s]] = mySingleSFP(coords = this.coords, values = this.values, img.grob = this.img.grob, my.pt.size = sample_pt_size[s]*pt.size.multiplier, zoom.out.factor = 0.05, col.min = min.val, col.max = max.val)
  }
  
  # Get color legend
  leg_p = ggplot(data.frame(a = 1, b = 1), aes(a, b, color = a)) + geom_point() + scale_color_gradientn(colors=pal(100), limits = c(min.val, max.val), name = NULL) + theme(legend.position = 'bottom', legend.background = element_blank())
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
      this.angle = switch(s, "c2a" = 155, "c2b" = 145, "c2c" = -115, "c2d" = 155,
                             "b2a" =  90, "b2b" =  95, "b2c" =  100, "b2d" =  85,
                             "c1a" =  90, "c1b" =  95, "c1c" =   98, "c1d" =  90,
                             "b1a" =   0, "b1b" =   0, "b1c" = -115, "b1d" =   0)
      
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


