# Load libraries
library("stringr")
library("ggplot2")
library("Seurat")
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
# library("grobblR")
library("png")
library("grid")

# Functions
SpatialTransform <- function(data, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) {
  # Quick input argument checking
  if (!all(sapply(X = list(xlim, ylim), FUN = length) == 2)) {
    stop("'xlim' and 'ylim' must be two-length numeric vectors", call. = FALSE)
  }
  # Save original names
  df.names <- colnames(x = data)
  colnames(x = data)[1:2] <- c('x', 'y')
  # Rescale the X and Y values
  data <- transform_position(
    df = data,
    trans_x = function(df) {
      return(rescale(x = df, from = xlim))
    },
    trans_y = function(df) {
      return(rescale(x = df, from = ylim))
    }
  )
  # Something that ggplot2 does
  data <- transform_position(
    df = data,
    trans_x = squish_infinite,
    trans_y = squish_infinite
  )
  # Restore original names
  colnames(x = data) <- df.names
  return(data)
}

mySingleSFP = function(obj, feature, assay, slot) {
  # TODO: assays
  # TODO: metadata
  coords = GetTissueCoordinates(object = obj)
  coords$value = FetchData(object = obj, vars = feature, slot = slot)[,1]
  img.grob = GetImage(obj)
  # img.grob = grobblR::grob_image(paste0(data_dir, "images/small/V11M15-295_A1_for_brianna.png"))
  # img.grob = png::readPNG(source = paste0(data_dir, "images/small/V11M15-295_A1_for_brianna.png"))
  # raster.img.grob = grid::rasterGrob(img.grob)
  image = obj[["slice1"]]
  
  my.imagerow.min = min(image@coordinates$imagerow)
  my.imagecol.min = min(image@coordinates$imagecol)
  my.imagerow.max = max(image@coordinates$imagerow)
  my.imagecol.max = max(image@coordinates$imagecol)
  
  my.imagerow.min = min(coords$imagerow)
  my.imagecol.min = min(coords$imagecol)
  my.imagerow.max = max(coords$imagerow)
  my.imagecol.max = max(coords$imagecol)
  
  vp <- viewport(
    x = image@coordinates$imagerow,
    y = image@coordinates$imagecol,
    width = unit(x = , units = "npc"),
    height = unit(x = 1, units = "npc"),
    just = c("left", "bottom")
  )  
  vp <- viewport(
    # x = unit(x = 0.5, units = "npc"),
    # y = unit(x = 0.5, units = "npc"),
    width = unit(x = 0.5, units = "npc"),
    height = unit(x = 0.5, units = "npc"),
    just = c("left", "bottom")
  )
  
  # # z = data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image)))
  # z = SpatialTransform(data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))), c())
  # z$y <- -rev(z$y) + 1
  # wdth <- z$x[2] - z$x[1]
  # hgth <- z$y[2] - z$y[1]
  # vp <- viewport(
  #   x = unit(x = z$x[1], units = "npc"),
  #   y = unit(x = z$y[1], units = "npc"),
  #   width = unit(x = wdth, units = "npc"),
  #   height = unit(x = hgth, units = "npc"),
  #   just = c("left", "bottom")
  # )
  zoom.out.factor = 0.01
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
  
  coords2 = coords
  coords2$px = coords2$imagecol
  coords2$py = -coords2$imagerow
  
  # ggplot(coords, aes(x=imagecol, y=-imagerow, color = value, fill = value)) + annotation_custom(img.grob, xmin = min(coords$imagecol)-my_left, ymin = max(-coords$imagerow)+my_bot, xmax = max(coords$imagecol)+my_right, ymax = min(-coords$imagerow)-my_top) + geom_point() + scale_fill_gradientn(colors=viridis(100)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  # ggplot(coords, aes(x=imagecol, y=-imagerow, color = value)) + annotation_custom(img.grob.test.grob) + geom_point() + scale_color_gradientn(colors=viridis(100)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  ggplot(coords, aes(x=imagecol, y=-imagerow, color = value)) + annotation_custom(img.grob.test.grob) + geom_point() + scale_color_gradientn(colors=viridis(100)) + scale_x_continuous(expand=c(0,0), limits = c(my.y.min, my.y.max)) + scale_y_continuous(expand=c(0,0), limits = c(-my.x.max, -my.x.min))
  # ggplot(coords, aes(x=imagecol, y=-imagerow, color = value, fill = value)) + annotation_custom(img.grob, xmin = my.imagerow.min, xmax = my.imagerow.max, ymin = my.imagecol.min, ymax = my.imagecol.max) + geom_point() + scale_fill_gradientn(colors=viridis(100)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
}
