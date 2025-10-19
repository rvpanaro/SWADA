ggsavepdf <- function(pngpath, ...){
  
  library(ggplot2)
  ggsave(pngpath, ...)
  
  library(magick)
  img <- image_read(pngpath)
  
  out <- sub("\\.png$", ".pdf", pngpath)
  image_write(img, path = out, format = "pdf")
  
  file.remove(pngpath)  
}
