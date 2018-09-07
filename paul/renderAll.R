## Script for convenient rendering of both html slides and notes from
## a single Rmd file.

library(rmarkdown)
library(revealjs)

formats <- c("html_document",
             #"revealjs::revealjs_presentation")
             "ioslides_presentation")
extensions <- c(".html","_slides.html")
renderAll <- function(filename, output_format="all", ...) {
  slidename <- sub(".Rmd$","_slides.Rmd",filename)
  if (!file.exists(slidename)) {
    system(sprintf("ln -s %s %s",filename, slidename))
  }
  if (output_format=="all") {
    #for (i in 1:length(formats)) {
    #  fmt  <- formats[i]
    #  ext <- extensions[i]
    slides <- FALSE
    render(filename,output_format="html_document")
    slides <- TRUE
    render(slidename,output_format="ioslides_presentation")
    }
  else {
      render(filename, output_format=output_format, ...)
  }
}
