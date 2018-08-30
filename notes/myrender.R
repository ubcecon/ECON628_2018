## Script for convenient rendering of both html slides and notes from
## a single Rmd file.

library(rmarkdown)
library(revealjs)

formats <- c("html_document",
             "revealjs::revealjs_presentation")
extensions <- c(".html","_slides.html")
myrender <- function(filename, output_format="all", ...) {
  if (output_format=="all") {
    for (i in 1:length(formats)) {
      fmt  <- formats[i]
      ext <- extensions[i]
      render(filename,output_format=fmt,
             output_file=sub(".Rmd$",ext,filename), ...)
    }
  } else {
      render(filename, output_format=output_format, ...)
  }
}
