renderVignette = function(){
  dh = readLines("Rmarkdown/dev_header.Rmd")
  ch = readLines("Rmarkdown/cran_header.Rmd")
  b = readLines("Rmarkdown/body.Rmd")
  
  da = append(dh, b)
  ca = append(ch, b)
  
  writeLines(da, "Rmarkdown/bound/assignR.Rmd")
  writeLines(ca, "vignettes/assignR.Rmd")
  
  rmarkdown::render(
    "Rmarkdown/bound/assignR.Rmd",
    output_file = "../../docs/index.html"
  )
  
  return(0)
}
