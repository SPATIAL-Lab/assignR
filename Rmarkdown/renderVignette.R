rmarkdown::render(
  "Rmarkdown/assignR.Rmd",
  params = list(date = format(Sys.Date(), "%B %d, %Y")),
  output_file = "../docs/index.html"
)
