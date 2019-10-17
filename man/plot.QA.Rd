\name{plot.QA}
\alias{plot.QA}
\title{
Plot quality assessment metrics for geographic assignment.
}
\description{
Plot the output from QA results. Spatial precision, accuracy, sensitivity and odds ratio of known locations for validation samples.
}
\usage{
  \method{plot}{QA}(x, ..., savePNG = FALSE)
}
\arguments{
  \item{x}{
    Either a single QA object or a list of QA objects
  }
  \item{...}{
    Other arguments to be passed to plot
  }
  \item{savePNG}{
    logical. Save the plots as PNG file (TRUE) or not (FALSE [default]).
  }
}
\references{
Ma et al., 2019

Vander Zanden, H.B., Wunder, M.B., Hobson, K.A., Van Wilgenburg, S.L., Wassenaar, L.I., Welker, J.M. and Bowen, G.J., 2014. Contrasting assignment of migratory organisms to geographic origins using long-term versus year-specific precipitation isotope maps. Methods in Ecology and Evolution, 5(9), pp.891-900.
}

\seealso{
\code{\link[ggplot2]{ggplot}}
\code{\link[assignR]{QA}}
}

\examples{
# load data
data("naMap") # North America 
data("d2h_world") # precipitation hydrogen isotope of the world
data("knownOrig") # isotopes of known-origin samples

# extract some known-origin data
bird_d2h = subOrigData(taxon = "Lanius ludovicianus")

\dontrun{# run quality assessment based hydrogen isotope from precipitation and known-origin bird
d2h_QA = QA(isoscape = d2h_world, known = bird_d2h, valiStation = 2, 
                    valiTime = 5, mask = naMap)
                    
# plot the QA result
plot(d2h_QA)}
}
