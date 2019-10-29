\name{plot.QA}

\alias{plot.QA}

\title{
Plot quality metrics for geographic assignments
}

\description{
Plot the output from \code{QA}, including spatial precision, bias, sensitivity and odds ratio of known locations for validation samples.
}

\usage{
  \method{plot}{QA}(x, ..., outDir = NULL)
}

\arguments{
  \item{x}{
    One or more QA objects
  }
  \item{...}{
    Other arguments to be passed to plot
  }
  \item{outDir}{
  character string. Directory path to which output will be saved. If NULL no files are written.
  }
}

\references{
Ma et al. (in review) Does transpiration matter? Comparing geographic assignment with precipitation- and plant-based isoscapes using IsoMAP and assignR software. \emph{Movement Ecology}.

Vander Zanden, H. B. et. al (2014) Contrasting assignment of migratory organisms to geographic origins using long-term versus year-specific precipitation isotope maps. \emph{Methods in Ecology and Evolution} \strong{5} 891--900.
}

\seealso{
\code{\link[assignR]{QA}}
}

\examples{
# load data
data("naMap") # North America 
data("d2h_world") # precipitation hydrogen isotope of the world
data("knownOrig") # isotopes of known-origin samples

# extract some known-origin data
d1 = subOrigData(taxon = "Charadrius montanus")
d2 = subOrigData(taxon = "Buteo lagopus")

# run quality assessment based on precipitation hydrogen isotopes and 
# known-origin birds; small values of valiStation and valiTime and aggregated
# grid used in example to reduce run time

# first with one example
qa1 = QA(isoscape = raster::aggregate(d2h_world, 6), known = d1, 
          valiStation = 1, valiTime = 2, mask = naMap, name = "Charadrius")
                    
# plot the qa result
plot(qa1)

# now compare
\donttest{qa2 = QA(isoscape = raster::aggregate(d2h_world, 6), known = d2, 
          valiStation = 1, valiTime = 2, mask = naMap, name = "Buteo")
plot(qa1, qa2)}
}
