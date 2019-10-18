\name{plot.QA}

\alias{plot.QA}

\title{
Plot quality metrics for geographic assignments
}

\description{
Plot the output from \code{QA}, including spatial precision, bias, sensitivity and odds ratio of known locations for validation samples.
}

\usage{
  \method{plot}{QA}(x, ..., savePNG = FALSE)
}

\arguments{
  \item{x}{
    One or more QA objects
  }
  \item{...}{
    Other arguments to be passed to plot
  }
  \item{savePNG}{
    logical. Save the plots as PNG file.
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
d1 = subOrigData(taxon = "Lanius ludovicianus")
d2 = subOrigData(taxon = "Buteo lagopus")

\dontrun{# run quality assessment based hydrogen isotope from precipitation and known-origin bird
qa1 = QA(isoscape = d2h_world, known = d1, valiStation = 2, 
                    valiTime = 5, mask = naMap, name = "Lanius")
                    
qa2 = QA(isoscape = d2h_world, known = d2, valiStation = 2, 
                    valiTime = 5, mask = naMap, name = "Buteo")

# plot the QA result
plot(qa1, qa2)}
}
