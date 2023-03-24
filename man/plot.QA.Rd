\name{plot.QA}

\alias{plot.QA}

\title{
Plot quality metrics for geographic assignments
}

\description{
Plot the output from \code{\link{QA}}, including spatial precision, bias, sensitivity and odds ratio of known locations for validation samples.
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
Ma, C. et al. (2020) assignR : An R package for isotope-based geographic assignment. \emph{Methods in Ecology and Evolution} \strong{11} 996--1001. \doi{10.1111/2041-210X.13426}.

Vander Zanden, H. B. et al. (2014) Contrasting assignment of migratory organisms to geographic origins using long-term versus year-specific precipitation isotope maps. \emph{Methods in Ecology and Evolution} \strong{5} 891--900. \doi{10.1111/2041-210X.12229}
}

\seealso{
\code{\link[assignR]{QA}}
}

\examples{
library(terra)

# load North America boundary and global isoscape
data("naMap")
d2h_lrNA = rast(system.file("extdata/d2h_lrNA.tif", package = "assignR"))

# extract some known-origin data
d1 = subOrigData(taxon = "Buteo lagopus")

# run quality assessment based on precipitation hydrogen isotopes and 
# known-origin samples; small values of valiStation and valiTime 
# are used in example to reduce run time

# first with one example
# gives warning because a small number of samples are available
qa1 = QA(isoscape = d2h_lrNA, known = d1, valiStation = 1, 
          valiTime = 2, by = 10, mask = naMap, name = "Buteo")
                    
# plot the qa result
plot(qa1)

# now compare with a second data set
\donttest{d2 = subOrigData(taxon = "Charadrius montanus")
qa2 = QA(isoscape = d2h_lrNA, known = d2, valiStation = 1, 
          valiTime = 2, by = 10, mask = naMap, name = "Charadrius")
plot(qa1, qa2)}
}
