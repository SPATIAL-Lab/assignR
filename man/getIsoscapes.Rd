\name{getIsoscapes}

\alias{getIsoscapes}

\title{
Download and unpack isoscapes from the web
}

\description{
This function retrieves gridded isotope maps from waterisotopes.org, unpacks the zip archives, and bundles the map layers as a RasterStack.
}

\usage{
getIsoscapes(isoType = "GlobalPrecipGS", timeout = 1200)
}

\arguments{
  \item{isoType}{
character string indicating which isoscapes are requested: see 'Details'.
}
  \item{timeout}{
numeric. Maximum alowed file download time, in seconds. Some isoscape archives exceed 2 GB in size and may require long download times on slow connections.
}
}

\details{
Accepted \code{isoType} values are:
\describe{
\item{"GlobalPrecipGS"}{Global growing-season precipitation H and O isotope values}
\item{"GlobalPrecipMA"}{Global mean-annual precipitation H and O isotope values}
\item{"GlobalPrecipMO"}{Global monthly precipitation H and O isotope values}
\item{"GlobalPrecipMA"}{Global mean-annual and monthly precipitation H and O isotope values}
\item{"USPrecipMA"}{High-resolution contiguous USA mean-annual precipitation H and O isotope values}
\item{"USPrecipMO"}{High-resolution contiguous USA monthly precipitation H and O isotope values}
\item{"USPrecipMA"}{High-resolution contiguous USA mean-annual and monthly precipitation H and O isotope values}
\item{"USSurf"}{High-resolution contiguous USA surface water H and O isotope values}
\item{"USTap"}{High-resolution contiguous USA surface water H and O isotope values}
\item{"USSr"}{High-resolution contiguous USA Sr isotope ratios}
\item{"CaribSr"}{High-resolution Sr isotope ratios for the circum-Caribbean region}
}
}

\value{
RasterStack containing the requested isoscape layers.
}

\references{
\url{https://wateriso.utah.edu/waterisotopes/pages/data_access/ArcGrids.html}
}

\examples{
\dontrun{iso = getIsoscapes("CaribSr")
}
}