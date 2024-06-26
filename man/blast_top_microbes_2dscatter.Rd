% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blastanalyse.R
\name{blast_top_microbes_2dscatter}
\alias{blast_top_microbes_2dscatter}
\title{Top Hits 2d-scatter}
\usage{
blast_top_microbes_2dscatter(
  blast,
  annotate_strong_hits = TRUE,
  zoomed_out = FALSE,
  legend_hide = FALSE,
  size = 1,
  svg_width = 6,
  svg_height = 5
)
}
\arguments{
\item{blast}{blast object. See \code{\link[=blast_parse]{blast_parse()}}}

\item{annotate_strong_hits}{Draw ggforce annotation box around strong hits (boolean)}
}
\value{
ggiraph object
}
\description{
Top Hits 2d-scatter
}
