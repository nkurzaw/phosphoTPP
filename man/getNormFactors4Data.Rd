% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{getNormFactors4Data}
\alias{getNormFactors4Data}
\title{Derive normalization factors for replicates}
\usage{
getNormFactors4Data(
  pep_tab,
  temperature_anno,
  n_rep,
  filter_criteria = list(rel_fc_127H_low = 0, rel_fc_129H_low = 0.4, rel_fc_129H_high =
    0.6, rel_fc_130H_low = 0, rel_fc_130H_high = 0.3, rel_fc_131L_low = 0,
    rel_fc_131L_high = 0.2)
)
}
\arguments{
\item{pep_tab}{data frame with information on qunatified peptides}

\item{temperature_anno}{data frame indication which TMT channel 
corresponds to which temperature}

\item{n_rep}{number of replicates peptides used for normalization should
be found in}

\item{filter_criteria}{list of filtering criteria that should be applied
to retrived a overlapping high quality dataset for retriving normalization
factors}
}
\description{
Derive normalization factors for replicates
}
