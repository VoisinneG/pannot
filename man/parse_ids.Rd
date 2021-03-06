% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pannot.R
\name{parse_ids}
\alias{parse_ids}
\title{Parse protein ids.
This is useful to remove additionnal information such as isoform or entry name.}
\usage{
parse_ids(x, sep_split = ";", sep_secondary = c("|", "-"), sep_collapse = ";")
}
\arguments{
\item{x}{Character string with protein ids}

\item{sep_split}{Character separating different protein ids}

\item{sep_secondary}{Character separating UniProt entry from other infromation in protein ids}

\item{sep_collapse}{Character used to separate different protein ids after parsing}
}
\description{
Parse protein ids.
This is useful to remove additionnal information such as isoform or entry name.
}
\examples{
ids <- "A2AMW0|A2AMW0_MOUSE; P47757-2|CAPZB_MOUSE; P47757-4|CAPZB_MOUSE; Q3TVK4|Q3TVK4_MOUSE"
parse_ids(ids, sep = "; ", sep_secondary=c("|", "-"), sep_collapse = ";")

ids <- c("Q5SWU9|ACACA_MOUSE", "Q9ES52-2|SHIP1_MOUSE; Q9ES52|SHIP1_MOUSE", "Q8VDD5|MYH9_MOUSE")
parse_ids(ids, sep = "; ", sep_secondary=c("|", "-"), sep_collapse = ";")

}
