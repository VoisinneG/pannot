# import HPRD database

THPRD <- read.table(paste("~/ownCloud/++Work/++Research/Resources-Databases/HPRD/",
                          "HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt",
                          sep=""),
                    sep="\t",header = TRUE,quote="\"")



# save in ./R/sysdata.rda  ---------------------------------------------------------------------------------------

devtools::use_data(
  THPRD,
  pkg=".",
  internal = TRUE,
  overwrite = TRUE)

