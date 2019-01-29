# R Package : pannot

The pannot package provides tools to annotate and analyze annotation enrichment within sets of proteins or genes

Installation
---
Install the package from the github repository using:
```
devtools::install_github("VoisinneG/pannot")
```

Usage
---
Retrieve annotations from the enrichR database "GO_Biological_Process_2018" for a set of genes.

```
library(pannot)
genes <- c("Itsn2","Eps15l1","Cbl","Cblb","Cltc1","Cd5","Cd6")
df <- get_annotations_enrichr(data = genes, dbs = "GO_Biological_Process_2018")
print(df)
```

Perform annotation enrichment analysis on a subset of the genes. 
```
idx_subset = which(df$names %in% c("Eps15l1", "Cblb","Cbl"))
res <- annotation_enrichment_analysis(df, idx_subset = idx_subset, sep=";")
print(res)
```

