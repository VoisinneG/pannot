---
output:
  html_document: default
  pdf_document: default
---
# pannot

The pannot package provides tools to annotate and analyze annotation enrichment within sets of proteins or genes

Installation
---
Install the package from the github repository using:
```
devtools::install_github("VoisinneG/pannot")
```

Usage
---
Retrieve annotations from the enrichR database "GO_Biological_Process_2018" for genes Itsn2 and Eps15l1

```
library(pannot)
df <- get_annotations_enrichr(data = c("Itsn2","Eps15l1"), dbs = "GO_Biological_Process_2018")
print(df)
```

