[![Travis-CI Build
Status](https://travis-ci.org/VoisinneG/pannot.svg?branch=master)](https://travis-ci.org/VoisinneG/pannot)

R Package : pannot
==================

The pannot package provides tools to annotate and analyze annotation
enrichment within sets of proteins or genes

Install
-------

Install the package from the github repository using:

    devtools::install_github("VoisinneG/pannot")
    library(pannot)

    ## Warning: replacing previous import 'IRanges::desc' by 'plyr::desc' when
    ## loading 'PSICQUIC'

Usage
-----

Retrieve annotations from the enrichR database
"GO\_Biological\_Process\_2018" for a set of genes:

    genes <- c("Itsn2","Eps15l1","Cbl","Cblb","Cltc1","Cd5","Cd6")
    df <- get_annotations_enrichr(data = genes, dbs = "GO_Biological_Process_2018")

    ## Uploading data to Enrichr... Done.
    ##   Querying GO_Biological_Process_2018... Done.
    ## Parsing results... Done.

    print(df)

    ##     names
    ## 1     Cbl
    ## 2    Cblb
    ## 3     Cd5
    ## 4     Cd6
    ## 5   Cltc1
    ## 6 Eps15l1
    ## 7   Itsn2
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 GO_Biological_Process_2018
    ## 1 negative regulation of epidermal growth factor-activated receptor activity (GO:0007175);negative regulation of epidermal growth factor receptor signaling pathway (GO:0042059);negative regulation of protein tyrosine kinase activity (GO:0061099);regulation of epidermal growth factor-activated receptor activity (GO:0007176);epidermal growth factor receptor signaling pathway (GO:0007173);negative regulation of ERBB signaling pathway (GO:1901185);negative regulation of receptor activity (GO:2000272);entry of bacterium into host cell (GO:0035635);regulation of epidermal growth factor receptor signaling pathway (GO:0042058);ERBB signaling pathway (GO:0038127);positive regulation of epidermal growth factor receptor signaling pathway (GO:0045742);cellular response to interleukin-6 (GO:0071354);interleukin-6-mediated signaling pathway (GO:0070102);positive regulation of receptor-mediated endocytosis (GO:0048260);positive regulation of ERBB signaling pathway (GO:1901186);regulation of receptor-mediated endocytosis (GO:0048259);entry into host cell (GO:0030260);cellular response to transforming growth factor beta stimulus (GO:0071560);transmembrane receptor protein serine/threonine kinase signaling pathway (GO:0007178);cellular response to fibroblast growth factor stimulus (GO:0044344);positive regulation of phosphatidylinositol 3-kinase signaling (GO:0014068);transforming growth factor beta receptor signaling pathway (GO:0007179);fibroblast growth factor receptor signaling pathway (GO:0008543);positive regulation of endocytosis (GO:0045807);regulation of phosphatidylinositol 3-kinase signaling (GO:0014066);modification-dependent protein catabolic process (GO:0019941);transmembrane receptor protein tyrosine kinase signaling pathway (GO:0007169);positive regulation of intracellular signal transduction (GO:1902533);negative regulation of programmed cell death (GO:0043069);regulation of apoptotic process (GO:0042981);ubiquitin-dependent protein catabolic process (GO:0006511);protein modification by small protein conjugation (GO:0032446);cellular response to cytokine stimulus (GO:0071345);cytokine-mediated signaling pathway (GO:0019221);protein ubiquitination (GO:0016567);negative regulation of apoptotic process (GO:0043066)
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      negative regulation of epidermal growth factor-activated receptor activity (GO:0007175);negative regulation of epidermal growth factor receptor signaling pathway (GO:0042059);negative regulation of protein tyrosine kinase activity (GO:0061099);regulation of epidermal growth factor-activated receptor activity (GO:0007176);epidermal growth factor receptor signaling pathway (GO:0007173);negative regulation of receptor activity (GO:2000272);ERBB signaling pathway (GO:0038127);NLS-bearing protein import into nucleus (GO:0006607);protein import into nucleus (GO:0006606);proteolysis (GO:0006508)
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       immunological synapse formation (GO:0001771);lipopolysaccharide-mediated signaling pathway (GO:0031663);positive regulation of cytokine production involved in inflammatory response (GO:1900017);regulation of cytokine production involved in inflammatory response (GO:1900015);heterophilic cell-cell adhesion via plasma membrane cell adhesion molecules (GO:0007157);acute inflammatory response (GO:0002526);regulation of T cell proliferation (GO:0042129);positive regulation of inflammatory response (GO:0050729);cellular response to lipopolysaccharide (GO:0071222);positive regulation of T cell proliferation (GO:0042102);positive regulation of lymphocyte proliferation (GO:0050671);positive regulation of cytokine production (GO:0001819);response to lipopolysaccharide (GO:0032496);positive regulation of T cell activation (GO:0050870);response to molecule of bacterial origin (GO:0002237);inflammatory response (GO:0006954);response to lipid (GO:0033993);cell-cell adhesion via plasma-membrane adhesion molecules (GO:0098742)
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          negative regulation of epidermal growth factor receptor signaling pathway (GO:0042059);negative regulation of ERBB signaling pathway (GO:1901185);regulation of epidermal growth factor receptor signaling pathway (GO:0042058)
    ## 7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          positive regulation of dendrite extension (GO:1903861);regulation of dendrite extension (GO:1903859);positive regulation of developmental growth (GO:0048639);positive regulation of cell growth (GO:0030307);positive regulation of neuron projection development (GO:0010976)

Perform annotation enrichment analysis on a subset of the genes:

    idx_subset = which(df$names %in% c("Eps15l1", "Cblb","Cbl"))
    res <- annotation_enrichment_analysis(df, idx_subset = idx_subset, sep=";", showProgress = FALSE)

    print(res)

    ##                                                                                  annot_terms
    ## 2     negative regulation of epidermal growth factor receptor signaling pathway _GO:0042059!
    ## 1    negative regulation of epidermal growth factor-activated receptor activity _GO:0007175!
    ## 3                       negative regulation of protein tyrosine kinase activity _GO:0061099!
    ## 4             regulation of epidermal growth factor-activated receptor activity _GO:0007176!
    ## 5                            epidermal growth factor receptor signaling pathway _GO:0007173!
    ## 6                                 negative regulation of ERBB signaling pathway _GO:1901185!
    ## 7                                      negative regulation of receptor activity _GO:2000272!
    ## 9              regulation of epidermal growth factor receptor signaling pathway _GO:0042058!
    ## 10                                                       ERBB signaling pathway _GO:0038127!
    ## 8                                             entry of bacterium into host cell _GO:0035635!
    ## 11    positive regulation of epidermal growth factor receptor signaling pathway _GO:0045742!
    ## 12                                           cellular response to interleukin-6 _GO:0071354!
    ## 13                                     interleukin-6-mediated signaling pathway _GO:0070102!
    ## 14                         positive regulation of receptor-mediated endocytosis _GO:0048260!
    ## 15                                positive regulation of ERBB signaling pathway _GO:1901186!
    ## 16                                  regulation of receptor-mediated endocytosis _GO:0048259!
    ## 17                                                         entry into host cell _GO:0030260!
    ## 18                cellular response to transforming growth factor beta stimulus _GO:0071560!
    ## 19     transmembrane receptor protein serine/threonine kinase signaling pathway _GO:0007178!
    ## 20                       cellular response to fibroblast growth factor stimulus _GO:0044344!
    ## 21               positive regulation of phosphatidylinositol 3-kinase signaling _GO:0014068!
    ## 22                   transforming growth factor beta receptor signaling pathway _GO:0007179!
    ## 23                          fibroblast growth factor receptor signaling pathway _GO:0008543!
    ## 24                                           positive regulation of endocytosis _GO:0045807!
    ## 25                        regulation of phosphatidylinositol 3-kinase signaling _GO:0014066!
    ## 26                             modification-dependent protein catabolic process _GO:0019941!
    ## 27             transmembrane receptor protein tyrosine kinase signaling pathway _GO:0007169!
    ## 28                     positive regulation of intracellular signal transduction _GO:1902533!
    ## 29                                 negative regulation of programmed cell death _GO:0043069!
    ## 30                                              regulation of apoptotic process _GO:0042981!
    ## 31                                ubiquitin-dependent protein catabolic process _GO:0006511!
    ## 32                            protein modification by small protein conjugation _GO:0032446!
    ## 33                                       cellular response to cytokine stimulus _GO:0071345!
    ## 34                                          cytokine-mediated signaling pathway _GO:0019221!
    ## 35                                                       protein ubiquitination _GO:0016567!
    ## 36                                     negative regulation of apoptotic process _GO:0043066!
    ## 37                                      NLS-bearing protein import into nucleus _GO:0006607!
    ## 38                                                  protein import into nucleus _GO:0006606!
    ## 39                                                                  proteolysis _GO:0006508!
    ## 41                                              immunological synapse formation _GO:0001771!
    ## 42                                lipopolysaccharide-mediated signaling pathway _GO:0031663!
    ## 43 positive regulation of cytokine production involved in inflammatory response _GO:1900017!
    ## 44          regulation of cytokine production involved in inflammatory response _GO:1900015!
    ## 45  heterophilic cell-cell adhesion via plasma membrane cell adhesion molecules _GO:0007157!
    ## 46                                                  acute inflammatory response _GO:0002526!
    ## 47                                           regulation of T cell proliferation _GO:0042129!
    ## 48                                 positive regulation of inflammatory response _GO:0050729!
    ## 49                                      cellular response to lipopolysaccharide _GO:0071222!
    ## 50                                  positive regulation of T cell proliferation _GO:0042102!
    ## 51                              positive regulation of lymphocyte proliferation _GO:0050671!
    ## 52                                   positive regulation of cytokine production _GO:0001819!
    ## 53                                               response to lipopolysaccharide _GO:0032496!
    ## 54                                     positive regulation of T cell activation _GO:0050870!
    ## 55                                     response to molecule of bacterial origin _GO:0002237!
    ## 56                                                        inflammatory response _GO:0006954!
    ## 57                                                            response to lipid _GO:0033993!
    ## 58                    cell-cell adhesion via plasma-membrane adhesion molecules _GO:0098742!
    ## 59                                    positive regulation of dendrite extension _GO:1903861!
    ## 60                                             regulation of dendrite extension _GO:1903859!
    ## 61                                  positive regulation of developmental growth _GO:0048639!
    ## 62                                           positive regulation of cell growth _GO:0030307!
    ## 63                         positive regulation of neuron projection development _GO:0010976!
    ##                    annot_type
    ## 2  GO_Biological_Process_2018
    ## 1  GO_Biological_Process_2018
    ## 3  GO_Biological_Process_2018
    ## 4  GO_Biological_Process_2018
    ## 5  GO_Biological_Process_2018
    ## 6  GO_Biological_Process_2018
    ## 7  GO_Biological_Process_2018
    ## 9  GO_Biological_Process_2018
    ## 10 GO_Biological_Process_2018
    ## 8  GO_Biological_Process_2018
    ## 11 GO_Biological_Process_2018
    ## 12 GO_Biological_Process_2018
    ## 13 GO_Biological_Process_2018
    ## 14 GO_Biological_Process_2018
    ## 15 GO_Biological_Process_2018
    ## 16 GO_Biological_Process_2018
    ## 17 GO_Biological_Process_2018
    ## 18 GO_Biological_Process_2018
    ## 19 GO_Biological_Process_2018
    ## 20 GO_Biological_Process_2018
    ## 21 GO_Biological_Process_2018
    ## 22 GO_Biological_Process_2018
    ## 23 GO_Biological_Process_2018
    ## 24 GO_Biological_Process_2018
    ## 25 GO_Biological_Process_2018
    ## 26 GO_Biological_Process_2018
    ## 27 GO_Biological_Process_2018
    ## 28 GO_Biological_Process_2018
    ## 29 GO_Biological_Process_2018
    ## 30 GO_Biological_Process_2018
    ## 31 GO_Biological_Process_2018
    ## 32 GO_Biological_Process_2018
    ## 33 GO_Biological_Process_2018
    ## 34 GO_Biological_Process_2018
    ## 35 GO_Biological_Process_2018
    ## 36 GO_Biological_Process_2018
    ## 37 GO_Biological_Process_2018
    ## 38 GO_Biological_Process_2018
    ## 39 GO_Biological_Process_2018
    ## 41 GO_Biological_Process_2018
    ## 42 GO_Biological_Process_2018
    ## 43 GO_Biological_Process_2018
    ## 44 GO_Biological_Process_2018
    ## 45 GO_Biological_Process_2018
    ## 46 GO_Biological_Process_2018
    ## 47 GO_Biological_Process_2018
    ## 48 GO_Biological_Process_2018
    ## 49 GO_Biological_Process_2018
    ## 50 GO_Biological_Process_2018
    ## 51 GO_Biological_Process_2018
    ## 52 GO_Biological_Process_2018
    ## 53 GO_Biological_Process_2018
    ## 54 GO_Biological_Process_2018
    ## 55 GO_Biological_Process_2018
    ## 56 GO_Biological_Process_2018
    ## 57 GO_Biological_Process_2018
    ## 58 GO_Biological_Process_2018
    ## 59 GO_Biological_Process_2018
    ## 60 GO_Biological_Process_2018
    ## 61 GO_Biological_Process_2018
    ## 62 GO_Biological_Process_2018
    ## 63 GO_Biological_Process_2018
    ##                                                                                  annot_names
    ## 2     negative regulation of epidermal growth factor receptor signaling pathway (GO:0042059)
    ## 1    negative regulation of epidermal growth factor-activated receptor activity (GO:0007175)
    ## 3                       negative regulation of protein tyrosine kinase activity (GO:0061099)
    ## 4             regulation of epidermal growth factor-activated receptor activity (GO:0007176)
    ## 5                            epidermal growth factor receptor signaling pathway (GO:0007173)
    ## 6                                 negative regulation of ERBB signaling pathway (GO:1901185)
    ## 7                                      negative regulation of receptor activity (GO:2000272)
    ## 9              regulation of epidermal growth factor receptor signaling pathway (GO:0042058)
    ## 10                                                       ERBB signaling pathway (GO:0038127)
    ## 8                                             entry of bacterium into host cell (GO:0035635)
    ## 11    positive regulation of epidermal growth factor receptor signaling pathway (GO:0045742)
    ## 12                                           cellular response to interleukin-6 (GO:0071354)
    ## 13                                     interleukin-6-mediated signaling pathway (GO:0070102)
    ## 14                         positive regulation of receptor-mediated endocytosis (GO:0048260)
    ## 15                                positive regulation of ERBB signaling pathway (GO:1901186)
    ## 16                                  regulation of receptor-mediated endocytosis (GO:0048259)
    ## 17                                                         entry into host cell (GO:0030260)
    ## 18                cellular response to transforming growth factor beta stimulus (GO:0071560)
    ## 19     transmembrane receptor protein serine/threonine kinase signaling pathway (GO:0007178)
    ## 20                       cellular response to fibroblast growth factor stimulus (GO:0044344)
    ## 21               positive regulation of phosphatidylinositol 3-kinase signaling (GO:0014068)
    ## 22                   transforming growth factor beta receptor signaling pathway (GO:0007179)
    ## 23                          fibroblast growth factor receptor signaling pathway (GO:0008543)
    ## 24                                           positive regulation of endocytosis (GO:0045807)
    ## 25                        regulation of phosphatidylinositol 3-kinase signaling (GO:0014066)
    ## 26                             modification-dependent protein catabolic process (GO:0019941)
    ## 27             transmembrane receptor protein tyrosine kinase signaling pathway (GO:0007169)
    ## 28                     positive regulation of intracellular signal transduction (GO:1902533)
    ## 29                                 negative regulation of programmed cell death (GO:0043069)
    ## 30                                              regulation of apoptotic process (GO:0042981)
    ## 31                                ubiquitin-dependent protein catabolic process (GO:0006511)
    ## 32                            protein modification by small protein conjugation (GO:0032446)
    ## 33                                       cellular response to cytokine stimulus (GO:0071345)
    ## 34                                          cytokine-mediated signaling pathway (GO:0019221)
    ## 35                                                       protein ubiquitination (GO:0016567)
    ## 36                                     negative regulation of apoptotic process (GO:0043066)
    ## 37                                      NLS-bearing protein import into nucleus (GO:0006607)
    ## 38                                                  protein import into nucleus (GO:0006606)
    ## 39                                                                  proteolysis (GO:0006508)
    ## 41                                              immunological synapse formation (GO:0001771)
    ## 42                                lipopolysaccharide-mediated signaling pathway (GO:0031663)
    ## 43 positive regulation of cytokine production involved in inflammatory response (GO:1900017)
    ## 44          regulation of cytokine production involved in inflammatory response (GO:1900015)
    ## 45  heterophilic cell-cell adhesion via plasma membrane cell adhesion molecules (GO:0007157)
    ## 46                                                  acute inflammatory response (GO:0002526)
    ## 47                                           regulation of T cell proliferation (GO:0042129)
    ## 48                                 positive regulation of inflammatory response (GO:0050729)
    ## 49                                      cellular response to lipopolysaccharide (GO:0071222)
    ## 50                                  positive regulation of T cell proliferation (GO:0042102)
    ## 51                              positive regulation of lymphocyte proliferation (GO:0050671)
    ## 52                                   positive regulation of cytokine production (GO:0001819)
    ## 53                                               response to lipopolysaccharide (GO:0032496)
    ## 54                                     positive regulation of T cell activation (GO:0050870)
    ## 55                                     response to molecule of bacterial origin (GO:0002237)
    ## 56                                                        inflammatory response (GO:0006954)
    ## 57                                                            response to lipid (GO:0033993)
    ## 58                    cell-cell adhesion via plasma-membrane adhesion molecules (GO:0098742)
    ## 59                                    positive regulation of dendrite extension (GO:1903861)
    ## 60                                             regulation of dendrite extension (GO:1903859)
    ## 61                                  positive regulation of developmental growth (GO:0048639)
    ## 62                                           positive regulation of cell growth (GO:0030307)
    ## 63                         positive regulation of neuron projection development (GO:0010976)
    ##    N_annot freq_annot fold_change    p_value p_value_adjust_fdr
    ## 2        3  1.0000000    2.333333 0.02857143          0.4285714
    ## 1        2  0.6666667    2.333333 0.14285714          0.4285714
    ## 3        2  0.6666667    2.333333 0.14285714          0.4285714
    ## 4        2  0.6666667    2.333333 0.14285714          0.4285714
    ## 5        2  0.6666667    2.333333 0.14285714          0.4285714
    ## 6        2  0.6666667    2.333333 0.14285714          0.4285714
    ## 7        2  0.6666667    2.333333 0.14285714          0.4285714
    ## 9        2  0.6666667    2.333333 0.14285714          0.4285714
    ## 10       2  0.6666667    2.333333 0.14285714          0.4285714
    ## 8        1  0.3333333    2.333333 0.42857143          0.4285714
    ## 11       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 12       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 13       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 14       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 15       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 16       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 17       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 18       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 19       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 20       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 21       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 22       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 23       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 24       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 25       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 26       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 27       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 28       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 29       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 30       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 31       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 32       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 33       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 34       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 35       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 36       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 37       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 38       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 39       1  0.3333333    2.333333 0.42857143          0.4285714
    ## 41       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 42       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 43       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 44       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 45       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 46       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 47       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 48       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 49       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 50       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 51       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 52       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 53       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 54       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 55       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 56       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 57       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 58       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 59       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 60       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 61       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 62       0  0.0000000    0.000000 1.00000000          1.0000000
    ## 63       0  0.0000000    0.000000 1.00000000          1.0000000
    ##         nodes_annot p_value_adjust_bonferroni N_annot_background
    ## 2  Cbl;Cblb;Eps15l1                         1                  3
    ## 1          Cbl;Cblb                         1                  2
    ## 3          Cbl;Cblb                         1                  2
    ## 4          Cbl;Cblb                         1                  2
    ## 5          Cbl;Cblb                         1                  2
    ## 6       Cbl;Eps15l1                         1                  2
    ## 7          Cbl;Cblb                         1                  2
    ## 9       Cbl;Eps15l1                         1                  2
    ## 10         Cbl;Cblb                         1                  2
    ## 8               Cbl                         1                  1
    ## 11              Cbl                         1                  1
    ## 12              Cbl                         1                  1
    ## 13              Cbl                         1                  1
    ## 14              Cbl                         1                  1
    ## 15              Cbl                         1                  1
    ## 16              Cbl                         1                  1
    ## 17              Cbl                         1                  1
    ## 18              Cbl                         1                  1
    ## 19              Cbl                         1                  1
    ## 20              Cbl                         1                  1
    ## 21              Cbl                         1                  1
    ## 22              Cbl                         1                  1
    ## 23              Cbl                         1                  1
    ## 24              Cbl                         1                  1
    ## 25              Cbl                         1                  1
    ## 26              Cbl                         1                  1
    ## 27              Cbl                         1                  1
    ## 28              Cbl                         1                  1
    ## 29              Cbl                         1                  1
    ## 30              Cbl                         1                  1
    ## 31              Cbl                         1                  1
    ## 32              Cbl                         1                  1
    ## 33              Cbl                         1                  1
    ## 34              Cbl                         1                  1
    ## 35              Cbl                         1                  1
    ## 36              Cbl                         1                  1
    ## 37             Cblb                         1                  1
    ## 38             Cblb                         1                  1
    ## 39             Cblb                         1                  1
    ## 41                                          1                  1
    ## 42                                          1                  1
    ## 43                                          1                  1
    ## 44                                          1                  1
    ## 45                                          1                  1
    ## 46                                          1                  1
    ## 47                                          1                  1
    ## 48                                          1                  1
    ## 49                                          1                  1
    ## 50                                          1                  1
    ## 51                                          1                  1
    ## 52                                          1                  1
    ## 53                                          1                  1
    ## 54                                          1                  1
    ## 55                                          1                  1
    ## 56                                          1                  1
    ## 57                                          1                  1
    ## 58                                          1                  1
    ## 59                                          1                  1
    ## 60                                          1                  1
    ## 61                                          1                  1
    ## 62                                          1                  1
    ## 63                                          1                  1
    ##    freq_annot_background nodes_annot_background
    ## 2              0.4285714       Cbl;Cblb;Eps15l1
    ## 1              0.2857143               Cbl;Cblb
    ## 3              0.2857143               Cbl;Cblb
    ## 4              0.2857143               Cbl;Cblb
    ## 5              0.2857143               Cbl;Cblb
    ## 6              0.2857143            Cbl;Eps15l1
    ## 7              0.2857143               Cbl;Cblb
    ## 9              0.2857143            Cbl;Eps15l1
    ## 10             0.2857143               Cbl;Cblb
    ## 8              0.1428571                    Cbl
    ## 11             0.1428571                    Cbl
    ## 12             0.1428571                    Cbl
    ## 13             0.1428571                    Cbl
    ## 14             0.1428571                    Cbl
    ## 15             0.1428571                    Cbl
    ## 16             0.1428571                    Cbl
    ## 17             0.1428571                    Cbl
    ## 18             0.1428571                    Cbl
    ## 19             0.1428571                    Cbl
    ## 20             0.1428571                    Cbl
    ## 21             0.1428571                    Cbl
    ## 22             0.1428571                    Cbl
    ## 23             0.1428571                    Cbl
    ## 24             0.1428571                    Cbl
    ## 25             0.1428571                    Cbl
    ## 26             0.1428571                    Cbl
    ## 27             0.1428571                    Cbl
    ## 28             0.1428571                    Cbl
    ## 29             0.1428571                    Cbl
    ## 30             0.1428571                    Cbl
    ## 31             0.1428571                    Cbl
    ## 32             0.1428571                    Cbl
    ## 33             0.1428571                    Cbl
    ## 34             0.1428571                    Cbl
    ## 35             0.1428571                    Cbl
    ## 36             0.1428571                    Cbl
    ## 37             0.1428571                   Cblb
    ## 38             0.1428571                   Cblb
    ## 39             0.1428571                   Cblb
    ## 41             0.1428571                    Cd6
    ## 42             0.1428571                    Cd6
    ## 43             0.1428571                    Cd6
    ## 44             0.1428571                    Cd6
    ## 45             0.1428571                    Cd6
    ## 46             0.1428571                    Cd6
    ## 47             0.1428571                    Cd6
    ## 48             0.1428571                    Cd6
    ## 49             0.1428571                    Cd6
    ## 50             0.1428571                    Cd6
    ## 51             0.1428571                    Cd6
    ## 52             0.1428571                    Cd6
    ## 53             0.1428571                    Cd6
    ## 54             0.1428571                    Cd6
    ## 55             0.1428571                    Cd6
    ## 56             0.1428571                    Cd6
    ## 57             0.1428571                    Cd6
    ## 58             0.1428571                    Cd6
    ## 59             0.1428571                  Itsn2
    ## 60             0.1428571                  Itsn2
    ## 61             0.1428571                  Itsn2
    ## 62             0.1428571                  Itsn2
    ## 63             0.1428571                  Itsn2
