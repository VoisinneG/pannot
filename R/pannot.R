utils::globalVariables(c("Organism.ID", "nb"))

#' Parse protein ids.
#' This is useful to remove additionnal information such as isoform or entry name.
#' 
#' @param x Character string with protein ids
#' @param sep_split Character separating different protein ids
#' @param sep_secondary Character separating UniProt entry from other infromation in protein ids
#' @param sep_collapse Character used to separate different protein ids after parsing
#' @examples
#' ids <- "A2AMW0|A2AMW0_MOUSE; P47757-2|CAPZB_MOUSE; P47757-4|CAPZB_MOUSE; Q3TVK4|Q3TVK4_MOUSE"
#' parse_ids(ids, sep = "; ", sep_secondary=c("|", "-"), sep_collapse = ";")
#' 
#' ids <- c("Q5SWU9|ACACA_MOUSE", "Q9ES52-2|SHIP1_MOUSE; Q9ES52|SHIP1_MOUSE", "Q8VDD5|MYH9_MOUSE")
#' parse_ids(ids, sep = "; ", sep_secondary=c("|", "-"), sep_collapse = ";")
#' 
#' @export
parse_ids <- function(x, sep_split = ";", sep_secondary = c("|", "-"), sep_collapse=";"){
  
  
  if(typeof(x)!="character"){
    x <- as.character(x)
    if(typeof(x)!="character"){
      stop("Could not convert x into a character")
    }
  }
  
  if(length(x)>1){
    return(sapply(x, 
                  parse_ids, 
                  sep_split = sep_split, 
                  sep_secondary = sep_secondary,  
                  sep_collapse = sep_collapse, 
                  USE.NAMES = FALSE))
  }
  else{
    prot_ids_all <- NULL
    
    prot_ids <- strsplit(x, split = sep_split, fixed = TRUE)[[1]]
    
    if(length(prot_ids) > 0){
      for(j in 1:length(prot_ids)){
        
        prot_id_int <- prot_ids[j]
        if(length(sep_secondary)>0){
          for(k in 1:length(sep_secondary)){
            prot_id_int <- strsplit(prot_id_int, split = sep_secondary[k], fixed = TRUE)[[1]][1]
          }
        }
        prot_ids_all <- c(prot_ids_all, prot_id_int)
        
      }
    }else{
      prot_ids_all <- NA
    }
    
    
    
    return( paste(unique(prot_ids_all), collapse = sep_collapse) )
  }
}
  
#' Identify which protein IDs correspond to SwissProt reviewed entries
#' @param ids vector of protein ids
#' @param sep character separating different protein ids
#' @param organism taxon id (10090 for Mus musculus, 9606 for human). If null, the most represented organism for ten first protein ids is chosen.
#' @import queryup
#' @import dplyr
#' @return a vector of reviwed protein ids
#' @export
identify_reviewed_proteins_ids <- function(ids, sep = ";", organism = NULL){
  unique_ids <- unique( strsplit( paste(ids, collapse = sep), split = sep)[[1]] )
  
  #identify most represented organism for 10 first protein ids
  if(is.null(organism)){
    cat("Guessing taxon id from 10 first protein IDs\n")
    df <- queryup::query_uniprot(query = list("id" = unique_ids[1:10]), 
                                 columns = c("id", "organism-id"))
    dfgroup <- df %>% 
      group_by(`Organism ID`) %>% 
      summarise(nb = n()) %>% 
      filter(nb == max(nb))
    organism <- dfgroup[["Organism ID"]]
    cat(paste("taxon : ", organism, "\n"))
  }
  
  #get all reviewed protein ids for selected organisms
  cat("Getting SwissProt reviewed entries\n")
  
  df <- queryup::query_uniprot(query = list("reviewed" = "yes", "organism" = organism ),
                               columns = c("id", "organism", "reviewed"))
  
  reviewed_ids <- rep(NA, length(ids))
  is_reviewed <- rep(FALSE, length(ids))
  for(i in 1:length(ids)){
    protein_ids <- strsplit(as.character(ids[i]), split = sep)[[1]]
    idx_match <- match(protein_ids, df[["Entry"]])
    idx_match <- idx_match[!is.na(idx_match)]
    if(length(idx_match)>0){
      reviewed_ids[i] <- paste(df[["Entry"]][idx_match], collapse = sep)
      is_reviewed[i] <- TRUE
    }else{
      reviewed_ids[i] <- paste(protein_ids, collapse = sep)
    }
  }
  
  return(data.frame(id = reviewed_ids, reviewed = is_reviewed))
}

#' Create a data.frame with UniProt annotations corrresponding to a set of UniProt IDs
#' @param id Character vector with UniProt IDs
#' @param sep Character separating different protein ids
#' @param columns names of uniprot data columns to retrieve. Examples include "id",
#' "genes", "keywords", "sequence", "go" (use \code{list_data_columns()} to see the full list)
#' @param max_keys maximum number of field items submitted
#' @param updateProgress used to display progress in shiny apps
#' @param show_progress Show progress bar
#' @importFrom queryup query_uniprot
#' @return a data.frame
#' @examples
#' id <- c("P26450", "O00459")
#' df <-  get_annotations_uniprot(id = id)
#' @export
get_annotations_uniprot <- function(id,
                                    sep = ";",
                                    columns = c("genes", "organism", "reviewed", "keywords", "families", "go") ,
                                    max_keys = 400,
                                    updateProgress = NULL,
                                    show_progress = TRUE){
                            
  
  idx <- which(!is.na(id))
  unique_ids <- unique( strsplit( paste(id[idx], collapse = sep), split = sep)[[1]] )
  
  query <- list("id" = unique_ids)
  columns <- union("id", columns)
  
  cat("Getting annotations from UniProt... \n")
  df_annot <- tryCatch({
    
    query_uniprot(query = query,
                  columns = columns,
                  max_keys = max_keys,
                  updateProgress = updateProgress,
                  show_progress = show_progress)
    
  }, error = function(err){
    warning("Query failed. Please retry later.")
    NULL
  })
  
  if(is.null(df_annot)) return(NULL)
  
  
  df_annot_merge <- list()
  
  for(i in 1:length(id)){
    df_annot_merge[[i]] <- rep(NA, dim(df_annot)[2])
    names(df_annot_merge[[i]]) <- names(df_annot)
    protein_ids <- strsplit(as.character(id[i]), split = sep)[[1]]
    idx_match <- match(protein_ids, df_annot[["Entry"]])
    idx_match <- idx_match[!is.na(idx_match)]
    if(length(idx_match)>0){
      for(j in 1:dim(df_annot)[2]){
        df_annot_merge[[i]][[j]] <- paste(df_annot[[j]][idx_match], collapse = "|") 
      }
    }
  }
  
  df_annot_merge <- as.data.frame( do.call(rbind, df_annot_merge) )
  
  df <- data.frame(query_id = id, df_annot_merge)
  
  return(df)
}

#' Get annotations using enrichR
#' @description Get annotations from an enrichR database for a set of genes.
#' @param data a vector of gene names or a data.frame with gene names in column \code{name_id}
#' @param name_id column name used to map gene names
#' @param dbs name of the enrichR database. Use \code{enrichR::listEnrichrDbs()} to see available databases.
#' @param append_to_data logical, append annotations as a new column
#' @return an annotated data.frame
#' @examples
#' df <- get_annotations_enrichr(c("Itsn2","Eps15l1"))
#' print(df)
#' @import enrichR
#' @export
get_annotations_enrichr <- function(data, name_id = "names", dbs = "GO_Biological_Process_2018", append_to_data = TRUE){
  #library(enrichR)
  #enrichR::listEnrichrDbs()
  #dbs<-"GO_Biological_Process_2017"
  
  df <- data
  name_id_0 <- name_id
  
  if(typeof(data) == "character"){
    df <- list("names" = data)
    name_id_0 <- "names"
  }
  
  if( length(setdiff(dbs, names(df)))==0 ){
    warning("Annotations already loaded")
    #return(df)
  }
  
  dbs_int <- setdiff(dbs, names(df))
  enriched <- enrichR::enrichr(as.character(df[[name_id_0]]), dbs_int)
  if(is.null(enriched)){
    return(NULL)
  }
  
  annot <- vector("list", length(dbs_int))
  names(annot) <- dbs_int
  
  
  
  for(i in 1:length(dbs_int)){
    
    annot[[i]] <- rep("", length(df[[name_id_0]]) )
    
    for(j in 1:length(enriched[[dbs_int[i]]]$Term)){
      genes <- strsplit(enriched[[dbs_int[i]]]$Genes[j], split = ";")[[1]]
      idx_match <- match(genes, toupper(df[[name_id]]))
      if(length(idx_match)>0){
        for(k in 1:length(idx_match)){
          if(nchar(annot[[i]][idx_match[k]])>0){
            annot[[i]][idx_match[k]] <- paste( c(annot[[i]][idx_match[k]], enriched[[dbs_int[i]]]$Term[j]), collapse = ";")
          }else{
            annot[[i]][idx_match[k]] <- enriched[[dbs_int[i]]]$Term[j]
          }
          
        }
      }
      
    }
    
    
  }
  
  annot <- as.data.frame(annot)
  annot[[name_id_0]] <- df[[name_id_0]]
  
  if(append_to_data){
    df <- merge(df, annot, by = name_id_0)
    return(df)
  }else{
    return(annot)
  }
  
}

#' Retrieve protein-protein interaction information using PSICQUIC
#' @param gene_name the gene name for which to retrieve PPI
#' @param taxon_ID taxon ID for which to retrieve PPI
#' @param provider database from which to retrieve PPI
#' @return a data.frame PPI information
#' @importFrom S4Vectors DataFrame
#' @import PSICQUIC
#' @export
get_PPI_from_psicquic <- function( gene_name, taxon_ID = c(9606,10090) , provider = c("IntAct","MINT") ){
  
  psicquic <- PSICQUIC::PSICQUIC()
  
  for (k in 1:length(taxon_ID) ){
    
    tbl <- PSICQUIC::interactions(psicquic, 
                                  gene_name, 
                                  species = taxon_ID[k] , 
                                  provider = provider )
    
    
    s<-strsplit(tbl$aliasA, split="|", fixed = TRUE);
    gene_name_A <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        ign<-grep("(gene name)",s[[i]],fixed=TRUE)
        if(length(ign)>0){
          gene_name_A[i] <- strsplit( strsplit(s[[i]][ign],split=":")[[1]][2], 
                                      split="(" ,fixed=TRUE )[[1]][1]
          
        }
      }
    }
    
    s<-strsplit(tbl$aliasB,split="|",fixed = TRUE);
    gene_name_B <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        ign<-grep("(gene name)",s[[i]],fixed=TRUE)
        if(length(ign)>0){
          gene_name_B[i] <- strsplit( strsplit(s[[i]][ign],split=":")[[1]][2], 
                                      split="(" ,fixed=TRUE )[[1]][1]
        }
      }
    }
    
    s<-strsplit(tbl$A,split = ":")
    uniprot_A <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        uniprot_A[i] <- s[[i]][2]
      }
    }
    
    s<-strsplit(tbl$B,split = ":")
    uniprot_B <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        uniprot_B[i] <- s[[i]][2]
      }
    }
    
    s<-strsplit(tbl$publicationID,split = "|",fixed=TRUE)
    Pubmed_ID <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        ip<-grep("pubmed",s[[i]],fixed=TRUE);
        Pubmed_ID[i] <- strsplit(s[[i]][ip],split=":")[[1]][2];
      }
    }
    
    s<-strsplit(tbl$type,split = "(",fixed=TRUE)
    Int_type <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        Int_type[i] <- strsplit(s[[i]][2],split=")",fixed=TRUE)[[1]][1];
      }
    }
    
    s<-strsplit(tbl$detectionMethod,split = "(",fixed=TRUE)
    Detection_method <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        Detection_method[i] <- strsplit(s[[i]][2],split=")",fixed=TRUE)[[1]][1];
      }
    }
    
    s<-strsplit(tbl$firstAuthor,split = " ",fixed=TRUE)
    Author <- rep("",length(s))
    if(length(s)>0){
      for (i in 1:length(s) ){
        Author[i] <- paste(s[[i]][1], s[[i]][length(s[[i]])],sep=" ");
      }
    }
    
    Encoding( Author ) <- "latin1"
    taxon<-rep(taxon_ID[k],length(s) );
    
    if(k>1){
      #df2<-data.frame(gene_name_A, gene_name_B, uniprot_A, uniprot_B, taxon, Int_type, Detection_method, Author=Author, pubmed_ID=pubmed_ID, Database=tbl$provider)
      df2<-data.frame(gene_name_A, gene_name_B, taxon, Int_type, 
                      Detection_method, 
                      Author=iconv(Author, "latin1", "ASCII", sub="_"), 
                      Pubmed_ID=Pubmed_ID, 
                      Database=tbl$provider)
      df1<-rbind(df1,df2);
      
    }
    else{
      #df1<-data.frame(gene_name_A, gene_name_B, uniprot_A, uniprot_B, taxon, Int_type, Detection_method, Author=Author, pubmed_ID=pubmed_ID, Database=tbl$provider)
      df1<-data.frame(gene_name_A, gene_name_B, taxon, Int_type, 
                      Detection_method, 
                      Author=iconv(Author, "latin1", "ASCII", sub="_"), 
                      Pubmed_ID=Pubmed_ID, 
                      Database=tbl$provider)
      
    }
    
    df1 <- df1[nchar(as.character(df1$gene_name_A))>0 & nchar(as.character(df1$gene_name_B))>0 , ]
    
  }
  
  return(df1)
}

#' Retrieve protein-protein interaction information from BioGRID
#' @param gene_name the gene name for which to retrieve PPI
#' @param taxon_ID taxon ID for which to retrieve PPI
#' @return a data.frame PPI information
#' @importFrom utils read.table
#' @export
get_PPI_from_BioGRID <- function( gene_name, taxon_ID = c(9606,10090) ){
  
  access_key <- "7ad36061b7644111aa9f5b3948429fb2"
  
  for (k in 1:length(taxon_ID) ){
    
    url_adress <- paste("http://webservice.thebiogrid.org/interactions?searchNames=true&geneList=",
                        gene_name,"&includeInteractors=true&format=tab2&includeHeader=true&taxId=",
                        taxon_ID[k],"&accesskey=",
                        access_key,sep="");
    
    Tbiogrid <- utils::read.table(url_adress, header=TRUE, fill=TRUE, sep="\t", comment.char="", quote="\"")
    
    Tbiogrid <- Tbiogrid[Tbiogrid$Organism.Interactor.A == Tbiogrid$Organism.Interactor.B,]
    
    taxon_biogrid <- rep(taxon_ID[k],dim(Tbiogrid)[1] );  
    
    s<-strsplit(as.character(Tbiogrid$Author),split = " ",fixed=TRUE)
    Author_Biogrid <- rep("",length(s))
    
    if(length(s)>0){
      for (i in 1:length(s) ){
        Author_Biogrid[i] <- paste(s[[i]][1], s[[i]][length(s[[i]])],sep=" ");
      }
    }
    
    
    Encoding( Author_Biogrid ) <- "latin1"
    
    if(k>1){
      df_biogrid_2 <- data.frame(gene_name_A=Tbiogrid$Official.Symbol.Interactor.A, 
                                 gene_name_B = Tbiogrid$Official.Symbol.Interactor.B, 
                                 taxon=taxon_biogrid, 
                                 Int_type=Tbiogrid$Experimental.System.Type, 
                                 Detection_method=Tbiogrid$Experimental.System, 
                                 Author=iconv(Author_Biogrid, "latin1", "ASCII", sub="_"), 
                                 Pubmed_ID=Tbiogrid$Pubmed.ID,  
                                 Database=Tbiogrid$Source.Database )
      df_biogrid_1 <- rbind(df_biogrid_1,df_biogrid_2);
    }
    else{
      df_biogrid_1 <- data.frame(gene_name_A=Tbiogrid$Official.Symbol.Interactor.A, 
                                 gene_name_B = Tbiogrid$Official.Symbol.Interactor.B, 
                                 taxon=taxon_biogrid, 
                                 Int_type=Tbiogrid$Experimental.System.Type, 
                                 Detection_method=Tbiogrid$Experimental.System,
                                 Author=iconv(Author_Biogrid, "latin1", "ASCII", sub="_"), 
                                 Pubmed_ID=Tbiogrid$Pubmed.ID,  
                                 Database=Tbiogrid$Source.Database )
      
    }
    
  }
  
  return(df_biogrid_1)
  
}

#' Retrieve protein-protein interaction information from HPRD
#' @param gene_name the gene name for which to retrieve PPI
#' @export
get_PPI_from_HPRD <- function( gene_name ){
  
  THPRD <- THPRD[which(THPRD$Gene_symbol_1 == toupper(gene_name) | THPRD$Gene_symbol_2 == toupper(gene_name)), ]
  
  df_HPRD <- data.frame(gene_name_A = THPRD$Gene_symbol_1, 
                        gene_name_B = THPRD$Gene_symbol_2, 
                        taxon = rep(9606, dim(THPRD)[1] ), 
                        Int_type = rep("NA", dim(THPRD)[1] ), 
                        Detection_method = THPRD$Experiment_type, 
                        Author = rep("NA", dim(THPRD)[1] ), 
                        Pubmed_ID = THPRD$Pubmed_id, 
                        Database = rep("HPRD", dim(THPRD)[1] ));
  return(df_HPRD)
  
}

#' Retrieve protein-protein interaction information from databses 
#' IntAct, MINT, BioGRID and HPRD
#' @param gene_name the gene name for which to retrieve PPI
#' @return a data.frame PPI information
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
create_summary_table_PPI <- function(gene_name){
  
  cat("Fetching PPI from databases...\n")
  pb <- utils::txtProgressBar(min = 0, max = 3, style = 3)
  
  df_psicquic <- try(get_PPI_from_psicquic(gene_name = gene_name), silent = FALSE)
  utils::setTxtProgressBar(pb, 1)
  df_biogrid <- get_PPI_from_BioGRID(gene_name = gene_name)
  utils::setTxtProgressBar(pb, 2)
  df_HPRD <- get_PPI_from_HPRD(gene_name = gene_name)
  utils::setTxtProgressBar(pb, 3)
  close(pb)
  
  
  
  df_tot <- rbind(df_biogrid, df_psicquic, df_HPRD)
  
  uInteractors <- sort(unique(toupper(c(as.character(df_tot$gene_name_A), as.character(df_tot$gene_name_B) ) )))
  uInteractors <- uInteractors[ uInteractors != toupper(gene_name) ];
  
  
  
  
  if(length(uInteractors)>0){
    
    N_pub_0 <- rep(0,length(uInteractors));
    Authors_0 <- rep("",length(uInteractors));
    Pubmed_ID_0 <- rep("",length(uInteractors));
    Detection_method_0 <- rep("",length(uInteractors));
    Int_type_0 <- rep("",length(uInteractors));
    Database_0 <- rep("",length(uInteractors));
    
    for (i in 1:length(uInteractors) ){
      
      i_int <- which( toupper(as.character(df_tot$gene_name_A)) == uInteractors[i] | toupper(as.character(df_tot$gene_name_B)) == uInteractors[i]  )
      
      if(length(i_int)>0){
        Authors_0[i] <- paste(as.character(unique(df_tot$Author[i_int])), collapse="|")
        
        Pubmed_ID_0[i] <- paste(as.character(unique(df_tot$Pubmed_ID[i_int])), collapse=",")
        spl <- strsplit(Pubmed_ID_0[i], split=",");
        Pubmed_ID_0[i] <- paste(spl[[1]], collapse="|");
        N_pub_0[i] <- length(unique(spl[[1]]));
        
        #N_pub[idx_int] <- length(unique(df_tot$pubmed_ID[i_int]));
        
        Detection_method_0[i] <- paste(as.character(unique(df_tot$Detection_method[i_int])), collapse="|")
        Int_type_0[i] <- paste(as.character(unique(df_tot$Int_type[i_int])), collapse="|")
        Database_0[i] <- paste(as.character(unique(df_tot$Database[i_int])), collapse="|")
      }
      
      
    }
    
    df_summary <- data.frame(gene_name_A=rep(toupper(gene_name),length(uInteractors)), 
                            gene_name_B=uInteractors, 
                            N_pub = N_pub_0, 
                            Authors = Authors_0, 
                            Pubmed_ID = Pubmed_ID_0,
                            Detection_method = Detection_method_0,
                            Int_type = Int_type_0,
                            Database= Database_0)
  }else{
    message(paste("No interactions involving", gene_name, "found in databases"))
    df_summary <- NULL
  }
  
  return(df_summary)
}


#' Perform enrichment analysis
#' @description Perform enrichment analysis using a hypergeometric test for protein annotations stored in a formatted data.frame
#' @param df a data.frame with annotations corresponding to each row. Types of annotations are organized by columns. 
#' For a given type of annotations, annotations are separated by \code{sep}.
#' @param sep_split Character string separating different annotation groups.
#' @param sep Character string separating different annotations.
#' @param idx_subset indexes of the foreground set.
#' @param annotation_selected set of annotations on which to perform the analysis. Annotations selectd must be a subset of df's names.
#' @param col_names df's column name containing gene names.
#' @param two_sided logical, perform a two-sided hypergeometric test
#' @param updateProgress logical, function to show progress in shiny app
#' @param showProgress logical, show progress in console
#' @param orderOutput logical, order annotations by enrichment p-values in the output data.frame
#' @return a data.frame
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats phyper p.adjust
#' @export
annotation_enrichment_analysis <- function( df,
                                            sep_split="|",
                                            sep = NULL,
                                            idx_subset, 
                                            annotation_selected = names(df)[2], 
                                            col_names = names(df)[1], 
                                            two_sided = FALSE,
                                            updateProgress = NULL, 
                                            showProgress = TRUE,
                                            orderOutput = TRUE){
  
  
  if( is.null(df) | (sum(annotation_selected %in% names(df)) != length(annotation_selected)) ){
    stop("Annotations not available. Import annotations first.")
  }else if ( length(annotation_selected) == 0) {
    stop("No annotations selected. Change selected annotations")
  }
  
  if (showProgress) cat("Perform annotation enrichment analysis...\n")
  
  #list annotation terms found in the dataset ------------------------------------------------
  
  df_int <- df
  
  for(i in length(annotation_selected)){
    df_int[[annotation_selected[i]]] <- as.character(df_int[[annotation_selected[i]]])
  }
  
  annot_terms <- NULL
  annot_type <- NULL
  annot_names <- NULL
  
  if(is.null(sep)){
    warning("No separation character provided. Using ';' by default")
    collapse_sep <- ";"
  }else{
    collapse_sep <- sep
  }
  
  for (annot_type_sel in annotation_selected){
    
    df_int[[annot_type_sel]] <- gsub(sep_split, collapse_sep, df_int[[annot_type_sel]], fixed = TRUE)
    
    df_int[[annot_type_sel]] <- gsub("(", "_", df_int[[annot_type_sel]], fixed = TRUE)
    df_int[[annot_type_sel]] <- gsub(")", "!", df_int[[annot_type_sel]], fixed = TRUE)
    df_int[[annot_type_sel]] <- gsub("[", "_", df_int[[annot_type_sel]], fixed = TRUE)
    df_int[[annot_type_sel]] <- gsub("]", "!", df_int[[annot_type_sel]], fixed = TRUE)
    
    # if( annot_type_sel %in% c("Protein.families", "Keywords") ) collapse_sep <- "; "
    
    
    u_annot<-paste(unique(df_int[[annot_type_sel]]), collapse = collapse_sep)
    terms <- unique(strsplit(u_annot, split = collapse_sep)[[1]])
    
    annot_names_int <- terms
    
    
    annot_terms <- c(annot_terms,  terms)
    annot_type <- c(annot_type, rep(annot_type_sel, length(terms)))
    annot_names <- c(annot_names, annot_names_int)
    
  }
  
  df.annot <- data.frame(annot_terms = annot_terms, annot_type = annot_type, annot_names = annot_names)
  df.annot <- df.annot[which(df.annot$annot_terms != ""), ]
  
  df.annot$annot_names <- gsub("_", "(", df.annot$annot_names, fixed = TRUE)
  df.annot$annot_names <- gsub("!", ")", df.annot$annot_names, fixed = TRUE)
  
  
  # Compute Background -------------------------------------------------------------------------------------
  
  nodes_tot <- as.character(df[[col_names]]);
  
  u_annot_nodes_collapse <- rep("", length(nodes_tot));
  idx_tot <- rep(0, length(nodes_tot));
  
  for ( i in 1:length(nodes_tot) ){
    s <- NULL
    for (annot_type in annotation_selected) {
      s <- c(s, as.character(df_int[[annot_type]][ i ]))
    }
    u_annot_nodes_collapse[i] <- paste(s, collapse = collapse_sep)
  }
  
  N_background = length(nodes_tot);
  
  n_annot <- dim(df.annot)[1]
  
  N_annot_background <- rep(0, n_annot);
  freq_annot_background <- rep(0, n_annot);
  nodes_annot_background <- rep("", n_annot);
  
  if (showProgress & typeof(idx_subset)!="list") pb <- utils::txtProgressBar(min = 0, max = 2*n_annot, style = 3)
  count<-0
  
  for ( k in 1:dim(df.annot)[1] ){
    
    annot <- df.annot$annot_terms[k]
    
    idx_annot <- grep(paste("(",collapse_sep,"|^)", annot, "($|",collapse_sep,")",sep=""), 
                      u_annot_nodes_collapse, fixed=FALSE)
    
    #idx_annot <- grep(df.annot$annot_terms[k], u_annot_nodes_collapse, fixed=TRUE)
    
    N_annot_background[k] = length(idx_annot);
    nodes_annot_background[k] = paste(nodes_tot[idx_annot], collapse=";")
    freq_annot_background[k] = N_annot_background[k]/N_background;
    
    count <- count +1
    if (showProgress & typeof(idx_subset)!="list") utils::setTxtProgressBar(pb, count)
    # progress bar
    if (is.function(updateProgress)) {
      text <- paste0( round(count/(2*n_annot)*100, 0), " %")
      updateProgress(value = count/(2*n_annot)*100, detail = text)
    }
  }
  
  N_annotation_test <- length(which(N_annot_background>0))
  
  # Perform enrichment test for each annotation and each subset of indices ---------------------------------------
  
  if (typeof(idx_subset)=="list"){
    n_sets <- length(idx_subset)
  } else {
    n_sets = 1
  }
  
  df.annot.tot <- list()
  
  if (showProgress  & typeof(idx_subset)=="list") {
    pb <- utils::txtProgressBar(min = 0, max = n_sets, style = 3)
  }
  
  for (i in 1:n_sets){
    
    if (showProgress & typeof(idx_subset)=="list"){
      utils::setTxtProgressBar(pb, i)
    } 
    
    if (typeof(idx_subset)=="list"){
      idx_d <- idx_subset[[i]]
    } else {
      idx_d <- idx_subset
    }
    
    N_annot <- rep(0, n_annot);
    freq_annot <- rep(0, n_annot);
    nodes_annot <- rep("", n_annot);
    p_value <- rep(0, n_annot);
    fold_change <- rep(0, n_annot);
    p_value_adjust <- rep(0, n_annot);
    
    for( k in 1:n_annot ){
      
      annot <- df.annot$annot_terms[k]
      
      idx_annot <- idx_d[ grep(paste("(",collapse_sep,"|^)", annot, "($|",collapse_sep,")",sep=""), 
                               u_annot_nodes_collapse[idx_d], fixed=FALSE) ]
      
      #idx_annot <- idx_d[ grep(df.annot$annot_terms[k], u_annot_nodes_collapse[idx_d],fixed=TRUE) ]
      
      N_annot[k]=length(idx_annot);
      N_sample = length(idx_d);
      
      freq_annot[k] = N_annot[k]/N_sample;
      nodes_annot[k]=paste(nodes_tot[idx_annot], collapse=";")
      
      inclusive_upper_tail <- 1-phyper(N_annot[k]-1,
                                       N_annot_background[k],
                                       N_background-N_annot_background[k],  
                                       N_sample)
      
      inclusive_lower_tail <- phyper(N_annot[k],
                                     N_annot_background[k],  
                                     N_background-N_annot_background[k],  
                                     N_sample)
      
      two_sided_hypergeometric_p_value <- 2.0*min(c(inclusive_upper_tail, inclusive_lower_tail))
      
      
      if(two_sided){
        p_value[k] = two_sided_hypergeometric_p_value
      }else{
        p_value[k] = inclusive_upper_tail
      }
      
      
      fold_change[k] = freq_annot[k]/freq_annot_background[k];
      
      count <- count +1
      if (showProgress & typeof(idx_subset)!="list"){
        utils::setTxtProgressBar(pb, count)
      }
      # progress bar
      if (is.function(updateProgress)) {
        text <- paste0( round(count/(2*n_annot)*100, 0), " %")
        updateProgress(value = count/(2*n_annot)*100, detail = text)
      }
    }
    
    if (showProgress & typeof(idx_subset)!="list"){
      close(pb)
      cat("Done.\n")
    }
    
    if(two_sided){
      idx_annot_exist <-  which(N_annot_background>0)
    }else{
      idx_annot_exist <-  which(N_annot>0)
    }
    
    p_value_adjust_fdr <- rep( 1,length(p_value) );
    p_value_adjust_bonferroni <- rep( 1,length(p_value) );
    p_value_adjust_bonferroni[idx_annot_exist] <- p.adjust(p_value[idx_annot_exist], method = "bonferroni");
    p_value_adjust_fdr[idx_annot_exist] <- p.adjust(p_value[idx_annot_exist], method = "fdr");
    
    df.annot.set <- data.frame(
      N_annot,
      freq_annot,
      fold_change, 
      p_value, 
      p_value_adjust_fdr,
      nodes_annot,
      p_value_adjust_bonferroni,
      N_annot_background, 
      freq_annot_background,
      nodes_annot_background)
    
    df.annot.set <- cbind(df.annot, df.annot.set)
    
    if (orderOutput) df.annot.set <- df.annot.set[ order(df.annot.set$p_value, decreasing = FALSE), ]
    
    df.annot.tot[[i]] <- df.annot.set
    
  }
  
  if (showProgress & typeof(idx_subset)=="list") close(pb)
  
  if (typeof(idx_subset)=="list"){
    return(df.annot.tot)
  } else {
    return(df.annot.tot[[1]])
  }
  
}

#' Plot the result of the annotation enrichment analysis
#' @param df a formatted data.frame obtained by the function \code{annotation_enrichment_analysis()}
#' @param p_val_max threshold for the enrichment p-value
#' @param method_adjust_p_val method to adjust p-value for multiple comparisons
#' @param fold_change_min threshold for the enrichment fold-change
#' @param N_annot_min minimum number of elements that are annotated in the foreground set
#' @param test_depletion logical, test for annotation depletion as well as enrichment
#' @return a data.frame
#' @export
filter_annotation_results <- function(df, 
                                      p_val_max=0.05, 
                                      method_adjust_p_val = "fdr", 
                                      fold_change_min =2,
                                      N_annot_min=2, 
                                      test_depletion = FALSE
                                      ){
  
  if(length(df) == 0 ){
    warning("Empty input...")
  }else if( dim(df)[1] == 0){
    warning("Empty input...")
  }
  
  name_p_val <- switch(method_adjust_p_val,
                       "none" = "p_value",
                       "fdr" = "p_value_adjust_fdr",
                       "bonferroni" = "p_value_adjust_bonferroni")
  
  df$p_value <- df[[name_p_val]]
  
  
  if(test_depletion){
    idx_filter <-  which(df$p_value <= p_val_max & 
                           (df$fold_change >= fold_change_min | df$fold_change <= 1/fold_change_min) & 
                           df$N_annot >= N_annot_min)
  } else {
    idx_filter <-  which(df$p_value <= p_val_max & 
                           df$fold_change >= fold_change_min & 
                           df$N_annot >= N_annot_min)
  }
  
  
  
  if(length(idx_filter) == 0){
    warning("No annotation left after filtering. You might want to change input parameters")
    return(NULL)
  }
  df_filter <- df[ idx_filter, ]
  
  return(df_filter)
}


#' Plot the result of the annotation enrichment analysis
#' @param df a formatted data.frame obtained by the function \code{annotation_enrichment_analysis()}
#' @param method_adjust_p_val name of the p-value variable. Possible choices : "none","fdr" or "bonferroni".
#' By default : "fdr".
#' @param fold_change_limits limits for fold-change scale. By default : c(1,4).
#' @param save_file path where the plot will be saved
#' @param filter Apply filtering to annotation results (using \code{filter_annotation_results()})
#' @param type Palette type (one of seq (sequential), div (diverging)) 
#' @param flip flip x-y coordinates
#' @param angle angle of x tick labels
#' @param font_size font size
#' @param trans trans object applied to fill scale
#' @param scale_factor_width factor to adjust plot width
#' @param scale_factor_height factor to adjust plot height
#' @param ... parameters passed to function \code{filter_annotation_results()}
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @import scales
#' @import RColorBrewer
#' @return a plot
#' @export
plot_annotation_results <- function(df,
                                    method_adjust_p_val = "fdr",
                                    fold_change_limits = c(1,4),
                                    save_file = NULL,
                                    filter = FALSE,
                                    type = "seq",
                                    flip = FALSE,
                                    angle = 90,
                                    font_size = 10,
                                    trans = NULL,
                                    scale_factor_width = 1,
                                    scale_factor_height = 1,
                                    ...){
  
  if(length(df) == 0 ){
    warning("Empty input...")
  }else if( dim(df)[1] == 0){
    warning("Empty input...")
  }
  
  if(is.character(trans)){
    trans <- switch(trans,
                    "linear" = identity_trans(),
                    "log2" = log2_trans())
  }
  if(is.null(trans)) trans <- identity_trans()
  
  name_p_val <- switch(method_adjust_p_val,
                       "none" = "p_value",
                       "fdr" = "p_value_adjust_fdr",
                       "bonferroni" = "p_value_adjust_bonferroni")
  
  df_filter <- df
  if(filter){
    df_filter <- filter_annotation_results(df = df, 
                                           method_adjust_p_val = method_adjust_p_val, 
                                           ...)
  }
  
  if(length(df_filter$fold_change) == 0){
    warning("No annotation left after filtering. You might want to change input parameters")
    return(NULL)
  }
  
  df_filter$fold_change_sign <- df_filter$fold_change
  df_filter$fold_change_sign[df_filter$fold_change>=1] <- 1
  df_filter$fold_change_sign[df_filter$fold_change<1] <- -1
  
  df_filter$annot_names[df_filter$p_value == 0] <- paste("*", df_filter$annot_names[df_filter$p_value == 0])
  df_filter$p_value[df_filter$p_value == 0] <- min(df_filter$p_value[df_filter$p_value != 0])*0.1
  
  df_filter <- df_filter[ order(df_filter$fold_change_sign * (-log10(df_filter$p_value)), decreasing = FALSE), ]
  #df_filter <- df_filter[ order(df_filter$p_value, decreasing = TRUE), ]
  df_filter$order <- 1:dim(df_filter)[1]
  df_filter$fold_change[df_filter$fold_change >= fold_change_limits[2]] <- fold_change_limits[2]
  df_filter$fold_change[df_filter$fold_change <= fold_change_limits[1]] <- fold_change_limits[1]
  
  df_filter$minus_log10_p_value <- -log10(df_filter$p_value)
  df_filter$log2_fold_change = log2(df_filter$fold_change)
  
  p <- ggplot( df_filter, aes_string(x="order", 
                                     y="minus_log10_p_value",
                                     fill = "fold_change")) + 
    theme(
      axis.text.y = element_text(size=font_size),
      axis.text.x = element_text(size=font_size, angle = angle, hjust = 1, vjust = 1),
      axis.title.x = element_text(size=font_size)
    ) +
    scale_x_continuous(name = NULL, breaks=df_filter$order, labels=df_filter$annot_names) +
    scale_y_continuous(name = paste("-log10(",name_p_val,")",sep="")) 
    #scale_fill_distiller(palette = palette, direction = direction, trans = trans, limits = c(fold_change_limits[1], fold_change_limits[2])) + 
    
  if(type == "seq"){
    p <- p + scale_fill_gradientn(colors = RColorBrewer::brewer.pal(name = "RdBu", n = 11)[6:1],
                         trans = trans, limits = fold_change_limits)
  }else if(type == "div"){
    p <- p + scale_fill_gradientn(colors = RColorBrewer::brewer.pal(name = "RdBu", n = 11)[11:1],
                                  trans = trans, limits = fold_change_limits)
  }
    
  p <- p +  geom_col()
  
  if(flip){
    p <- p + coord_flip() + ylab("")
  }else{
    p <- p + xlab("")
  }
  
  if(!is.null(save_file)){
    plot_height <- 0.1*( 0.5*max( sapply(as.character(df_filter$annot_terms), nchar) ) + 20 )
    plot_width <- 0.13*(1.5*length(unique(df_filter$annot_terms)) + 35)
    if(flip){
      plot_width <- 0.13*( 0.5*max( sapply(as.character(df_filter$annot_terms), nchar) ) + 35 )
      plot_height <- 0.1*(1.5*length(unique(df_filter$annot_terms)) + 20)
    }
    pdf(save_file, plot_width*scale_factor_width, plot_height*scale_factor_height)
    print(p)
    dev.off()
  }
  
  return(p)
  
}

