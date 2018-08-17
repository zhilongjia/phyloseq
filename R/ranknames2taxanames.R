#' convert ranknames to taxanames
#' 
#' convert ranknames, like Genus, to taxanames, like otu1.
#' 
#' @param taxa A character vector of the taxa in object x that you 
#' want to keep. See \code{\link{prune_taxa}} for details.
#' @param x A phylogenetic object
#' @param rankby one of rank_names(x)
#' @export
#' 
#' @return a character vector of taxanames
#' data("GlobalPatterns")
#' rankby = "Genus"
#' TopNrank <- names(sort(taxa_sums(GlobalPatterns, rankby), TRUE)[1:5])
#' TopNOTUs <- ranknames2taxanames(TopNrank, GlobalPatterns, rankby)
#' rankby_prune  <- prune_taxa(TopNOTUs, GlobalPatterns)

ranknames2taxanames <- function(taxa, x, rankby="Genus") {
    
    rankby <- match.arg(rankby, rank_names(x) )
    
    tax_tab_raw <- as.data.frame(tax_table(x)@.Data)
    tax_tab_filtered <- tax_tab_raw[as.character(tax_tab_raw[[rankby]]) %in% taxa, ]
    rownames(tax_tab_filtered)
}


