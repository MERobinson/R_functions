# Biomart
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("genes", dataset = "hsapiens_gene_ensembl",  
                     host = "https://dec2021.archive.ensembl.org", mirror = "www")
  mouse = useEnsembl("genes", dataset = "mmusculus_gene_ensembl",
                     host = "https://dec2021.archive.ensembl.org", mirror = "www")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , 
                   mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

