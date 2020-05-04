#' Blast Sequence
#'
#' @param refseq Reference sequence to be BLAST searched.
#' @param blastDB The database to be searched against.
#' @param hits Sets limit to length of hits list.
#' @param targetName Name of the target species/taxon.
#'
#' @return Returns dataframe of Blast results.
#' @export
#' @import rBLAST
#' @import Biostrings
#' @import plyr
#'
#' @examples
blast_seq <- function(refseq ,blastDB,hits,targetName){
  BR <- rBLAST::predict(blastDB, refseq) #Blast sequence against 16S database, returns up to 500 hits
  BRhitLen <- length(BR$SubjectID) #return number of blast hits
  if(BRhitLen>hits) {
    BRhitLen <- hits #Limit hit length
  }
  trun_BR <- as.vector(BR[1:BRhitLen,])
  write.table(trun_BR$SubjectID, file = "./GBACC.txt", row.names = F, col.names = F, quote = F)
  system('./retrieveSeq.sh')
  BRSeqs <- Biostrings::readDNAStringSet("seqs.txt")
  trun_BR$Taxa <- paste(gsub('16S .*', '',BRSeqs@ranges@NAMES))
  trun_BR$Seqs <- substr(paste(BRSeqs),BR$S.end,BR$S.start)
  targetDF <- data.frame(Taxa=targetName, Seqs=reverseComplement(refseq)) #Create target DF
  GBDF <- rbind.fill(trun_BR,targetDF) #Add target DF to main DF
  return(GBDF)
}

#' Add New Target
#'
#' @param SeqDF Existing sequence dataframe to be appended.
#' @param newTar Sequence of new target.
#' @param targetName Name of new target species/taxon.
#'
#' @return Returns new sequence dataframe with new target blast data
#' @export
#'
#' @examples
add_target <- function(SeqDF,newTar,targetName){
  tar2 <- data.frame(Taxa=targetName, Seqs=Biostrings::reverseComplement(newTar))
  tar2DF <-merge (SeqDF,tar2, all.y = T)
  finalDF <- rbind (SeqDF,tar2DF)
}


#' Write Fasta
#'
#' @param SeqDF Dataframe of sequences
#' @param path Path where fasta file will be saved
#'
#' @return None
#' @export
#' @import seqRFLP
#'
#' @examples
to_fasta <- function(SeqDF,path){
  fasDF <- data.frame('names' = paste(SeqDF$Taxa), 'sequences' = SeqDF$Seqs)
  fasDF.fasta <- seqRFLP::dataframe2fas(fasDF, file=path)
}


#' Align Sequences
#'
#' @param SeqDF Dataframe of sequences to be aligned
#'
#' @return Returns aligned sequence in the form of a DNAbin object.
#' @export
#' @import ape
#'
#' @examples
align_seqs <- function(SeqDF){
  GbkDNA <- sapply(paste(SeqDF$Seqs),strsplit,split="")
  names(GbkDNA) <- paste(SeqDF$Taxa)
  GbkDNA <- ape::as.DNAbin(GbkDNA)
  GbkAlign <- ape::muscle(GbkDNA,quiet=F)
  ape::checkAlignment(GbkAlign, what=1)
  return(GbkAlign)
}


#' Convert to phyDat format
#'
#' @param aligned_seqs DNAbin object of aligned sequences.
#'
#' @return Returns PhyloSeq Object
#' @export
#' @import phangorn
#'
#' @examples
convert_phyobject<-function(aligned_seqs){
  phy_object<-phangorn::phyDat(aligned_seqs, type = "DNA")
  return(phy_object)
}


#' Optimize parsimony for distance based methods: NJ, UPGMA
#'
#' @param phy_object phyDat format object.
#'
#' @return Returns pml object
#' @export
#' @import phangorn
#'
#' @examples
optimize_db_parsimony<-function(phy_object) {
  dna_dist <- phangorn::dist.ml(phy_object)
  treeNJ <- phangorn::NJ(dna_dist)
  treeUPGMA <- phangorn::upgma(dna_dist)
  if (phangorn::parsimony(treeNJ,phy_object) < phangorn::parsimony(treeUPGMA,phy_object)) {
    fit = phangorn::pml(treeNJ,phy_object)
  } else {
    fit = phangorn::pml(treeUPGMA,phy_object)
  }
  return(fit)
  
}


#' Model selection for tree optimization
#'
#' @param phy_object phyDat format object.
#' @param pml_object PML object.
#'
#' @return Returns pml object
#' @export
#' @import phangorn
#'
#' @examples
optimize_likelihood<-function(phy_object, pml_object) {
  model_test <- phangorn::modelTest(phy_object)
  model_test <- model_test %>% plyr::arrange(BIC)
  print("----------------------------------------------------------------------")
  print(model_test[1,])
  params <- strsplit(paste(model_test[1,]$Model),"\\+")
  model <- params[[1]][1]
  Gamma <- TRUE
  Inv <- TRUE
  if (length(params[[1]]) == 2 & params[[1]][2] == "I") {
    Gamma = FALSE
  } else {
    Inv = FALSE
  }
  print(paste0("Gamma = ",Gamma))
  print(paste0("Inv = ", Inv))
  print("----------------------------------------------------------------------")
  fit <- phangorn::optim.pml(pml_object, model = model, optInv = Inv, optGamma = Gamma, rearrangement = "stochastic")
  print("----------------------------------------------------------------------")
  print(pml_object)
  print("----------------------------------------------------------------------")
  print(fit)
  print("----------------------------------------------------------------------")
  return(fit)
}

#' Bootstrap tree and save as PDF
#'
#' @param fitted_model PML object.
#' @param bs_iterations Integer for bootstrap iterations.
#' @param path Path where PDF will be saved.
#'
#' @return None
#' @export
#' @import phangorn
#' @import ape
#'
#' @examples
bootstrap_tree<-function(fitted_model,bs_iterations,path){ 
  bs <- phangorn::bootstrap.pml(fitted_model, bs=bs_iterations, optNni=TRUE, multicore=F, control = pml.control(trace=0))
  pdf(path,width=16,height=12)
  phangorn::plotBS(midpoint(fitted_model$tree),bs, p = 0, type="p")
  ape::add.scale.bar(length = 0.01, cex = 0.9, font = 2)
  dev.off()
}

