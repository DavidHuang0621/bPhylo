#Iterative Phylogeny Builder V2 by David Huang 04/26/2020
#Required Libraries:
library(Biostrings) 
library(rBLAST) 
library(genbankr) 
library(plyr) 
library(ape) 
library(phangorn) 
library(ggtree)
library(seqRFLP)

#Input files
seqdir="~/ColauttiLabScratch/QIIME2/02:04:2020/UniqSeqs.fasta" # seq 8 interesting check BS tree
dbdir="~/ColauttiLabScratch/16SMicrobialDB/16SMicrobial" #https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

#Initialize
seqs <- readDNAStringSet(seqdir) #load sequences
blastDB<- blast(dbdir) #load 16S database

#Blast refseq, create DF containing Blast results
blast_seq<-function(refseq ,blastDB,hits,targetName){ 
  BR<-predict(blastDB, refseq) #Blast sequence against 16S database, returns up to 500 hits
  BRhitLen<-length(BR$SubjectID) #return number of blast hits
  if(BRhitLen>hits) BRhitLen=hits #Limit hit length
  BRhits<-as.vector(BR$SubjectID[1:BRhitLen]) #Create truncated hit vector
  GBAlist<-lapply(BRhits, GBAccession) #convert BRsublist objects to GenBankAccession objects
  GbRs<<-lapply(GBAlist, readGenBank) #Obtain sequence data from GenBank
  BrSeqDF<-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("ID", "Taxa", "Seqs"))
  for (j in 1:length(GbRs)) {
    tempGBR <- GbRs[[j]]
    tempDF <- data.frame(ID=paste(GBAlist[j]), Taxa=paste(gsub('16S .*', '', tempGBR@definition), BR$SubjectID[j]), Seqs=substr(paste(tempGBR@sequence),BR$S.end[j],BR$S.start[j]))
    BrSeqDF <- rbind(BrSeqDF, tempDF)
  } #Extract sequence data, append to list
  targetDF<-data.frame(Taxa=targetName, Seqs=reverseComplement(refseq)) #Create target DF
  GBDF<-cbind(BR[1:BRhitLen,],BrSeqDF) #Add blast result information to DF
  GBDF<-rbind.fill(GBDF,targetDF) #Add target DF to main DF
  return(GBDF)
}
SeqDF<-blast_seq(seqs[2],blastDB,50,'Borrelia')

#Add more targets
add_target<-function(SeqDF,newTar,targetName){
  tar2 <-data.frame(Taxa=targetName, Seqs=reverseComplement(newTar))
  tar2DF <-merge(SeqDF,tar2, all.y = T)
  finalDF<-rbind(SeqDF,tar2DF)
}
SeqDF<-add_target(SeqDF,seqs[3],'Borreliella')

#Write Fasta
to_fasta<-function(SeqDF,path){
  fasDF <- data.frame('names' = paste(SeqDF$Taxa), 'sequences' = SeqDF$Seqs)
  fasDF.fasta = dataframe2fas(fasDF, file=path)
}
to_fasta(SeqDF,"~/ColauttiLabScratch/FastasForClaire/BorreliaTrimmed")

#Align sequences 
align_seqs<-function(SeqDF){
  GbkDNA<-sapply(paste(SeqDF$Seqs),strsplit,split="")
  names(GbkDNA)<-paste(SeqDF$Taxa)
  GbkDNA<-as.DNAbin(GbkDNA)
  GbkAlign<-muscle(GbkDNA,quiet=F)
  checkAlignment(GbkAlign, what=1)
  return(GbkAlign)
} 
AlignedSeqs <- align_seqs(SeqDF)

#Build phylogenic trees
bootstrap_tree<-function(AlignedSeqs){ 
  phyObject<-phyDat(AlignedSeqs, type = "DNA")
  dna_dist <- dist.ml(phyObject, model="F81")
  dna_UPGMA <- upgma(dna_dist)
  fit <- pml(dna_UPGMA, phyObject)
  fitJC <- optim.pml(fit, model = "F81", rearrangement = "stochastic")
  bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=F, control = pml.control(trace=0))
  plotBS(midpoint(fitJC$tree),bs, p = 0, type="p", main = "P = 0")
}
Bs<-bootstrap_tree(AlignedSeqs)
plot(Bs)

NJ_Tree<-function(AlignedSeqs){
  SeqDNA<-AlignedSeqs
  GbkDM<-dist.dna(SeqDNA, model="K80")
  Tree<-nj(GbkDM)
  PhyloTree<-ggtree(Tree,layout="rectangular") + geom_tiplab() + xlim(NA,0.3)
  return(PhyloTree)
}
NJ<-NJ_Tree(AlignedSeqs)
NJ
