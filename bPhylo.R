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
#seqdir="~/ColauttiLabScratch/QIIME2/02:04:2020/UniqSeqs.fasta" # seq 8 interesting check BS tree
seqdir="~/ColauttiLabScratch/20200229_top18.fasta"
dbdir="~/ColauttiLabScratch/16SMicrobialDB/16SMicrobial" #https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

#Initialize
seqs <- readDNAStringSet(seqdir) #load sequences
blastDB<- blast(dbdir) #load 16S database

#Blast refseq, create DF containing Blast results
blast_seq<-function(refseq ,blastDB,hits,targetName){ 
  BR<-predict(blastDB, refseq) #Blast sequence against 16S database, returns up to 500 hits
  BRhitLen<-length(BR$SubjectID) #return number of blast hits
  if(BRhitLen>hits) {
    BRhitLen=hits #Limit hit length
  }
  trun_BR<-as.vector(BR[1:BRhitLen,])
  write.table(trun_BR$SubjectID, file = "./GBACC.txt", row.names = F, col.names = F, quote = F)
  system('./retrieveSeq.sh')
  BRSeqs<-readDNAStringSet("seqs.txt")
  BR$Taxa<-paste(gsub('16S .*', '',BRSeqs@ranges@NAMES))
  BR$Seqs<-substr(paste(BRSeqs),BR$S.end,BR$S.start)
  targetDF<-data.frame(Taxa=targetName, Seqs=reverseComplement(refseq)) #Create target DF
  GBDF<-rbind.fill(BR,targetDF) #Add target DF to main DF
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
