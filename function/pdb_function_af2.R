source("function/translator.R")

#Writes the contact PDB output file
write_contacts_pdb<-function(shared_cons_contacts,name,tax,pdb_ID,contact_thresh,conservation_thresh,contact_number,mainchain) {
#read representative pdb for projection

file<-paste0("data/pdb/all_structures/","pdb",pdb_ID,".pdb")
pdb<-read.pdb(file)

#give matching alignment entry
path<-paste0("data_output/",tax,"_",pdb_ID,"_CN_",contact_thresh,"_CONS_",conservation_thresh,"_contacts_",contact_number,"_contact_network.pdb")
#Remove everything but the CA
ca.inds <- atom.select(pdb, "calpha", chain = mainchain, inverse=F)
only_ca <- trim.pdb(pdb, ca.inds)
#write as pdb (with resnumber for atomid)
write.pdb(only_ca, file=path,eleno=only_ca$atom$resno)

#translate the shared & conserved contacts back into structure-specific (e.g. Ago2) residue positions
shared_cons_contacts_trans<-c()
shared_cons_contacts_trans$pos1<-translate2_top(name, shared_cons_contacts$from)
shared_cons_contacts_trans$pos2<-translate2_top(name, shared_cons_contacts$to)
shared_cons_contacts_trans<-as.data.frame(shared_cons_contacts_trans)


#Remove END from pdb file
data <- readLines(path,warn = FALSE)
data<-data[-length(data)]
file<-file(path)
writeLines(data, file)
close(file)

#Format for pymol "conect" output, append to CA-only structure for visualization
sink(path, append = T)
cat("TER")
cat("\n")
for (i in 1:dim(shared_cons_contacts_trans)[1]){
  line<-c(shared_cons_contacts_trans[i,1], shared_cons_contacts_trans[i,2])
  cat("CONECT",format(line, width = 5, justify = "right"),sep="")
  cat("\n")
}
cat("END")
cat("\n")
sink()
}

#Writes the conservation PDB output file
write_conservation_pdb<-function(scores,name,tax,pdb_ID,contact_thresh,conservation_thresh,contact_number,mainchain,IDs) {
  #read representative pdb for projection
  file<-paste0("data/pdb/all_structures/","pdb",pdb_ID,".pdb")
  
  suppressWarnings(pdb<-read.pdb(file))
  
  #give matching alignment entry
  path<-paste0("data_output/",tax,"_",pdb_ID,"_CN_",contact_thresh,"_CONS_",conservation_thresh,"_contacts_",contact_number,"_conservation.pdb")
  pdb<-trim.pdb(pdb,chain="A")
  suppressWarnings(pdb<-clean.pdb(pdb, rm.wat = TRUE,rm.lig=TRUE))
  
  scores<-get_cons(IDs)
  keys<-translate2_top(name,pdb$atom$resno)
  
  cons<-score_cons(keys,"cons",scores)
  cons[is.na(cons)]<-0
  pdb$atom$b<-cons
  write.pdb(pdb, file=path)
  
}
