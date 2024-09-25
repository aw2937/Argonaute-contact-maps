
#Read data
DF<-read.csv("data/reference_alignment_position_table.csv",row.names = 1)
phylo<-read.csv("data/phylo_all_sequences.csv",row.names = 1)
ref_structures<-read.csv("./data/Ago_structures.csv",row.names = 1)
align<-bio3d::read.fasta(file="data/reference_alignment.fasta")
scores_ori<-read.csv("data/reference_alignment_annotated_table.csv",row.names = 1,header=T)

file_names <- c("universal_sig.txt", "PIWI_sig.txt", "AGO_sig.txt", "eAGO_sig.txt", "pAGO_sig.txt","fAGO_sig.txt","plAGO_sig.txt","mAGO_sig.txt","vAGO_sig.txt")
df_list <- list()
for (file in file_names) {
  data <- read.table(paste("filter/", file, sep = ""), row.names = 1, header = TRUE)
  data <- as.data.frame(t(data))
  top <- rownames(data)
  data$top <- top
  data <- data.table(data)
  colnames(data) <- c("V1", "top")
  setkey(data, "top")
  df_list[[file]] <- data
}
suppressMessages(attach(df_list))

#Translates between topology and residues positions
translate2_top<-function(ali,keys){
  input<-keys
  if(nchar(ali)<5){ali<-ref_structures[grep(ali, rownames(ref_structures)), 1]}
  DT<-DF[grep(ali, rownames(DF)), ]
  
  if(dim(DT)[1]>1){print("More than 1 one match from ID");return(NA)}
  if(dim(DT)[1]==0){print("No match from ID");return(NA)}
  DT<-as.data.frame(t(DT))
  top<-rownames(DT)
  #Replace gaps with NAs for data.table
  suppressWarnings(DT<-as.data.frame(lapply(DT, function(x) as.numeric(as.character(x)))))
  DT$top<-top
  DT<-data.table(DT)
  colnames(DT)<-c("V1","top")
  
  if (length(keys[grep("A:", keys)])!=0){
    new_keys<-c()
      for (i in keys){
      new_keys<-append(new_keys,strsplit(i,split=':', fixed=TRUE)[[1]][[3]])
      }
    keys<-as.numeric(new_keys)
    setkey(DT,"V1")
    keys<-as.data.frame(keys)
    keys<-data.table(keys)
    translation<-DT[keys]$top
    if(any(is.na(translation))){print("Some values could not be translated")}
    return(translation)
  }
  
  if (is.numeric(keys)){
    setkey(DT,"V1")
    keys<-as.data.frame(keys)
    keys<-data.table(keys)
    translation<-DT[keys]$top
    if(any(is.na(translation))){print("Some values could not be translated")}
    return(translation)
  }
  
  if (is.character(keys)){
    setkey(DT,"top")
    keys<-as.data.frame(keys)
    keys<-data.table(keys)
    translation<-DT[keys]$V1
    if(any(is.na(translation))){print("Some values could not be translated")}
    return(translation)
  }
  
}

#Calculates conservation
get_cons<-function(IDs){
  align$ali<-align$ali[grep(paste(IDs,collapse="|"), align$id),]
  align$id<-align$id[grep(paste(IDs,collapse="|"), align$id)]
  score_bio3d<-conserv(align,method="similarity",sub.matrix = "blosum62",normalize.matrix = T)
  score_bio3d_df<-as.data.frame(list(score_bio3d,names(scores_ori)))
  score_bio3d_df<-data.table(score_bio3d_df)
  colnames(score_bio3d_df)<-c("V1","top")
  setkey(score_bio3d_df,"top")
  return(score_bio3d_df)
}

#Scores conservation or applies signature filter
score_cons<-function(keys,filter,scores){
  if (filter=="cons"){scores<-scores}
  #set filter for signature
  if (filter=="universal"){scores<-universal_sig.txt}
  if (filter=="PIWI"){scores<-PIWI_sig.txt}
  if (filter=="AGO"){scores<-AGO_sig.txt}
  if (filter=="eAgo"){scores<-eAGO_sig.txt}
  if (filter=="pAgo"){scores<-pAGO_sig.txt}
  if (filter=="animal AGO"){scores<-mAGO_sig.txt}
  if (filter=="vertebrate AGO"){scores<-vAGO_sig.txt}
  if (filter=="plant AGO"){scores<-plAGO_sig.txt}
  if (filter=="fungi AGO"){scores<-fAGO_sig.txt}
  
  
  keys<-as.data.frame(keys)
  keys<-data.table(keys)
  translation<-scores[keys]$V1
  if(any(is.na(translation))){print("Some values could not be translated")}
  return(translation)
}


