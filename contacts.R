#Load packages
packages<-c("shiny", "tidyverse","tidygraph", "data.table", "circlize", "network", "igraph", "networkD3", "chorddiag", "NGLVieweR", "shinyalert", "rotl", "ggtree", "ggnewscale","ggplot2","ggraph", "htmlwidgets", "ape","reshape2","bio3d")
suppressWarnings(suppressMessages(lapply(packages, require, character.only = TRUE)))

#Source functions
source("function/translator.R")
source("function/pdb_function_af2.R")
source("function/overlap_pos_contacts.R")

#Load data

#List of structures from the PDB
ref_structures<-read.csv("./data/Ago_structures.csv",row.names = 1)
ref_structures_dt<-ref_structures
ref_structures_dt$ID<-rownames(ref_structures_dt)

#Phylogenetic information of all Argonautes
phylo<-read.csv("./data/phylo_all_sequences.csv",row.names = 1)
#phylogenetic tree of Argonautes
tree <- read.tree("./data/tree.nwk")
#AF2 plDDT score of all Argonautes by residue
af2_score<-read.csv("./data/reference_alignment_af2_score_table.csv",row.names = 1)
#List of all available contacts (no plDDT threshold applied)
all<-read.csv("data/all_contacts.csv")

#List only taxa with at least 10 reps
all_phylo<-unlist(strsplit(phylo$Taxonomic.lineage..ALL.,", "))
phylo_freq<-table(all_phylo)
phylo_sel<-names(phylo_freq[phylo_freq>=10])

#Define UI for application
ui <- fluidPage(
 
    #Title
    titlePanel("Argonaute contacts"),

    #Sidebar architecture
    sidebarLayout(
        sidebarPanel(
            radioButtons(inputId="cons_type", label="Sequence conservation ", 
                       choices=c("signature positions","custom threshold")),
            uiOutput("cons_slider"),
            
            radioButtons(inputId="contact_type", label="Contact conservation", 
                         choices=c("SCN threshold","custom threshold")),
            uiOutput("contact_slider"),
            
            radioButtons(inputId="tax_sele", label="Argonaute subset selection", 
                         choices=c("signature-associated set","custom")),
            uiOutput("tax_pop"),
            uiOutput("filter_slider"),
            
            checkboxInput("show_all_contacts", "Show all contacts", value = TRUE),
            
            checkboxInput("only_interelement", "Include only inter-topology element contacts.", value = FALSE),
            
            checkboxInput("af2", "Include AlphaFold2 models", value = FALSE),
            uiOutput("af2slider"),
            
            checkboxInput("representative", "Include only one PDB representative per species", value = TRUE),
           
            selectInput(inputId = "pdb",
                        label = "Select structure (PDB code) for contact visualization",
                        choices = c('4OLB','5G5T','4N41','3DLB','5UX0','7KX7','7R8K','3DLH')),
            actionButton("click", "Let's go."),
        ),

        #Plotting
        mainPanel(
            tabsetPanel(
            tabPanel("Network",forceNetworkOutput("simple", width = "100%", height = "600px")),
            suppressMessages(tabPanel("Tree",plotOutput("top",height = "200px"),plotOutput("tree",height = "800px"))),
            tabPanel("Chordplot",chorddiagOutput("chord", height = 600)), 
            tabPanel("Structure",NGLVieweROutput("structure",height = 600)),
            tabPanel("Contacts",NGLVieweROutput("network",height = 600))))

    )
)

#Define server logic
server <- function(input, output, session) {
    
  #Set variables for all sessions/events
    dim<-NULL
    name<-NULL
    

    #Render UI depending on choices
    output$af2slider = renderUI({
      if (input$af2 == 0) {
        return(NULL)
      }
      if (input$af2 == 1) {
        list(
          sliderInput("af2thresh",
                      "AF2 plDDT score threshold",
                      min = 0,
                      max = 100,
                      value = 80)
        )
      }
      
    })
    
    output$cons_slider = renderUI({
      if (input$cons_type == "signature positions") {
        selectInput(inputId = "filter",
                    label = "Select sequence signature",
                    choices = c("universal","pAgo","eAgo","PIWI","AGO","fungi AGO","plant AGO","animal AGO","vertebrate AGO"),
                    selected = "universal")
      }
      
      else {
        sliderInput("cons",
                    "Conservation",
                    min = 0,
                    max = 1,
                    value = 0.0)
        
      }
    })
    
    output$contact_slider = renderUI({
      if (input$contact_type == "SCN threshold") {
          return(NULL)
      }
      
      else {
        sliderInput("contact",
                    "Contact",
                    min = 0,
                    max = 1,
                    value = 0.5)
        
      }
    })
    
    
    output$tax_pop = renderUI({
      
      if (input$tax_sele == "signature-associated set") {
      return(NULL)
      
      }
      else {
        radioButtons(inputId="filter_type", label="Select filter type", 
                     choices=c("Phylogeny","Clades"),
                     selected = "Phylogeny")
      }})
    
    
    output$filter_slider = renderUI({
      if (input$tax_sele == "signature-associated set") {
        return(NULL)
      }
      
      
      if (input$tax_sele == "custom") {
      
      if (input$filter_type == "Clades") {
        selectInput(inputId = "branch",
                    label = "Clade",
                    choices = c("PIWI","AGO","animal AGO","animal PIWI","long-B pAgo"),
                    selected = "PIWI")
      }
      else {
        selectInput(inputId = "tax",
                    label = "Phylogeny",
                    choices = c("Prokaryota",phylo_sel),
                    selected = "cellular organisms")
        
      }
      }

    })
    
    #Calculation
    observeEvent(input$click, {

    IDs<-c()
    #Set type of conservational filter (e.g. cons,eAgo)
    if(input$cons_type == "custom threshold"){conservation_thresh<-input$cons;filter<-"cons"}
    if(input$cons_type == "signature positions"){conservation_thresh<-0.1;filter=input$filter}
    
    #Set 1xSD contact conservation for respective group (SCN)
    if(input$contact_type == "SCN threshold"){
      if(input$filter=="eAgo"){contact_thresh<-0.34}
      else if(input$filter=="pAgo"){contact_thresh<-0.2}
      else if(input$filter=="PIWI"){contact_thresh<-0.58}
      else if(input$filter=="AGO"){contact_thresh<-0.4}
      else if(input$filter=="plant AGO"){contact_thresh<-0.78}
      else if(input$filter=="fungi AGO"){contact_thresh<-0.49}
      else if(input$filter=="animal AGO"){contact_thresh<-0.99}
      else if(input$filter=="vertebrate AGO"){contact_thresh<-0.99}
      else if(input$filter=="universal"){contact_thresh<-0.34} 
    }
    
    #Set custom contact conservation
    if(input$contact_type != "SCN threshold"){
      contact_thresh<-input$contact}
    
    #Check if tax-sele SCN without signature
    if(!exists("contact_thresh")){shinyalert("Input Error", "Cannot calculate SCN, because no signature was selected.", type = "error")}
    validate(need(exists("contact_thresh"), "Please select a data set"))
    
    #PDB structure for contact visualization
    pdb_ID<-input$pdb
    input_pdb <- ref_structures %>% filter(row.names(ref_structures) %in% c(pdb_ID))
    name<<-input_pdb$uniprot
    if (str_length(input_pdb$mainchain)>1){mainchain<-strsplit(input_pdb$mainchain,split=',', fixed=TRUE)[[1]][1];guidechain<-strsplit(input_pdb$chain_guide,split=',', fixed=TRUE)[[1]][1]}
    else{mainchain<-input_pdb$mainchain;guidechain<-input_pdb$chain_guide}
    
    #Parse af2 structures
    af2<-input$af2
    af2_thresh<-0
    if(af2==TRUE){af2_thresh<-input$af2thresh}
    only_interelement<-input$only_interelement
    
    #Select set of proteins by signature
    if (input$tax_sele == "signature-associated set") {
      tax<-input$filter
      if(input$filter=="universal"){IDs<-rownames(phylo[grep("cellular organism", phylo$Taxonomic.lineage..ALL.),])}
      if(input$filter=="eAgo"){IDs<-rownames(phylo[grep("Eukaryota", phylo$Taxonomic.lineage..ALL.),])}
      if(input$filter=="pAgo"){IDs<-rownames(phylo[-grep("Eukaryota", phylo$Taxonomic.lineage..ALL.),])}
      if(input$filter=="PIWI"){subtree<-extract.clade(tree,node = 526); subtree$tip.label<-sub("\\_.*", "", subtree$tip.label);IDs<-subtree$tip.label;IDs[str_length(IDs)<5]<-setDT(ref_structures_dt, key = 'ID')[J(toupper(IDs[str_length(IDs)<5]))]$uniprot}
      if(input$filter=="AGO"){subtree<-extract.clade(tree,node = 380);  subtree$tip.label<-sub("\\_.*", "", subtree$tip.label);IDs<-subtree$tip.label;IDs[str_length(IDs)<5]<-setDT(ref_structures_dt, key = 'ID')[J(toupper(IDs[str_length(IDs)<5]))]$uniprot}
      if(input$filter=="fungi AGO"){IDs<-rownames(phylo[grep("Fungi", phylo$Taxonomic.lineage..ALL.),])}
      if(input$filter=="plant AGO"){IDs<-rownames(phylo[grep("Viridiplantae", phylo$Taxonomic.lineage..ALL.),])}
      if(input$filter=="animal AGO"){subtree<-extract.clade(tree,node = 392);  subtree$tip.label<-sub("\\_.*", "", subtree$tip.label);IDs<-subtree$tip.label;IDs[str_length(IDs)<5]<-setDT(ref_structures_dt, key = 'ID')[J(toupper(IDs[str_length(IDs)<5]))]$uniprot}
      if(input$filter=="vertebrate AGO"){
      subtree<-extract.clade(tree,node = 392)
      subtree$tip.label<-sub("\\_.*", "", subtree$tip.label)
      IDs<-subtree$tip.label
      IDs[str_length(IDs)<5]<-setDT(ref_structures_dt, key = 'ID')[J(toupper(IDs[str_length(IDs)<5]))]$uniprot
      IDs2<-rownames(phylo[grep("Vertebrata", phylo$Taxonomic.lineage..ALL.),])
      IDs<-IDs2[IDs2%in% IDs]
      }
    }
    if (input$tax_sele == "custom") {
      
    #Select Argonautes based on Phylogeny
    if (input$filter_type == "Phylogeny"){
    tax<-input$tax
    if(tax=="Prokaryota"){IDs<-rownames(phylo[-grep("Eukaryota", phylo$Taxonomic.lineage..ALL.),])}
    if(tax!="Prokaryota"){IDs<-rownames(phylo[grep(tax, phylo$Taxonomic.lineage..ALL.),])}
    }

    #Select Argonautes based on branch in tree (Clades)
    if (input$filter_type == "Clades"){
      tax<-input$branch
      if(tax=="PIWI"){subtree<-extract.clade(tree,node = 526)}
      if(tax=="AGO"){subtree<-extract.clade(tree,node = 380)}
      if(tax=="long-B pAgo"){subtree<-extract.clade(tree,node = 570)}
      if(tax=="animal AGO"){subtree<-extract.clade(tree,node = 392)}
      if(tax=="animal PIWI"){subtree<-extract.clade(tree,node = 536)}
      
      subtree$tip.label<-sub("\\_.*", "", subtree$tip.label)
      IDs<-subtree$tip.label
      IDs[str_length(IDs)<5]<-setDT(ref_structures_dt, key = 'ID')[J(toupper(IDs[str_length(IDs)<5]))]$uniprot


      }
    }
    
    #If no IDs have been selected calculation cannot start
    if(length(IDs) == 0){shinyalert("Input Error", "Cannot calculate on signature-associated set, because no signature was selected.", type = "error")}
    validate(need(length(IDs) != 0, "Please select a data set"))
    
    #Selection annotation for tree
    matches_sel <- unique(grep(paste(IDs,collapse="|"),tree$tip.label, value=TRUE))
    Taxon<-rep(tax,length(matches_sel))
    selected_tree <- data.frame(Taxon)
    rownames(selected_tree) <- sub("\\_.*", "", matches_sel)

    #Tree annotation
    kingdom_anno<-c()
    matches_all<-c()
    for(king in c("Archaea","Eukaryota","Bacteria")){
      anno<-rownames(phylo[grep(king, phylo$Taxonomic.lineage..ALL.),])
      anno[str_length(anno)<5]<-setDT(ref_structures_dt, key = 'ID')[J(toupper(anno[str_length(anno)<5]))]$uniprot
      matches <- unique(grep(paste(anno,collapse="|"),tree$tip.label, value=TRUE))
      kingdom_anno <- rbind(kingdom_anno,data.frame(tax=rep(king,length(matches))))
      matches_all<-append(matches_all,matches)
    }
    rownames(kingdom_anno) <- sub("\\_.*", "", matches_all)
    tree$tip.label<-sub("\\_.*", "", tree$tip.label)
    circ<-ggtree(tree, layout="circular",branch.length = "none")+ geom_tiplab(size=2)
    

    #Get PDB structures included in selection
    ref_structures_sub<-ref_structures[ref_structures$uniprot %in% IDs, ]
    
    #Get conservation scores
    if (conservation_thresh!=0){scores<-get_cons(IDs)}
    
    chains<-substr(ref_structures_sub$mainchain, 1,1)
    Ago_list_red_sub<-list()
    #Prepare contacts_df
    contacts_df<-as.data.frame(matrix(ncol = 3))
    colnames(contacts_df)<-c("from","to","value")
    
    #Remove redundant PDB structures and retain one per species
    if (input$representative ==TRUE){ 
      original_structures<-c("4F1N","4F3T","6OON","5VM9","4KRE", "5GUH","6KR6","7KX7","7SWF","3DLB","6QZK","5AWH","5I4A","1YVU", "7R8K","5G5T","1Z26")
      ref_structures_sub<-subset(ref_structures_sub,rownames(ref_structures_sub) %in% original_structures)
    }
    
    pdbIDs<-rownames(ref_structures_sub)
    print("Parsing PDB structures")
    
    #Parse contacts from PDBs 
    i<-0
    for (y in pdbIDs){
        i<-i+1
        chain<-chains[i]
        if (file.exists(paste0("data/pdb/red_CN/",y,".txt"))){
        Ago_list_red_sub[[y]]<-read.delim(paste0("data/pdb/red_CN/",y,".txt"), comment.char = "#",header=T)
        Ago_list_red_sub[[y]]$struc_name<-rep(y)
        }
    }
    #Read contacts from AF2 models 
    if (af2==TRUE){
      print("Parsing AlphaFold2 models and filtering against threshold")
        for (i in IDs){
            path<-paste0("data/af2/red_CN/",i,".txt")
            if (file.exists(path)){Ago_list_red_sub[[i]]<-read.delim(path, comment.char = "#",header=T)
            Ago_list_red_sub[[i]]$struc_name<-rep(i)
            }
            #Filter residues by af2 score
            if(af2_thresh!=0){
            af2_score_filt<-af2_score[grep(i, rownames(af2_score)), ]
            above_thresh<-colnames(af2_score_filt[,af2_score_filt>af2_thresh])
            Ago_list_red_sub[[i]]<-Ago_list_red_sub[[i]][Ago_list_red_sub[[i]]$pos1%in%above_thresh & Ago_list_red_sub[[i]]$pos2%in%above_thresh,]
            }
        }
    }
    
    #Check if enough structures/models for calculation
    if(length(Ago_list_red_sub)==0){shinyalert("Input Error", "No structures found. Try to include predicted AlphaFold2 models to increase sample size.", type = "error")}
    validate(need(length(Ago_list_red_sub)!=0, "Please select a data set"))
    
    #Define minimal representation for contact threshold
    min_rep<-round(contact_thresh*length(Ago_list_red_sub))
    
    #Merge lists of contacts
    fuse<-c()
    for (i in 1:length(Ago_list_red_sub)){
        fuse<-rbind(fuse,Ago_list_red_sub[[i]])
    }
    #Contact accounting and summary
    red = fuse %>% 
    mutate( V1 = pmin(pos1, pos2), 
            V2 = pmax(pos1, pos2) ) %>%
    group_by(V1,V2) %>% 
    summarise(num_structures=n(),avg_contacts=sum(num_contacts)/n(),.groups = 'drop',source_structures=paste0(struc_name,",", collapse = ""))
    red<-unique(red)
    colnames(red)[c(1,2)]<-c("pos1","pos2")
    
    #Apply contact conservation threshold
    cons<-as.data.frame(subset(red,num_structures>=min_rep))
    
    #Check if contacts exist at contact conservation threshold.
    if(length(cons$pos1)==0){shinyalert("Input Error", "Contact threshold too high. No contacts found at this threshold.", type = "error")}
    validate(need(length(cons$pos1)!=0, "Please select a data set"))
    
    #Apply sequence conservation threshold
    retain<-c()
    cons<-na.omit(cons)
    if (conservation_thresh!=0){
      print("Filtering by conservation")
      for (i in 1:dim(cons)[1]){
          if(input$show_all_contacts==FALSE){
              if(score_cons(cons[i,1],filter,scores)>conservation_thresh & score_cons(cons[i,2],filter,scores)>conservation_thresh){retain<-rbind(retain, cons[i,])}
          }
          if(input$show_all_contacts==TRUE){
              if(score_cons(cons[i,1],filter,scores)>conservation_thresh | score_cons(cons[i,2],filter,scores)>conservation_thresh){retain<-rbind(retain, cons[i,])}
          }
      }
    }
    if (conservation_thresh==0){retain<-cons}
    
    #Check if contacts exist at sequence conservation threshold.
    if(length(retain$pos1)==0){shinyalert("Input Error", "Conservation threshold too high. No contacts found at this threshold.", type = "error")}
    validate(need(length(retain$pos1)!=0, "Please select a data set"))
    
    
    shared_cons_contacts<-as.data.frame(retain)
    shared_cons_contacts<-na.omit(shared_cons_contacts)
      
    #Filter Contacts for inter-SSE contacts
      if(only_interelement==TRUE){
        filt_shared_cons_contacts<-c()
        for (i in 1:length(shared_cons_contacts$pos1)){
          group_pos1<-paste0(strsplit(shared_cons_contacts[i,1],split='.', fixed=TRUE)[[1]][1],".",strsplit(shared_cons_contacts[i,1],split='.', fixed=TRUE)[[1]][2])
          group_pos2<-paste0(strsplit(shared_cons_contacts[i,2],split='.', fixed=TRUE)[[1]][1],".",strsplit(shared_cons_contacts[i,2],split='.', fixed=TRUE)[[1]][2])
          if(group_pos1==group_pos2){next}
          filt_shared_cons_contacts<-rbind(filt_shared_cons_contacts,shared_cons_contacts[i,])
        }
      shared_cons_contacts<-filt_shared_cons_contacts
      #Check if inter-SSE contacts exist at sequence conservation threshold.
      if(length(shared_cons_contacts$pos1)==0){shinyalert("Input Error", "No inter-SSE element contacts found at this threshold.", type = "error")}
      validate(need(length(shared_cons_contacts$pos1)!=0, "Please select a data set"))
    }
    
    contact_number<-dim(shared_cons_contacts)[1]
    contacts_df<-shared_cons_contacts
    colnames(contacts_df)<-c("from","to","value","avg_contacts","source_structures")
    contacts_df <- contacts_df[order(contacts_df$value),]
    contacts_df$num_struc<-length(Ago_list_red_sub)
    contacts_df$min_rep<-min_rep
    
    #Create matrix for chordplot
    ma_contacts_df<-contacts_df
    nameVals <- sort(unique(unlist(ma_contacts_df[1:2])))
    myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
    myMat[as.matrix(ma_contacts_df[c("from", "to")])] <- ma_contacts_df[["value"]]
    myMat[as.matrix(ma_contacts_df[c("to", "from")])] <- ma_contacts_df[["value"]]
    vals<-as.vector(colSums(myMat))
    normalized = (vals-min(vals))/(max(vals)-min(vals))
    #Define colour palette
    pal <- colorRamp(c("grey", "red"))
    
    
    #Precalculated universal SCN (see paper)
    if(input$filter=="universal" & input$contact_type == "SCN threshold"){contacts_df<-read.csv("data/uni_SCN.csv")} 
    if(input$contact_type == "SCN threshold"){contact_thresh<-"SCN"}
    if(input$cons_type == "signature positions"){conservation_thresh<-"SIG"}
    contact_number<-dim(contacts_df)[1]
    contact_network_csv<-paste0("data_output/",tax,"_",pdb_ID,"_CN_",contact_thresh,"_CONS_",conservation_thresh,"_contacts_",contact_number,"_contact_network.csv")
    write_csv(contacts_df,file = contact_network_csv)
    
    write_conservation_pdb(scores,name,tax,pdb_ID,contact_thresh,conservation_thresh,contact_number,mainchain,IDs)
    
  
    #Write output files   
    file<-paste0("data/pdb/all_structures/","pdb",pdb_ID,".pdb")
    pos<-translate2_top(name, unique(c(contacts_df$from, contacts_df$to)))
    sele<-paste(as.character(pos), collapse=" or ")
    contact_network_pdb<-paste0("data_output/",tax,"_",pdb_ID,"_CN_",contact_thresh,"_CONS_",conservation_thresh,"_contacts_",contact_number,"_contact_network.pdb")
    contact_network_html<-paste0("data_output/",tax,"_",pdb_ID,"_CN_",contact_thresh,"_CONS_",conservation_thresh,"_contacts_",contact_number,"_contact_network.html")
    write_contacts_pdb(contacts_df,name,tax,pdb_ID,contact_thresh,conservation_thresh,contact_number,mainchain)
    conservation_pdb<-paste0("data_output/",tax,"_",pdb_ID,"_CN_",contact_thresh,"_CONS_",conservation_thresh,"_contacts_",contact_number,"_conservation.pdb")
    
    print(paste("File written: ", contact_network_pdb))
    print(paste("File written: ", contact_network_csv))
    print(paste("File written: ", contact_network_html))
    print(paste("File written: ", conservation_pdb))
    
    #Save old contacts_df
    contacts_df_o<-read.csv(contact_network_csv)
    
    print("DONE!")
    
    
    
    #Plotting
    
    data<-list(all=all,
               select=contacts_df_o
    )
    result<-overlap_contacts(data,FALSE,FALSE)
    
    #Get all structures for tree heatmap
    select<-result[[1]]$all_select_overlap
    all_contacts<-table(unlist(str_split(select$source_structures,pattern=",")))/dim(select)[1]
    all_contacts_df<-as.data.frame(all_contacts,row.names = T)
    rownames(all_contacts_df)[str_length(rownames(all_contacts_df))<5]<-tolower(rownames(all_contacts_df)[str_length(rownames(all_contacts_df))<5])
    
    
    #Create force network visualization
    nodes<-unique(c(contacts_df$from,contacts_df$to))
    nodes_df<-data.frame(nodes)
    group<-c()
    #Create node table
    for (i in nodes){
        group<-append(group,paste0(strsplit(i,split='.', fixed=TRUE)[[1]][1],".",strsplit(i,split='.', fixed=TRUE)[[1]][2]))
    }
    nodes_df<-cbind(nodes_df,group)
    nodes_df$node_size<-rep(10)
    
    nodes_df<-data.table(nodes_df)
    setkey(nodes_df,"nodes")
    nodes_df$index<-c(0:(length(nodes)-1))
    
    contacts_df$from_index<-nodes_df[contacts_df$from]$index
    contacts_df$to_index<-nodes_df[contacts_df$to]$index
    #Weigh connector size by average contacts conservation
    contacts_df$con_viz<-(contacts_df$value/contacts_df$num_struc)*5
    MyClickScript <- "Shiny.onInputChange('nodename', d.name);"
    fn<-forceNetwork(Links = contacts_df, Nodes = nodes_df,
                 Source = "from_index", Target = "to_index",
                 Value = "con_viz", NodeID = "nodes",
                 Group = "group", opacity = 0.8,
                 opacityNoHover = TRUE, 
                 fontSize = 16,
                 linkDistance = 50, 
                 Nodesize = "node_size",
                 zoom = T,
                 charge = -30,
                 clickAction = MyClickScript
    )
    output$simple <- renderForceNetwork({fn})
    saveWidget(fn, file=contact_network_html,selfcontained=T)
    
    #Create chord diagram
    output$chord <- renderChorddiag({
            chorddiag(myMat, groupColors = rgb(pal(normalized)/255), showTicks = F, clickAction = "Shiny.onInputChange('sourceIndex', d.source.index+1);
                                 Shiny.onInputChange('targetIndex', d.target.index+1);", clickGroupAction="Shiny.onInputChange('groupIndex', d.index+1)")})
    
    #Save dim for chord selection
    dim<<-nameVals
    output$shiny_return <- renderPrint({paste0("Clicked chord: ", dim[input$sourceIndex], " <-> ",dim[input$targetIndex])})
    
    #Render structure
    output$structure <- renderNGLVieweR({
            NGLVieweR(conservation_pdb) %>%
            stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
            addRepresentation("cartoon",param = list(name = "cartoon", colorScheme= "bfactor", sele = paste0(":",mainchain))
            )%>%
            zoomMove(zoom = paste0(":",mainchain), center = paste0(":",mainchain)
            )%>%
            addRepresentation("ball+stick",param = list(name = "cartoon", color = "orange", sele = paste0(":",guidechain))
            )%>%
            addRepresentation("ball+stick",
                              param = list(
                                  name = "interacting",
                                  color = "grey",
                                  sele = paste0(":",mainchain," and (",sele, ")" ))
            
            ) %>%
            addRepresentation("ball+stick",
                              param = list(
                                  name = "selection",
                                  color = "green",
                                  sele = "none" )
                              
            ) %>%
            addRepresentation("contact",
                              param = list(
                                  name = "contact",
                                  sele = paste0(":",mainchain," and (",sele, ")" ),
                                  filterSele = list(pos,pos),
                                  labelVisible = TRUE,
                                  hydrophobic=TRUE,
                                  maxHydrophobicDist= 5.0)
            )
    })
    
    #Draw tree
    output$tree <- renderPlot({
      p1 <- gheatmap(circ, kingdom_anno, offset=5, width=.2,colnames = F) +
      scale_fill_manual(name="Kingdom",na.value = "grey90",values = c("#E3998C", "#C84630", "#456990"))
      p2 <- p1 + new_scale_fill()
      p3<-gheatmap(p2, selected_tree, offset=13, width=.1,colnames = F) +
      scale_fill_viridis_d(option="A", name="Selected Kingdom")
      p4<- p3 + new_scale_fill()
      gheatmap(p4, all_contacts_df, offset=15, width=.2,colnames = F) +
      scale_fill_gradient(low = "white", high = "#3FA261",name="Contact conservation",limits = c(0, 1))})
    output$top <- renderPlot({
      plot_top(contacts_df_o)
    })
    
    #Render network
    output$network <- renderNGLVieweR({
        NGLVieweR(contact_network_pdb) %>%
        stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
        addRepresentation("cartoon",param = list(name = "cartoon", color = "grey",sele = paste0(":",mainchain))
        )%>%
        zoomMove(zoom = paste0(":",mainchain), center = paste0(":",mainchain)
        )%>%
        addRepresentation("licorice",
                          param = list(
                              name = "interacting",
                              colorValue = "red",
                              sele = paste0(":",mainchain," and (",sele, ")" ))
        )%>%
        addRepresentation("licorice",
                            param = list(
                                name = "selection",
                                colorValue = "green",
                                sele = "none")
                          
        )
      })
    })
    
    #Update zoom when clicked on chord
    observeEvent(input$sourceIndex, {
        select<-c(translate2_top(name,dim[input$sourceIndex]),translate2_top(name,dim[input$targetIndex]))
        select<-paste(as.character(select), collapse=" or ")
        NGLVieweR_proxy("structure") %>%
            updateZoomMove(
                center = toString(select),
                zoom = toString(select),
                z_offSet = -20
                             )
        NGLVieweR_proxy("structure") %>%
            updateSelection(name = "selection", sele = select)
        
        NGLVieweR_proxy("network") %>%
            updateZoomMove(
                center = select,
                zoom = select,
                z_offSet = -20
            )
        NGLVieweR_proxy("network") %>%
            updateSelection(name = "selection", sele = select)
            
    })
    observeEvent(input$groupIndex, {
        select<-c(translate2_top(name,dim[input$groupIndex]))
        select<-paste(as.character(select), collapse=" or ")
        NGLVieweR_proxy("structure") %>%
            updateZoomMove(
                center = toString(select),
                zoom = toString(select),
                z_offSet = -20
            )
        NGLVieweR_proxy("structure") %>%
            updateSelection(name = "selection", sele = select)
        
        NGLVieweR_proxy("network") %>%
            updateZoomMove(
                center = select,
                zoom = select,
                z_offSet = -20
            )
        NGLVieweR_proxy("network") %>%
            updateSelection(name = "selection", sele = select)
        
    })
    observeEvent(input$nodename, {
        select<-c(translate2_top(name,input$nodename),translate2_top(name,input$nodename))
        select<-paste(as.character(select), collapse=" or ")
        NGLVieweR_proxy("structure") %>%
            updateZoomMove(
                center = toString(select),
                zoom = toString(select),
                z_offSet = -20
            )
        NGLVieweR_proxy("structure") %>%
            updateSelection(name = "selection", sele = select)
        
        NGLVieweR_proxy("network") %>%
            updateZoomMove(
                center = select,
                zoom = select,
                z_offSet = -20
            )
        NGLVieweR_proxy("network") %>%
            updateSelection(name = "selection", sele = select)
        
    })
    }

# Run the application 
shinyApp(ui = ui, server = server)
