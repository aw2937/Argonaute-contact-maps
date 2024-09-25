DF<-read.csv("data/reference_alignment_position_table.csv",row.names = 1)

#Overlaps two contact lists
overlap_contacts<-function(contacts,plot,retain_values) {
  
  set1<-contacts[[1]]
  set2<-contacts[[2]]
  
  overlap<-inner_join(set1, set2, by = c('from'= 'from','to'='to'))
  overlap2<-inner_join(set1, set2, by = c('from'= 'to','to'='from'))
  overlap_final<-rbind(overlap2,overlap)
  overlap_switch<-overlap_final
  overlap_switch$to<-overlap_final$from
  overlap_switch$from<-overlap_final$to
  
  overlap_comb<-rbind(overlap_final,overlap_switch)
  
  set1_spec<-anti_join(set1, overlap_comb, by = c('from'= 'from','to'='to'))
  set2_spec<-anti_join(set2, overlap_comb, by = c('from'= 'from','to'='to'))
  
  
  source_structures<-paste0(overlap_final$source_structures.x,overlap_final$source_structures.y)
  source_structures<-lapply(str_split(source_structures,pattern=","), unique)
  source_structures<-unlist(lapply(source_structures, function(x) paste0(x, sep=",",collapse = "")))
  
  avg_contacts<-(overlap_final$avg_contacts.x+overlap_final$avg_contacts.y)/2
  value<-overlap_final$value.x+overlap_final$value.y
  num_struc<-overlap_final$num_struc.x+overlap_final$num_struc.y
  to<-overlap_final$to
  from<-overlap_final$from
  min_rep<-overlap_final$min_rep.y

  if(retain_values==TRUE){overlap_final_sorted<-data.frame(from,to,value,overlap_final$value.x,overlap_final$value.y,avg_contacts,source_structures,num_struc,min_rep)}
  else{overlap_final_sorted<-data.frame(from,to,value,num_struc,avg_contacts,source_structures)}
  
  result<-list(overlap=overlap_final_sorted,set1_specific=set1_spec,set2_specific=set2_spec)
  names(result)<-c(paste0(names(contacts)[1],"_", names(contacts)[2],"_overlap"),paste0(names(contacts),"_specific"))
        
  if(plot==TRUE){
  venn<-list(set1=set1, set2 = set2)
  names(venn)<-names(contacts)
  venn<-ggvenn(venn,fill_color = c("#B5B8BA","#B5B8BA"),stroke_size = 0.5, set_name_size = 4,fill_alpha = 0.3)+
    ggtitle("Overlap")+
    theme(plot.title = element_text(hjust = 0.5))
  
  par(mfrow=c(3,1))

  barplot(sort(table(c(overlap_final$pos1,overlap_final$pos2))),las=2,main = names(result)[1])
  barplot(sort(table(c(set1_final$pos1,set1_final$pos2))),las=2,main = names(result)[2])
  barplot(sort(table(c(set2_final$pos1,set2_final$pos2))),las=2,main = names(result)[3])
  }
  
  return(list(result))
}


#Stripchart plot
plot_top<-function(contacts) {
  
  #hago2 domains: N 53-139, L1 139-229, PAZ 229-347, L2 347-445, MID 445-580, PIWI 580-859
  pos<-translate2_top("Q9UKV8",c(1,139,229,347,445,580,859))
  pos<-match(pos,colnames(DF))
  
  annotations <- data.frame(
    X = pos[-7],
    Y =  c(-0, -0,-0,-0,-0, -0),
    text = c("N","L1","PAZ","L2","MID","PIWI"),
    x_adjust = c(0,0,0,0,0,0),
    y_adjust = c(2.5,2.5,2.5,2.5,2.5,2.5))
  
  if(is.data.frame(contacts)==F){
    positions<-match(contacts,colnames(DF))
    df<-data.frame(X=sort(positions),Y=rep(0))
    p1<-ggplot(df) +
      geom_segment(aes(x = 17, y = 0, xend = 721, yend = 0), color = "#D3D3D3", size=10)+  
      geom_segment(aes(x = 722, y = 0, xend = 1186, yend = 0), color = "#F08080", size=10)+  
      geom_segment(aes(x = 1187, y = 0, xend = 1866, yend = 0), color = "#9CBFA7", size=10)+  
      geom_segment(aes(x = 1867, y = 0, xend = 2491, yend = 0), color = "#FFB238", size=10)+  
      geom_segment(aes(x = 2492, y = 0, xend = 3251, yend = 0), color = "#ADD8E6", size=10)+  
      geom_segment(aes(x = 3252, y = 0, xend = 4369, yend = 0), color = "#C9F299", size=10)+  
      geom_text(data=annotations, aes(x=X,y=Y,hjust=x_adjust,vjust=y_adjust,label=text))+
      geom_segment(aes(x = X-5, y = 0, xend = X+5, yend = 0), color = "black", size=10)+
      coord_cartesian(clip = "off") + 
      theme_graph()+
      theme(legend.position = "none")
    return(p1)
    }
    
  
  if(typeof(contacts)=="list"){

  empty<-data.frame(from=colnames(DF),to=colnames(DF),value=rep(0),avg_contacts=rep(0),source_structures=rep(0),num_struc=rep(0),min_rep=rep(0))
  all_df<-rbind(contacts,empty)
  all_df<-all_df[order(match(all_df$from,colnames(DF))),]
  net <- graph.data.frame(all_df, directed=FALSE)
  
  edgelist <- get.edgelist(net)
  # get vertex labels
  label <- igraph::get.vertex.attribute(net, "name")
  # get vertex degree
  degrees <- igraph::degree(net)
  
  # data frame with groups, degree, labels and id
  nodes <- data.frame(degrees,label, id=1:vcount(net))
  nodes$label <- factor(nodes$label, levels = unique(nodes$label))
  nodes <- as_tibble(nodes)

  # prepare data for edges
  edges <- as_tibble(edgelist)
  
  net.tidy <- tbl_graph(nodes = nodes, edges = edges, directed = F, node_key = "label")
  
  ggraph(net.tidy, layout = "linear") + 
    geom_segment(aes(x = 17, y = 0, xend = 721, yend = 0), color = "#D3D3D3", size=10)+  
    geom_segment(aes(x = 722, y = 0, xend = 1186, yend = 0), color = "#F08080", size=10)+  
    geom_segment(aes(x = 1187, y = 0, xend = 1866, yend = 0), color = "#9CBFA7", size=10)+  
    geom_segment(aes(x = 1867, y = 0, xend = 2491, yend = 0), color = "#FFB238", size=10)+  
    geom_segment(aes(x = 2492, y = 0, xend = 3251, yend = 0), color = "#ADD8E6", size=10)+  
    geom_segment(aes(x = 3252, y = 0, xend = 4369, yend = 0), color = "#C9F299", size=10)+  
    geom_text(data=annotations, aes(
      x=X,y=Y,hjust=x_adjust,vjust=y_adjust,label=text))+
    geom_edge_arc(alpha = 0.5) + 
    scale_edge_width(range = c(0.2, 2)) +
    #scale_colour_manual(values= vrtxc) +
    geom_node_point(aes(size = degrees)) +
    #geom_node_text(aes(label = label), angle = 90, hjust = 1, nudge_y = -0.2, size = 4) +   
    coord_cartesian(clip = "off") + 
    theme_graph()+
    theme(legend.position = "none")
  }

  
}



