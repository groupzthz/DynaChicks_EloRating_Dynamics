# Social Network Analysis #
nodes <- data.frame(sort(unique(PenD$Loser)))
colnames(nodes)[1] <- "BirdID"
nodes$ID <- 1:nrow(nodes)
#change order
nodes = select(nodes, "ID", "BirdID")

rm(links)
src <- data.frame(PenD$Winner)
colnames(src)[1] <- "src"
src.1 <- merge(nodes, src, by.x = "BirdID", by.y = "src",all.y=T)
src.1$ID <- src.1$ID %>% replace_na(20)
colnames(src.1)[2] <- "src"

dt_add_row <- data.table::copy(nodes)                          # Replicate dt_orig
new_row    <- data.table("ID" = 20, "BirdID" = "GN")
dt_add_row <- rbindlist(list(dt_add_row, new_row))               # Row-Wisely combine dt_add_row and new_row        
nodes <- data.frame(dt_add_row)


target <- data.frame(PenD$Loser)
colnames(target)[1] <- "target"
target.1 <- merge(nodes, target, by.x = "BirdID", by.y = "target",all=T)
#target.1$ID <- target.1$ID %>% replace_na(20)
colnames(target.1)[2] <- "target"

links <- data.frame(src.1$src, target.1$target)

# get rid of duplicates and names connected with themselves
links <- data.frame(src, target) %>%  
  filter(!(src == target)) 
links <- unique(links[c("src", "target")])   
library(dplyr)
library(data.table)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(igraph)
library(CINNA)
str(nodes)
str(links)
set.seed(11)
social_net_tbls <- tbl_graph(nodes = nodes, 
                             edges = links, 
                             directed = T)
V(social_net_tbls)$size=degree(social_net_tbls, mode = "all")

social_net <- ggraph(social_net_tbls, layout = "graphopt") +                                                                                                         
  geom_node_point(size = 3) +                                         
  geom_node_text(aes(label = BirdID), nudge_y = 1.5, nudge_x = 1.5)+ 
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm')) +
  theme_void()
show(social_net)

degree.cent <- centr_degree(social_net_tbls, mode = "all")
degree.cent$res
# bird 9 highest centrality

closeness.cent <- closeness(social_net_tbls, mode="all")
closeness.cent
# again bird 9 highest

layout1 <- layout.fruchterman.reingold(social_net_tbls)
plot(social_net_tbls, edge.arrow.size=0.25,edge.arrow.mode = "-")

set.seed(123)
social_net_tbls %>%
  activate(nodes) %>%
  mutate(community = as.factor(group_infomap())) %>% 
  ggraph(layout = "graphopt") + 
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm'),width = 1, colour = "lightgray") +
  geom_node_point(aes(colour = community), size = 4) +
  geom_node_text(aes(label = BirdID), repel = TRUE)+
  theme_graph()


socwiki_OutDegree <- degree(social_net_tbls, mode = "out")
socwiki_OutDegree <- as.data.frame(socwiki_OutDegree)
degree(social_net_tbls, mode = "in")
socwiki_InDegree <- degree(social_net_tbls, mode = "in")
socwiki_InDegree <- as.data.frame(socwiki_InDegree)
set.seed(3952)

layout1 <- layout.fruchterman.reingold(social_net_tbls,niter=500)

V(social_net_tbls)$size=degree(social_net_tbls, mode = "all")/5
#V(social_net_tbls)$color <- ifelse(social_net_tbls[V(BirdID), 2] == "BirdID", "blue", "red")
E(social_net_tbls)$color <- "grey"
plot(social_net_tbls, edge.arrow.size=0.25,edge.arrow.mode = "-", vertex.label = NA)

socwiki_Betweeness <- betweenness(social_net_tbls)
socwiki_Betweeness <- as.data.frame(socwiki_Betweeness)

set.seed(3952)
soc_graph <- graph.data.frame(social_net_tbls, directed=TRUE)
layout1 <- layout.fruchterman.reingold(soc_graph,niter=500)

V(social_net_tbls)$size=betweenness(social_net_tbls) 

plot(social_net_tbls, edge.arrow.size=0.25,edge.arrow.mode = "-", vertex.label = NA)

closeness(social_net_tbls, mode="in")
socwiki_InCloseness <- closeness(social_net_tbls, mode="in")
socwiki_InCloseness <- as.data.frame(socwiki_InCloseness)
socwiki_OutCloseness <- closeness(social_net_tbls, mode="out")
socwiki_OutCloseness <- as.data.frame(socwiki_OutCloseness)

closeness(social_net_tbls, mode="out")
closeness(social_net_tbls, mode="all")
set.seed(3952)
layout1 <- layout.fruchterman.reingold(social_net_tbls,niter=500)
V(social_net_tbls)$size=closeness(socwiki_Graph, mode = "out")/.05
V(social_net_tbls)$color <- ifelse(nodes[V(social_net_tbls), 2] == "BirdID", "blue", "red")
E(social_net_tbls)$color <- "grey"

plot(social_net_tbls, edge.arrow.size=0.25,edge.arrow.mode = "-", vertex.label = NA)

