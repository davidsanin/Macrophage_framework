# ==============================================================================
# =                                                                            =
# =                   Predictive Framework of Macrophage                       =
# =                                Activation                                  =
# =                                                                            =
# =                            David E Sanin 2021                              =
# =                             Data processing                                =
# ==============================================================================

# Code includes
# 1. Pseudotime calculation
# 2. Gene expression over pseudotime
# 3. Benchmarking label transfer anchors
# 4. Transfer data from reference dataset
# 5. Examine label assignment
# 6. Differential gene expression analysis in query datasets
# 7. Indexed data dimensionality reduction and clustering
# 8. Gene expression network
# 9. Join all data with no supervision

# ============= 1. Pseudotime calculation ============= ####
require(Slingshot)
data #Reference fully processed Seurat object
#this takes the dimensionality reduction and clustering as an input
cl <- Idents(data)
rd <- Embeddings(object = data, reduction = "umap")[,1:2]

lin <- getLineages(rd, cl, start.clus = '4') #retrieve lineage breaking points
lin <- getCurves(lin) #calculate curves and pseudotime

data@meta.data <- cbind(data@meta.data,slingPseudotime(lin)) #Add pseudotime calculation to Seurat object

# ============= 2. Gene expression over pseudotime ============= ####
require(gam)

data #Reference fully processed Seurat object
#the pseudotime calculations are in columns labelled as "Path_#"
pst <- data@meta.data[,grep('Path',colnames(data@meta.data))]
sets <- colnames(pst)
pval <- list()
for(i in sets){
  cells <- rownames(pst[!is.na(pst[,i]),]) #get cells for each path
  sub <-  subset(data, cells = cells) #subset object to retain onle relevant cells
  df <- GetAssayData(object = sub, assay = 'integrated', slot = 'scale.data') #retrieve data

  var2k <- names(sort(apply(df, 1, var),decreasing = TRUE))[1:2000] #recover top 2,000 most variable genes within Path cells
  df <- df[var2k, ]  # subset counts to only highly variable

  # Fit GAM for each gene using pseudotime as independent variable, retrieving coefficient p value
  t <- FetchData(object = sub, vars = i)
  gam.pval <- apply(df, 1, function(z){
    d <- data.frame(z=z, t=t)
    tmp <- gam(z ~ lo(t), data=d)
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  pval[[i]] <- gam.pval # Store p values per gene and path
}

for(i in sets){
  pval[[i]] <- p.adjust(pval[[i]], method = 'fdr') # Correct p values for multiple comparissons
}

#  ============= 3. Benchmarking label transfer parameters ============= ####
n1 #Fully processed dataset of mixed macrophages and other immune cells
data #Reference fully processed Seurat object

n1 <- AddModuleScore(object = n1, name = 'Mac',
                     features = list('Mac' = c('H2-Ab1','Lyz2',
                                               'Csf1r','Adgre1',
                                               'Mertk','Cd164',
                                               'Cd68','Itgam'))) #identify macrophages

#Benchmarking anchors
#k.filter
an <- list()
for(i in c(5,25,50,100,150,200,250,300)){
  t <- FindTransferAnchors(reference = data,
                           query = n1,
                           dims = 1:50,
                           normalization.method = 'SCT',
                           npcs = 50, k.filter = i)
  an[[i]] <- t@anchors
}

o <- data.frame(an[[300]]) %>%
  left_join(data.frame(an[[250]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[200]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[150]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[100]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[50]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[25]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[5]]), by = c('cell1','cell2'))

colnames(o) <- c('cell1','cell2',300,250,200,150,100,50,25,5)
t <- melt(o, id.vars = c('cell1','cell2'))

ggplot(t, aes(x = value, color = variable)) +
  geom_density() +
  theme_classic() +
  scale_color_manual(values = color_hue(8)) +
  guides(color = guide_legend(title = 'k.filter')) +
  labs(x = 'Anchor score', y = 'Density') +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10,
                                 colour = 'black'))

#max.features
an <- list()
for(i in c(50,100,200,400,800)){
  t <- FindTransferAnchors(reference = data,
                           query = n1,
                           dims = 1:50,
                           normalization.method = 'SCT',
                           npcs = 50, k.filter = 5,
                           max.features = i)
  an[[i]] <- t@anchors
}

o <- data.frame(an[[800]]) %>%
  left_join(data.frame(an[[400]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[200]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[100]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[50]]), by = c('cell1','cell2'))

colnames(o) <- c('cell1','cell2',800,400,200,100,50)
t <- melt(o, id.vars = c('cell1','cell2'))

ggplot(t, aes(x = value, color = variable)) +
  geom_density() +
  theme_classic() +
  scale_color_manual(values = color_hue(5)) +
  guides(color = guide_legend(title = 'max.features')) +
  labs(x = 'Anchor score', y = 'Density') +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10,
                                 colour = 'black'))

#k.anchors
an <- list()
for(i in c(5,10,15,20,25)){
  t <- FindTransferAnchors(reference = data,
                           query = n1,
                           dims = 1:50,
                           normalization.method = 'SCT',
                           npcs = 50, k.filter = 5,
                           max.features = 100, k.anchor = i)
  an[[i]] <- t@anchors
}

o <- data.frame(an[[25]]) %>%
  left_join(data.frame(an[[20]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[15]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[10]]), by = c('cell1','cell2')) %>%
  left_join(data.frame(an[[5]]), by = c('cell1','cell2'))

colnames(o) <- c('cell1','cell2',25,20,15,10,5)
t <- melt(o, id.vars = c('cell1','cell2'))

ggplot(t, aes(x = value, color = variable)) +
  geom_density() +
  theme_classic() +
  scale_color_manual(values = color_hue(5)) +
  guides(color = guide_legend(title = 'k.anchor')) +
  labs(x = 'Anchor score', y = 'Density') +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10,
                                 colour = 'black'))

#Based on these, selected less features (100) with smaller neighborhoods (5)
#Benchmarking labels
a <- FindTransferAnchors(reference = data,
                         query = n1,
                         dims = 1:50,
                         normalization.method = 'SCT',
                         npcs = 50, k.filter = 5,
                         max.features = 100, k.anchor = 5)

#k.weight
an <- list()
for(i in c(5,10,25,50,75,100)){
  p <- TransferData(anchorset = a,
                    refdata = data$stage,
                    dims = 1:30, k.weight = i,
                    sd.weight = 1)
  p <- p[,c(1,ncol(p))]
  colnames(p) <- c('stage','stage.score')
  an[[i]] <- p
}

o <- data.frame(rownames_to_column(an[[5]])) %>%
  left_join(data.frame(rownames_to_column(an[[10]])),
            by = c("rowname",'stage')) %>%
  left_join(data.frame(rownames_to_column(an[[25]])),
            by = c("rowname",'stage')) %>%
  left_join(data.frame(rownames_to_column(an[[50]])),
            by = c("rowname",'stage')) %>%
  left_join(data.frame(rownames_to_column(an[[75]])),
            by = c("rowname",'stage')) %>%
  left_join(data.frame(rownames_to_column(an[[100]])),
            by = c("rowname",'stage'))

colnames(o) <- c('Cell','Cell Type',5,10,25,50,75,100)
t <- melt(o)

ggplot(t, aes(x = value, color = variable)) +
  geom_density() +
  theme_classic() +
  guides(color = guide_legend(title = 'k.weight')) +
  scale_color_manual(values = color_hue(6)) +
  labs(x = 'Label probability', y = 'Density') +
  geom_vline(xintercept = 0.8,
             color = 'blue', linetype = 'dashed') +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10,
                                 colour = 'black'))

#sd.weight
an <- list()
for(i in c(1,2,3,4)){
  p <- TransferData(anchorset = a,
                    refdata = data$stage,
                    dims = 1:30, k.weight = 75,
                    sd.weight = i)
  p <- p[,c(1,ncol(p))]
  colnames(p) <- c('stage','stage.score')
  an[[i]] <- p
}

o <- data.frame(rownames_to_column(an[[1]])) %>%
  left_join(data.frame(rownames_to_column(an[[2]])),
            by = c("rowname",'stage')) %>%
  left_join(data.frame(rownames_to_column(an[[3]])),
            by = c("rowname",'stage')) %>%
  left_join(data.frame(rownames_to_column(an[[4]])),
            by = c("rowname",'stage'))

colnames(o) <- c('Cell','Cell Type',1,2,3,4)
t <- melt(o)

ggplot(t, aes(x = value, color = variable)) +
  geom_density() +
  theme_classic() +
  guides(color = guide_legend(title = 'sd.weight')) +
  scale_color_manual(values = color_hue(6)) +
  labs(x = 'Label probability', y = 'Density') +
  geom_vline(xintercept = 0.8,
             color = 'blue', linetype = 'dashed') +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10,
                                 colour = 'black'))

#Based on these, selected neighborhood of size 75

# ============= 4. Transfer data from reference dataset ============= ####
data #Reference fully processed Seurat object
q.data #Fully processed query dataset Seurat object

#Find anchors
anchors <- FindTransferAnchors(reference = data,
                             query = q.data,
                             dims = 1:50,
                             normalization.method = 'SCT',
                             npcs = 50, k.filter = 5,
                             max.features = 100, k.anchor = 5)

#Transfer cell scores
p <- TransferData(anchorset = anchors,
                  refdata = data$stage,
                  dims = 1:30, k.weight = 25,
                  sd.weight = 1)

p <- p[,c(1,ncol(p))]
colnames(p) <- c('stage.o','stage.score')

p$stage <- ifelse(p[,2] < 0.8, 'Not classified', p[,1]) #Set threshold for labeling

q.data <- AddMetaData(q.data, metadata = p) #Add labeling results

#Imputate gene expression
q.data[['transfer']] <- TransferData(anchorset = anchors,
                                    refdata = GetAssayData(data[['integrated']]),
                                    dims = 1:30, k.weight = 25,
                                    sd.weight = 1)

DefaultAssay(q.data) <- 'transfer'
q.data@assays$transfer@var.features <- data@assays$integrated@var.features
q.data<- ScaleData(q.data, assay = 'transfer', do.scale = F)

# ============= 5. Examine label assignment ============= ####
df <- FetchData(object = q.data,
                vars = c('ident','stage.o','stage.score'))

ggplot(data = df, aes(x = factor(stage.o), fill = factor(stage.o),
                            y = stage.score))+
  geom_violin() +
  geom_jitter(size = 0.1, show.legend = F) +
  scale_fill_manual(values = colors[-11], drop = FALSE) +
  facet_wrap(~res_0.7, ncol = 5) +
  geom_hline(yintercept = 0.8,
             color = 'blue',
             linetype = 'dashed') +
  labs(x = 'Cell type', y = "Label probability") +
  theme_classic() +
  guides(fill = guide_legend(title = 'Cell type')) +
  theme(legend.position = 'bottom',
        axis.text.y = element_text(color = 'black',
                                   size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = NA,
                                        color = "black",
                                        linetype = 'dashed')) #plot distribution

#Examine dominant labels per cluster
q.data$confidence <- ifelse(q.data$stage.score > 0.8, 'HC','LC')
i <- table(Idents(q.data),q.data$stage.o)
i <- i/rowSums(i)*100
i <- melt(i)
j <- table(Idents(q.data),q.data$confidence)
j <- j/rowSums(j)*100
j <- melt(j)
i <- i %>%
  group_by(Var1) %>%
  top_n(1) %>%
  left_join(j[j$Var2 == 'HC',], by = 'Var1') %>%
  arrange(Var1) %>%
  select_('cluster' = 'Var1','cell.type' = 'Var2.x',
          'perc.dominant.label' = 'value.x',
          'perc.HC.label' = 'value.y')

# ============= 6. Differential gene expression analysis in query datasets ============= ####
q.data #Fully processed query dataset Seurat object with transferred labels

Idents(object=q.data) <- q.data$stage
DefaultAssay(object = q.data) <- 'SCT'
markers <- FindMarkers(object = q.data, only.pos = FALSE,
                        min.pct = 0.4,
                        logfc.threshold = 0.5,
                        ident.1 = 'Condition_1', #biological condition of interest within dataset
                        subset.ident = "Late.P1", #activation stage of interest within dataset
                        return.thresh = 0.01,
                        group.by = 'condition')

# ============= 7. Indexed data dimensionality reduction and clustering ============= ####
q.data #Fully processed query dataset Seurat object with transferred labels
x #Indexed flow cytometry data with barcoded cell annotation

fc.umap <- uwot::umap(x,
                      verbose = TRUE)
rownames(fc.umap) <- rownames(x)

#UMAP clusterering based on K means
inertia <- c()
for(i in c(1:20)){
  inertia[i] <- kmeans(x = fc.umap, centers = i, iter.max = 50)$tot.withinss
}

#Visualize inertia
ggplot(data.frame(inertia, k = seq(1:20)),
             aes(x=k, y=inertia, group =1))+
  geom_point(size = 2)+
  geom_line()+
  labs(x = "Number of clusters", y = "Total whithin cluster SS")+
  theme_classic()+
  theme(axis.text = element_text(size = 8))

#Add data to object
q.data <- AddMetaData(q.data, metadata = fc, col.name = c('FC.UMAP_1', 'FC.UMAP_2', 'Flow.Cluster'))

#score mean fluoresence
out <- x %>%
  mutate(cluster = fc$flow.cluster)

df <- melt(out)
colnames(df) <- c('cluster','marker', 'value')
res <- aov(value ~ marker*cluster, data = df)

out <- x %>%
  mutate(cluster = fc$flow.cluster) %>%
  group_by(cluster) %>%
  dplyr::summarise(
    MHCII = mean(MHCII),
    FSC = mean(FSC.A),
    SSC = mean(SSC.A),
    CD11b = mean(CD11b),
    CD45 = mean(CD45),
    CD301b = mean(CD301b),
    DAPI = mean(DAPI),
    F4_80 = mean(F4_80),
    Ly6c = mean(Ly6c)
  )

out[,2:10] <- scale(out[,2:10])
marker.order <- colnames(out[,2:10])[hclust(dist(t(out[,2:10])))$order]
out <- melt(out)
out$variable <- factor(out$variable, levels = marker.order)

ggplot(out, aes(x = cluster, y = variable, fill = value))+
  geom_tile()+
  labs(fill = "Marker (MFI)")+
  scale_fill_gradientn(colours = c("#191970","#FFFFFF","#FF8C00") )+
  theme(text = element_text(size = 8),
        legend.key.size = unit(4, 'mm'),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 8, color = 'black'),
        legend.position = 'bottom',
        panel.background = element_blank())

# ============= 8. Gene expression network ============= ####
require(igraph)
require(ggraph)

#load and prepare data
pval #object with p values for each gene in each path
edges #Load String network of genes in pval

#calculate aggregate p values
df1 <- data.frame(id = names(pval$Path_1), L1 = pval$Path_1)
df2 <- data.frame(id = names(pval$Path_2), L2 = pval$Path_2)
df3 <- data.frame(id = names(pval$Path_3), L3 = pval$Path_3)
df4 <- data.frame(id = names(pval$Path_4), L4 = pval$Path_4)

pval.df <- df1 %>%
  full_join(df2, by = 'id') %>%
  full_join(df3, by = 'id') %>%
  full_join(df4, by = 'id')

pval.df[pval.df == 0] <- 1e-320 #set lowest p val
pval.df[is.na(pval.df)] <- 1 #fill NAs with 1
rownames(pval.df) <- pval.df$ids #set gene names as rownames
pval.df[,2:5] <- -1*log(pval.df[,2:5])
pval.df[,2:5] <- scale(pval.df[,2:5],center=F) #scale per path
pval.df$agg.pval <- rowSums(pval.df[,2:5]) #aggregate

#filter string network
edges <- edges[,c(1:2,13)]
colnames(edges) <- c('from','to','score')

#Build node data.frame with genes with p value < 10^-9
nodes <- c()
for(i in pval){
  genes <- names(sort(i[i < 1e-9], decreasing = FALSE)) # Identify genes with the most significant time-dependent model fit.
  genes <- genes[-grep('Rpl',genes)] #remove ribosomal genes
  genes <- genes[-grep('Rps',genes)] #remove ribosomal genes
  nodes <- unique(c(nodes,genes))
}

nodes <- pval.df[pval.df$id %in% nodes,]

#Calculate edge weight
edges <- edges %>%
  left_join(nodes[,-c(2:5)], by = c("from" = "id")) %>%
  left_join(nodes[,-c(2:5)], by = c("to" = "id"), suffix = c("_from","_to")) %>%
  mutate(weight = (agg.pval_from+agg.pval_to)*score)
edges <- edges[,-c(3:5)]

#Filter out non-overlapping nodes
nodes <- nodes[nodes$id %in% unique(c(edges$to, edges$from)),]
ex <- setdiff(unique(c(edges$to, edges$from)),nodes$id)
edges <- edges[!edges$to %in% ex, ]
edges <- edges[!edges$from %in% ex, ]
nodes <- nodes[nodes$id %in% unique(c(edges$to, edges$from)),]

#Build network
g <- graph_from_data_frame(d = edges,
                           vertices = nodes,
                           directed = FALSE) # Create an igraph object with attributes directly from dataframes

#filter disconnected nodes
gd <- induced_subgraph(g,
                       V(g)[!components(g)$membership %in% which.max(components(g)$csize)])
g <- induced_subgraph(g,
                      V(g)[components(g)$membership == which.max(components(g)$csize)])

l <- layout_with_fr(graph = g, dim = 2,
                    niter = vcount(g)*50, grid = 'nogrid',
                    weights = E(g)$weight) #Calculate layout

#Add node properties
V(g)$degree <- degree(g)
V(g)$bw <- betweenness(g, directed = FALSE, weights = 1/E(g)$weight)
V(g)$eg <- eigen_centrality(g, directed = FALSE)$vector
V(g)$strength <- strength(g)
V(g)$cls <- closeness(g)
E(g)$ebw <- edge_betweenness(g, directed = FALSE, weights = E(g)$weight)

keep <- as_ids(V(g)[components(g)$membership == which.max(components(g)$csize)])
df <- nodes[nodes$id %in% keep, ] %>%
  mutate(
    degree = V(g)$degree,
    betweenness = V(g)$bw,
    eigen_centrality = V(g)$eg,
    strength = V(g)$strength,
    closeness = V(g)$cls
  ) %>%
  dplyr::arrange(desc(betweenness))#annotate df with properties

#Cluster network
lc <- cluster_louvain(g)
V(g)$cluster <- membership(lc)

# ============= 9. Join all data with no supervision ============= ####
require(Seurat)
datasets #list of Datasets from all tissues, fully processed and with label assignments but no other meta.data

d <- list()
n = 1
max.f = 100000
for(i in datasets2){
  DefaultAssay(i) <- 'RNA'
  i <- DietSeurat(i, assays = 'RNA')
  if('stage.score' %in% colnames(i@meta.data)){
    t <- subset(i, subset = stage.score > 0.8, slot = 'counts') #retain high score cells from query datasets
  }else{
    t <- subset(i, cells = sample(colnames(i),500)) #retain 500 cells from reference data
  }
  t <- SCTransform(object = t, assay = "RNA") #SCTransform for normalization
  if(nrow(t) < max.f){
    max.f <- nrow(t)
  }
  d[[n]] <- t
  n = n+1
}

#Prep objects for integration and select features
features <- SelectIntegrationFeatures(object.list = d, nfeatures = max.f)
max.size = 350000*max.f #To correct the globals
options(future.globals.maxSize= max.size)
d <- PrepSCTIntegration(object.list = d, anchor.features = features)

#Run PCA in each object with integration features
d <- lapply(X = d, FUN = function(x) {
  x <- RunPCA(x, verbose = FALSE, npcs = 50, features = features)
})

#Find anchors
anchors <- FindIntegrationAnchors(object.list = d,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = "rpca",
                                  dims = 1:50)

#integrate dataset
d <- IntegrateData(anchorset = anchors,
                   new.assay.name = 'integrated',
                   normalization.method = "SCT",
                   dims = 1:50)

#Batch correct with Hamony
d <- d %>%
  RunHarmony(group.by.vars = 'tissue',
             assay.use = 'SCT', reduction = "pca") %>%
  RunUMAP(reduction = "harmony", dims = 1:20,
          reduction.name = 'harmony.umap',
          assay = 'SCT', reduction.key = 'hUMAP_')
