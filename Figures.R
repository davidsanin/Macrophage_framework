# ==============================================================================
# =                                                                            =
# =                   Predictive Framework of Macrophage                       =
# =                                Activation                                  =
# =                                                                            =
# =                            David E Sanin 2021                              =
# =                                  Figures                                   =
# ==============================================================================

# Code to generate figures

#load libraries
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(clustree)
  library(ggplot2)
  library(plyr)
  library(scales)
  library(ggrepel)
  library(biomaRt)
  library(tidyverse)
  library(patchwork)
  library(future)
  library(reticulate)
  library(SeuratWrappers)
  library(destiny)
  library(gam)
  library(slingshot)
  library(TSP)
  library(pheatmap)
  library(biomaRt)
  library(reshape2)
  library(goseq)
  library(RColorBrewer)
  library(igraph)
  library(ggraph)
  library(ggtree)
  library(ggcorrplot)
  library(rgl)
  library(gasper)
  library(gganimate)
  library(directlabels)
  library(RcisTarget)
  library(wordcloud2)
  library(KEGGREST)
  library(GO.db)
})

#Functions
color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# ===================== Data ===================== ####
pval #genes annotated with p values
tg # differentially regulated genes - pseudotime analysis with lower than 10^-9 p values
n1 # negative control dataset1 H.poly SVF
n2 # negative control dataset2 Lamina propria large intestine
n3 # negative control dataset2 L.mono SVF
markers # differentially regulated genes
lin # pseudotime analysis results
c0 #fully processed reference dataset
c1 #fully processed lamina propria dataset
c2 #fully processed sciatic nerve dataset
c3 #fully processed breast cancer dataset
c4 #fully processed artery dataset
c5 #fully processed lung dataset
c6 #fully processed liver dataset
c7 #fully processed heart dataset
c8 #fully processed retina dataset
c9 #fully processed muscle dataset
cs #fully processed skin dataset
cm #fully processed microglia dataset
c.all #fully processed integrated datasets
node_master #annotated network nodes
g_old #network
l_old #network layout
fc.umap #flow cytometry based UMAP
datasets <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9)

# ===================== Color scheme and factors ===================== ####
colors <- c(brewer.pal(10,'Paired'),'#d4d4d4')[c(4,2,1,7,8,5,6,10,9,3,11)]
wound <- c('#ecc19c', '#1e847f')
col.paths <- c("#FF7F00","#E31A1C","#CAB2D6","#6A3D9A")
col.pst1 <- colorRampPalette(c("#33A02C","#1F78B4","#A6CEE3",'#FDBF6F','#FF7F00'))(100)
col.pst2 <- colorRampPalette(c("#33A02C","#1F78B4","#A6CEE3",'#FB9A99','#E31A1C'))(100)
col.pst3 <- colorRampPalette(c("#33A02C","#1F78B4","#A6CEE3",'#CAB2D6'))(100)
col.pst4 <- colorRampPalette(c("#33A02C","#1F78B4","#A6CEE3",'#6A3D9A'))(100)
col.pst <- colorRampPalette(c('#191970'	,'#008B8B',	'#00FF00'))(100)

colors.vec <- list(
  "Path_1" = col.pst1,
  "Path_2" = col.pst2,
  "Path_3" = col.pst3,
  "Path_4" = col.pst4
)

order <- c('Initial','Early','Intermediate',
           'Late.P1', 'Final.P1','Late.P2',
           'Final.P2','Final.P4','Final.P3',
           'Cycling')

c0$stage <- factor(x = c0$stage,
                       levels = order)

cs[['stage']] <- factor(x = cs$stage,
                       levels = c(order,'Not classified'))

path.names <- c("Phagocytic", "Oxidative stress",
                "Inflammatory", "Remodelling")

titles <- c(
  "Path_1" = "Phagocytic path",
  "Path_2" = "Oxidative stress path",
  "Path_3" = "Inflammatory path",
  "Path_4" = "Remodelling path"
)

titles2 <- c(
  'Lamina propria\nHigh fat vs. control diet',
  'Sciatic Nerve\nNaive vs. injured',
  'Breast tumor\nWT vs. Dab2 KO',
  'Atherosclerotic plaque\nProgressing vs. Regressing\nlesion',
  'Lung\nNaive vs. Cryptococcus\nneoformans infection',
  'Liver\nHealthy vs. Fibrotic tissue',
  'Heart\nHealthy vs. Infarcted tissue',
  'Retina\nHealthy vs. Neurodegenerative\ncondition',
  'Skeletal muscle\nNaive vs. T. gondii\ninfected'
)

titles3 <- c(
  'Reference (infection)',
  'Lamina propria (diet)',
  'Sciatic Nerve (injury)',
  'Breast (cancer)',
  'Artery (lesion)',
  'Lung (infection)',
  'Liver (fibrosis)',
  'Heart (infarction)',
  'Retina (neurodegeneration)',
  'Skeletal muscle (infection)',
  'Skin (wound)'
)

n = 1
for(i in datasets){
  i[['stage']] <- factor(x = i$stage,
                         levels = c(order,'Not classified'))
  datasets[[n]] <- i
  n=n+1
}

colors.datasets <- list('colors.c1' = colors[-c(6:8,10)],
                        'colors.c2' = colors[-c(6:8)],
                        'colors.c3' = colors[-c(6:8)],
                        'colors.c4' = colors[-c(6,8)],
                        'colors.c5' = colors[-c(6:8)],
                        'colors.c6' = colors[-c(2,6:7,8)],
                        'colors.c7' = colors[-c(3,6,8,10)],
                        'colors.c8' = colors[-c(6,8,10)],
                        'colors.c9' = colors[-c(6:8)])

colors.cs <- colors[-c(3,6:8)]

colors.all <- list('colors.c0' = colors,
                   'colors.c1' = colors[-c(6:8,10)],
                   'colors.c2' = colors[-c(6:8)],
                   'colors.c3' = colors[-c(6:8)],
                   'colors.c4' = colors[-c(6,8)],
                   'colors.c5' = colors[-c(6:8)],
                   'colors.c6' = colors[-c(2,6:7,8)],
                   'colors.c7' = colors[-c(3,6,8,10)],
                   'colors.c8' = colors[-c(6,8,10)],
                   'colors.c9' = colors[-c(6:8)],
                   'colors.c10' = colors[-c(3,6:8)])

# ===================== Figure 1 and S1 ===================== ####

#Fig. 1B cartoon of experimental set up then UMAP####
f1b <- DimPlot(object = c0, reduction = "umap",
               pt.size = 0.1, label = T,
               label.size = 5) +
  theme(plot.title = element_text(face = 'plain', size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  NoLegend()

#Fig. S1A original datasets - H. poly####
fs1a <- c(DimPlot(object = n1, reduction = "umap",
                  pt.size = 0.01, label.size = 3,
                  label = T, combine = F),
          FeaturePlot(object = n1,
                      features = c("Cd68",'Adgre1','H2-Ab1',
                                   "Lyz2","Itgam",'Csf1r'),
                      reduction = "umap",
                      min.cutoff = 'q40', max.cutoff = 'q80',
                      pt.size = 0.01, combine =F))

fs1a <- lapply(fs1a,
             function(x){x + theme(axis.title=element_blank(),
                                   axis.text=element_blank(),
                                   axis.ticks=element_blank(),
                                   axis.line = element_blank(),
                                   plot.title = element_text(face = 'italic',
                                                             size = 8),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank())+
                 NoLegend()
             })

#Fig. S1B original datasets - L. mono####
fs1b <- c(DimPlot(object = n3, reduction = "umap",
                  pt.size = 0.01, label.size = 3,
                  label = T, combine = F),
          FeaturePlot(object = n3,
                      features = c("Cd68",'Adgre1','H2-Ab1',
                                   "Lyz2","Itgam",'Csf1r'),
                      reduction = "umap",
                      min.cutoff = 'q40', max.cutoff = 'q80',
                      pt.size = 0.01, combine =F))

fs1b <- lapply(fs1b,
               function(x){x + theme(axis.title=element_blank(),
                                     axis.text=element_blank(),
                                     axis.ticks=element_blank(),
                                     axis.line = element_blank(),
                                     plot.title = element_text(face = 'italic',
                                                               size = 8),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.background = element_blank())+
                   NoLegend()
               })

#Fig. S1C - heatmap####
top20 <- markers %>% group_by(cluster) %>% top_n(15, avg_logFC)

fs1c <- DoHeatmap(object = c0, features = top20$gene,
                  angle = 0, size = 2) +
  NoLegend() +
  scale_fill_gradientn(colors = c("darkblue", "white", "darkorange")) +
  theme(axis.text = element_text(face='italic', color = 'black', size = 6))

#Fig. 1C Scores and cluster annotation####
ident.scores <- list("Complement & Phagocytosis" = c('C1qc','F13a1','C1qa','C4b','Cfh','C5ar1',
                                                     'Snx2','Tgfbr2','Dab2','Folr2','Cltc','Wwp1',
                                                     'Cd209d','Mrc1','Cd209f','Cd209g','Cd36',
                                                     'Ctsb','Lgmn','Cltc','Cd63'),
                     'ECM & actin regulation' = c('Cd44','Sdc1','Fn1',
                                                  'Pfn1','Fn1','Actg1','Tmsb4x'),
                     "Antigen presentation" = c('H2-Ab1','H2-Aa','H2-Eb1'),
                     "Antigen processing & presentation" =  c('H2-Oa','H2-Ab1','H2-DMb2',
                                                              'H2-Aa','H2-Eb1','H2-Ob',
                                                              'Cd74','Klrd1','H2-DMb1'),
                     "Phagosome" = c('Ctss','Cyba','Msr1',
                                     'Fcgr1','Coro1a','Thbs1',
                                     'Ncf4','Fcgr3'),
                     "Oxidative stress" = c('Prdx5','Txn1','Gsr','Ptgs2',
                                            'Ccs','Prdx6','Gpx4','Sesn1',
                                            'Sod3','Ltc4s'),
                     'ECM organization' = c('Col1a1','Nid1','Dpt','B4galt1',
                                            'Lum','Col3a1','Ccdc80',
                                            'Ramp2','Serpinh1','Ddr2'),
                     'Cycling' = c('Racgap1','Cks1b','Stmn1','Ran',
                                   'Cep57','Smc4','Top2a','Cks2',
                                   'Ube2s','Ube2c','Cenpw','Smc2',
                                   'Anp32b','Ranbp1','Cenpa')
                     )

c0 <- AddModuleScore(object = c0, features = ident.scores ,
                     assay = 'integrated',
                     name = 'Enrichment')

f1c <- FeaturePlot(object = c0,
                   features = paste0(rep('Enrichment',8),c(1:8)),
                   reduction = "umap",
                   min.cutoff = 'q40', max.cutoff = 'q80',
                   pt.size = 0.01, combine =F)

f1c <- lapply(f1c,
              function(x){x + theme(axis.title=element_blank(),
                                    axis.text=element_blank(),
                                    axis.ticks=element_blank(),
                                    axis.line = element_blank(),
                                    plot.title = element_text(face = 'plain', size = 10),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank())+
                  NoLegend()
              })

for(i in 1:8){
  f1c[[i]] <- f1c[[i]] +
    labs(title = names(ident.scores)[[i]])
}

#Fig. 1D Finding trajectories####
cl <- Idents(c0)
rd <- Embeddings(object = c0, reduction = "umap")[,1:2]

rd <- data.frame(rd)
idx <- match(rownames(rd),names(cl))
rd$cluster <- as.vector(cl)[idx]
tib <- rd %>%
  group_by(cluster) %>%
  dplyr::summarise(
    x = mean(UMAP_1),
    y = mean(UMAP_2)
  )
p1 <- tib[c(5,10,3,2,4,1),]
p2 <- tib[c(5,10,3,2,6,8),]
p3 <- tib[c(5,10,3,2,7),]
p4 <- tib[c(5,10,3,2,9),]

f1d <- DimPlot(object = c0, reduction = "umap",
               pt.size = 0.1, label = F) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  NoLegend() +
  geom_point(data = tib, aes(x = x, y = y),
             size = 2, color = 'black') +
  geom_path(data = p1, aes(x = x, y = y),
            color = 'black', size = 0.5) +
  geom_path(data = p2, aes(x = x, y = y),
            color = 'black', size = 0.5) +
  geom_path(data = p3, aes(x = x, y = y),
            color = 'black', size = 0.5) +
  geom_path(data = p4, aes(x = x, y = y),
            color = 'black', size = 0.5)

#Fig. 1E - plotting paths####
f1e <- FeaturePlot(object = c0, features = c('Path_1','Path_2',
                                             'Path_3','Path_4'),
                   reduction = "umap", order = T,
                   pt.size = 0.01, combine =F)

f1e <- lapply(f1e,
               function(x){x + theme(axis.text = element_blank(),
                                     axis.ticks = element_blank(),
                                     axis.line = element_blank(),
                                     axis.title = element_blank(),
                                     plot.title = element_text(face = 'plain',
                                                               size = 10),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.background = element_blank())+
                   NoLegend()
               })

set.seed(1)
df <- data.frame(lin@curves$curve1$s)
df_tsp <- ETSP(df)
o <- solve_TSP(df_tsp, control = list(start = which.min(df$UMAP_1)))
df <- df[o,]
df <- rbind(df[(which.max(df$UMAP_1)+10):nrow(df),], df[1:(which.max(df$UMAP_1)-10),])
f1e[[1]] <- f1e[[1]] +
  geom_path(data = df[1:(nrow(df)-10),],
            aes(x = UMAP_1, y = UMAP_2),
            arrow = arrow(ends = 'last', type = 'closed',
                          length = unit(0.1, 'inches')),
            size = 0.5) +
  scale_color_gradientn(na.value = 'lightgrey',
                        colours = col.pst1)

df <- data.frame(lin@curves$curve2$s)
df <- df[order(df$UMAP_1, decreasing = T),]
f1e[[2]] <- f1e[[2]] +
  geom_path(data = df[1:(nrow(df)-50),],
            aes(x = UMAP_1, y = UMAP_2),
            arrow = arrow(ends = 'last', type = 'closed',
                          length = unit(0.1, 'inches')),
            size = 0.5) +
  scale_color_gradientn(na.value = 'lightgrey',
                        colours = col.pst2)

df <- data.frame(lin@curves$curve3$s)
df <- df[order(df$UMAP_1, decreasing = T),]
f1e[[3]] <- f1e[[3]] +
  geom_path(data = df[1:(nrow(df)-80),],
            aes(x = UMAP_1, y = UMAP_2),
            arrow = arrow(ends = 'last', type = 'closed',
                          length = unit(0.1, 'inches')),
            size = 0.5) +
  scale_color_gradientn(na.value = 'lightgrey',
                        colours = col.pst3)

set.seed(1)
df <- data.frame(lin@curves$curve4$s)
df_tsp <- ETSP(df)
o <- solve_TSP(df_tsp, control = list(start = which.max(df$UMAP_1)))
df <- df[o,]
df <- rbind(df[which.max(df$UMAP_1):nrow(df),], df[1:(which.max(df$UMAP_1)+10),])
f1e[[4]] <- f1e[[4]] +
  geom_path(data = df[10:1950,],
            aes(x = UMAP_1, y = UMAP_2),
            arrow = arrow(ends = 'last', type = 'closed',
                          length = unit(0.1, 'inches')),
            size = 0.5) +
  scale_color_gradientn(na.value = 'lightgrey',
                        colours = col.pst4)

for(i in 1:4){
  f1e[[i]] <- f1e[[i]] +
    labs(title = paste0(i,'. ',path.names[[i]],' path'))
}

#Fig. 1F - Stages####
f1f <- DimPlot(object = c0, reduction = "umap",
               pt.size = 1, label = T, label.size = 4,
               group.by = 'stage', cols = colors) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())

#Fig. 1G - conditions distributions####
tb <- table(c0$stage, c0$condition, c0$tissue)
tb <- data.frame(tb)
colnames(tb) <- c('stage','condition','tissue','percentage')
tb$percentage <- tb$percentage/5
tb$condition <- factor(x = tb$condition,
                       levels = c("Naive", "H.poly", "L.mono"))
tb$tissue <- factor(x = tb$tissue,
                       levels = c("Adipose (m)", "Adipose (p)"))

tb <- tb[-c(11:20,31:40),]

f1g <- ggplot() +
  geom_bar(aes(y = percentage, x = condition, fill = stage),
           color ='black',
           data = tb[tb$tissue == 'Adipose (m)',],
           stat="identity")+
  labs(y="Label distribution (%)") +
  scale_fill_discrete(type = colors) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.text.y = element_text(colour="black", size = 8),
        axis.text.x = element_text(colour="black", size = 8, angle = 30,
                                   hjust = 1, vjust = 1),
        axis.title = element_text(colour="black", size = 10),
        axis.title.x = element_blank()) +
  ggplot() +
  geom_bar(aes(y = percentage, x = condition, fill = stage),
           color ='black',
           data = tb[tb$tissue == 'Adipose (p)',],
           stat="identity")+
  scale_fill_discrete(type = colors) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour="black", size = 8, angle = 30,
                                   hjust = 1, vjust = 1),
        axis.title = element_blank())

# ===================== Figure 2 and S2 ===================== ####

#Fig. 2A - Gene regulation####
f2a1 <- DimPlot(object = c0, reduction = "umap",
              pt.size = 1, label = T, label.size = 3,
              group.by = 'stage', cols = colors) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())

cells <- rownames(c0@meta.data[!is.na(c0@meta.data[,'Path_1']),])
cells <- sample(x = cells, size = 300)
sub <-  subset(c0, cells = cells)
sub$bin <- cut(sub[['Path_1']][,1],100)
f2a2 <- FeatureScatter(object = sub, combine = F,
                      feature1 = 'Path_1',
                      feature2 = 'Retnla',
                      cells = cells, group.by = 'bin',
                      cols = col.pst, pt.size = 0.2) +
  theme(axis.text=element_blank(),
        axis.title = element_text(size = 6),
        axis.ticks=element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  geom_smooth(method = gam::gam, colour = 'black',
              se = F, formula = y ~ lo(x), size = 0.5) +
  labs(x = "Pseudotime - t", y = 'Gene expression - G') +
  NoLegend() +
  coord_cartesian(clip = 'off')

#Fig. S2A - p value distribution####
fs2a <- list()
for(i in names(pval)){
  t <- data.frame(pval[[i]],
                  y = names(pval[[i]]))[order(pval[[i]]),]
  t$y <- factor(t$y, levels = t$y)
  colnames(t) <- c('pval','gene')
  t[t == 0] <- 1e-320

  fs2a[[i]] <- ggplot(t[-c(1:10),],
                      aes(x = pval, y = gene, group = 1)) +
    geom_line() +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 1e-9, color = 'blue',
                linetype = 'dashed') +
    theme(axis.text.y = element_blank(),
          plot.title = element_text(size = 8),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 6,
                                     color = 'black'),
          axis.title = element_text(size = 8)) +
    labs(title = paste(titles[[i]]),
         x = 'p value', y = 'Genes')
}

#Fig. 2B - Top regulated genes in each path####
f2b <- list()

for(i in names(tg)){
  cells <- rownames(c0@meta.data[!is.na(c0@meta.data[,i]),])
  sub <-  subset(c0, cells = cells)
  sub$bin <- cut(sub[[i]][,1],100)
  f2b[[i]] <- FeatureScatter(object = sub, combine = F,
                             feature1 = i, feature2 = tg[[i]][1],
                             cells = cells, group.by = 'bin',
                             cols = colors.vec[[i]], pt.size = 0.2) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_text(face = 'italic',
                                      size = 10),
          axis.title.x = element_text(size = 8),
          plot.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 6))+
    geom_smooth(method = gam::gam, colour = 'black', se = F,
                formula = y ~ lo(x))+
    NoLegend() +
    labs(x = titles[[i]]) +
    annotate('text', size = 2,
             label = paste0('p.val: ',
                            format(pval[[i]][tg[[i]][1]],
                                   digits = 3)),
             x = 5, y = Inf) +
    coord_cartesian(clip = 'off')
}

#Fig. 2C - Well-known gene expression####
f2c <- list()
genes <- c('Il4ra','Clec10a','Mrc1',
           'Il1b','H2-Ab1','Il6',
           'Hif1a','mt-Co1')
for(i in genes){
  t <- FetchData(object = c0,
                 vars = c(names(tg),i),
                 slot = 'scale.data')
  t <- t[!is.na(t[,1]),]
  t <- t[order(t[,1],decreasing = F),] #only counts for regulated genes
  t <- melt(t, id.vars = i)
  colnames(t) <- c('Y','Path','X')

  f2c[[i]] <- ggplot(data = t, aes(x = X, y = Y,
                                 color = Path)) +
    scale_color_manual(values = col.paths) +
    theme(axis.text = element_blank(),
          strip.background = element_blank(),
          strip.placement = 'inside',
          strip.text = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8,
                                      face = 'italic'),
          axis.line = element_line(),
          axis.ticks = element_blank(),
          legend.key = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    geom_smooth(method = gam::gam, se = F,
                formula = y ~ lo(x), size = 0.5) +
    labs(y = paste(i)) +
    facet_wrap(~Path, nrow = 1, ncol = 4, scales = 'free_x')
}

#Fig. S2B - Violins of genes####
fs2b <- VlnPlot(c0, features = genes, cols = colors,
                combine = F, pt.size = 0.01,
                group.by = 'stage')

fs2b[[1]] <- fs2b[[1]] +
  theme(axis.title =element_blank(),
        plot.title = element_text(size = 10,
                                  face = 'italic'),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_blank()) +
  coord_flip() +
  NoLegend()

fs2b[2:9] <- lapply(fs2b[2:9], function(x){
  x +
    theme(axis.title =element_blank(),
          plot.title = element_text(size = 10,
                                    face = 'italic'),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank()) +
    coord_flip() +
    NoLegend()})

#Fig. 2D - Apoptosis####
t <- FetchData(object = c0,
               vars = c(names(tg),'Apoptosis1'),
               slot = 'scale.data')

f2d <- list()

f2d[[1]] <- ggplot()+
  theme(axis.text = element_blank(),
        axis.title = element_text(size = 8),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_smooth(method = gam::gam, se = F,
              formula = y ~ lo(x), size = 0.5,
              data = t[!is.na(t[,1]),] %>% arrange(Path_1),
              aes(x = Path_1, y = Apoptosis1), color = "#FF7F00") +
  geom_smooth(method = gam::gam, se = F,
              formula = y ~ lo(x),  size = 0.5,
              data = t[!is.na(t[,2]),] %>% arrange(Path_2),
              aes(x = Path_2, y = Apoptosis1), color = "#E31A1C") +
  geom_smooth(method = gam::gam, se = F,
              formula = y ~ lo(x),  size = 0.5,
              data = t[!is.na(t[,3]),] %>% arrange(Path_3),
              aes(x = Path_3, y = Apoptosis1), color = "#CAB2D6") +
  geom_smooth(method = gam::gam, se = F,
              formula = y ~ lo(x),  size = 0.5,
              data = t[!is.na(t[,4]),] %>% arrange(Path_4),
              aes(x = Path_4, y = Apoptosis1), color = "#6A3D9A") +
  labs(x = 'Paths',
       y = 'Apoptosis score')

genes <- c('Acin1','Acvr1c','Aifm1','Aifm3','Akt1',
           'Ano6','Apaf1','Bax','Bbc3','Bcl2l1',
           'Bcl10','Blcap','Bok','Casp3','Casp7',
           'Casp8','Casp14','Casp16','Cd24a','Cdk5rap3',
           'Cdkn2a','Cecr2','Cidea','Cideb','Cidec','Cxcr3',
           'Dedd2','Dffa','Dffb','Dicer1','Dlc1','Dnase1l3',
           'Dnase2a','Dnase2b','Endog','Ern2','Exog','Fap',
           'Foxl2','Fzd3','Gcg','Gm20594','Gper1','Hsf1','Igfbp3','Il6',
           'Madd','Nmnat1','Pak2','Pam16','Plscr1',
           'Plscr2','Ptgis','Rffl','Rnf34','Rps3','Sharpin','Sirt2',
           'Stk24','Taok1','Tnf','Top2a','Trp53','Trp53bp2','Trpc5',
           'Xkr4','Xkr5','Xkr6','Xkr7','Xkr8','Xkr9','Zc3h12a')

genes <- genes[genes %in% rownames(c0)]
genes <- genes[genes %in% names(pval$Path_2)]
cells <- rownames(c0@meta.data[!is.na(c0$Path_2),])
df <- FetchData(object = c0, cells = cells, vars = genes)
o <- hclust(dist(t(df)))
genes <- genes[o$order][c(2:6,11:12)]
f2d[[2]] <- DoHeatmap(object = c0, features = genes, raster = F,
                      cells = cells, group.by = "L2_bin", group.bar = T,
                      group.colors = col.pst2, label = F, draw.lines = F) +
  NoLegend() +
  scale_fill_gradientn(colors = c("darkblue", "white", "darkorange")) +
  theme(axis.text = element_text(face='italic', color = 'black', size = 6))

#Fig. 2E - Retnla####
genes <- c('Retnla','Ear2')

f2e <- VlnPlot(c0, features = genes, cols = colors,
                combine = F, pt.size = 0.01,
                group.by = 'stage')

f2e[[1]] <- f2e[[1]] +
  theme(axis.title =element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  NoLegend()

f2e[[2]] <- f2e[[2]] +
    theme(axis.title =element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8)) +
  NoLegend()

t <- FetchData(object = c0,
               vars = c('Path_1',genes),
               slot = 'scale.data')
t <- t[!is.na(t[,1]),]
t <- t[order(t[,1],decreasing = F),] #only counts for regulated genes
t <- melt(t, id.vars = 'Path_1')
colnames(t) <- c('X','Gene','Y')

f2e[[3]] <- ggplot(data = t, aes(x = X, y = Y,
                                 color = Gene,
                                 fill = Gene))+
  theme(axis.text = element_blank(),
        axis.title = element_text(size =8),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  geom_smooth(method = gam::gam, se = F,
              formula = y ~ lo(x), size = 0.5) +
  labs(x = titles['Path_1'],
       y = 'Expression Level') +
  annotate('text', color = color_hue(2)[1], size = 2,
           label = paste0('Retnla p.val: ',
                          format(pval$Path_1[['Retnla']],digits = 3)),
           x = 10, y = 2) +
  annotate('text', color = color_hue(2)[2], size = 2,
           label = paste0('Ear2 p.val: ',
                          format(pval$Path_1[['Ear2']],digits = 3)),
           x = 4, y = Inf) +
  coord_cartesian(clip = 'off')

#Fig. 2G - Quantification Retnla####
df <-data.frame('2' = c(4.58,11.07,16.5,1.52,
                        4.89,1.73,6.74,11.8,
                        4.37,7.9),
                '4' = c(12.4,10.1,38.2,42,
                        50.8,35.1,37.4,
                        66.4,24.2,50.9),
                '8' = c(53.7,33.3,79.1,72.3,
                        66.7,37,63.5,
                        76.7,77.2,66.1))

df <- melt(df)
colnames(df) <- c('day','percentage')
df$day <- factor(as.numeric(gsub('X','',df$day)), levels = c(2,4,8))
res <- aov(percentage ~ day, data = df)
summary(res)
TukeyHSD(res)

f2g <- ggplot(data = df, aes(x = day,
                             y = percentage,
                             fill = day)) +
  geom_jitter(size = 2, height = 0, width = 0.1,
              color = 'black', pch = 21) +
  ylab(expression(paste(Retnla^'+',' RELM',alpha^'+',' (%)'))) +
  labs(x = 'Days after adoptive transfer') +
  scale_fill_manual(values = c('#ef9d10', '#3b4d61', '#6b7b8c')) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.text.y = element_text(colour="black", size = 8),
        axis.text.x = element_text(colour="black", size = 8,
                                   angle = 0,
                                   hjust = 1, vjust = 1),
        axis.title = element_text(colour="black", size = 8)) +
  annotate('text',
           label = '****',
           x = 2, y = Inf) +
  annotate('text',
           label = '****',
           x = 3, y = Inf) +
  coord_cartesian(clip = 'off')

#Fig. 2h####
df1 <- data.frame('Host' = c(4.740,
                             4.980,
                             9.990,
                             4.680),
                  'Transferred' = c(14.30,
                                    2.13,
                                    4.81,
                                    3.77))
df2 <- data.frame('Host' = c(10.530,
                             20.520,
                             17.780,
                             19.330,
                             4.960,
                             3.970,
                             2.500,
                             5.180),
                  'Transferred_WT' = c(33.8500,
                                       57.7200,
                                       15.0400,
                                       38.4600,
                                       NA,NA,NA,NA),
                  'Transferred_KO' = c(40.500000,
                                       31.760000,
                                       24.230000,
                                       NA,NA,NA,NA,NA))
df1 <- melt(df1)
df2 <- melt(df2)
colnames(df1) <- c('condition','percentage')
colnames(df2) <- c('condition','percentage')
res1 <- aov(percentage ~ condition, data = df1)
summary(res1)
TukeyHSD(res1)
res2 <- aov(percentage ~ condition, data = df2)
summary(res2)
TukeyHSD(res2)

f2h1 <- ggplot(data = df1, aes(x = condition,
                             y = percentage,
                             fill = condition)) +
  geom_jitter(size = 2, height = 0, width = 0.1,
              pch =21, color = 'black') +
  xlab('Peritoneal macrophages\n(Day 4 post-adoptive transfer)') +
  ylab(expression(paste(TIM4^'+',' RELM',alpha^'+',' (%)'))) +
  scale_fill_manual(values = c('#e2d810','#322e2f')) +
  scale_y_continuous(limits = c(0,80)) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.text.y = element_text(colour="black", size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(colour="black", size = 8)) +
  annotate('text',
           label = 'ns',
           x = 2, y = Inf) +
  coord_cartesian(clip = 'off')

f2h2 <- ggplot(data = df2, aes(x = condition,
                              y = percentage,
                              fill = condition)) +
  geom_jitter(size = 2, height = 0, width = 0.1,
              pch =21, color = 'black') +
  xlab('Monocytes\n(Day 4 post-adoptive transfer)') +
  ylab(expression(paste(TIM4^'-',' RELM',alpha^'+',' (%)'))) +
  scale_fill_manual(values = c('#e2d810','#d9138a','#12a4d9')) +
  scale_y_continuous(limits = c(0,80)) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.text.y = element_text(colour="black", size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(colour="black", size = 8))  +
  annotate('text',
           label = '**',
           x = 2, y = 65) +
  annotate('text',
           label = '*',
           x = 3, y = 65) +
  coord_cartesian(clip = 'off') +
  geom_linerange(xmin=2,xmax=3,y =70)+
  annotate('text',
           label = 'ns',
           x = 2.5, y = 75)

# ===================== Figure 3 and S3-4 ===================== ####

#Fig. S3A-H####
n1 <- AddModuleScore(object = n1, name = 'Mac',
                     features = list('Mac' = c('H2-Ab1','Lyz2',
                                               'Csf1r','Adgre1',
                                               'Mertk','Cd164',
                                               'Cd68','Itgam')))
n2 <- AddModuleScore(object = n2, name = 'Mac',
                     features = list('Mac' = c('H2-Ab1','Lyz2',
                                               'Csf1r','Adgre1',
                                               'Mertk','Cd164',
                                               'Cd68','Itgam')))
#Transfer data from reference dataset
an_n1 <- FindTransferAnchors(reference = c0,
                             query = n1,
                             dims = 1:50,
                             normalization.method = 'SCT',
                             npcs = 50, k.filter = 5,
                             max.features = 100, k.anchor = 5)

an_n2 <- FindTransferAnchors(reference = c0,
                             query = n2,
                             dims = 1:50,
                             normalization.method = 'SCT',
                             npcs = 50, k.filter = 5,
                             max.features = 100, k.anchor = 5)

#Transfer cell scores
p <- TransferData(anchorset = an_n1,
                  refdata = c0$stage,
                  dims = 1:30,  k.weight = 25,
                  sd.weight = 1)
p <- p[,c(1,ncol(p))]
colnames(p) <- c('stage','stage.score')
p_n1 <- p

p <- TransferData(anchorset = an_n2,
                  refdata = c0$stage,
                  dims = 1:30,  k.weight = 25,
                  sd.weight = 1)
p <- p[,c(1,ncol(p))]
colnames(p) <- c('stage','stage.score')
p_n2 <- p

p_n1[,1] <- ifelse(p_n1[,2] < 0.8, 'Not classified', p_n1[,1])
p_n2[,1] <- ifelse(p_n2[,2] < 0.8, 'Not classified', p_n2[,1])

n1 <- AddMetaData(n1, metadata = p_n1)
n2 <- AddMetaData(n2, metadata = p_n2)

n1$stage <- factor(x = n1$stage,
                   levels = c(order,'Not classified'))
n2$stage <- factor(x = n2$stage,
                   levels = c(order,'Not classified'))

#Plots
d <- c(
  area(1,2,2,5),
  area(1,1),
  area(2,1),
  area(3,2),
  area(3,3),
  area(3,4),
  area(3,5)
)

#Hpoly
p1 <- DimPlot(object = n1, reduction = "umap",
              pt.size = 0.1, group.by = 'stage',
              cols = colors[-6]) + NoLegend() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10,
                                  face = 'italic'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  labs(title = 'H. polygyrus')

p2 <- ggplot(p_n1, aes(x = stage.score)) +
  geom_density() +
  theme_classic() +
  labs(x = 'Label probability', y = 'Density') +
  geom_vline(xintercept = 0.8,
             color = 'blue', linetype = 'dashed') +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8,
                                 colour = 'black'))

p3 <- FeatureScatter(object = n1,
                     feature1 = 'Mac1',
                     feature2 = 'stage.score',
                     group.by = 'stage',
                     pt.size = 0.1, cols = colors[-6],
                     combine = F) +
  NoLegend() +
  labs(x = 'Macrophage score',
       y = 'Label probability') +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))  +
  geom_hline(yintercept = 0.8,
             color = 'blue', linetype = 'dashed')

p4 <- FeaturePlot(object = n1,
                  features = c('H2-Ab1','Lyz2','Csf1r','Adgre1'),
                  min.cutoff = 'q40', max.cutoff = 'q80',
                  combine = F, pt.size = 0.1, order = T)
p4 <- lapply(p4,
             function(x){x + theme(plot.title = element_text(face = 'italic',
                                                             size = 10),
                                   axis.title=element_blank(),
                                   axis.text=element_blank(),
                                   axis.ticks=element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank(),
                                   axis.line = element_blank()) +
                 NoLegend()
             })

fs3abcd <- p1+p2+p3+p4+
  plot_layout(design = d, widths = c(1.5,0.5,0.5,0.5,0.5))

#LP
p1 <- DimPlot(object = n2, reduction = "umap",
              pt.size = 0.1, group.by = 'stage',
              cols = colors[-c(2,6:10)]) + NoLegend() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10, face = 'plain'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  labs(title = 'Lamina Propria')

p2 <- ggplot(p_n2, aes(x = stage.score)) +
  geom_density() +
  theme_classic() +
  labs(x = 'Label probability', y = 'Density') +
  geom_vline(xintercept = 0.8,
             color = 'blue', linetype = 'dashed') +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8,
                                 colour = 'black'))

p3 <- FeatureScatter(object = n2,
                     feature1 = 'Mac1',
                     feature2 = 'stage.score',
                     group.by = 'stage',
                     pt.size = 0.1, cols = colors[-c(2,6:10)],
                     combine = F) +
  NoLegend() +
  labs(x = 'Macrophage score',
       y = 'Label probability') +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))  +
  geom_hline(yintercept = 0.8,
             color = 'blue', linetype = 'dashed')

p4 <- FeaturePlot(object = n2,
                  features = c('H2-Ab1','Lyz2','Csf1r','Adgre1'),
                  min.cutoff = 'q40', max.cutoff = 'q80',
                  combine = F, pt.size = 0.1, order = T)
p4 <- lapply(p4,
             function(x){x + theme(plot.title = element_text(face = 'italic',
                                                             size = 10),
                                   axis.title=element_blank(),
                                   axis.text=element_blank(),
                                   axis.ticks=element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank(),
                                   axis.line = element_blank()) +
                 NoLegend()
             })

fs3efdh <- p1+p2+p3+p4+
  plot_layout(design = d, widths = c(1.5,0.5,0.5,0.5,0.5))

#Fig. 3 UMAPs####
f3_1 <- list()
n = 1
for(i in datasets){
  f3a[[n]] <- DimPlot(object = i, reduction = "umap",
                 pt.size = 1, label = T, label.size = 2,
                 group.by = 'stage', cols = colors.datasets[[n]], repel = T) +
    theme(legend.position = 'none',
          plot.title = element_text(size = 8, face = 'plain',
                                    lineheight = 0.9),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank()) +
    labs(title = titles2[[n]])
  n = n+1
}

#Fig. 3 Probability distributions####
f3_2 <- list()
n = 1
for(i in datasets){
  f3b[[n]] <- ggplot(i[['stage.score']],
                     aes(x = stage.score)) +
    geom_density() +
    theme_classic() +
    scale_x_continuous(breaks = c(0.4,0.8)) +
    labs(x = 'Label probability', y = 'Density') +
    geom_vline(xintercept = 0.8,
               color = 'blue', linetype = 'dashed') +
    theme(plot.title = element_blank(),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6,
                                   colour = 'black'))
  n = n+1
}

#Fig. 3 Bars cell type conditon####
f3_3 <- list()
n = 1
for(i in datasets){
  tb <- table(i$stage, i$condition)
  tb <- t(t(tb)/colSums(tb)*100)
  tb <- data.frame(tb)
  colnames(tb) <- c('stage','condition','percentage')

  f3c[[n]] <- ggplot() +
    geom_bar(aes(y = percentage, x = condition, fill = stage),
             color ='black',
             data = tb,
             stat="identity")+
    labs(y="Label distribution (%)") +
    scale_fill_discrete(type = colors) +
    theme(axis.line = element_line(size=0.5, colour = "black"),
          panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_blank(),
          legend.position = 'none',
          text=element_text(family="Arial"),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour="black", size = 6),
          axis.text.x = element_blank(),
          axis.title = element_text(colour="black", size = 6),
          axis.title.x = element_blank())
  n = n+1
}

#Fig. 3 RELMa and Ear2####
f3_4 <- list()
n = 1
for(i in datasets){
  DefaultAssay(i) <- 'transfer'
  f3d[[n]] <- FeaturePlot(object = i, features = 'Retnla',
                          reduction = "umap",
                          min.cutoff = 'q40', max.cutoff = 'q80',
                          pt.size = 0.01, combine = F)
  f3d[[n]] <- lapply(f3d[[n]], function(x){
    x +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks= element_blank(),
            axis.line = element_blank(),
            plot.title = element_text(face = 'italic',
                                      size = 6),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      NoLegend()})
  n=n+1
}

f3_5 <- list()
n = 1
for(i in datasets){
  DefaultAssay(i) <- 'transfer'
  f3e[[n]] <- FeaturePlot(object = i, features = 'Ear2',
                          reduction = "umap",
                          min.cutoff = 'q40', max.cutoff = 'q80',
                          pt.size = 0.01, combine = F)
  f3e[[n]] <- lapply(f3e[[n]], function(x){
    x +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks= element_blank(),
            axis.line = element_blank(),
            plot.title = element_text(face = 'italic',
                                      size = 6),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      NoLegend()})
  n=n+1
}

#Fig. S3I-Q Data imputation####
fs3i_q <- list()
n = 1
genes <- list(c('Mrc1','H2-Ab1'),c('Mrc1','H2-Ab1'),
              c('Mrc1','H2-Ab1'),c('Mrc1','H2-Ab1'),
              c('Mrc1','H2-Ab1'),c('Mrc1','H2-Ab1'),
              c('Mrc1','Il1b'),c('Il1b','H2-Ab1'),
              c('Mrc1','H2-Ab1'))
for(i in datasets){
  DefaultAssay(i) <- 'transfer'
  p1 <- FeaturePlot(object = i, features = genes[[n]],
                    reduction = "umap",
                    min.cutoff = 'q40', max.cutoff = 'q80',
                    pt.size = 0.01, combine =F)
  DefaultAssay(i) <- 'integrated'
  p2 <- FeaturePlot(object = i, features = genes[[n]],
                    reduction = "umap",
                    min.cutoff = 'q40', max.cutoff = 'q80',
                    pt.size = 0.01, combine =F)
  fs3[[n]] <- c(p1,p2)
  fs3[[n]] <- lapply(fs3[[n]],
                     function(x){x + theme(axis.title=element_blank(),
                                           axis.text=element_blank(),
                                           axis.ticks=element_blank(),
                                           axis.line = element_blank(),
                                           plot.title = element_text(face = 'italic',
                                                                     size = 8),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.background = element_blank())+
                         NoLegend()
                     })
  n=n+1
}

#Figure S4####
fs4 <- list()
n = 1
for(i in datasets){
  p1 <- DimPlot(object = i, reduction = "umap",
                pt.size = 0.5, label.size = 3,
                label = T) +
    theme(legend.position = 'none',
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank()) +
    NoLegend()

  i[['stage.o']] <- factor(x = i$stage.o,
                           levels = order)
  df <- FetchData(object = i,
                  vars = c('ident','stage.o','stage.score'))

  p2 <- ggplot(data = df, aes(x = stage.o, fill = stage.o,
                              y = stage.score))+
    geom_violin(size = 0.1) +
    geom_jitter(size = 0.2, show.legend = F, stroke = 0, shape = 16) +
    scale_fill_manual(values = colors[-11], drop = FALSE) +
    facet_wrap(~ident, ncol = 6) +
    geom_hline(yintercept = 0.8,
               color = 'blue',
               linetype = 'dashed') +
    labs(y = "Label probability") +
    theme_classic() +
    guides(fill = guide_legend(title = 'Cell type', ncol = 2)) +
    theme(legend.key.size = unit(0.3,'cm'),
          legend.text = element_text(size = 6),
          axis.text.y = element_text(color = 'black',
                                     size = 8),
          legend.title = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          panel.background = element_rect(fill = NA,
                                          color = "black",
                                          linetype = 'dashed'))
  a <- ceiling(length(unique(df$ident))/6)
  d <- c(
    area(1,1,1,1),
    area(1,2,a,2),
    area(2,1,a,1)
  )

  fs4[[n]] <- p1+p2+guide_area()+
    plot_layout(design = d, widths = c(1,3),guides = 'collect')

  n = n+1
}

# ===================== Figure 4 ===================== ####

#Fig. 4A - Microglia####
f4a <- ggplot() +
  geom_density(data = cm[["E14"]]@meta.data,
               aes(x = stage.score),
               color = color_hue(4)[1]) +
  geom_density(data = cm[["P4_5"]]@meta.data,
               aes(x = stage.score),
               color = color_hue(4)[2]) +
  geom_density(data = cm[["P30"]]@meta.data,
               aes(x = stage.score),
               color = color_hue(4)[3]) +
  geom_density(data = cm[["P100"]]@meta.data,
               aes(x = stage.score),
               color = color_hue(4)[4]) +
  theme_classic() +
  theme(text = element_text(face = 'plain', size = 8),
        axis.text = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 8, color = 'black')) +
  geom_vline(xintercept = 0.8,
             color = 'blue',
             linetype = 'dashed') +
  labs(x = 'Label Probability')

#Fig. 4B-E - Atherosclerotic plaque####
f4b <- DimPlot(object = datasets[[4]], reduction = "umap",
               pt.size = 0.5, label = T, label.size = 3,
               group.by = 'stage', cols = colors.datasets[[4]], repel = T) +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10, face = 'plain',
                                  lineheight = 0.9),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  labs(title = 'Atherosclerotic plaque\nProgressing vs. Regressing lession')

tb <- table(datasets[[4]]$stage, datasets[[4]]$condition)
tb <- t(t(tb)/colSums(tb)*100)
tb <- data.frame(tb)
colnames(tb) <- c('stage','condition','percentage')
tb <- tb %>%
  filter(stage == 'Late.P1')

f4c <- ggplot() +
  geom_bar(aes(y = percentage, x = condition, fill = condition),
           color ='black',
           data = tb,
           stat="identity")+
  labs(y="Late.P1 cells (%)") +
  scale_fill_manual(values = c('#e52165','#0d1137')) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black", size = 6),
        axis.text.x = element_blank(),
        axis.title = element_text(colour="black", size = 6),
        axis.title.x = element_blank())

m # differential gene expression analysis data

top <- rownames_to_column(m, var = 'gene') %>%
  filter(p_val_adj < 0.01) %>%
  dplyr::arrange(desc(avg_logFC))

q <- top$gene[top$gene %in% node_master$id]

datasets[[4]]$seurat_clusters <- Idents(datasets[[4]])
Idents(object=datasets[[4]]) <- datasets[[4]]$stage
f4d <- VlnPlot(datasets[[4]], features = q,
               assay = 'SCT',
               idents = "Late.P1",
               cols = c('#e52165','#0d1137'), combine = F,
               group.by = 'condition', pt.size = 0)
Idents(object=datasets[[4]]) <- datasets[[4]]$seurat_clusters

n = 1
for(i in f4d){
  f4d[[n]] <- f4d[[n]] +
    coord_flip() +
    xlab(paste(q[[n]])) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8,
                                      face = 'italic'),
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size = 6),
          plot.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank())
  n = n+1
}

genes <- q
t <- FetchData(object = c0,
               vars = c('Path_1',genes),
               slot = 'scale.data')
t <- t[!is.na(t[,1]),]
t <- t[order(t[,1],decreasing = F),] #only counts for regulated genes
t <- melt(t, id.vars = 'Path_1')
colnames(t) <- c('X','Gene','Y')

f4e <- ggplot(data = t, aes(x = X, y = Y,
                            color = Gene,
                            fill = Gene))+
  theme(axis.text = element_blank(),
        axis.title = element_text(size =8),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(size = 8,
                                  face = 'italic', hjust = 0),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  guides(color = guide_legend(ncol = 3)) +
  geom_smooth(method = gam::gam, se = F,
              formula = y ~ lo(x), size = 0.5) +
  labs(x = "Phagocytic path",
       y = 'Expression Level') +
  geom_vline(xintercept = quantile(c0@meta.data[c0$stage == 'Late.P1',]$Path_1)[['25%']],
             color = 'black', size = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = quantile(c0@meta.data[c0$stage == 'Late.P1',]$Path_1)[['75%']],
             color = 'black', size = 0.5, linetype = 'dashed')  +
  facet_wrap(~Gene, nrow = 3, ncol = 3, scales = 'free')

#Fig. 4f-j - Dab2 - breast tumor####
f4f <- DimPlot(object = datasets[[3]], reduction = "umap",
               pt.size = 0.5, label = T, label.size = 3,
               group.by = 'stage', cols = colors.datasets[[3]], repel = T) +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10, face = 'plain',
                                  lineheight = 0.9),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  labs(title = 'Breast tumor\nWT vs. Dab2 KO')

f4g <- VlnPlot(datasets[[3]], features = 'Dab2', group.by = 'stage',
               cols = colors.datasets[[3]], combine = F,
               pt.size = 0)

f4g <- lapply(f4g, function(x){
  x +
    xlab('Dab2') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8,
                                      face = 'italic'),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank()) +
    coord_flip() +
    NoLegend()})

tb <- table(datasets[[3]]$stage, datasets[[3]]$condition)
tb <- t(t(tb)/colSums(tb)*100)
tb <- data.frame(tb)
colnames(tb) <- c('stage','condition','percentage')
tb <- tb %>%
  filter(stage == 'Late.P1')


f4h <- ggplot() +
  geom_bar(aes(y = percentage, x = condition, fill = condition),
           color ='black',
           data = tb,
           stat="identity")+
  labs(y="Late.P1 cells (%)") +
  scale_fill_manual(values = c("#E69F00","#56B4E9")) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black", size = 6),
        axis.text.x = element_blank(),
        axis.title = element_text(colour="black", size = 6),
        axis.title.x = element_blank())

m # differential gene expression analysis data

top <- rownames_to_column(m, var = 'gene') %>%
  filter(p_val_adj < 0.01) %>%
  dplyr::arrange(desc(avg_logFC))

q <- top$gene[top$gene %in% node_master$id]

datasets[[3]]$seurat_clusters <- Idents(datasets[[3]])
Idents(object=datasets[[3]]) <- datasets[[3]]$stage
f4i <- VlnPlot(datasets[[3]], features = q,
               assay = 'SCT',
               idents = "Late.P1",
               cols = c("#E69F00","#56B4E9"), combine = F,
               group.by = 'condition', pt.size = 0)
Idents(object=datasets[[3]]) <- datasets[[3]]$seurat_clusters

n = 1
for(i in f4i){
  f4i[[n]] <- f4i[[n]] +
    coord_flip() +
    xlab(paste(q[[n]])) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8,
                                      face = 'italic'),
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size = 6),
          plot.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank())
  n = n+1
}

genes <- q
t <- FetchData(object = c0,
               vars = c('Path_1',genes),
               slot = 'scale.data')
t <- t[!is.na(t[,1]),]
t <- t[order(t[,1],decreasing = F),] #only counts for regulated genes
t <- melt(t, id.vars = 'Path_1')
colnames(t) <- c('X','Gene','Y')

f4j <- ggplot(data = t, aes(x = X, y = Y,
                            color = Gene,
                            fill = Gene))+
  theme(axis.text = element_blank(),
        axis.title = element_text(size =8),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = 'inside',
        strip.text = element_text(size = 8,
                                  face = 'italic', hjust = 0),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  guides(color = guide_legend(ncol = 2)) +
  geom_smooth(method = gam::gam, se = F,
              formula = y ~ lo(x), size = 0.5) +
  labs(x = "Phagocytic path",
       y = 'Expression Level') +
  geom_vline(xintercept = quantile(c0@meta.data[c0$stage == 'Late.P1',]$Path_1)[['25%']],
             color = 'black', size = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = quantile(c0@meta.data[c0$stage == 'Late.P1',]$Path_1)[['75%']],
             color = 'black', size = 0.5, linetype = 'dashed')  +
  facet_wrap(~Gene, nrow = 3, ncol = 5, scales = 'free')

# ===================== Figure 5 and S5 ===================== ####

#Fig. 5A-D - UMAP and label transfer####
f5a <- DimPlot(object = cs, reduction = "umap",
               pt.size = 0.5, label = T, label.size = 3,
               group.by = 'stage', cols = colors.cs) +
  theme(legend.position = 'none',
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  annotate('text', label = 'Skin wound - d4 vs. d14',
           size = 4, color = 'black', x = 0, y = Inf) +
  coord_cartesian(clip = 'off')

f5c <- ggplot(cs[['stage.score']],
                    aes(x = stage.score)) +
  geom_density() +
  theme_classic() +
  labs(x = 'Label probability', y = 'Density') +
  geom_vline(xintercept = 0.8,
             color = 'blue', linetype = 'dashed') +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8,
                                 colour = 'black'))

tb <- table(cs$stage, cs$condition)
tb <- t(t(tb)/colSums(tb)*100)
tb <- data.frame(tb)
colnames(tb) <- c('stage','condition','percentage')

f5d <- ggplot() +
  geom_bar(aes(y = percentage, x = condition, fill = stage),
           color ='black',
           data = tb,
           stat="identity")+
  labs(y="Label distribution (%)") +
  scale_x_discrete(labels = c('4 dpw', '14 dpw')) +
  scale_fill_discrete(type = colors) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black", size = 8),
        axis.text.x = element_text(colour="black", size = 8),
        axis.title = element_text(colour="black", size = 8),
        axis.title.x = element_blank())

#Fig. S5A####
p1 <- DimPlot(object = cs, reduction = "umap",
              pt.size = 0.5, label.size = 5,
              label = T) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  NoLegend()

cs[['stage.o']] <- factor(x = cs$stage.o,
                         levels = order)
df <- FetchData(object = cs,
                vars = c('ident','stage.o','stage.score'))

p2 <- ggplot(data = df, aes(x = stage.o, fill = stage.o,
                            y = stage.score))+
  geom_violin(size = 0.1) +
  geom_jitter(size = 0.4, show.legend = F, stroke = 0, shape = 16) +
  scale_fill_manual(values = colors[-11], drop = FALSE) +
  facet_wrap(~ident, ncol = 5) +
  geom_hline(yintercept = 0.8,
             color = 'blue',
             linetype = 'dashed') +
  labs(y = "Label probability") +
  theme_classic() +
  guides(fill = guide_legend(title = 'Cell type', ncol = 2)) +
  theme(legend.key.size = unit(0.3,'cm'),
        legend.text = element_text(size = 8),
        axis.text.y = element_text(color = 'black',
                                   size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = NA,
                                        color = "black",
                                        linetype = 'dashed'))

a <- ceiling(length(unique(df$ident))/5)
d <- c(
  area(1,1,1,1),
  area(1,2,a,2),
  area(2,1,a,1)
)

fs5ab <- p1+p2+guide_area()+
  plot_layout(design = d, widths = c(1,3),guides = 'collect')

#Fig. 5E - Monocytes UMAP####
cs@meta.data[,c(8:22,26)] <- as.data.frame(apply(cs@meta.data[,c(8:22,26)], 2, as.numeric))
cs$TxMo <- ifelse(cs$tdRFP >= 5000, 'Yes', 'No')

f5e <- DimPlot(object = cs, reduction = "umap",
               pt.size = 0.5, label = F, na.value = 'lightgrey',
               group.by = 'TxMo', cols = c('lightgrey','red')) +
  theme(legend.position = 'none',
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  annotate('text', label = 'Transferred Monocytes',
           size = 4, color = 'red', x = 0, y = Inf) +
  coord_cartesian(clip = 'off')

#Fig. S5B-C - Monocyte clusters####
cells <- rownames(cs@meta.data[cs$exp.date %in% c('200904','200914'),])
css <- subset(cs, cells = cells)

css$tdRFP <- css$tdRFP/10000
p1 <- VlnPlot(object = css, features = 'tdRFP',
              pt.size = 0.5, group.by = 'condition',
              combine = F, cols = wound)
p2 <- VlnPlot(object = css, features = 'tdRFP',
              pt.size = 0.5, group.by = 'monocytes.day',
              combine = F, cols = c('#79cbb8', '#500472'))

fs5c <- c(p1,p2)

fs5d <- lapply(fs5b,
               function(x){x + theme(axis.title.x = element_blank(),
                                     plot.title = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.text.y = element_text(size = 8),
                                     axis.title.y = element_text(size = 8),
                                     axis.ticks.x = element_blank(),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.background = element_blank()) +
                   NoLegend() +
                   ylab(expression(paste('tdTomato (MFI x ',10^4,')'))) +
                   geom_hline(yintercept = 0.5, linetype= 'dashed', color = 'black')
               })

#Fig. 5F - Monocyte clusters####
cells <- rownames(cs@meta.data[cs$exp.date %in% c('200904','200914'),])
css <- subset(cs, cells = cells)
css$tdRFP <- css$tdRFP/10000
f5f <- VlnPlot(object = css, features = 'tdRFP',
               pt.size = 0.2, group.by = 'stage',
               combine = F, cols = colors.cs)

f5f <- lapply(f5f,
              function(x){x + theme(axis.title.x = element_blank(),
                                    plot.title = element_blank(),
                                    axis.text.x = element_text(angle = 30,
                                                               size = 8),
                                    axis.text.y = element_text(size = 8),
                                    axis.title.y = element_text(size = 8),
                                    axis.ticks.x = element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank()) +
                  NoLegend() +
                  ylab(expression(paste('tdTomato (MFI x ',10^4,')'))) +
                  geom_hline(yintercept = 0.5, linetype= 'dashed', color = 'black')
              })

#Fig. 5G - Bar####
cells <- rownames(cs@meta.data[cs$TxMo == 'Yes',])
css <- subset(cs, cells = cells)

tb <- table(css$stage, css$condition)
tb <- t(t(tb)/colSums(tb)*100)
tb <- data.frame(tb)
colnames(tb) <- c('stage','condition','percentage')

f5g <- ggplot() +
  geom_bar(aes(y = percentage, x = condition, fill = stage),
           color ='black',
           data = tb,
           stat="identity")+
  labs(y="Label distribution (%)") +
  scale_x_discrete(labels = c('4 dpw', '14 dpw')) +
  scale_fill_discrete(type = colors) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black", size = 8),
        axis.text.x = element_text(colour="black", size = 8),
        axis.title = element_text(colour="black", size = 8),
        axis.title.x = element_blank())

#Fig. 5H-I - Retnla and Ear2 in TxMono####
cells <- rownames(cs@meta.data[cs$TxMo == 'Yes',])
css <- subset(cs, cells = cells)
f5h1 <- FeaturePlot(object = css, features = c('Retnla', 'Ear2'),
                   reduction = "umap",
                   min.cutoff = 'q40', max.cutoff = 'q80',
                   pt.size = 0.01, combine =F)
f5h1 <- lapply(f5h1,
              function(x){x + theme(axis.title = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks= element_blank(),
                                    axis.line = element_blank(),
                                    plot.title = element_text(face = 'italic',
                                                              size = 8),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank())+
                  NoLegend()
              })


f5h2 <- VlnPlot(object = css, features = c('Retnla','Ear2'),
              pt.size = 0.2, group.by = 'stage',
              combine = F, cols = colors[-c(3,6:8,10)])

f5h2[[1]] <- f5h2[[1]] +
  theme(axis.title.y = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  NoLegend()

f5h2[[2]] <- f5h2[[2]] +
  theme(axis.title.y = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  NoLegend()

#Fig. 5I-M - Surface markers####
x <- cs@meta.data[,c(8,11,16:19,21:22)] #removed TREM2 as it is absent on 2nd data

#UMAP clusters
inertia <- c()
for(i in c(1:20)){
  inertia[i] <- kmeans(x = fc.umap, centers = i, iter.max = 50)$tot.withinss
}

fs5e <- ggplot(data.frame(inertia, k = seq(1:20)),
             aes(x=k, y=inertia, group =1))+
  geom_point(size = 2)+
  geom_line()+
  labs(x = "Number of clusters", y = "Total whithin cluster SS")+
  theme_classic()+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))

k <- 6
kc <- kmeans(x = fc.umap, centers = k, iter.max = 50)
fc <- data.frame(fc.umap)
fc$flow.cluster <- factor(kc$cluster)
fs5d <- fs5d +
  geom_vline(xintercept = k, linetype = 'dashed', color = 'red') +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

cs <- AddMetaData(cs, metadata = fc,
                  col.name = c('FC.UMAP_1', 'FC.UMAP_2', 'Flow.Cluster'))
f5ijk <- list()
f5ijk[[1]] <- FeatureScatter(cs, feature1 = 'FC.UMAP_1', feature2 = 'FC.UMAP_2',
                           group.by = 'stage', cols = colors.cs, pt.size = 0.5)
f5ijk[[2]] <- FeatureScatter(cs, feature1 = 'FC.UMAP_1', feature2 = 'FC.UMAP_2',
                           group.by = 'Flow.Cluster', pt.size = 0.1)
f5ijk[[3]] <- FeatureScatter(cs, feature1 = 'FC.UMAP_1', feature2 = 'FC.UMAP_2',
                           group.by = 'TxMo', pt.size = 0.1)

f5ijk <- lapply(f5ijk, function(x){
  x +
    theme(legend.position = 'bottom',
          legend.key.size = unit(6,'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          plot.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank()) +
    coord_cartesian(clip = 'off')
})

f5ijk[[3]] <- f5ijk[[3]]+
  scale_color_manual(na.value = 'lightgrey',
                     values = c('lightgrey','red'))


fs5f <- FeatureScatter(cs, feature1 = 'FC.UMAP_1', feature2 = 'FC.UMAP_2',
                       group.by = 'condition', cols = wound,
                       pt.size = 0.5)+
  theme(legend.position = 'bottom',
        legend.key.size = unit(6,'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  coord_cartesian(clip = 'off')

tb <- table(cs$stage, cs$Flow.Cluster)
tb <- t(t(tb)/colSums(tb)*100)
tb <- data.frame(tb)
colnames(tb) <- c('stage','flow.cluster','percentage')

f5l <- ggplot() +
  geom_bar(aes(y = percentage, x = flow.cluster, fill = stage),
           color ='black',
           data = tb,
           stat="identity")+
  labs(y="Label distribution (%)") +
  scale_fill_discrete(type = colors) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black", size = 8),
        axis.text.x = element_blank(),
        axis.title = element_text(colour="black", size = 8),
        axis.title.x = element_blank())

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
    F4_80 = mean(F4_80),
    Ly6c = mean(Ly6c)
  )

out[,2:9] <- scale(out[,2:9])
marker.order <- colnames(out[,2:9])[hclust(dist(t(out[,2:9])))$order]
out <- melt(out)
out$variable <- factor(out$variable, levels = marker.order)

f5m <- ggplot(out, aes(x = cluster, y = variable, fill = value))+
  geom_tile()+
  labs(fill = "Marker (MFI)")+
  scale_fill_gradientn(colours = c("#191970","#FFFFFF","#FF8C00") )+
  theme(text = element_text(size = 6),
        legend.key.size = unit(4, 'mm'),
        legend.text = element_text(size = 6, color = 'black'),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 6, color = 'black'),
        legend.position = 'bottom',
        panel.background = element_blank())

# ===================== Figure 6 ===================== ####

#Fig. 6A Stage distribution####
datasets2 <- c(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,cs)

order2 <- c('Initial','Early','Intermediate',
            'Late.P1', 'Final.P1','Late.P2',
            'Final.P2','Final.P4','Final.P3',
            'Cycling','Not classified')

n = 1
for(i in datasets2){
  i[['stage']] <- factor(x = i$stage,
                         levels = order2)
  datasets2[[n]] <- i
  n=n+1
}

m <- data.frame()
for(i in datasets2){
  tb <- table(i$stage, i$condition)
  tb <- t(t(tb)/colSums(tb)*100)
  tb <- data.frame(tb)
  colnames(tb) <- c('stage','condition','percentage')
  tb$tissue <- unique(i$tissue)[1]
  tb$condition <- paste(tb$condition,tb$tissue,sep='_')
  m <- rbind(m,tb[,-4])
}


id <- data.frame('from' = c("Naive_Adipose (m)","H.poly_Adipose (m)","L.mono_Adipose (m)",
                            "Naive_Lamina propria","HFD_Lamina propria",
                            "Naive_Sciatic nerve","Injury (d1)_Sciatic nerve","Injury (d5)_Sciatic nerve",
                            "Tumor associated macrophage (WT)_Tumor","Tumor associated macrophage (KO)_Tumor",
                            "Lession progression_Atherosclerotic plaque","Lession regression_Atherosclerotic plaque",
                            "naive_Lung","infected_Lung",
                            "ctrl_Liver","w2_Liver","w4_Liver",
                            "naive_Heart","d4_Heart",
                            "dark_Retina","light_Retina",
                            "naive_Skeletal muscle","inf_Skeletal muscle",
                            "Wound (d4)_Skin","Wound (d14)_Skin"),
                 'to' = c("Naive (fat)","H.poly","L.mono",
                          "Control diet","HFD",
                          "Naive (sn)","1 dpw","5 dpw",
                          "WT","Dab2 KO",
                          "Prog. lesion","Reg. lesion",
                          "Naive (lung)","C.neoformans",
                          "Healty (liver)","Week 2","Week 4",
                          "Healthy (heart)","Infarcted",
                          "Healthy (dark)","Degenerating (light)",
                          "Naive (sm)","T.gondii",
                          "4 dpw","14 dpw"))

idx <- match(m$condition,id$from)
m$condition <- id$to[idx]
m$condition <- factor(m$condition, levels = id$to)
grad <- colorRampPalette(c('#191970'	,'#FFFFFF', '#FF8C00'))(255)

f6a <- ggplot(data = m, aes(y = stage,
                            x = condition,
                            fill = percentage)) +
  geom_tile() +
  guides(fill = guide_colorbar(title = 'Percentage (%)',
                               barwidth = unit(25, 'mm'),
                               barheight = unit(2, 'mm'))) +
  scale_fill_gradientn(colours = grad, trans = 'sqrt',
                       na.value = "white") +
  scale_x_discrete(position = 'top') +
  theme(axis.text.x.top = element_text(angle = 90,
                                       size = 6,
                                       hjust = 0, vjust = 0,
                                       color = 'black'),
        axis.text.y = element_text(size = 6,
                                   color = 'black'),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6))

#Fig. 6C-D New Maps####
f6c <- DimPlot(object = c.all, reduction = "umap",
               group.by = 'res_0.5',
               pt.size = 0.5, label = T, label.size = 4) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())



f6d <- DimPlot(object = c.all, reduction = "umap",
               pt.size = 0.5, label = T, label.size = 4,
               group.by = 'stage', cols = colors) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())

#Fig. 6E - Surface markers####
m #differential gene expression results
csg #cell surface genes

surf <- m %>%
  filter(gene %in% csg, p_val_adj < 0.05)

tops <- surf %>%
  group_by(cluster) %>%
  top_n(4, wt = avg_logFC)

f6e <- DotPlot(object = c.all, dot.scale = 2,
               features = unique(tops$gene),
               group.by = 'stage') +
  guides(color = guide_colorbar(ticks = F, title = NULL,
                                direction = 'horizontal',
                                barwidth = unit(10, 'mm'),
                                barheight = unit(2, 'mm')),
         size = guide_legend(ncol = 2, title = NULL,
                             keywidth = unit(2, 'mm'),
                             keyheight = unit(2, 'mm'))) +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, face = 'italic',
                                   angle = 90,
                                   hjust = 1, vjust = 0.5),
        legend.text = element_text(size = 6),
        panel.grid = element_line(color = 'black', size = 0.05,
                                  linetype = 'dotted'))

#Fig. 6F - Cluster markers####
top5 <- m %>% group_by(cluster) %>% top_n(5, avg_logFC)

f6f <- DoHeatmap(object = c.all, features = top5$gene,
                 draw.lines = T,
                 angle = 90, size = 3, group.colors = colors) +
  NoLegend() +
  scale_fill_gradientn(colors = c("darkblue", "white", "darkorange")) +
  theme(axis.text = element_text(face='italic', color = 'black', size = 6))

# ===================== Figure 7 and S6 ===================== ####

#Fig. S6A - Subsetting network####
fs6a <- ggraph(g_old, layout = l_old) +
  geom_edge_link(aes(alpha = weight)) +
  geom_node_point(size = 0.5, color = 'black',
                  stroke = 0, shape = 16) +
  guides(edge_alpha = FALSE) +
  theme(legend.key = element_blank(),
        legend.position = 'bottom',
        panel.background = element_blank())

tb <- data.frame(E(g_old)$weight)
fs6b <- ggplot(tb,
               aes(x = E.g_old..weight)) +
  geom_density() +
  theme_classic() +
  labs(x = 'Edge weight', y = 'Density') +
  geom_vline(xintercept = 10,
             color = 'blue', linetype = 'dashed') +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8,
                                 colour = 'black'))

g <- subgraph.edges(g,E(g)[[weight>=10]])
g <- induced_subgraph(g,
                      V(g)[components(g)$membership == which.max(components(g)$csize)])
V(g)$degree <- degree(g)
V(g)$bw <- betweenness(g, directed = FALSE, weights = 1/E(g)$weight)
V(g)$eg <- eigen_centrality(g, directed = FALSE)$vector
V(g)$strength <- strength(g)
V(g)$cls <- closeness(g)
E(g)$ebw <- edge_betweenness(g, directed = FALSE, weights = E(g)$weight)
l <- layout_with_fr(graph = g, dim = 2,
                    niter = vcount(g)*50, grid = 'nogrid',
                    weights = E(g)$weight)
lc <- cluster_louvain(g)
V(g)$cluster <- membership(lc)

#Fig. S6C GO of clusters####
vc <- as_data_frame(g, what = 'vertices')
net <- list()
for(i in unique(vc$cluster)){
  net[[i]] <- vc[vc$cluster == i,'name']
}

genome <- 'mm10'
goResults <- list()
geneU <- rownames(c0)
ensembl = useMart( "ensembl",
                   dataset = "mmusculus_gene_ensembl" )
gs <- getBM(attributes = c("ensembl_gene_id",
                           "transcript_length",
                           "mgi_symbol"),
            filters = "mgi_symbol",
            values = geneU,
            mart = ensembl)

TxLen <- list()
for(gene in levels(as.factor(gs$ensembl_gene_id))){
  out <- gs[gs$ensembl_gene_id %in% gene,]$transcript_length
  TxLen[[gene]] <- median(out)
}

goResults <- list()
for(i in 1:length(net)){
  qs <- getBM(attributes = "ensembl_gene_id",
              filters = "mgi_symbol",
              values = net[[i]],
              mart = ensembl)$ensembl_gene_id
  qs <- qs[!is.na(qs)]

  genes <- as.integer(names(TxLen) %in% qs)
  names(genes) <- names(TxLen)

  pwf <- nullp(DEgenes = genes,
               genome = genome,
               id = "ensGene",
               bias.data = as.numeric(TxLen))

  goResults[[i]] <- goseq(pwf,
                          genome = genome,
                          id = "ensGene",
                          test.cats=c("GO:BP")) %>%
    mutate(hitsPerc=numDEInCat*100/numInCat) %>%
    filter(over_represented_pvalue < 0.01, numDEInCat > 2)  %>%
    dplyr::select(category,over_represented_pvalue,
                  hitsPerc,numDEInCat,term)
}

go2 <- list(c(1:3),
            c(1:3),
            c(1,3,6),
            c(2,30,42),
            c(1,8,36),
            c(6,11,17),
            c(1:2),
            c(1,3,10),
            c(1,4,26),
            c(1,2,37))

fs6c <- list()
leg <- list()
for(i in 1:length(goResults)){
  fs6c[[i]] <- goResults[[i]][go2[[i]],] %>%
    mutate(pval = -log(over_represented_pvalue+1e-320)) %>%
    dplyr::arrange(pval) %>%
    ggplot(aes(x=pval,
               y=term,
               size=numDEInCat,
               colour=hitsPerc)) +
    geom_point() +
    expand_limits(x=0) +
    coord_cartesian(clip = 'off') +
    scale_size_continuous(range = c(1,3)) +
    xlab(expression(paste(-log[10],'p value'))) +
    labs(colour="Hits (%)",
         size="Count",
         title = paste0('Cluster ',i)) +
    theme(panel.background = element_blank(),
          legend.position = 'bottom',
          plot.title = element_text(size = 8,
                                    color = 'black'),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 8,
                                      color = 'black'),
          axis.text.y = element_text(size = 6,
                                     color = 'black'),
          axis.text.x = element_text(size = 8,
                                     color = 'black'),
          axis.line = element_line(color = 'black',
                                   linetype = 'solid'),
          legend.key.size = unit(3, 'mm'),
          legend.key = element_blank())
  leg[[i]] <- cowplot::get_legend(fs6c[[i]])
  fs6c[[i]] <- fs6b[[i]]+theme(legend.position = 'none')
}

#Fig. 7A - Network####
go <- data.frame('cluster' = seq(1:10),
                 'Path' = c('Chemotaxis',
                            'Cytoskeleton\norganization',
                            'Lipid transport',
                            'Myeloid\ndifferentiation',
                            'Antigen processing\nand presentation',
                            'Response to\nexternal stimulus',
                            'protein complex\noligomerization',
                            'Leukotriene\nsynthesis',
                            'Superoxide\nmetabolic process',
                            'Lipid synthesis'))
cent <- data.frame(l)
cent$cluster <- V(g)$cluster
cent <- cent %>%
  group_by(cluster) %>%
  dplyr::summarise(x = mean(X1), y = mean(X2), ) %>%
  left_join(go, by = 'cluster')


f7a <- ggraph(g, layout = l) +
  geom_edge_link(aes(alpha = weight)) +
  geom_node_point(aes(color = factor(cluster),
                      size = strength)) +
  scale_size_continuous(range = c(0.2,2)) +
  guides(edge_alpha = guide_legend(title = "Weight"),
         size = guide_legend(title = "Strenght"),
         color = FALSE) +
  theme(legend.key = element_blank(),
        legend.position = 'bottom',
        legend.key.size = unit(4,'mm'),
        panel.background = element_blank(),
        text = element_text(size = 6)) +
  geom_text(data = cent, aes(x = x, y = y,
                             label = Path,
                             color = factor(cluster)),
            size = 3)

#Fig. 7B Network Per path####
g1 <- induced_subgraph(g, tg$Path_1[tg$Path_1 %in% V(g)$name])
V(g1)$degree <- degree(g1)
g1 <- induced_subgraph(g1, V(g1)[[degree>0]])
V(g1)$bw <- betweenness(g1, directed = FALSE, weights = 1/E(g1)$weight)
V(g1)$strength <- strength(g1)
l1 <- layout_with_fr(graph = g1, dim = 2,
                     niter = vcount(g1)*50, grid = 'nogrid',
                     weights = E(g1)$weight)

f7b1 <- ggraph(g1, layout = l1) +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = factor(cluster)),
                  size = 0.4) +
  scale_color_manual(values = color_hue(10)) +
  theme(legend.key = element_blank(),
        legend.position = 'none',
        panel.background = element_blank())

g2 <- induced_subgraph(g, tg$Path_2[tg$Path_2 %in% V(g)$name])
V(g2)$degree <- degree(g2)
g2 <- induced_subgraph(g2, V(g2)[[degree>0]])
V(g2)$bw <- betweenness(g2, directed = FALSE, weights = 1/E(g2)$weight)
V(g2)$strength <- strength(g2)
l2 <- layout_with_fr(graph = g2, dim = 2,
                     niter = vcount(g2)*50, grid = 'nogrid',
                     weights = E(g2)$weight)

f7b2 <- ggraph(g2, layout = l2) +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = factor(cluster)),
                  size = 0.4) +
  scale_color_manual(values = color_hue(10)[-8]) +
  theme(legend.key = element_blank(),
        legend.position = 'none',
        panel.background = element_blank())

g3 <- induced_subgraph(g, tg$Path_3[tg$Path_3 %in% V(g)$name])
V(g3)$degree <- degree(g3)
g3 <- induced_subgraph(g3, V(g3)[[degree>0]])
V(g3)$bw <- betweenness(g3, directed = FALSE, weights = 1/E(g3)$weight)
V(g3)$strength <- strength(g3)
l3 <- layout_with_fr(graph = g3, dim = 2,
                     niter = vcount(g3)*50, grid = 'nogrid',
                     weights = E(g3)$weight)

f7b3 <- ggraph(g3, layout = l3) +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = factor(cluster)),
                  size = 0.4) +
  scale_color_manual(values = color_hue(10)[-c(3,8)]) +
  theme(legend.key = element_blank(),
        legend.position = 'none',
        panel.background = element_blank())

g4 <- induced_subgraph(g, tg$Path_4[tg$Path_4 %in% V(g)$name])
V(g4)$degree <- degree(g4)
g4 <- induced_subgraph(g4, V(g4)[[degree>0]])
V(g4)$bw <- betweenness(g4, directed = FALSE, weights = 1/E(g4)$weight)
V(g4)$strength <- strength(g4)
l4 <- layout_with_fr(graph = g4, dim = 2,
                     niter = vcount(g4)*50, grid = 'nogrid',
                     weights = E(g4)$weight)

f7b4 <- ggraph(g4, layout = l4) +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = factor(cluster)),
                  size = 0.4) +
  scale_color_manual(values = color_hue(10)) +
  theme(legend.key = element_blank(),
        legend.position = 'none',
        panel.background = element_blank())

#Fig. 7C - Crossings####
cr <- crossing(graph = g, communities = lc)
gx <- delete_edges(g, names(cr[!cr])) # delete edges
bw_q75 <- quantile(V(g)$bw, 0.75)
gx <- induced_subgraph(gx, V(gx)[degree(gx)>0 & betweenness(gx) > as.numeric(bw_q75)])
V(gx)$degree <- degree(gx)

f7c <- ggraph(gx, layout = 'auto') +
  geom_edge_link(aes(alpha = weight)) +
  geom_node_point(aes(color = factor(cluster),
                      size = strength)) +
  scale_size_continuous(range = c(1,3)) +
  guides(edge_alpha = guide_legend(title = "Weight"),
         size = guide_legend(title = "Strenght"),
         color = FALSE) +
  theme(legend.key = element_blank(),
        legend.position = 'bottom',
        panel.background = element_blank(),
        legend.key.size = unit(4,'mm')) +
  geom_node_text(aes(label = name), color = 'red',
                 repel = TRUE, size = 1.5,
                 segment.size = 0.1)

#Fig. 7D TF enrichment####
vc <- as_data_frame(g, what = 'vertices')
net <- list()
for(i in unique(vc$cluster)){
  net[[i]] <- vc[vc$cluster == i,'name']
}

names(net) <- paste0('Cluster_',seq(1:10))
geneLists <- net
motifRankings # annotation mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
data(motifAnnotations_mgi)

motifEnrichmentTable_wGenes <- cisTarget(geneLists,
                                         motifRankings,
                                         motifAnnot=motifAnnotations_mgi)

mET <- motifEnrichmentTable_wGenes %>%
  filter(nEnrGenes > 0)

anotatedTfs <- lapply(split(mET$TF_highConf,
                            mET$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                      })

tf.s <- list()
for(i in names(tg)){
  tf.s[[i]] <- list()
  for(j in names(anotatedTfs)){
    tf.s[[i]][[j]] <- tg[[i]][tg[[i]] %in% anotatedTfs[[j]]]
  }
}

wc <- data.frame(unlist(tf.s)) %>%
  group_by(unlist.tf.s.) %>%
  dplyr::summarise(freq = n())

colnames(wc) <- c('word','freq')

dev.new(width = 2000, height = 2000, unit = "px")
set.seed(1) # for reproducibility
wordcloud2(data = wc, fontFamily = 'SansSerif',
           color = rep(brewer.pal(8, "Dark2"),nrow(wc)/8))

wc.p <- wc %>%
  arrange(dplyr::desc(freq)) %>%
  top_n(6, freq)

f7d <- list()
genes <- wc.p$word
for(i in genes){
  t <- FetchData(object = c0,
                 vars = c(names(tg),i),
                 slot = 'scale.data')
  t <- t[!is.na(t[,1]),]
  t <- t[order(t[,1],decreasing = F),] #only counts for regulated genes
  t <- melt(t, id.vars = i)
  colnames(t) <- c('Y','Path','X')

  f7d[[i]] <- ggplot(data = t, aes(x = X, y = Y,
                                   color = Path)) +
    scale_color_manual(values = col.paths) +
    theme(axis.text = element_blank(),
          strip.background = element_blank(),
          strip.placement = 'inside',
          strip.text = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8,
                                      face = 'italic'),
          axis.line = element_line(),
          axis.ticks = element_blank(),
          legend.key = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    geom_smooth(method = gam::gam, se = F,
                formula = y ~ lo(x), size = 0.5) +
    labs(y = paste(i)) +
    facet_wrap(~Path, nrow = 1, ncol = 4, scales = 'free_x')
}

#Fig. 7E-F Process variation####
GOtop <- c()
for(i in 1:10){
  tmp <- goResults[[i]][go2[[i]],'category']
  GOtop <- c(GOtop, tmp)
}

keggResults <- list()
for(i in 1:length(net)){
  qs <- getBM(attributes = "ensembl_gene_id",
              filters = "mgi_symbol",
              values = net[[i]],
              mart = ensembl)$ensembl_gene_id
  qs <- qs[!is.na(qs)]

  genes <- as.integer(names(TxLen) %in% qs)
  names(genes) <- names(TxLen)

  pwf <- nullp(DEgenes = genes,
               genome = genome,
               id = "ensGene",
               bias.data = as.numeric(TxLen))

  keggResults[[i]] <- goseq(pwf,
                            genome = genome,
                            id = "ensGene",
                            test.cats=c("KEGG")) %>%
    mutate(hitsPerc=numDEInCat*100/numInCat) %>%
    filter(over_represented_pvalue < 0.01, numDEInCat > 2)  %>%
    dplyr::select(category,over_represented_pvalue,
                  hitsPerc,numDEInCat)
}

KEGGtop <- c()
for(i in 1:10){
  tmp <- keggResults[[i]][,'category']
  KEGGtop <- c(KEGGtop, tmp)
}

KEGGtop <- unique(KEGGtop)
kg <- list()
for(i in KEGGtop){
  acc <- paste0("mmu",i)
  input <- keggGet(acc)
  term <- input[[1]]$PATHWAY_MAP
  genes <- gsub('([A-Za-z0-9]+)(;)(.+)','\\1',input[[1]]$GENE)[c(F,T)]
  genes <- genes[genes %in% vc$name]
  kg[[term]] <- genes
}

gg <- list()
for(i in GOtop){
  genes <- unlist(mget(get(i, org.Mm.egGO2ALLEGS),org.Mm.egSYMBOL))
  genes <- genes[genes %in% vc$name]
  gg[[i]] <- genes
}

GOnames <- c()
for(i in 1:10){
  tmp <- goResults[[i]][go2[[i]],'term']
  GOnames <- c(GOnames, tmp)
}

names(gg) <- GOnames
gg <- lapply(gg, unname)
gg <- lapply(gg, unique)

kg <- kg[c(1:5,7,9,12,15,
           16,20,21,35,
           39,40)]

c0 <- AddModuleScore(object = c0, features = gg,
                     assay = 'integrated',
                     name = 'GO')

c0 <- AddModuleScore(object = c0, features = kg,
                     assay = 'integrated',
                     name = 'KEGG')

pst <- c0@meta.data[,grep('Path',colnames(c0@meta.data))]
sets <- colnames(pst)
pval.go <- list()
var.go <- list()
for(i in sets){
  cells <- rownames(pst[!is.na(pst[,i]),])
  sub <-  subset(c0, cells = cells)
  df <- t(sub@meta.data[,grep('GO',colnames(sub@meta.data))])
  # Fit GAM for each gene using pseudotime as independent variable.
  t <- FetchData(object = sub, vars = i)
  gam.pval <- apply(df, 1, function(z){
    d <- data.frame(z=z, t=t)
    tmp <- gam(z ~ lo(t), data=d)
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  pval.go[[i]] <- p.adjust(gam.pval, method = 'fdr')
  var.go[[i]] <- apply(df, 1, var)
  var.go[[i]][which(pval.go[[i]] > 0.01)] <- 0
}

ph <- as.matrix(data.frame(var.go, row.names = GOnames))
ann_colors <- list(Paths = c('Phagocytic path'="#FF7F00",
                             'Oxidative stress path'="#E31A1C",
                             'Inflammatory path'="#CAB2D6",
                             'Remodelling path'="#6A3D9A"))

f7e <- pheatmap(ph, show_colnames = F,
                treeheight_row = 10 , cluster_cols = T,
                annotation_col = data.frame('Paths' = titles),
                annotation_colors = ann_colors,
                color = colorRampPalette(c('#191970'	,'#FFFFFF', '#FF8C00'))(255),
                border_color = 'grey', scale = 'row',
                fontsize_row = 6, fontsize = 8,
                cellwidth = 6, cellheight = 6)

#KEGG
pval.kg <- list()
var.kg <- list()
for(i in sets){
  cells <- rownames(pst[!is.na(pst[,i]),])
  sub <-  subset(c0, cells = cells)
  df <- t(sub@meta.data[,grep('KEGG',colnames(sub@meta.data))])
  # Fit GAM for each gene using pseudotime as independent variable.
  t <- FetchData(object = sub, vars = i)
  gam.pval <- apply(df, 1, function(z){
    d <- data.frame(z=z, t=t)
    tmp <- gam(z ~ lo(t), data=d)
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  pval.kg[[i]] <- p.adjust(gam.pval, method = 'fdr')
  var.kg[[i]] <- apply(df, 1, var)
  var.kg[[i]][which(pval.kg[[i]] > 0.01)] <- 0
}

ph <- as.matrix(data.frame(var.kg, row.names = names(kg)))
ann_colors <- list(Paths = c('Phagocytic path'="#FF7F00",
                             'Oxidative stress path'="#E31A1C",
                             'Inflammatory path'="#CAB2D6",
                             'Remodelling path'="#6A3D9A"))

f7f <- pheatmap(ph, show_colnames = F,
                treeheight_row = 10 , cluster_cols = T,
                annotation_col = data.frame('Paths' = titles),
                annotation_colors = ann_colors,
                color = colorRampPalette(c('#191970'	,'#FFFFFF', '#FF8C00'))(255),
                border_color = 'grey', scale = 'row',
                fontsize_row = 6, fontsize = 8,
                cellwidth = 6, cellheight = 6)
