---
title: "Dual organism analysis - integrated libraries - malaria"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 3
    toc_float: true
    number_sections: true
    code_folding: hide
  pdf_document:
    keep_md: true
    toc: true
    toc_depth: 3
    number_sections: true
    df_print: kable
always_allow_html: true
---


```r
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
knitr::opts_chunk$set(fig.width = 10, fig.height = 8) 
```


```r
library(Seurat)
library(patchwork)
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(broom)
library(stringr)
library(ggplot2)
library(SeuratData)
#library(celldex)
#install.packages("cowplot")
library(cowplot)
#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
#BiocManager::install('multtest', force = T)
library(multtest)
library(RColorBrewer)
# install.packages('metap')
#library(metap)
#remotes::install_github("satijalab/seurat-data")
```


```r
# DotPlot function

myDotPlot <- function(
  object,
  assay = NULL,
  features,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA,
  text.size = 10
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot() +
    theme(axis.text.y=element_text(size=text.size)) +
    theme(axis.text.x=element_text(size=text.size, angle = 45, vjust = 1, hjust=1)) +
    coord_flip()
  
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  return(plot)
}
```



```r
# EXP1_l2 <- readRDS("Jan27analysis/EXP1_l2_res0.9.rds")
# EXP2_l1 <- readRDS("Jan27analysis/EXP2_l1_res0.9.rds")
# EXP2_l2 <- readRDS("Jan27analysis/EXP2_l2_res0.9.rds")
# EXP2_l3 <- readRDS("Jan27analysis/EXP2_l3_res0.9.rds")
file_1_l1 = list.files("robjects/",pattern = "^dual_EXP1_l1.*\\_cellbender_plasmo.bg.removed.rds$")
file_1_l2 = list.files("robjects/",pattern = "^dual_EXP1_l2.*\\_cellbender_plasmo.bg.removed.rds$")
file_2_l1 = list.files("robjects/",pattern = "^dual_EXP2_l1.*\\_cellbender_plasmo.bg.removed.rds$")
file_2_l2 = list.files("robjects/",pattern = "^dual_EXP2_l2_reseq.*\\_cellbender_plasmo.bg.removed.rds$")
file_2_l3 = list.files("robjects/",pattern = "^dual_EXP2_l3.*\\_cellbender_plasmo.bg.removed.rds$")

EXP1_l1 <- readRDS(paste0("robjects/",file_1_l1))
EXP1_l2 <- readRDS(paste0("robjects/",file_1_l2))
EXP2_l1 <- readRDS(paste0("robjects/",file_2_l1))
EXP2_l2 <- readRDS(paste0("robjects/",file_2_l2))
EXP2_l3 <- readRDS(paste0("robjects/",file_2_l3))
```

# Before integration {.tabset}

## Included cells

```r
DimPlot(EXP1_l1, reduction = "umap", group.by = "donor_id") + plot_annotation(title = "EXP1_l1")
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
ggplot(EXP1_l1@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(binwidth = 50) +
  facet_wrap( ~ donor_id)
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```r
DimPlot(EXP1_l2, reduction = "umap", group.by = "donor_id") + plot_annotation(title = "EXP1_l2")
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

```r
ggplot(EXP1_l2@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(binwidth = 50) +
  facet_wrap( ~ donor_id)
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-4.png)<!-- -->

```r
DimPlot(EXP2_l1, reduction = "umap", group.by = "donor_id") + plot_annotation(title = "EXP2_l1")
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-5.png)<!-- -->

```r
ggplot(EXP2_l1@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(binwidth = 50) +
  facet_wrap( ~ donor_id)
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-6.png)<!-- -->

```r
DimPlot(EXP2_l2, reduction = "umap", group.by = "donor_id") + plot_annotation(title = "EXP2_l2")
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-7.png)<!-- -->

```r
ggplot(EXP2_l2@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(binwidth = 50) +
  facet_wrap( ~ donor_id)
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-8.png)<!-- -->

```r
DimPlot(EXP2_l3, reduction = "umap", group.by = "donor_id") + plot_annotation(title = "EXP2_l3")
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-9.png)<!-- -->

```r
ggplot(EXP2_l3@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(binwidth = 50) +
  facet_wrap( ~ donor_id)
```

![](integrate_dual_files/figure-html/unnamed-chunk-4-10.png)<!-- -->


```r
exp.list <- list(EXP1_l1 = EXP1_l1,
  EXP1_l2 = EXP1_l2,
  EXP2_l1 = EXP2_l1,
  EXP2_l2 = EXP2_l2,
  EXP2_l3 = EXP2_l3
            )

# normalize and identify variable features for each dataset independently
exp.list <- lapply(X = exp.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = exp.list)
# retrieve anchors

anchors <- FindIntegrationAnchors(object.list = exp.list, anchor.features = features)

# this command creates an 'integrated' data assay
dual.integrated <- IntegrateData(anchorset = anchors)
```


```r
# Let’s set the assay to RNA and visualize the datasets before integration.
dual.integrated.RNA <- dual.integrated
DefaultAssay(dual.integrated.RNA) <- "RNA"

#Let’s do normalization, HVG finding, scaling, PCA, and UMAP on the un-integrated (RNA) assay:

dual.integrated.RNA <- NormalizeData(dual.integrated.RNA, verbose = F)
dual.integrated.RNA <- FindVariableFeatures(dual.integrated.RNA, selection.method = "vst", nfeatures = 2000, verbose = F)
dual.integrated.RNA <- ScaleData(dual.integrated.RNA, verbose = F)
dual.integrated.RNA <- RunPCA(dual.integrated.RNA, npcs = 30, verbose = F)
dual.integrated.RNA <- RunUMAP(dual.integrated.RNA, reduction = "pca", dims = 1:30, verbose = F)
```


```r
## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## This message will be shown once per session

# UMAP plot of the datasets before integration shows clear separation. Note that we can use patchwork syntax with Seurat plotting functions:

DimPlot(dual.integrated.RNA,reduction = "umap", group.by = "orig.ident") + ggtitle("Before Integration")
```

![](integrate_dual_files/figure-html/unnamed-chunk-7-1.png)<!-- -->
# After integration {.tabset}


```r
DefaultAssay(dual.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
dual.integrated <- ScaleData(dual.integrated, verbose = FALSE)
dual.integrated <- RunPCA(dual.integrated, npcs = 30, verbose = FALSE)
dual.integrated <- RunUMAP(dual.integrated, reduction = "pca", dims = 1:30)
dual.integrated <- FindNeighbors(dual.integrated, reduction = "pca", dims = 1:30)
dual.integrated <- FindClusters(dual.integrated, resolution = 0.5)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 15543
## Number of edges: 582785
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8872
## Number of communities: 12
## Elapsed time: 2 seconds
```

```r
dual.integrated@meta.data <- dual.integrated@meta.data %>% 
  mutate(Donorid = case_when(
    orig.ident %in% c("EXP1_l2", "EXP1_l1") & donor_id == "donor0" ~ "D1",
    orig.ident %in% c("EXP1_l2", "EXP1_l1") & donor_id == "donor1" ~ "D2",
    orig.ident %in% c("EXP2_l1", "EXP2_l2", "EXP2_l3") & donor_id == "donor0" ~ "D3",
    orig.ident %in% c("EXP2_l1", "EXP2_l2", "EXP2_l3") & donor_id == "donor1" ~ "D4",
    orig.ident %in% c("EXP2_l1", "EXP2_l2", "EXP2_l3") & donor_id == "donor2" ~ "D5"
  ))

#saveRDS(dual.integrated, "dual.integrated_all_libs_preprocessed.rds")

DimPlot(dual.integrated,reduction = "umap", group.by = "orig.ident") + ggtitle("After Integration")
```

![](integrate_dual_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


```r
#DefaultAssay(dual.integrated) <- "RNA"
DimPlot(dual.integrated, group.by = "Donorid",reduction = "umap", label = "F")
```

![](integrate_dual_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


```r
DefaultAssay(dual.integrated) <- "RNA"
FeaturePlot(dual.integrated, features = c("XIST"))
```

![](integrate_dual_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


```r
FeaturePlot(dual.integrated, features = c("CD14", "FCGR3A"))
```

![](integrate_dual_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


```r
DefaultAssay(dual.integrated) <- "integrated"
DimPlot(dual.integrated, group.by = "seurat_clusters",reduction = "umap", label = "T")
```

![](integrate_dual_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


```r
DimPlot(dual.integrated, group.by = "Group",reduction = "umap", label = "T")
```

![](integrate_dual_files/figure-html/unnamed-chunk-13-1.png)<!-- -->
 
## Plasmodium RNA in samples

```r
breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

ggplot(dual.integrated@meta.data, aes(x = percent.plasmo)) +
  geom_histogram(bins = 70, fill = "white", colour = "blue") +
  facet_wrap(~ factor(Group, levels = c("control", "sporozoite_bystander", "sporozoite_infected",
                      "RBC_control", "RBC_bystander", "RBC_infected"))) +
  theme_bw() +
  scale_y_log10(minor_breaks = c(5, 50, 500, 5000)) +
  theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "gray75"),
      panel.grid.minor = element_line(colour = "gray"),
      panel.border = element_blank(),
      panel.background = element_blank()) +
      theme(axis.text = element_text(size = 11)) +
      theme(strip.background =element_rect(fill="aliceblue", colour = "white"),
            strip.text.x = element_text(size = 11))
```

![](integrate_dual_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
#ggsave("plasmo.percent_by_conditions_cellbender_plasmo.bg.removed.pdf")
```


```r
plasmo.percent.table <- dual.integrated@meta.data %>% 
  tibble::rownames_to_column("Cell") %>% 
  select(Cell, Group, percent.plasmo, percent.mt) %>% 
  group_by(Group) %>% 
  summarise(plasmo_median = median(percent.plasmo, na.rm = T),
            plasmo_max = max(percent.plasmo, na.rm = T),
            plasmo_min = min(percent.plasmo, na.rm = T),
            mt_median = median(percent.mt, na.rm = T),
            number_of_cells = n())

DT::datatable(plasmo.percent.table)
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-0c2f40e546d48a362064" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0c2f40e546d48a362064">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","sporozoite_infected"],[0,0,0,0,0,0],[97.4459724950884,32.4724172517553,93.574297188755,36.2407467009977,32.199151257957,1.91997129014893],[0,0,0,0,0,0],[5.34048805417146,5.77695396972505,6.22156365281484,6.38019642336999,5.76335241482799,4.97453974733218],[1795,4035,2098,4898,2665,52]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Group<\/th>\n      <th>plasmo_median<\/th>\n      <th>plasmo_max<\/th>\n      <th>plasmo_min<\/th>\n      <th>mt_median<\/th>\n      <th>number_of_cells<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false},"selection":{"mode":"multiple","selected":null,"target":"row","selectable":null}},"evals":[],"jsHooks":[]}</script>
```
 
## Features in clusters

```r
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
lib_name = "integrated"
# percent.mt on umap

meta <- dual.integrated@meta.data %>% 
  rownames_to_column("Cell") %>% 
  select("Cell", "nFeature_RNA", "percent.mt", "nCount_RNA", "percent.plasmo")

umap_embed <- dual.integrated@reductions[["umap"]]@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column("Cell") %>%
  left_join(., meta) %>% 
  mutate(lib = lib_name)

# mt.percent on umap
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, max(umap_embed$percent.mt)))
ggplot(umap_embed, aes(UMAP_1, UMAP_2, colour = percent.mt)) +
  geom_point(size = 0.3, alpha = 0.8) +
  facet_wrap( ~ lib, ncol = 2, scales = "free") +
  sc +
  theme_bw()
```

![](integrate_dual_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

```r
# nF on umap
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(umap_embed$nFeature_RNA), max(umap_embed$nFeature_RNA)))
ggplot(umap_embed, aes(UMAP_1, UMAP_2, colour = nFeature_RNA)) +
  geom_point(size = 0.3, alpha = 0.8) +
  facet_wrap( ~ lib, ncol = 2, scales = "free") +
  sc +
  theme_bw()
```

![](integrate_dual_files/figure-html/unnamed-chunk-16-2.png)<!-- -->

```r
# nCount on umap
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(umap_embed$nCount_RNA), max(umap_embed$nCount_RNA)))
ggplot(umap_embed, aes(UMAP_1, UMAP_2, colour = nCount_RNA)) +
  geom_point(size = 0.3, alpha = 0.8) +
  facet_wrap( ~ lib, ncol = 2, scales = "free") +
  sc +
  theme_bw()
```

![](integrate_dual_files/figure-html/unnamed-chunk-16-3.png)<!-- -->

```r
# parasite on umap
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(umap_embed$percent.plasmo), max(umap_embed$percent.plasmo)))
ggplot(umap_embed, aes(UMAP_1, UMAP_2, colour = percent.plasmo)) +
  geom_point(size = 0.3, alpha = 0.8) +
  facet_wrap( ~ lib, ncol = 2, scales = "free") +
  sc +
  theme_bw()
```

![](integrate_dual_files/figure-html/unnamed-chunk-16-4.png)<!-- -->

```r
#ggsave("plasmo.percent_umap.pdf")
```
## Clusters

```r
DimPlot(dual.integrated, group.by = "seurat_clusters",reduction = "umap", label = "T")
```

![](integrate_dual_files/figure-html/unnamed-chunk-17-1.png)<!-- -->
## Plasmo percent in clusters

```r
plasmo.percent.clusters.table <- dual.integrated@meta.data %>% 
  tibble::rownames_to_column("Cell") %>% 
  select(Cell, Group, percent.plasmo, percent.mt, seurat_clusters) %>% 
  group_by(seurat_clusters, Group) %>% 
  summarise(plasmo_median = median(percent.plasmo, na.rm = T),
            plasmo_max = max(percent.plasmo, na.rm = T),
            plasmo_min = min(percent.plasmo, na.rm = T),
            mt_median = median(percent.mt, na.rm = T),
            number_of_cells = n())

DT::datatable(plasmo.percent.clusters.table)
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-c27cb71094a3ca8b006a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-c27cb71094a3ca8b006a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60"],["0","0","0","0","0","0","1","1","1","1","1","1","2","2","2","2","2","3","3","3","3","3","4","4","4","4","4","4","5","5","5","5","6","6","6","6","6","7","7","7","7","7","8","8","8","8","8","8","9","9","9","9","9","10","10","10","10","10","11","11"],["RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","sporozoite_infected","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","sporozoite_infected","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","sporozoite_infected","RBC_bystander","RBC_infected","sporozoite_bystander","sporozoite_infected","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","sporozoite_infected","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","RBC_control","control"],[0,0,0,0,0,0,0,0,0,0,0,0.249478153603642,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8.88826025459689,12.1673416959066,12.4423137237331,10.5039309189329,8.95009077154923,0,0,0,0,0,0,0,0,0,0,0,0,0.0965178852663381,0,0,0,0,0],[0.00551724137931034,1.00637370010064,31.2552113398756,2.98455964806349,0,0,26.0606896551724,0.0068667170225915,48.9840086569677,2.9004329004329,3.84926598081147,1.91997129014893,0.00513267977210902,4.21443020903574,0.0511116790186558,0.502318392581144,0,0,0.537634408602151,0,11.8350898946063,0,0.559952619393744,0,19.9033266988911,0,32.199151257957,0.134394706916158,18.6173131304605,37.1052425358514,6.64556962025316,0.58008378988076,0,1.10584518167457,0,5.50284629981025,0,97.4459724950884,32.4724172517553,93.574297188755,36.2407467009977,25.6024096385542,0.0134093194770365,0.0667556742323097,0.410201533797039,0.104493207941484,0.0534045393858478,0,0.146929180135175,0.0179856115107914,8.01389481382168,0.148927720413026,6.40396307138032,1.23922413793103,15.0375939849624,25.4183354364404,12.1043201274139,10.0125156445557,0.0190876121397213,0.0617093489663684],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.27421555252387,5.28459046737621,2.94544794408641,2.81255670477227,4.5618689660255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[4.70045605710078,6.37712804096134,7.50856110526958,6.81172008823269,4.27535004516712,12.4819740426214,5.59609076912779,5.69520816967793,6.38020833333333,5.85499642260911,5.59008063073444,5.08212045429624,3.06468880890336,4.66462824838549,6.17713008658813,5.35894843276036,3.81610044144828,1.55317935803484,6.32485601477979,11.6231438812084,6.82268464388903,3.01887692282854,4.65805211095807,5.48703750314624,5.84484590860786,5.49533273110509,5.97145488029466,4.15789945223268,5.25644143920281,6.49260628465804,5.80368059757077,4.85610722018436,4.26026557392947,5.65966799423909,13.3136094674556,5.77681276530907,5.82511554059971,5.3674965893588,5.52040136411163,5.81057112451224,5.7997557997558,5.55655287660134,6.57220178424933,6.63839021503255,7.33050591312758,7.3814086498357,6.38005496662741,5.72139303482587,6.64943982081773,6.67379679144385,9.88470609123928,8.04959527174611,6.4804469273743,6.76810073452256,7.60841719574665,6.62608110989122,6.17764441119652,4.54751975055661,6.14954018392643,4.96244024584566],[4,1506,8,1952,38,1,811,7,947,11,854,14,6,1228,12,1339,4,2,752,1,1054,2,480,1,265,3,1047,12,338,349,512,22,4,298,1,282,2,17,24,447,33,18,32,89,24,76,59,3,54,59,14,59,71,47,50,30,40,58,21,49]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>seurat_clusters<\/th>\n      <th>Group<\/th>\n      <th>plasmo_median<\/th>\n      <th>plasmo_max<\/th>\n      <th>plasmo_min<\/th>\n      <th>mt_median<\/th>\n      <th>number_of_cells<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false},"selection":{"mode":"multiple","selected":null,"target":"row","selectable":null}},"evals":[],"jsHooks":[]}</script>
```

## Top 10 markers in each cluster

```r
imm.comb.markers <- FindAllMarkers(dual.integrated, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)

imm.comb.markers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top15

DT::datatable(top15)
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-9cfd0b120550a8bb9708" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-9cfd0b120550a8bb9708">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180"],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.02315030974325e-282,9.34030571662843e-98,0,0,0,0,0,0,0,0,0,0,0,0,0,4.77655054293364e-93,4.34356273422999e-89,0,0,0,0,0,0,0,0,0,0,0,1.57788973009247e-270,1.13405918797995e-250,2.05092132331414e-224,4.04684239884145e-54,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.57293304600581e-219,2.0047999766087e-184,1.39487415069459e-304,6.62344668675551e-219,2.70052569850533e-194,1.76390267778197e-189,1.25648022300792e-188,8.09444694739177e-184,2.66067416829749e-178,6.79267131706937e-160,1.57508486174634e-136,1.14667702629575e-89,4.83548890033289e-79,4.68790198493485e-72,4.73565956858582e-72,7.94420353423464e-66,1.77293496556957e-64,0,0,0,0,0,0,0,0,0,0,0,0,0,2.95834312861033e-288,1.01925610108114e-213,0,0,1.11031554325206e-222,5.96672570132027e-157,1.2387739639661e-156,1.64884947628201e-148,3.88934676548272e-133,8.79050703209621e-133,2.16957235503936e-126,5.15960942627117e-126,1.26232280169728e-119,3.54924613676527e-117,2.98631369216987e-105,5.18139781597379e-84,4.934566678198e-73,4.83673539531148e-141,2.9686060679564e-121,4.86594027057332e-114,7.86803621060564e-106,4.73271427208662e-100,1.85686682204531e-98,1.41132083205128e-95,3.54070288026722e-93,1.51486482642133e-85,3.33260213357913e-75,5.0236170106549e-72,1.59000946092959e-71,5.0884930092223e-65,2.83957834207059e-52,9.42641617455235e-12,5.2631381908886e-71,6.10169165308325e-50,2.50416186499542e-27,2.66330427978073e-27,1.64918321991121e-26,4.51308465504807e-25,2.11839822092305e-24,2.24526703816107e-23,1.22562758878398e-18,2.38995656550293e-17,5.62305388848679e-17,9.82338244387169e-17,9.48148881357151e-16,1.75699361701766e-11,2.74260796278311e-11,2.83159792792454e-45,9.59194168224323e-43,7.28798705065416e-41,4.64499498803637e-30,7.40532628315074e-29,5.0343571318721e-25,7.29488931016864e-20,1.34577486269165e-13,3.58374840255203e-13,3.67500488999417e-13,1.27795174750041e-12,4.0026432688998e-11,3.03980258331622e-10,1.93237292217606e-05,0.000554126803832089],[1.64635211182031,1.5921674921418,1.46958171433367,1.46544724575509,1.39048204735267,1.25203399918599,1.2145542192569,1.18564981491581,1.0906574365374,1.08810382868369,1.07885050306291,1.06405572160724,1.06243912635326,1.05792821652586,1.0494667726906,2.19475500382871,1.91994484189888,1.70307217652622,1.67383563545982,1.64905034256661,1.62249082613083,1.57234530902446,1.57198231170454,1.52791870655851,1.52757792555402,1.49220683106461,1.47323623417234,1.43073294015054,1.77952605050762,2.05911010285759,1.87308713347854,1.85332873971839,1.48552901483155,1.3489247331536,1.34295342455541,1.25379166669603,1.18364409973477,1.08712531738707,1.03279202141193,1.00917967378101,0.967866343810586,0.901950287566915,0.900732313718769,0.983785872096832,0.838468491030242,1.72752201707352,1.51489385270028,1.45186929605984,1.35583236374618,1.33495829161768,1.2454900591799,1.15615613934205,1.154645436696,1.08328794909652,1.02853023141871,0.998074281183812,1.00470508629359,1.52966952741692,0.999034732751798,1.09211788299605,3.7019773898761,3.44590575846871,2.71923871202051,2.68269220760423,2.63494482425092,2.59460650441769,2.57934622364437,2.37382807315048,2.30224975425133,2.29122363608439,2.17529113732593,2.0684816868545,2.06711513812586,2.05205432791596,1.96462334007841,3.16742984717104,3.1539195673053,3.08161280528728,3.04913053430372,2.96931388385359,2.96618312268353,2.70324583176744,2.48479125306132,2.41014070517554,2.27159740520602,2.26127273968735,2.19721239425389,1.96693249254843,2.25519893574326,1.98659008666436,2.70108368226607,1.38325912944802,1.28437064890897,1.49180331254198,2.98680600714069,1.56668710243237,1.22181902286537,1.57734551519695,1.02033738435895,1.04194639619125,1.09895252457641,1.16650301908364,1.06680097957014,1.00859382120363,1.03951117164663,4.92342376901689,4.41229699045906,4.38428057331431,4.31226393374733,4.26822019008158,4.22946070990769,3.42061993364208,3.36282866327461,3.10563609209866,3.02877520649897,2.97779315872418,2.94345298399783,2.73246257925687,2.74430815063008,2.72745991855998,2.56725011143928,2.45030866499946,2.56610855851818,2.48174865619425,2.87824430598644,3.94289131711501,2.69031994623528,4.95197610814512,2.46274821500405,2.9879445767065,2.79962968579053,3.13333079890579,3.88747754984005,2.52669155667444,2.52227793459982,2.30557787060835,2.46829845280456,4.51024760312886,2.89259085901258,2.62446037384393,4.02418606974318,3.56342971555545,3.07961693395484,2.6803400820722,2.36773312798264,2.56773523679008,2.29814466432275,3.34538197591987,2.46778186852417,2.37976780325922,1.99344565109117,1.16564233078827,1.05050607075897,1.32425449187032,1.07573770096506,0.837757810601856,0.80193742021375,1.69726723421269,0.812269393353164,0.698094284760726,0.760495569616814,1.02202931140677,0.664261818735551,0.636509657736839,0.807225357029982,3.67551787549602,2.38898472402759,4.03437033196216,2.24531768694351,2.54776446118306,1.79431033758562,1.67983777558452,1.02523809692599,1.04264427810981,0.959230869166949,1.02304673512751,1.12742601487271,0.992818152885984,0.903264987118981,1.1488689870041],[0.956,0.951,0.981,0.999,0.951,0.863,0.982,0.858,0.935,0.863,0.9,0.871,0.906,0.999,0.814,1,0.88,0.969,0.924,1,0.908,0.997,0.995,0.991,0.889,0.944,0.813,0.986,0.764,0.65,0.992,0.996,0.988,0.958,0.994,0.98,0.996,0.869,0.961,0.979,0.966,0.976,0.971,0.826,0.657,0.925,0.98,0.929,0.903,0.876,0.972,0.814,0.977,0.961,0.97,0.897,0.949,0.808,0.854,0.66,0.863,0.834,0.994,0.981,0.933,0.89,0.97,0.936,0.975,0.986,0.931,0.997,0.835,0.934,0.899,0.977,1,1,0.98,0.999,0.742,0.989,0.962,0.988,0.989,0.975,0.994,0.975,0.811,0.803,0.988,0.997,0.991,0.985,0.917,0.988,0.925,0.968,0.872,0.976,0.942,0.905,0.923,0.855,0.944,0.948,0.97,0.983,0.981,0.959,0.981,0.97,0.942,0.935,0.941,0.931,0.946,0.889,0.937,0.857,0.848,0.813,0.795,0.852,0.82,0.929,0.915,0.837,0.781,0.837,0.933,0.883,0.901,0.873,0.88,1,0.973,0.86,0.786,0.875,0.868,0.914,0.879,0.922,0.926,0.821,0.813,0.844,0.946,0.307,0.942,0.956,0.902,0.862,0.916,0.884,0.84,0.929,0.924,0.533,0.8,0.742,0.973,0.942,0.773,1,1,0.971,0.986,0.871,0.886,0.8,0.957,0.9,0.8,0.957,0.829,0.914,0.714,0.343],[0.658,0.819,0.897,0.959,0.714,0.632,0.905,0.566,0.864,0.646,0.656,0.687,0.731,0.976,0.569,0.866,0.637,0.763,0.612,0.969,0.655,0.891,0.855,0.861,0.627,0.712,0.606,0.86,0.624,0.578,0.855,0.895,0.84,0.792,0.922,0.802,0.95,0.641,0.841,0.821,0.736,0.876,0.815,0.783,0.608,0.762,0.845,0.751,0.598,0.591,0.836,0.509,0.789,0.805,0.812,0.677,0.787,0.6,0.635,0.583,0.42,0.402,0.714,0.674,0.573,0.589,0.779,0.672,0.704,0.704,0.689,0.881,0.581,0.661,0.567,0.714,0.8,0.791,0.751,0.83,0.226,0.782,0.598,0.464,0.762,0.626,0.763,0.761,0.622,0.543,0.782,0.98,0.975,0.835,0.755,0.86,0.572,0.814,0.553,0.9,0.827,0.737,0.844,0.632,0.875,0.099,0.052,0.111,0.414,0.134,0.419,0.378,0.096,0.035,0.089,0.185,0.021,0.118,0.388,0.321,0.095,0.025,0.085,0.187,0.123,0.526,0.593,0.187,0.203,0.181,0.558,0.519,0.667,0.685,0.667,0.978,0.779,0.246,0.089,0.48,0.366,0.736,0.51,0.676,0.804,0.539,0.521,0.66,0.738,0.212,0.709,0.755,0.719,0.789,0.785,0.694,0.673,0.739,0.825,0.252,0.693,0.597,0.936,0.892,0.703,0.789,0.721,0.614,0.894,0.652,0.668,0.662,0.922,0.858,0.634,0.86,0.623,0.828,0.746,0.498],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.04630061948651e-279,1.86806114332569e-94,0,0,0,0,0,0,0,0,0,0,0,0,0,9.55310108586729e-90,8.68712546845997e-86,0,0,0,0,0,0,0,0,0,0,0,3.15577946018493e-267,2.26811837595991e-247,4.10184264662829e-221,8.09368479768291e-51,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.14586609201162e-216,4.0095999532174e-181,2.78974830138917e-301,1.3246893373511e-215,5.40105139701067e-191,3.52780535556394e-186,2.51296044601584e-185,1.61888938947835e-180,5.32134833659498e-175,1.35853426341387e-156,3.15016972349268e-133,2.29335405259149e-86,9.67097780066578e-76,9.37580396986969e-69,9.47131913717164e-69,1.58884070684693e-62,3.54586993113914e-61,0,0,0,0,0,0,0,0,0,0,0,0,0,5.91668625722067e-285,2.03851220216227e-210,0,0,2.22063108650412e-219,1.19334514026405e-153,2.47754792793221e-153,3.29769895256402e-145,7.77869353096544e-130,1.75810140641924e-129,4.33914471007873e-123,1.03192188525423e-122,2.52464560339457e-116,7.09849227353054e-114,5.97262738433974e-102,1.03627956319476e-80,9.86913335639599e-70,9.67347079062296e-138,5.9372121359128e-118,9.73188054114664e-111,1.57360724212113e-102,9.46542854417324e-97,3.71373364409062e-95,2.82264166410256e-92,7.08140576053443e-90,3.02972965284267e-82,6.66520426715826e-72,1.00472340213098e-68,3.18001892185917e-68,1.01769860184446e-61,5.67915668414117e-49,1.88528323491047e-08,1.05262763817772e-67,1.22033833061665e-46,5.00832372999084e-24,5.32660855956147e-24,3.29836643982242e-23,9.02616931009614e-22,4.2367964418461e-21,4.49053407632214e-20,2.45125517756797e-15,4.77991313100585e-14,1.12461077769736e-13,1.96467648877434e-13,1.8962977627143e-12,3.51398723403532e-08,5.48521592556623e-08,5.66319585584908e-42,1.91838833644865e-39,1.45759741013083e-37,9.28998997607274e-27,1.48106525663015e-25,1.00687142637442e-21,1.45897786203373e-16,2.69154972538331e-10,7.16749680510406e-10,7.35000977998834e-10,2.55590349500082e-09,8.0052865377996e-08,6.07960516663243e-07,0.0386474584435211,1],["0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","4","4","4","4","4","4","4","4","4","4","4","4","4","4","4","5","5","5","5","5","5","5","5","5","5","5","5","5","5","5","6","6","6","6","6","6","6","6","6","6","6","6","6","6","6","7","7","7","7","7","7","7","7","7","7","7","7","7","7","7","8","8","8","8","8","8","8","8","8","8","8","8","8","8","8","9","9","9","9","9","9","9","9","9","9","9","9","9","9","9","10","10","10","10","10","10","10","10","10","10","10","10","10","10","10","11","11","11","11","11","11","11","11","11","11","11","11","11","11","11"],["FOS","GAPT","IFNGR1","S100A4","MNDA","NFXL1","CD14","CCR2","ZFP36L2","JAML","PECAM1","SORL1","OTULINL","CALM2","VCAN","CXCL8","CXCL3","IL1B","SGPP2","SOD2","CXCL2","IER3","SLC7A11","WTAP","TNIP3","CCL2","IL7R","ATP13A3","CXCL1","SERPINB2","LGALS2","STAT1","FGL2","VAMP5","PSME2","HLA-DQA1","GLUL","LRRK2","LAMP1","HLA-DQB1","TNFSF13B","FYB1","PSMB9","MT1X","MT1G","GPNMB","FUCA1","RNASE1","CD163","MRC1","ACP5","SLC7A8","CD36","DAB2","CD84","WWP1","CD9","LGMN","ARL4C","APOC1","CXCL10","CCL8","GBP1","IFIT3","RSAD2","IFIT2","ISG15","APOBEC3A","MX1","IDO1","TNFSF10","MT2A","CXCL9","HERC5","IFIT1","TNF","CCL4L2","CCL4","CCL3L3","CCL3","IL12B","IL1B","G0S2","DNAAF1","CD40","IL1RN","BCL2A1","CD300E","PTX3","INHBA","FCGR3A","HLA-DPA1","HLA-DPB1","PSMB9","VMO1","FGL2","GIMAP4","VAMP5","GIMAP7","MS4A7","HLA-DQA1","WARS1","HLA-DQB1","GBP5","LGALS2","gene-PF3D7-0207600","gene-PF3D7-1105100","gene-PF3D7-0617800","gene-PF3D7-0726000","gene-PF3D7-0930300","gene-PF3D7-0532000","gene-PF3D7-0725600","gene-PF3D7-1228600","gene-PF3D7-1133200","gene-PF3D7-0610400","gene-PF3D7-0617900","gene-PF3D7-0613800","gene-PF3D7-0818900","gene-PF3D7-0531600","gene-PF3D7-1016300","TRBC2","CD3D","KLRD1","SYNE2","TRBC1","NKG7","IFITM1","GNLY","GZMB","CD2","GIMAP7","CST7","CCL5","CD52","IL32","SEC61B","CRIP1","CCL22","RAMP1","GADD45A","CCR7","BIRC3","FSCN1","RAB9A","NUB1","LAMP3","CST7","LTB","FABP5","GZMB","CCL24","RHBDD2","RGCC","FBP1","CLEC5A","TSC22D1","EGR2","FABP5","MFSD12","CCL22","IL1R2","EREG","MMP9","BTG2","BHLHE41","MT1X","MT1F","MT1G","MT2A","MT1H","MT1E","MT1M","CD14","ADAP2","IER5L","FUCA1","VCAN","H1-2","CXCR4","CXCL5"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>p_val<\/th>\n      <th>avg_log2FC<\/th>\n      <th>pct.1<\/th>\n      <th>pct.2<\/th>\n      <th>p_val_adj<\/th>\n      <th>cluster<\/th>\n      <th>gene<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false},"selection":{"mode":"multiple","selected":null,"target":"row","selectable":null}},"evals":[],"jsHooks":[]}</script>
```

```r
#write.csv2(top10, "integrated_data_noTBcells_cluster_markers.csv", row.names = F, quote = F)
DefaultAssay(dual.integrated) <- "RNA"
#pdf("Heatmap_gene_expr_integrated_data_noTBcells.pdf")
#DoHeatmap(dual.integrated, label = F) + NoLegend() #+ theme(text = element_text(size = 6.5))
#dev.off()

top15 %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top5

DefaultAssay(dual.integrated) <- "RNA"
myDotPlot(dual.integrated, features = unique(top5$gene),  dot.scale = 2) + RotatedAxis()
```

![](integrate_dual_files/figure-html/unnamed-chunk-19-2.png)<!-- -->



```r
#dual.integrated <- subset(dual.integrated, subset = nFeature_RNA >= 1500)
```

# Remove T, B cells

```r
# remove T and B cells

# rename clusters
# dual.integrated <- RenameIdents(dual.integrated, 
#                                 `0` = "Classical monocytes I (FOS)", 
#                                 `1` = "Activated monocytes I (IL1B CXCL2/1)", 
#                                 `2` = "Classical monocytes II (STAT1)",
#                                 `3` = "Classical monocytes III (LGMN)",
#                                 `4` = "Activated monocytes II (IFIT2/3, ISG15)", 
#                                 `5` = "Activated monocytes III (IL12 CD40 TNF)",
#                                 `6` = "Plasmo",
#                                 `7` = "Non-classical monocytes (FCGR3A !CD14)",
#                                 `8` = "Classical monocytes IV (CCL24)",
#                                 `9` = "CD8 T-cells",
#                                 `10` = "cDC1", 
#                                 `11` = "B-cells")

# sept 11, 2023:
# rename clusters
dual.integrated <- RenameIdents(dual.integrated,
                                `0` = "Classical monocytes I (FOS)",
                                `1` = "Activated monocytes I (IL1B CXCL2/1)",
                                `2` = "Classical monocytes II (STAT1)",
                                `3` = "Classical monocytes III (LGMN)",
                                `4` = "Activated monocytes II (IFIT2/3, ISG15)",
                                `5` = "Activated monocytes III (IL12 CD40 TNF)",
                                `6` = "Non-classical monocytes (FCGR3A !CD14)",
                                `7` = "Plasmo",
                                `8` = "B/T",
                                `9` = "cDC1",
                                `10` = "Classical monocytes IV (CCL24)",
                                `11` = "B-cells")

dual.integrated$predicted.celltype <- Idents(dual.integrated)

pdf("dimplot_dual.integrated_withTB_cellbender_plasmo.bg.removed.pdf", height = 5, width = 10)
DimPlot(dual.integrated, label = F)
dev.off()
```

```
## png 
##   2
```

```r
top15 %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top5

#pdf("dotplot_dual.integrated_withTB_cellbender_plasmo.bg.removed.pdf", height = 5, width = 10)
myDotPlot(dual.integrated, features = unique(top5$gene),  dot.scale = 2, text.size = 5) #+ RotatedAxis() #, text.size = 10
```

![](integrate_dual_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

```r
#dev.off()

# saveRDS(dual.integrated, "dual_integrated_withTB_cellbender.rds")
```


```r
# plot plasmo.percent by celltype

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

ggplot(dual.integrated@meta.data, aes(x = percent.plasmo)) +
  geom_histogram(bins = 70, fill = "white", colour = "blue") +
  facet_wrap(~ factor(predicted.celltype)) +
  theme_bw() +
  scale_y_log10(minor_breaks = c(5, 50, 500, 5000)) +
  theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "gray75"),
      panel.grid.minor = element_line(colour = "gray"),
      panel.border = element_blank(),
      panel.background = element_blank()) +
      theme(axis.text = element_text(size = 11)) +
      theme(strip.background =element_rect(fill="aliceblue", colour = "white"),
            strip.text.x = element_text(size = 11))
```

![](integrate_dual_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
#ggsave("plasmo.percent_by_celltype_withTB_cellbender_plasmo.bg.removed.pdf", height = 10, width = 15)
```


```r
# remove T and B
#dual.integrated.TBremoved <- subset(dual.integrated, subset = predicted.celltype == "CD8 T-cells", invert = T)
dual.integrated.TBremoved <- subset(dual.integrated, subset = predicted.celltype == "B-cells", invert = T)
dual.integrated.TBremoved <- subset(dual.integrated.TBremoved, subset = predicted.celltype == "B/T", invert = T)

#pdf("dimplot_dual.integrated_noTB_cellbender_plasmo.bg.removed.pdf", height = 5, width = 10)
DimPlot(dual.integrated.TBremoved, label = F)
```

![](integrate_dual_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
#dev.off()

#saveRDS(dual.integrated.TBremoved, "dual.integrated.TBremoved.cellbender_plasmo.bg.removed.rds")
```


```r
#dual.TBremoved <- readRDS("robjects/dual.integrated.TBremoved.cellbender.rds")
dual.TBremoved <- dual.integrated.TBremoved
# remove T and B cells and only cells from human analysis
#dual.TBremoved <- dual.integrated[,colnames(dual.integrated) %in% colnames(immune.combined.TBremoved)]
```

<!-- ```{r} -->
<!-- # copy cell type and cell from only-human analysis (immune.combined.TBremoved) -->
<!-- # add in cell type -->
<!-- human.meta <- immune.combined.TBremoved@meta.data %>%  -->
<!--   rownames_to_column("CellHash") %>%  -->
<!--   select(CellHash, predicted.celltype) -->

<!-- dual.TBremoved@meta.data <- dual.TBremoved@meta.data %>%  -->
<!--   rownames_to_column("CellHash") %>% -->
<!--   left_join(., human.meta) %>%  -->
<!--   column_to_rownames("CellHash") -->
<!-- ``` -->

# Re-analyse after removing T, B cells {.tabset}

```r
DefaultAssay(dual.TBremoved) <- "RNA"

dual.TBremoved.list <- SplitObject(dual.TBremoved, split.by = "orig.ident")
dual.TBremoved.list <- lapply(X = dual.TBremoved.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = dual.TBremoved.list)
anchors <- FindIntegrationAnchors(object.list = dual.TBremoved.list, anchor.features = features)

# set k.weight so that it is not larger than the number of cells in the smallest dataset
# https://github.com/satijalab/seurat/issues/3936
dual.TBremoved <- IntegrateData(anchorset = anchors, k.weight = 46)
DefaultAssay(dual.TBremoved) <- "integrated"

dual.TBremoved <- ScaleData(dual.TBremoved, verbose = FALSE)
dual.TBremoved <- RunPCA(dual.TBremoved, npcs = 30, verbose = FALSE)
dual.TBremoved <- RunUMAP(dual.TBremoved, reduction = "pca", dims = 1:10)
dual.TBremoved <- FindNeighbors(dual.TBremoved, reduction = "pca", dims = 1:10)
dual.TBremoved <- FindClusters(dual.TBremoved, resolution = 0.5)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 15190
## Number of edges: 468542
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8897
## Number of communities: 12
## Elapsed time: 1 seconds
```


```r
#DefaultAssay(dual.TBremoved) <- "RNA"
DimPlot(dual.TBremoved, group.by = "Donorid",reduction = "umap", label = "F")
```

![](integrate_dual_files/figure-html/unnamed-chunk-26-1.png)<!-- -->


```r
DefaultAssay(dual.TBremoved) <- "RNA"
FeaturePlot(dual.TBremoved, features = c("XIST"))
```

![](integrate_dual_files/figure-html/unnamed-chunk-27-1.png)<!-- -->


```r
FeaturePlot(dual.TBremoved, features = c("CD14", "FCGR3A"))
```

![](integrate_dual_files/figure-html/unnamed-chunk-28-1.png)<!-- -->


```r
DefaultAssay(dual.TBremoved) <- "integrated"
DimPlot(dual.TBremoved, group.by = "seurat_clusters",reduction = "umap", label = "T")
```

![](integrate_dual_files/figure-html/unnamed-chunk-29-1.png)<!-- -->


```r
DimPlot(dual.TBremoved, group.by = "Group",reduction = "umap", label = "T")
```

![](integrate_dual_files/figure-html/unnamed-chunk-30-1.png)<!-- -->


```r
#pdf("dimplot_dual.integrated_noTB_cellbender_plasmo.bg.removed_reintegration.pdf", height = 5, width = 10)
DimPlot(dual.TBremoved, group.by = "predicted.celltype",reduction = "umap", label = F) + ggtitle("After re-integration")
```

![](integrate_dual_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

```r
#dev.off()
```


```r
dual.plasmo <- subset(dual.TBremoved, subset = percent.plasmo >= 1) # 7014 cells

#pdf("dimplot_dual.integrated_noTB_cellbender_plasmo.bg.removed_reintegration_plasmo_greater1.pdf", height = 5, width = 10)
DimPlot(dual.plasmo, group.by = "predicted.celltype", reduction = "umap") + ggtitle("percent.plasmo > 1")
```

![](integrate_dual_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

```r
#dev.off()
```

## Plasmodium RNA in samples

```r
breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

# remove cells that have plasmo.percent 0


ggplot(dual.plasmo@meta.data, aes(x = percent.plasmo)) +
  geom_histogram(bins = 70, fill = "white", colour = "blue") +
  facet_wrap(~ factor(Group, levels = c("control", "sporozoite_bystander", "sporozoite_infected",
                      "RBC_control", "RBC_bystander", "RBC_infected"))) +
  theme_bw() +
  #scale_y_log10(minor_breaks = c(5, 50, 500, 5000)) +
  theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "gray75"),
      panel.grid.minor = element_line(colour = "gray"),
      panel.border = element_blank(),
      panel.background = element_blank()) +
      theme(axis.text = element_text(size = 11)) +
      theme(strip.background =element_rect(fill="aliceblue", colour = "white"),
            strip.text.x = element_text(size = 11))
```

![](integrate_dual_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

```r
ggplot(dual.TBremoved@meta.data, aes(x = percent.plasmo)) +
  geom_histogram(bins = 70, fill = "white", colour = "blue") +
  facet_wrap(~ factor(Group, levels = c("control", "sporozoite_bystander", "sporozoite_infected",
                      "RBC_control", "RBC_bystander", "RBC_infected"))) +
  theme_bw() +
  scale_y_log10(minor_breaks = c(5, 50, 500, 5000)) +
  theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "gray75"),
      panel.grid.minor = element_line(colour = "gray"),
      panel.border = element_blank(),
      panel.background = element_blank()) +
      theme(axis.text = element_text(size = 11)) +
      theme(strip.background =element_rect(fill="aliceblue", colour = "white"),
            strip.text.x = element_text(size = 11))
```

![](integrate_dual_files/figure-html/unnamed-chunk-33-2.png)<!-- -->


```r
plasmo.percent.table <- dual.TBremoved@meta.data %>% 
  tibble::rownames_to_column("Cell") %>% 
  select(Cell, Group, percent.plasmo, percent.mt) %>% 
  group_by(Group) %>% 
  summarise(plasmo_median = median(percent.plasmo, na.rm = T),
            plasmo_max = max(percent.plasmo, na.rm = T),
            plasmo_min = min(percent.plasmo, na.rm = T),
            mt_median = median(percent.mt, na.rm = T),
            number_of_cells = n())

DT::datatable(plasmo.percent.table)
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7302bcbf94eaf1b49686" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7302bcbf94eaf1b49686">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["RBC_bystander","RBC_control","RBC_infected","control","sporozoite_bystander","sporozoite_infected"],[0,0,0,0,0,0],[97.4459724950884,32.4724172517553,93.574297188755,36.2407467009977,32.199151257957,1.91997129014893],[0,0,0,0,0,0],[5.30428387571245,5.75250836120401,6.21573153735122,6.37994873255483,5.7527359835327,4.85871812543074],[1763,3925,2074,4773,2606,49]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Group<\/th>\n      <th>plasmo_median<\/th>\n      <th>plasmo_max<\/th>\n      <th>plasmo_min<\/th>\n      <th>mt_median<\/th>\n      <th>number_of_cells<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false},"selection":{"mode":"multiple","selected":null,"target":"row","selectable":null}},"evals":[],"jsHooks":[]}</script>
```

## Features in clusters

```r
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
lib_name = "integrated"
# percent.mt on umap

meta <- dual.TBremoved@meta.data %>% 
  rownames_to_column("Cell") %>% 
  select("Cell", "nFeature_RNA", "percent.mt", "nCount_RNA", "percent.plasmo")

umap_embed <- dual.TBremoved@reductions[["umap"]]@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column("Cell") %>%
  left_join(., meta) %>% 
  mutate(lib = lib_name)

# mt.percent on umap
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, max(umap_embed$percent.mt)))
ggplot(umap_embed, aes(UMAP_1, UMAP_2, colour = percent.mt)) +
  geom_point(size = 0.3, alpha = 0.8) +
  facet_wrap( ~ lib, ncol = 2, scales = "free") +
  sc +
  theme_bw()
```

![](integrate_dual_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

```r
# nF on umap
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(umap_embed$nFeature_RNA), max(umap_embed$nFeature_RNA)))
ggplot(umap_embed, aes(UMAP_1, UMAP_2, colour = nFeature_RNA)) +
  geom_point(size = 0.3, alpha = 0.8) +
  facet_wrap( ~ lib, ncol = 2, scales = "free") +
  sc +
  theme_bw()
```

![](integrate_dual_files/figure-html/unnamed-chunk-35-2.png)<!-- -->

```r
# nCount on umap
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(umap_embed$nCount_RNA), max(umap_embed$nCount_RNA)))
ggplot(umap_embed, aes(UMAP_1, UMAP_2, colour = nCount_RNA)) +
  geom_point(size = 0.3, alpha = 0.8) +
  facet_wrap( ~ lib, ncol = 2, scales = "free") +
  sc +
  theme_bw()
```

![](integrate_dual_files/figure-html/unnamed-chunk-35-3.png)<!-- -->

```r
# parasite on umap
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(umap_embed$percent.plasmo), max(umap_embed$percent.plasmo)))
ggplot(umap_embed, aes(UMAP_1, UMAP_2, colour = percent.plasmo)) +
  geom_point(size = 0.3, alpha = 0.8) +
  facet_wrap( ~ lib, ncol = 2, scales = "free") +
  sc +
  theme_bw()
```

![](integrate_dual_files/figure-html/unnamed-chunk-35-4.png)<!-- -->


```r
#saveRDS(dual.TBremoved, "dual.integrated.TBremoved.cellbender_plasmo.bg.removed_reintegration.rds")
```


<!-- # Which Plasmodium reads were found in Ctrl? {.tabset} -->

<!-- ## Which Plasmo genes are found in ctrl and stim conditions? -->
<!-- ```{r} -->
<!-- # Collect Plasmo genes in all the conditions -->
<!-- plasmo.genes.id <- grep(pattern = "gene-PF3D7", rownames(dual.TBremoved)) -->
<!-- plasmo.genes <- dual.TBremoved[plasmo.genes.id,] -->

<!-- plasmo.genes.cellcount<- data.frame() -->
<!-- a = 1 -->
<!-- for(i in 1:nrow(plasmo.genes)) -->
<!-- { -->
<!--   for(j in 1:length(unique(plasmo.genes@meta.data$Group))) -->
<!--   { -->
<!--     print(c(i,j)) -->
<!--     gene <- rownames(plasmo.genes[i]) -->
<!--     plasmo.genes.cellcount[i,1] <- gene -->
<!--     #plasmo.genes.cellcount[a,2] <- unique(plasmo.genes@meta.data$Group)[j] -->
<!--     choose.cells.in.group <- plasmo.genes@meta.data %>%  -->
<!--       filter(Group == unique(plasmo.genes@meta.data$Group)[j]) %>%  -->
<!--       tibble::rownames_to_column("Cell") %>%  -->
<!--       pull(Cell) -->
<!--     plasmo.genes.cellcount[i,unique(plasmo.genes@meta.data$Group)[j]] <- sum(plasmo.genes[['RNA']]@counts[gene,choose.cells.in.group] > 0) -->
<!--     #a = a+1 -->
<!--   }  -->
<!-- } -->
<!-- colnames(plasmo.genes.cellcount)[1] <- "gene" -->
<!-- saveRDS(plasmo.genes.cellcount, "plasmo.genes.cellcount.rds") -->

<!-- plasmo.genes.cellcount.long <- plasmo.genes.cellcount %>%  -->
<!--   pivot_longer( -->
<!--     cols = RBC_bystander:sporozoite_bystander, -->
<!--     names_to = "Group", -->
<!--     values_to = "Number_of_cells" -->
<!--   ) -->

<!-- ``` -->

```r
# plot plasmo cells in human analysis
host <- readRDS("robjects/immune.combined.TBremoved.rds")

# get cell names from plasmo cluster
plasmo.cluster.cells <- dual.TBremoved@meta.data %>% 
  rownames_to_column("Cell") %>% 
  filter(predicted.celltype == "Plasmo") # 539 cells

# make meta data column in host to say yes plasmo or no: # 513 cells
host@meta.data <- host@meta.data %>% 
  rownames_to_column("Cell") %>% 
  mutate(plasmo_cluster = case_when(
    Cell %in% plasmo.cluster.cells$Cell ~ "yes",
    .default = "no"
  )) %>% 
  column_to_rownames("Cell")

# dimplot umap of host with new column
DimPlot(host, reduction = "umap", group.by = "plasmo_cluster")
```

![](integrate_dual_files/figure-html/unnamed-chunk-37-1.png)<!-- -->

```r
# group_cols = c('control' = '#e6af2e', 'RBC_bystander' = '#f7a78d', 'RBC_control' = '#ef4f1a', 'RBC_infected' = '#6b0504', 'sporozoite_bystander' = '#c2fff0', 'sporozoite_infected' = '#00a39b')

group_cols = c('yes' = '#6b0504', 'no' = '#e6af2e')

#pdf("dimplot_plasmo.cluster_in_host_cellbender_plasmo.bg.removed_reintegration.pdf", width = 5, height = 5)
DimPlot(host, group.by = "plasmo_cluster",reduction = "umap", label = "F", cols = group_cols) + ggtitle("plasmo cluster")
```

![](integrate_dual_files/figure-html/unnamed-chunk-37-2.png)<!-- -->

```r
#dev.off()
```

