
#' Converts gene names within a Seurat object to lowercase
#'
#' @param object A Seurat object.
#' @param integration Boolean value indicating whether the object has an integrated assay.
#'
#' @return Returns a Seurat object with genes in lowercase.
LowerCase_genes = function(object, integration = FALSE){
  
  rownames(object@assays$RNA@counts)<- tolower(rownames(object@assays$RNA@counts))
  rownames(object@assays$RNA@data)<- tolower(rownames(object@assays$RNA@data))
  rownames(object@assays$RNA@scale.data)<- tolower(rownames(object@assays$RNA@scale.data))
  object@assays$RNA@var.features <- tolower(object@assays$RNA@var.features)

  if (integration){
    rownames(object@assays$integrated@counts)<- tolower(rownames(object@assays$integrated@counts))
    rownames(object@assays$integrated@data)<- tolower(rownames(object@assays$integrated@data))
    rownames(object@assays$integrated@scale.data)<- tolower(rownames(object@assays$integrated@scale.data))
    object@assays$integrated@var.features <- tolower(object@assays$integrated@var.features)
  }
  
  return(object)
}


#' Converts gene names withing a Seurat object to uppercase 
#'
#' @param object A Seurat object.
#' @param integration Boolean value indicating whether the object has an integrated assay.
#'
#' @return Returns a Seurat object with genes in uppercase.
UpperCase_genes = function(object, integration = FALSE){
  
  rownames(object@assays$RNA@counts)<- toupper(rownames(object@assays$RNA@counts))
  rownames(object@assays$RNA@data)<- toupper(rownames(object@assays$RNA@data))
  rownames(object@assays$RNA@scale.data)<- toupper(rownames(object@assays$RNA@scale.data))
  object@assays$RNA@var.features <- toupper(object@assays$RNA@var.features)
  #rownames(object@assays$RNA@meta.features)<- toupper(rownames(object@assays$RNA@meta.features))
  
  if (integration){
    rownames(object@assays$integrated@counts)<- toupper(rownames(object@assays$integrated@counts))
    rownames(object@assays$integrated@data)<- toupper(rownames(object@assays$integrated@data))
    rownames(object@assays$integrated@scale.data)<- toupper(rownames(object@assays$integrated@scale.data))
    object@assays$integrated@var.features <- toupper(object@assays$integrated@var.features)
    #rownames(object@assays$integration@meta.features)<- toupper(rownames(object@assays$integration@meta.features))
  }
  
  return(object)
}


#' This function calculates the percentage of cells within a cluster that express a certain gene
#' 
#' @param A Seurat object
#' @param clusID Cluster of interest
#' @param tf Gene of interest
#' @param threshold Threshold above which a gene must be expressed to be counted as expressed within that cell (default 0)
#'
#' @return The percentage of cells within a cluster that express the gene
#'
#' @examples
#' tfPercentExpression(object = pbmc, clusID = 1, tf = "LTB")
tfPercentExpression = function(object, clusID, tf, threshold = 0){
  # Determine which cells are in each cluster
  cluster_cells <- WhichCells(object, idents = clusID)
  
  # Determine which cells express transcription factors over a given threshold
  positiveCellMat <- object@assays$RNA@data[tf,cluster_cells] > threshold
  
  # Calculate the percentage of cells that express the transcription factors
  if(class(positiveCellMat) == "logical"){
    numPositiveCells <- sum(positiveCellMat)
  }else{
  numPositiveCells <- Matrix::rowSums(positiveCellMat)
  }
  numCells <- length(cluster_cells)
  percent_express <- numPositiveCells / numCells
  return(percent_express)
}


#' Calculate the Shannon Entropy for a vector
#'
#' @param vector 
#' @param threshold A threshold set to remove background. When the vector is scaled to sum to one, any elements below the threshold will be set to zero and the vector will be scaled.
#'
#' @return The Shannon Entropy for the vector
#'
#' @examples
#' ShannonEntropy(vector = c(1, 0, 0, 1, 0))
#' ShannonEntropy(vector = c(10, .5, 0, 10, 0), threshold = .05)
ShannonEntropy = function(vector, threshold = 0){
  # Scale the vector by removing background
  vector <- vector / sum(vector)
  vector[vector<threshold] <- 0
  # Renormalize
  vector <- vector / sum(vector)
  # Calculate Shannon entropy
  Hx <- -sum(vector[vector>0] * log2(vector[vector>0]))
  return(Hx)
}

SimpsonIndex = function(table, invert = FALSE){
  table <- as.array(table)
  denom = sum(table) * (sum(table) - 1)
  numer = sum(table * (table - 1))
  index = (numer / denom)
  if (invert){
    return(1 / index)
  }
  else{
    return(index)
  }
}

RaoDiversityIndex <- function(table){
  
}


#' Calculate the Occupation Number for a vector
#'
#' @param vector 
#' @param threshold A threshold set to remove background. When the vector is scaled to sum to one, any elements below the threshold will be set to zero and the vector will be scaled.
#'
#' @return The Occupation Number for the vector
#'
#' @examples
#' OccupationNumber(vector = c(1, 0, 0, 1, 0))
#' OccupationNumber(vector = c(10, .5, 0, 10, 0), threshold = .05)
OccupationNumber = function(vector, threshold = 0){
  # Scale the vector by removing background
  vector <- vector / sum(vector)
  vector[vector<threshold] <- 0
  vector <- vector / sum(vector)
  
  return(1/sum(vector^2))
  
}


#' Finds genes that are uniquely expressed in a single cluster from a list of differentially expressed genes
#'
#' @param object A Seurat object
#' @param markers A dataframe generated by FindAllMarkers
#' @param DEgenes_to_check Number of DE genes to check for uniqueness
#' @param min_percent_expression If a gene is not expressed in a cluster above this threshold, it is not tested for uniqueness
#'
#' @return A vector of genes that are uniquely expressed in each cluster
#'
#' @examples
#' FindUniqueMarkers(object = pbmc, markers = pbmc.markers)
FindUniqueMarkers = function(object, markers, DEgenes_to_check = 30, min_percent_expression = .3){
  object <- UpperCase_genes(object)
  markers$gene <- toupper(markers$gene)
  unique_markers <- NULL
  
  for(i in levels(Idents(object))){
    low_enrichment_score <- Inf
    top_mark <- NULL
    top_genes <- head(subset(markers, cluster == i), DEgenes_to_check)$gene
    
    for(gene in top_genes){
      per_express_i <- tfPercentExpression(object, clusID = i, tf = gene, threshold = 0)
      if(per_express_i < min_percent_expression){
        next()
      }
      per_express_vec <- unlist(lapply(levels(Idents(object)), function(x) tfPercentExpression(object, clusID = x, tf=gene, threshold = 1)))
      enrichment_score <- sum(per_express_vec / per_express_i)
      if(enrichment_score < low_enrichment_score){
        low_enrichment_score = enrichment_score
        top_mark = gene
        }
      }
    unique_markers <- c(unique_markers, top_mark)
    }
    
  return(unique_markers)
}

#' This function calculates the percentage of cells within all clusters that express a certain gene
#'
#' @param object A Seurat object
#' @param gene Gene (or genes) of interest
#' @param threshold Threshold above which a gene must be expressed to be counted as expressed within that cell (default 0)
#'
#' @return A matrix of genes by expression level in each cluster
#'
#' @examples
#' ExpressionByCluster(object = pbmc, gene = "LTB")
ExpressionByCluster = function(object, gene, threshold = 0){
  idents <- levels(Idents(object))
  num_idents = length(idents)
  per_express <- lapply(idents, function(x) tfPercentExpression(object, clusID = x, tf=gene, threshold = threshold))
  per_express_mat <- matrix(unlist(per_express), ncol=num_idents, dimnames = list(gene, idents))
  return(per_express_mat)
}

#' Simplified workflow for clustering a Seurat object
#'
#' @param object A Seurat object 
#' @param nfeatures Number of features to select as top variable features
#' @param numPCs Number of principal components for clustering analysis
#' @param normalization.method Method for normalization (select either LogNormalize, CLR, or RC). See Seurat::NormalizeData for more information.
#' @param scale.factor Scale factor for cell-level normalization
#' @param selection.method How to choose top variable features (select either vst, mean.var.plot, or dispersion). See Seurat::FindVariableFeatures for more information.
#' @param scale.all Boolean indicating whether to scale all features, or just variably expressed features
#' @param cluster_resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param elbow Boolean indicating whether to show an Elbow Plot for better selection of numPCs
#'Fovea_RGC$reclustered
#' @return A clustered Seurat object with tSNE and UMAP embeddings
ClusterSeurat = function(object, nfeatures = 2000, numPCs = 20, normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst",  scale.all = FALSE,  cluster_resolution = .5, elbow = FALSE, integrate.by = NULL, k.weight = 100, do.tsne = FALSE){
  if(!is.null(integrate.by)){
    smallest_ident <- min(table(object@meta.data[,integrate.by]))
    obj.list <- SplitObject(object, split.by = integrate.by)
    rm(object)
    for (i in 1:length(obj.list)) {
      obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
      obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = selection.method, nfeatures = nfeatures, verbose = FALSE)
    }
    obj.anchors <- FindIntegrationAnchors(object.list = obj.list)
    
    if(smallest_ident < k.weight){
      print("Smallest ident less than k.weight, setting k.weight to size of smallest ident.")
      k.weight = smallest_ident
    }
    object <- IntegrateData(anchorset = obj.anchors, k.weight = k.weight)
    DefaultAssay(object) <- "integrated"
  }
  else{
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
  }
  if(scale.all){
    all.genes <- rownames(object)
    object <- ScaleData(object, features = all.genes)
  }
  else{
    object <- ScaleData(object)
  }
  
  object <- RunPCA(object)
  if(elbow){
    ElbowPlot(object, ndims = 50)
    numPCs <- readline(prompt = "Enter the number of PCs desired for clustering analysis: ")
  }
  object <- FindNeighbors(object, dims = 1:numPCs)
  object <- FindClusters(object, resolution = cluster_resolution)
  if(do.tsne){
    object <- RunTSNE(object, dims = 1:numPCs)
  }
  object <- RunUMAP(object, dims = 1:numPCs)
  
  return(object)
}


#' Diagonlizes a set of genes for a Dot Plot
#'
#' @param genes A vector of genes to plot
#' @param object A Seurat object
#' @param increasing Boolean indicating the direction of the diagonal
#'
#' @return The same vector of genes but in diagonal order
DiagonalizeGenes <- function(genes, object, increasing = FALSE){
  num_idents <- length(levels(Idents(object)))
  per_express <- lapply(1:num_idents, function(x) tfPercentExpression(object, clusID = x, tf=genes, threshold = 0))
  per_express_mat <- matrix(unlist(per_express), ncol=num_idents, dimnames = list(genes, 1:num_idents))
  per_express_mat <- per_express_mat[,as.numeric(levels(Idents(object)))]
  gene_max_cluster <- apply(per_express_mat, 1, which.max)
  order <- names(sort(gene_max_cluster))
  if(increasing){
    order <- rev(order)
  }
  return(order)
}

#' Constructs a dendrogram and creates a new column in the metadata that is factored according to the dendrogram
#'
#' @param object A Seurat object with a new column "dendro_order" that is factored by a dendrogram
DendroOrder <- function(object, nfeatures = 500){
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)
  object <- BuildClusterTree(object)
  tree_obj = object@tools$BuildClusterTree
  
  if(0 %in% levels(Idents(object))){
    left_clusts = Seurat:::GetLeftDescendants(tree_obj, length(levels(Idents(object)))+1) - 1 
    right_clusts = Seurat:::GetRightDescendants(tree_obj, length(levels(Idents(object)))+1) - 1
  }
  else{
    left_clusts = Seurat:::GetLeftDescendants(tree_obj, length(levels(Idents(object)))+1)
    right_clusts = Seurat:::GetRightDescendants(tree_obj, length(levels(Idents(object)))+1)
  }
  
  clust_order = c(left_clusts, right_clusts)
  object@meta.data$dendro_order = NaN
  object@meta.data$dendro_order = factor(Idents(object), levels = clust_order)
  return(object)
}


ConvertBarcodes <- function(object){
  barcode <- function(string){
    code <- strsplit(string, c(":"))
    barcode <- strsplit(code[[1]][2], "x")
    return (paste0(barcode[[1]][1], "-1"))
  }
  new_names <- unlist(lapply(colnames(object), barcode))
  colnames(object) <- new_names
  rownames(object@meta.data) <- new_names
  return(object)
}

IntegrateObject <- function(object, split.by, dims = 40, nfeatures = 1500, selection.method = "vst"){
  # Split object
  object.list <- SplitObject(object, split.by = split.by)
  # Normalize each dataset and find variable features
  for (i in 1:length(object.list)) {
    object.list[[i]] <- NormalizeData(object.list[[i]], verbose = FALSE)
    object.list[[i]] <- FindVariableFeatures(object.list[[i]], selection.method = selection.method, 
                                            nfeatures = nfeatures, verbose = FALSE)
  }
  # Find Integration anchors
  object.anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:dims)
  # Integrate Data
  object.integrated <- IntegrateData(anchorset = object.anchors, dims = 1:dims)
  return(object.integrated)
}


# Define functions 
ReadCountMatrix = function(filename, cell_file, gene_file){
  library(reticulate)
  np <- import("numpy")
  cell_names <- read.csv(cell_file)
  cell_names <- cell_names$X
  
  gene_names <- read.csv(gene_file)
  gene_names <- rownames(gene_names)
  py_mat <- np$load(filename, allow_pickle = TRUE)
  gene_mat <- py_mat[[1]][0:(length(cell_names)-1)]
  
  rownames(gene_mat) <- cell_names
  colnames(gene_mat) <- gene_names
  gene_mat <- t(gene_mat)
  return(gene_mat)
}


#' Merges cluster identities. Saves new cluster identities in "seurat_clusters" metadata column.
#'
#' @param object A Seurat object 
#' @param idents A list of cluster identities to merge
#' @param meta.data The metadata column to pull the idents from
MergeClusters <- function(object, idents, meta.data = "seurat_clusters", refactor = FALSE){
  Idents(object) <- meta.data
  merge_cells <- c()
  for (i in idents){
    merge_cells <- append(merge_cells, WhichCells(object, idents = i))
  }
  object@meta.data[merge_cells,]$seurat_clusters <- head(idents, 1)
  if(refactor){
    object@meta.data$seurat_clusters <- droplevels(object@meta.data$seurat_clusters)
    levels(object@meta.data$seurat_clusters) = 1:length(levels(object@meta.data$seurat_clusters))
  }
  Idents(object) <- "seurat_clusters"
  return(object)
}

#' Removes cluster identities and refactors remaining clusters in "seurat_clusters" metadata column.
#'
#' @param object A Seurat object 
#' @param idents A list of cluster identities to remove from the data
#' @param meta.data The metadata column to pull the idents from
DropClusters <- function(object, idents, meta.data = "seurat_clusters", refactor = FALSE){
  Idents(object) <- meta.data
  cells.remove <- WhichCells(object, idents= idents)
  object <- subset(object, cells = setdiff(colnames(object), cells.remove))
  
  if(refactor){
    object@meta.data$seurat_clusters <- droplevels(object@meta.data$seurat_clusters)
    levels(object@meta.data$seurat_clusters) = 1:length(levels(object@meta.data$seurat_clusters))
  }
  Idents(object) <- "seurat_clusters"
  return(object)
}

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

InvestigateClusters <- function(object, ident.1, ident.2){
  markers <- FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2)
  DotPlot(object, features = head(rownames(subset(markers, avg_log2FC > 0)), 10), assay = "RNA", idents = c(ident.1, ident.2)) + RotatedAxis()
  DotPlot(object, features = head(rownames(subset(markers, avg_log2FC < 0)), 10), assay = "RNA", idents = c(ident.1, ident.2)) + RotatedAxis()
  # return(markers)
}
