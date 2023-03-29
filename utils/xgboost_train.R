# Training a multiclass classifier
# input
#' train_Data: a gene expression matrix (variable genes x cells)
#' train_labels = factor of labels
#' var.genes = features to use for learning
#' do.scale = whether we want to z-score the features
#' train.frac = fraction of cells in each ident we want to use for training
XGBoost_train = function(train_Data, train_labels = NULL,var.genes=NULL, do.scale=FALSE,scale.mean = NULL, scale.var = NULL,max.cells.per.ident = 400, train.frac = 0.6, nround = 200){
  library(xgboost)
  if (!is.factor(train_labels)){
    train_labels = factor(train_labels)
  }
  
  if (length(train_labels) != ncol(train_Data)) stop("Error: There must be as many training IDs as there are cells (columns)")
  
  # Subsetting and Scaling
  train_Data = as.matrix(train_Data[var.genes,])
  if (do.scale){
    if (is.null(scale.mean) & is.null(scale.var)){
      scale.mean = rowMeans(train_Data)
      scale.var = apply(train_Data, 1, sd)
      train_Data = t(scale(t(train_Data), center=scale.mean, scale=scale.var))
    } else {
      train_Data = t(scale(t(train_Data), center=scale.mean, scale=scale.var))
    }
  }
  
  # Training set vs. validation set
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using either ", train.frac*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training, whichever is smaller"))
  
  
  for (cc in c(1:length(levels(train_labels)))){
    i = levels(train_labels)[cc]
    
    cells.in.clust = names(train_labels)[train_labels == i];
    n = min(c(max.cells.per.ident, round(length(cells.in.clust)*train.frac)))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(cc-1,length(train.temp))); validation.label = c(validation.label, rep(cc-1, length(validation.temp)));
    # Consider upsampling
  }
  
  train_matrix <- xgb.DMatrix(data = t(train_Data[,training.set]), label=training.label)
  validation_matrix <- xgb.DMatrix(data = t(train_Data[,validation.set]), label=validation.label)
  
  numberOfClasses <- length(unique(training.label))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "eta" = 0.2,"max_depth"=6, subsample = 0.6)
  nround    <- nround # number of XGBoost rounds
  print(1)

  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = nround)
  
  print(2)
  # Predict hold-out validation set
  validation_pred <- predict(bst_model, newdata = validation_matrix)
  validation_prediction <- matrix(validation_pred, nrow = numberOfClasses,
                                  ncol=length(validation_pred)/numberOfClasses)
  
  valid_predlabels=apply(validation_prediction,2,which.max)-1
  A = table(validation.label, valid_predlabels)
  colnames(A) = levels(train_labels); rownames(A) = levels(train_labels)
  plotConfusionMatrix(A, order="Row", xlab.use = "True", ylab.use = "Predicted")
  
  to.return = list()
  to.return$bst_model = bst_model
  to.return$scale_mean = scale.mean
  to.return$scale_var = scale.var
  to.return$test_mat = A
  return(to.return)
}

#' Trains a xgboost classification model
#'
#' @param object A Seurat object to learn
#' @param training_genes A vector of training genes.
#' @param train_ident Which ident for the model to learn. 
#' @param do.scale Boolean value indicating whether to scale the data or not.
#'
#' @return An xgboost model
TrainModel = function(object, training_genes, train_ident = NULL, do.scale = TRUE, assay = "RNA", slot = "data", train.frac = 0.6){
  library(reshape2)
  train_data = slot(object@assays[[assay]], slot)[training_genes,]

  if(!is.null(train_ident)){
    Idents(object) <- train_ident
  }
  train_id = Idents(object)
  
  bst_model = XGBoost_train(train_data, train_labels = train_id, 
                            var.genes = training_genes, 
                            do.scale = do.scale)
  return(bst_model)
}


#' Builds a confusion matrix given an xgboost classification model and a Seurat object to apply the model to 
#'
#' @param train A Seurat object that was used to train the xgboost model
#' @param test A Seurat object to apply the xgboost model to
#' @param model An xgboost model
#' @param scale.by.model If TRUE, the test dataset is scaled by the model mean and variance. If FALSE, the test dataset is z-scored.
#'
#' @return A Confusion matrix
#'
#' @examples
#' adult_to_larva <- BuildConfusionMatrix(larva, adult, model = larva.model)
BuildConfusionMatrix = function(test, model, test_ident = NULL, scale = FALSE, scale.by.model = FALSE, assay = "RNA", slot = "data"){
  
  genes.use <- model$bst_model$feature_names
  train_id = factor(colnames(model$test_mat), levels=unique(colnames(model$test_mat)))
  
  test_data = as.matrix(slot(test@assays[[assay]], slot))
  
  if(!is.null(test_ident)){
    Idents(test) <- test_ident
  }
  test_id = Idents(test)
  
  
  # Use trained model to predict on test data
  numberOfClasses <- model$bst_model$params$num_class
  
  test_xgb = t(test_data[genes.use,])
  
  if(scale){
    if (scale.by.model){
      test_xgb = scale(test_xgb, center=model$scale_mean, scale = model$scale_var)
    }
    else {
      test_xgb = scale(test_xgb)
    }
  }
  
  test_xgb = xgboost::xgb.DMatrix(test_xgb)
  
  test_pred <- predict(model$bst_model, newdata = test_xgb)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses)
  
  # Find best class for each cell
  test_pred_margins = apply(test_prediction,2,max)
  test_predlabels = apply(test_prediction,2,which.max)
  names(test_predlabels) = colnames(test_data)
  test_pred_names = levels(train_id)[test_predlabels]
  names(test_pred_names) = names(test_predlabels)
  
  confusion_matrix = table(test_id, test_pred_names)
  return(confusion_matrix)
}

#' Predicts labels given an xgboost classification model and a Seurat object to apply the model to 
#'
#' @param train A Seurat object that was used to train the xgboost model
#' @param test A Seurat object to apply the xgboost model to
#' @param model An xgboost model
#' @param scale.by.model If TRUE, the test dataset is scaled by the model mean and variance. If FALSE, the test dataset is z-scored.
#'
#' @return Predicted cell labels
#'
#' @examples
#' adult_to_larva <- BuildConfusionMatrix(larva, adult, model = larva.model)
PredictLabels = function(test, model, test_ident = NULL, scale = FALSE, scale.by.model = FALSE, assay = "RNA", slot = "data"){
  genes.use <- model$bst_model$feature_names
  train_id = factor(colnames(model$test_mat), levels = unique(colnames(model$test_mat)))
  
  test_data = as.matrix(slot(test@assays[[assay]], slot))
  
  if(!is.null(test_ident)){
    Idents(test) <- test_ident
  }
  test_id = Idents(test)
  
  
  # Use trained model to predict on test data
  numberOfClasses <- model$bst_model$params$num_class
  
  test_xgb = t(test_data[genes.use,])
  
  if(scale){
    if (scale.by.model){
      test_xgb = scale(test_xgb, center=model$scale_mean, scale = model$scale_var)
    }
    else {
      test_xgb = scale(test_xgb)
    }
  }
  
  test_xgb = xgboost::xgb.DMatrix(test_xgb)
  
  test_pred <- predict(model$bst_model, newdata = test_xgb)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses)
  
  # Find best class for each cell
  test_pred_margins = apply(test_prediction,2,max)
  test_predlabels = apply(test_prediction,2,which.max)
  names(test_predlabels) = colnames(test_data)
  test_pred_names = levels(train_id)[test_predlabels]
  names(test_pred_names) = names(test_predlabels)
  
  return(test_pred_names)
}


#' Generates a diagonlized confusion matrix plot
#'
#' @param C A Confusion matrix
#' @param xlab.use X label to use
#' @param ylab.use Y label to use
#' @param title.use Main title
#' @param col.low Color to use for low values
#' @param col.high Color to use for high values
#' @param order Whether to return the order of clusters plotted on the x axis
#' @return A confusion matrix plot as well as the row order
MakePrettyConfusionMatrix = function(C, xlab.use = "Training clusters", ylab.use = "Test clusters", title.use = "Performance Confusion Matrix", col.high = "darkblue",col.low = "white", x.lab.rot=FALSE, order = FALSE){
  library(gplots)
  hc <- heatmap.2(C,hclustfun = function(x) hclust(x,method="single"))
  C1=C[hc$rowInd, hc$colInd]
  row.max = apply(C1,1,which.max)
  names(row.max) = rownames(C1)
  row.ord = names(sort(row.max))
  
  plotConfusionMatrix(C1[row.ord,], xlab.use = xlab.use, ylab.use = ylab.use, col.high = col.high, col.low = col.low, x.lab.rot = x.lab.rot)
  
  if(order){
    return(row.ord)
  }
}
