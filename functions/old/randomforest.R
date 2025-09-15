run_random_forest_by_taxonomic_rank <- function(comm.tax.list, stations, 
                                                ranks = c("Order", "Family", "Genus", "Species"), 
                                                ntree = 10000, 
                                                plot = TRUE) {
  
  rf.models <- list()
  importance.values <- list()
  
  for (rank in ranks) {
    
    # Prepare data
    rf.data <- comm.tax.list[[rank]]
    rf.data$habitat <- as.factor(stations$Type.1)
    
    # Clean names
    names(rf.data) <- make.names(names(rf.data))
    
    # Check that there's more than just the habitat column
    if (ncol(rf.data) <= 1) {
      message("Skipping rank ", rank, ": no taxa to model.")
      next
    }
    
    # Fit model
    rf.model <- randomForest(habitat ~ ., data = rf.data, importance = TRUE, ntree = ntree)
    rf.models[[rank]] <- rf.model
    
    # Importance
    imp <- as.data.frame(importance(rf.model))
    imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
    importance.values[[rank]] <- imp
    
    # Optional plot
    if (plot) {
      message("Attempting plot for ", rank, "...")
      tryCatch({
        varImpPlot(rf.model, n.var = 20, main = paste("Top 20 Taxa Defining Habitat Types at", rank, "level"), type = 1)
      }, error = function(e) {
        message("Skipping plot for ", rank, ": varImpPlot() failed with error: ", e$message)
      })
    }
  }
  
  return(list(models = rf.models, importance = importance.values))
}