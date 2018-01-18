

# Motivation analysis
# Section 3.4 The presence of correlated metrics in the set of metrics selected by a feature selection technique

# Import packages
library(DefectData)
library(Hmisc)
library(dendextend)
# List of functions

# VarClus (Variable Clustering - see section 2.1)
applyVarClus <- function(dataset, sw.metrics, VarClus.threshold) {
  print('Variable Clustering : START')
  output <- {
    
  }
  correlated <- {
  }
  # Apply variable clustering
  vc <-
    varclus(~ .,
            similarity = 'spearman',
            data = dataset[, sw.metrics],
            trans = "abs")
  
  # Apply cutoff threshold
  var.clustered <-
    cutree(vc$hclust, h = (1 - VarClus.threshold))
  
  # Get surviving metrics
  non.correlated.metrics <-
    names(var.clustered)[var.clustered %in% names(table(var.clustered)[table(var.clustered) == 1])]
  print(paste0(
    length(non.correlated.metrics),
    ' non-correlated metrics : ',
    paste0(non.correlated.metrics, collapse = ', ')
  ))
  
  # Get correlated clusters index
  correlated.clusters.index <-
    names(table(var.clustered)[table(var.clustered) > 1])
  
  print(paste0((length(sw.metrics) - length(non.correlated.metrics)),
               ' correlated metrics from ',
               length(correlated.clusters.index),
               ' clusters'
  ))
  
  # For each cluster of correlated metrics, print out for manual selection
  cluster.count <- 1
  for (cluster.index in correlated.clusters.index) {
    print(paste0(
      'Cluster ',
      cluster.count,
      ' : ',
      paste0(names(var.clustered)[var.clustered == cluster.index], collapse = ', ')
    ))
    output <- rbind(output,
                    data.frame(
                      CID = cluster.count,
                      CorrelatedMetrics = paste0(names(var.clustered)[var.clustered == cluster.index], collapse = ',')
                    ))
    correlated <-
      c(correlated, names(var.clustered)[var.clustered == cluster.index])
    cluster.count <- cluster.count + 1
  }
  print('Variable Clustering : END')
  return (list(
    correlated.cluster = output,
    correlated.metrics = correlated
  ))
}

# VIF (Variance inflation factor - see Section 2.1)
applyVIF <- function(dataset,
                     sw.metrics,
                     defect,
                     VIF.threshold) {
  print('Variance Inflation Factor : START')
  step <- 1
  
  repeat {
    glm.model <-
      glm(as.formula(paste(
        defect, '~', paste0(sw.metrics, collapse = '+')
      )),
      data = dataset,
      family = binomial())
    VIF <- rms::vif(glm.model)
    inter.correlated.metrics <- names(VIF[VIF >= VIF.threshold])
    if (length(inter.correlated.metrics) == 0){
      print(paste0('nNon-inter-correlated metrics : ', length(sw.metrics)))
      print(paste0(sw.metrics, collapse=', '))
      break
    }
    print(paste0('Step ', step, ' nMetrics : ', length(sw.metrics)))
    print(paste0('Step ', step, ' nInter-correlated metrics : ',
                 length(inter.correlated.metrics)))
    sw.metrics <- sw.metrics[!(sw.metrics %in% inter.correlated.metrics)]
    step <- step + 1
  }
  
  print('Variance Inflation Factor : END')
  return(sw.metrics)
  
}

# Plot a hierarchical cluster view using the results of 
# the variable clustering analysis and variance inflation factor analysis
meaningfulVCVIFPlot <-
  function(dataset, indep, vcIndep, vifIndep, label) {
    
    # Modify from plot.varclus in Hmisc package
    .colorVCPlot <- function (vc,
                              metricsLabel,
                              abbrev = FALSE,
                              legend. = FALSE,
                              loc,
                              maxlen = 20,
                              labels = NULL,
                              ...)
    {
      trans <- vc$trans
      s <- c(
        hoeffding = "30 * Hoeffding D",
        spearman = switch(
          trans,
          square = expression(paste(Spearman, ~
                                      rho ^ 2)),
          abs = expression(paste(Spearman,
                                 ~
                                   abs(rho))),
          none = expression(paste(Spearman,
                                  ~
                                    rho))
        ),
        pearson = switch(
          trans,
          square = expression(paste(Pearson,
                                    ~
                                      r ^ 2)),
          abs = expression(paste(Pearson, ~ abs(r))),
          none = expression(paste(Pearson, ~
                                    r))
        ),
        bothpos = "Proportion",
        ccbothpos = "Chance-Corrected Proportion"
      )[vc$similarity]
      if ((is.expression(s) &&
           as.character(s) == "NULL") ||
          (!is.expression(s) && (is.na(s) || s == "")))
        s <- vc$similarity
      ylab <- s
      
      if (legend.)
        abbrev <- TRUE
      if (!length(labels))
        labels <- dimnames(vc$sim)[[2]]
      olabels <- labels
      if (abbrev)
        labels <- abbreviate(labels)
      if (!length(vc$hclust))
        stop("clustering was not done on similarity=\"ccbothpos\"")
      
      # Modify from
      # https://cran.r-project.org/web/packages/dendextend/vignettes/FAQ.html
      dend <- as.dendrogram(vc$hclust, hang = -1)
      metricColor <- metricsLabel$rawMetrics
      metricColor <- metricColor[order.dendrogram(dend)]
      textColor <-
        ifelse(metricColor == 1, 'darkgreen', '') # Clean metrics
      textColor[metricColor == 2] <- 'red' # Correlated metrics
      textColor[metricColor == 3] <-
        'blue' # Inter-correlated metrics
      textColor[metricColor == 4] <-
        'purple' # Correlated and Inter-correlated metrics
      labels_colors(dend) <- textColor
      
      # max pixels need to plot metrics
      maxLength <- max(unlist(lapply(metricsLabel$textMetrics, function(x) strwidth(x, font = 12, units = 'in'))))
      bottomMargin <- (5 + (4.544286 * (maxLength - 0.574)))
      par(mar = c(bottomMargin, 5, 5, 2.5)) 
      plot(
        dend,
        main = "",
        ylab = '',
        yaxt = 'n',
        labels = labels,
        ann = FALSE,
        axes = FALSE,
        hang = -1
      )
      ya <- pretty(range(1 - vc$hclust$height))
      axis(2, at = 1 - ya, labels = format(ya))
      title(ylab = ylab)
      s <- labels != olabels
      if (legend. && any(s)) {
        if (missing(loc)) {
          cat("Click mouse at upper left corner of legend\n")
          loc <- locator(1)
        }
        olabels <- ifelse(nchar(olabels) > maxlen,
                          substring(olabels,
                                    1, maxlen),
                          olabels)
        text(loc, paste(paste(labels[s], ":", olabels[s], "\n"),
                        collapse = ""), adj = 0)
      }
      abline(h = 0.3)
      legend('top', legend = c('Clean metric', 'Correlated metric', 'Inter-correlated metric', 'Correlated and Inter-correlated metric'), fill = c('darkgreen', 'red', 'blue', 'purple'), xpd = TRUE, 
             inset = c(0, (-0.1 - (0.006 * (bottomMargin - 5)))))
      invisible()
      
    }
    
    metrics <- rep(1, length(indep)) +
      as.numeric(!(indep %in% vcIndep)) + # correlated metrics
      (2 * as.numeric(!(indep %in% vifIndep))) # inter-correlated metrics
    
    dataPlot <- dataset[, indep]
    tMetrics <- data.frame(rawMetrics = metrics,
                           textMetrics = indep)
    names(dataPlot) <- tMetrics$textMetrics
    # Drop one of two variables that have perfect collinearity
    vc <- varclus(~ .,
                  similarity = 'spearman',
                  data = dataPlot,
                  trans = "abs")
    
    pdf(
      paste0(label, '.pdf'),
      width = 10,
      height = 10,
      paper = 'special'
    )
    
    .colorVCPlot(vc, tMetrics)
    abline(h = 0.3)
    dev.off()
    
}

# Main

# Init thresthold values
VarClus.threshold <- 0.7
VIF.threshold <- 5

# Load eclipse-2.0 dataset
Data <- loadData('eclipse-2.0')
# Retrieve data, name of defect-proneness and software metrics
dataset <- Data$data
defect <- Data$dep
sw.metrics <- Data$indep

# We use the implementation as provided by Weka with:
# - CfsSubsetEval attribute evaluation option
# - BestFirst search method option
# The set of metrics selected by a correlation-based feature selection technique consists of:
# - pre, FOUT_max, MLOC_sum, NOM_avg, NSM_max, NSM_sum, VG_sum
fs.metrics <- sw.metrics[c(1, 4, 8, 16, 24, 28, 32)]

correlated.metrics <-
  applyVarClus(dataset, fs.metrics, VarClus.threshold)
inter.correlated.metrics <- fs.metrics[!(fs.metrics %in% applyVIF(dataset, fs.metrics, defect, 5))]

meaningfulVCVIFPlot(dataset, fs.metrics, fs.metrics[!(fs.metrics %in% correlated.metrics$correlated.metrics)], fs.metrics[!(fs.metrics %in% inter.correlated.metrics)], '3.4_plot')
