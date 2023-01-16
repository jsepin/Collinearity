# R: Function for Bootstrap Collinearity Diagnostics
#' @param m Model fitted
#' @param B Number of bootstrap samples to perform. Default is 100.
#' @param returndata Logical or data frame. If TRUE, the bootstrap data frame used for plotting is returned. This data frame can then be given as input so that bootstrapping, which may be time consuming does not have to be executed again. Default is FALSE.
#' @param ncon Number of condition indices plotted
#' @param main Title for the plot
#' @param xdi Amount to increase xlim. Default is 0.
#' @param ydi Amount to increase ylim. Default is 0.
#' @param xlim Numeric vectors of length 2, giving the x coordinates ranges.
#' @param ylim Numeric vectors of length 2, giving the x coordinates ranges.
#' @param alpha Significance level (Type I error probability). Default is 0.05.
#' @param cex.vd The magnification to be used for the variance decomposition proportions relative to the current setting of cex. Default is 1.
#' @param cex.main The magnification to be used for main titles relative to the current setting of cex. Default is 2.
#' @param cex.sub The magnification to be used for sub-titles relative to the current setting of cex. Default is 1.
#' @param cex.prop The magnification to be used for the proportions of significant results relative to the current setting of cex. Default is 1.
#' 
#' @return Returns plot

cofingerprint <- function(m,
                    B = 100,
                    returndata = FALSE,
                    ncon = NULL,
                    main = "Collinearity Fingerprint",
                    xdi = 0,
                    ydi = 0,
                    xlim = NULL,
                    ylim = NULL,
                    alpha = 0.05,
                    cex.vd = 1,
                    cex.main = 2,
                    cex.sub = 1,
                    cex.prop  = 1
){
  # par reset
  on.exit(par(par(...,no.readonly=TRUE)))
  
  # Create Bootstraps
  va_decomp <- Collinearity::Var_decom_mat.lm(m, equilibration=TRUE)
  colnames(va_decomp) <- gsub("[`\\]", "", colnames(va_decomp) )
  
  cn_original <- max(va_decomp[,"cond_ind"])
  wald_original <- coef(m)/sqrt(diag(vcov(m)))
  
  cond_ind <- round(va_decomp[,"cond_ind"],1)
  if(!is.null(ncon)){
    va_decomp <- tail(va_decomp, ncon)
  }
  p <- length(coef(m))
  
  if(!is.list(returndata)){
    results <- data.frame(matrix(NA, nrow = 1,ncol = 5))
    colnames(results) <- c("variable", "est", "se", "cond_nu", "i")
    results <- results[-1,]
    
    for(i in 1:B){
      data <- cbind(model.frame(m)[,1,drop = FALSE], model.matrix(m))
      data = data[sample(x = 1:nrow(model.matrix(m)),size = nrow(data), replace = TRUE),]
      modelcall <- unlist(strsplit(toString(m$call), ","))
      if(modelcall[[1]]=="lm"){
        formula <- paste0("model <- ",modelcall[[1]],"(`",colnames(data)[1],"`~.-1, data = data)")
      }else if(modelcall[[1]]=="glm"){
        formula <- paste0("model <- ",modelcall[[1]],"(`",colnames(data)[1],"`~.-1 ,family=",modelcall[[3]] , ", data = data)")
      }else{
        formula <- paste0("model <- ",modelcall[[1]],"(`",colnames(data)[1],"`~. , data = data)")
      }
      res <- data.frame("variable" = NA,
                        "est" = NA,
                        "se" = NA,
                        "cond_nu" = NA,
                        "id" = i)
      
      try({
        eval(parse(text = formula))
        #recover name
        variable <- gsub("[`\\]", "", names(coef(model)))
        cn <- max(Collinearity::Var_decom_mat.lm(model, equilibration=TRUE)[,"cond_ind"])
        res <- data.frame("variable" = variable,
                          "est" = coef(model),
                          "se" = sqrt(diag(vcov(model))),
                          "cond_nu" = cn,
                          "id" = i)
      },silent = TRUE
      )
      
      #results[(i*p - p + 1):(i*p), ] <- res
      results <- rbind(results, res)
    }
    
    sign <- 1
    if(!is.null(model$negative)){
      if(model$negative){
        sign <- -1
      }
    }
    results$wald <- sign* results$est/results$se
    results$sig  <- as.integer(abs(results$wald) >= qnorm(1- alpha/2) )
  }else{
    results <- returndata
  }
  
  variables <- gsub("[`\\]", "", unique(results$variable))
  # variables <- gsub("[`\\]", "", results$variable[1:length(unique(results$variable))])#########################################
  
  
  # brighter color
  t_col <- function(color, percent = 50, name = NULL) {
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
  }
  brightred <- t_col("red", perc = 80)
  
  # Variance proportions:
  width <- max(results$cond_nu,na.rm = T)-min(results$cond_nu,na.rm = T)
  ratio <- width/(2*sum(va_decomp[,"cond_ind"]))
  radius <- ratio*va_decomp[,"cond_ind"]
  color_node <- apply(va_decomp[,3:ncol(va_decomp),drop = F], 2, function(x) { t_col("blue", perc = (1-x)*100)})
  
  # drop na
  results <- na.omit(results)
  
  # plotting frame
  if(is.null(ylim)){
    ylim <- c(min(results$wald, na.rm = TRUE),max(results$wald, na.rm = TRUE) + ydi)
  }else{
    count <- nrow(results)
    results <- results[results$wald >= ylim[1] & results$wald <= ylim[2], ]
    print(paste0("ylim: Dropping ", count-nrow(results), " out of ", count, " observations"))
    ylim <- c(ylim[1], ylim[2] + ydi)
  }
  if(is.null(xlim)){
    xlim <- c(min(results$cond_nu, na.rm = TRUE),max(results$cond_nu, na.rm = TRUE) + xdi)
  }else{
    count <- nrow(results)
    results <- results[results$cond_nu >= xlim[1] & results$cond_nu <= xlim[2], ]
    print(paste0("xlim: Dropping ", count-nrow(results), " out of ", count, " observations"))
    xlim <- c(xlim[1], xlim[2] + xdi)
  }


  # plotting
  par(mfrow = c(1, p), oma=c(5,5,10,5), mar = rep(0, 4))
  for(i in 1:length(variables)){
    prop <- results$sig[results$variable == variables[i]]
    prop <- round(sum(prop)/length(prop),2)
    if(i ==1){
      # plot frame
      plot(NA,NA, ylim = ylim, xlim = xlim)
      
      # variance decomposition proportions
      lab <- data.frame(
        "id" = 1:nrow(va_decomp),
        "lab" = va_decomp[,variables[i]],
        "lab_color" = t_col("black", perc = 100*(1-va_decomp[,variables[i]]) ))
      lab <- lab[order(lab$id, decreasing = T),]
      for(k in 1:nrow(lab)){
        lab[k,"lab"] <- paste0(paste0(rep("\n",k-1),collapse = ""),"  ",k,": ",  sprintf("%.3f", as.numeric(lab[k,"lab"]) ), recycle0 = T)
        mtext(text = lab$lab[k], side = 3,outer = , line = -2.5, padj = 1, col = lab$lab_color[k], cex = cex.vd)
      }
      
      # non-significant area
      rect(-10, qnorm(alpha/2), max(results$cond_nu, na.rm = T)*1.5, qnorm(1- alpha/2), col = brightred)
      abline(h = qnorm(1- alpha/2)*c(-1,1))
      
      # points
      with(results[results$variable==variables[i],], lines(x = cond_nu, y = wald, pch = 19,type = "p" ))
      
      # name of explanatory variable
      mtext(text = variables[i],side = 3,outer = 0, line = 0.5)
      # original Wald statistics
      lines(x = cn_original, y = wald_original[i], pch = 19, col = "red",type = "p")
      
      # proportion of significant results
      mtext(text = paste0("Prop: ",prop),side = 3,outer = , line = -1.5,cex = cex.prop)


    }else{
      # plot frame
      plot(NA,NA, ylim = ylim, xlim = xlim, yaxt = "n", ylab ="")
      # variance decomposition proportions
      lab <- data.frame(
        "id" = 1:nrow(va_decomp),
        "lab" = va_decomp[,variables[i]],
        "lab_color" = t_col("black", perc = 100*(1-va_decomp[,variables[i]]) ))
      lab <- lab[order(lab$id, decreasing = T),]
      for(k in 1:nrow(lab)){
        lab[k,"lab"] <- paste0(paste0(rep("\n",k-1),collapse = ""),"  ",k,": ",  sprintf("%.3f", as.numeric(lab[k,"lab"]) ), recycle0 = T)
        mtext(text = lab$lab[k], side = 3,outer = , line = -2.5, padj = 1, col = lab$lab_color[k], cex = cex.vd)
      }
      
      # non-significant area
      rect(-10, qnorm(alpha/2), max(results$cond_nu, na.rm = T)*1.5, qnorm(1- alpha/2), col = brightred)
      abline(h = qnorm(1- alpha/2)*c(-1,1))
      
      # points
      with(results[results$variable==variables[i],], lines(x = cond_nu, y = wald, pch = 19,type = "p" ))
      
      # name of explanatory variable
      mtext(text = variables[i],side = 3,outer = 0, line = 0.5)
      # original Wald statistics
      lines(x = cn_original, y = wald_original[i],pch = 19, col = "red",type = "p")
      
      # proportion of significant results
      mtext(text = paste0("Prop: ",prop),side = 3,outer = , line = -1.5,cex = cex.prop)
      
    }}
  # descriptions
  mtext(text = bquote("Condition Number"~kappa(bold(E)) ),side = 1,outer = 3, line = 3)
  mtext(text = bquote(bold(.(main) )),side = 3,outer = 3, line = 6,cex = cex.main)
  mtext(text = bquote("Condition Indices: "~.(paste0(sort(cond_ind, decreasing = T), collapse = ", "))),side = 3,outer = 3, line = 4,cex = cex.sub)
  mtext(text = bquote("AIC: "~.(paste0(round(AIC(m),3), collapse = ", "))~", BIC: " ~.(paste0(round(BIC(m),3), collapse = ", "))   ),side = 3,outer = 3, line = 2,cex = cex.sub)
  mtext(text = bquote("Sample size: "~.(nrow(model.matrix(m)))  ),side = 3,outer = F, line = 2,cex = cex.sub)
  mtext(text = "Wald Statistics",side = 2,outer = 1, line = 3)
  
  # return dataframe if wanted
  if(!is.list(returndata)){
    if(returndata){
      return(results)
    }
  }
  
}
