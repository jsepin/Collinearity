# R: Function for Bootstrap Collinearity Diagnostics
#' @param m Model fitted
#' @param B Number of bootstrap samples to perform. Default is 100.
#' @param returndata Logical or data frame. If TRUE, the bootstrap data frame used for plotting is returned. This data frame can then be given as input so that bootstrapping, which may be time consuming does not have to be executed again. Default is FALSE.
#' @param ncon Number of condition indices plotted
#' @param main Title for the plot
#' @param xdi Horizontal position of the variance decomposition proportions. Default is 0.
#' @param ydi Amount to increase ylim. Default is 0.
#' @param alpha Significance level (Type I error probability). Default is 0.05.
#' @param cex.vd The magnification to be used for the variance decomposition proportions relative to the current setting of cex. Default is 1.
#' @param cex.main The magnification to be used for main titles relative to the current setting of cex. Default is 2.
#' @param cex.sub The magnification to be used for sub-titles relative to the current setting of cex. Default is 1.
#' @param cex.prop The magnification to be used for the proportions of significant results relative to the current setting of cex. Default is 1.
#' 
#' @return Returns plot

lmcopri <- function(m,
                    B = 100,
                    returndata = FALSE,
                    ncon = NULL,
                    main = "Collinearity Fingerprint",
                    xdi = 0,
                    ydi = 0,
                    alpha = 0.05,
                    cex.vd = 1,
                    cex.main = 2,
                    cex.sub = 1,
                    cex.prop  = 1
){
  
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
    results <- data.frame(matrix(NA, nrow = B*p,ncol = 5))
    colnames(results) <- c("variable", "est", "se", "cond_nu", "i")
    
    for(i in 1:B){
      data <- cbind(model.frame(m)[,1,drop = FALSE], model.matrix(m))
      data = data[sample(x = 1:nrow(model.matrix(m)),size = nrow(data), replace = TRUE),]
      modelcall <- unlist(strsplit(toString(m$call), ","))
      if(modelcall[[1]]=="lm"){
        formula <- paste0("model <- ",modelcall[[1]],"(`",colnames(data)[1],"`~.-1, data = data)")
      }else{
        formula <- paste0("model <- ",modelcall[[1]],"(`",colnames(data)[1],"`~. , data = data)")
      }
      eval(parse(text = formula))
      
      #recover name
      variable <- gsub("[`\\]", "", names(coef(model)))
      
      cn <- max(Collinearity::Var_decom_mat.lm(model, equilibration=TRUE)[,"cond_ind"])
      res <- data.frame("variable" = variable,
                        "est" = coef(model),
                        "se" = sqrt(diag(vcov(model))),
                        "cond_nu" = cn,
                        "id" = i)
      results[(i*p - p + 1):(i*p), ] <- res
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
  
  variables <- gsub("[`\\]", "", results$variable[1:length(unique(results$variable))])
  
  
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
  width <- max(results$cond_nu)-min(results$cond_nu)
  ratio <- width/(2*sum(va_decomp[,"cond_ind"]))
  radius <- ratio*va_decomp[,"cond_ind"]
  xcoor <- min(results$cond_nu) + cumsum(2*ratio*va_decomp[,"cond_ind"]) - ratio*va_decomp[,"cond_ind"]
  ycoor <- rep(max(results$wald) + max(radius), length(xcoor) )
  color_node <- apply(va_decomp[,3:ncol(va_decomp)], 2, function(x) { t_col("blue", perc = (1-x)*100)})
  
  # plotting
  par(mfrow = c(1, p), oma=c(5,5,10,5), mar = rep(0, 4))
  for(i in 1:length(variables)){
    prop <- results$sig[results$variable == variables[i]]
    prop <- round(sum(prop)/length(prop),2)
    if(i ==1){
      plot(NA,NA, ylim = c(min(results$wald),max(results$wald)+ydi), xlim = range(results$cond_nu))
      lab <- data.frame(
        "id" = 1:nrow(va_decomp),
        "lab" = va_decomp[,variables[i]],
        "lab_color" = t_col("black", perc = 100*(1-va_decomp[,variables[i]]) ))
      
      lab <- lab[order(lab$id, decreasing = T),]
      
      for(k in 1:nrow(lab)){
        lab[k,"lab"] <- paste0(paste0(rep("\n",k-1),collapse = ""),"  ",k,": ",  sprintf("%.3f", as.numeric(lab[k,"lab"]) ), recycle0 = T)
      }
      text(x = mean(xcoor)+xdi,y = max(results$wald)+ydi,labels = lab$lab, cex = cex.vd, pos = 1,col = lab$lab_color )
      
      rect(-10, qnorm(alpha/2), max(results$cond_nu,na.rm = T)*1.5, qnorm(1- alpha/2), col = brightred)
      with(results[results$variable==variables[i],],lines(x = cond_nu, y = wald, pch = 19,type = "p" ))
      
      abline(h = qnorm(1- alpha/2)*c(-1,1))
      mtext(text = variables[i],side = 3,outer = 0, line = 0.5)
      lines(x = cn_original, y = wald_original[i],pch = 19, col = "red",type = "p")
      mtext(text = paste0("Prop: ",prop),side = 3,outer = , line = -1.5,cex = cex.prop)
      
    }else{
      plot(NA,NA, ylim = c(min(results$wald),max(results$wald)+ydi), xlim = range(results$cond_nu),yaxt = "n", ylab ="")
      lab <- data.frame(
        "id" = 1:nrow(va_decomp),
        "lab" = va_decomp[,variables[i]],
        "lab_color" = t_col("black", perc = 100*(1-va_decomp[,variables[i]]) ))
      lab <- lab[order(lab$id, decreasing = T),]
      
      for(k in 1:nrow(lab)){
        lab[k,"lab"] <- paste0(paste0(rep("\n",k-1),collapse = ""),"  ",k,": ",  sprintf("%.3f", as.numeric(lab[k,"lab"]) ), recycle0 = T)
      }
      text(x = mean(xcoor)+xdi,y = max(results$wald)+ydi,labels =lab$lab, cex = cex.vd, pos = 1,col = lab$lab_color )
      
      rect(-10, qnorm(alpha/2), max(results$cond_nu,na.rm = T)*2, qnorm(1- alpha/2), col = brightred)
      with(results[results$variable==variables[i],],lines(x = cond_nu, y = wald, pch = 19,type = "p"
      ))
      abline(h = qnorm(1- alpha/2)*c(-1,1))
      mtext(text = variables[i],side = 3,outer = 0, line = 0.5)
      mtext(text = paste0("Prop: ",prop),side = 3,outer = , line = -1.5,cex = cex.prop)
      
      lines(x = cn_original, y = wald_original[i],pch = 19, col = "red",type = "p")
    }}
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
