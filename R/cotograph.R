# R: Function for graph as alternative Collinearity Diagnostics
#' @param m Model fitted
#' @param voi Variable of interest
#' @param equilibrate Logical. If TRUE, the design matrix is equilibrated first. Default value is FALSE
#' @param main Title for the plot.
#' @param cex.node Text size of nodes. Default is 1.
#' @param cex.tovoi Text size of edge to voi-node. Default is 1.
#' @param col.edge.line Color of edges to voi-node. Default is blue.
#' @param col.edge.text Color of edges text to voi-node. Default is blue.
#' @param col.node.voi Color of the voi node. Default is green.
#' @param col.node.nonvoi Color of non-voi node(s). Default is lightblue.
#' @param radius_circle Size of node. Default is 0.2.
#' @param subR2 Logical. If TRUE, the R2 values of the non-voi nodes are plotted as well. Default value is FALSE.
#' @param lwdtovoi Maximum linewidth of edges. Default is 5.
#' @param mar Margin of the plot. c(2.1, 2.1,2.1, 2.1)
#' 
#' @return Returns plot

cotograph <- function(m,voi,
                      equilibrate = FALSE,
                      main = "Graph: Relation of explanatory variables\n Multivariable fitted Model",
                      cex.node = 1,
                      cex.tovoi = 1,
                      cex.main = 1,
                      col.edge.line = "blue",
                      col.edge.text = "blue",
                      col.node.voi = "green",
                      col.node.nonvoi = "lightblue",
                      radius_circle = 0.2,
                      subR2 = FALSE,
                      lwdtovoi = 5, 
                      mar = c(2.1, 2.1,2.1, 2.1)
){
  # par reset
  on.exit(par(par(...,no.readonly=TRUE)))
  
  X <- model.matrix(m)
  
  # equilibrate if needed
  if(equilibrate){
    X <- Collinearity::equilibrate_matrix(X)
    colnames(X) <- colnames(model.matrix(m))
    main <- paste0(main, " (equilibrated)")
  }
  
  if(!(gsub("[`]","",voi) %in% gsub("[`]","",colnames(X)) ) ){
    stop(paste0("voi not found! \n Try: ", paste0(gsub("[`]","",colnames(X)),collapse = ", ") ) )
  }
  
  # function for circle
  loc <- function(angle,r,cx = 0,cy = 0){
    x <- round(cos(angle)*r,12) + cx
    y <- round(sin(angle)*r,12) + cy
    out <- cbind("x"=x, "y" = y)
    return(out)
  }
  # data frame with characteristics and location
  dd <- cbind(loc(angle = seq(from = 0 - pi/2 , to = 2*pi  - pi/2, length.out = ncol(X) ), r = 1),
              "angle"= seq(from = 0 - pi/2 , to = 2*pi  - pi/2, length.out = ncol(X) ))
  dd <- dd[-nrow(dd),]
  # bring voi into the middle
  dd <- rbind(c(0,0,NA), dd)
  row.names(dd) <- gsub("[`]","",colnames(X))
  index_voi <- which(row.names(dd)==voi)
  row.names(dd) <- c(row.names(dd)[index_voi], row.names(dd)[-index_voi])
  dd <- as.data.frame(dd)
  dd$color <-  c(col.node.voi, rep(col.node.nonvoi,nrow(dd)-1 ) )
  dd$label <- gsub("[`]","",row.names(dd))
  
  # Arrows
  dd_from <- dd[,c("x", "y")]
  colnames(dd_from) <- c("x0", "y0")
  dd_from$from <- row.names(dd_from)
  row.names(dd_from) <- NULL
  dd_from <- dd_from[rep(1:nrow(dd_from), each = nrow(dd_from)),]
  dd_to <- dd[,c("x", "y")]
  colnames(dd_to) <- c("x1", "y1")
  dd_to$to <- row.names(dd_to)
  row.names(dd_to) <- NULL
  dd_to <- dd_to[rep(1:nrow(dd_to), times = nrow(dd_to)),]
  dd_arr <- cbind(dd_from,dd_to)
  dd_arr <- dd_arr[!(dd_arr$from == dd_arr$to),]# remove self-arrows
  dd_arr <- dd_arr[!(dd_arr$from==voi),]
  dd_arr$from <- gsub("[`]","",dd_arr$from)
  dd_arr$to <- gsub("[`]","",dd_arr$to)
  
  
  # Trim tips of arrow
  trimmer <- function(x0,y0,x1=0,y1 = 0,r = 1){
    xsection <- x0-x1
    ysection <- y0-y1
    if(xsection >= 0){
      xmult <- 1
    }else{
      xmult <- -1
    }
    if(ysection >= 0){
      ymult <- 1
    }else{
      ymult <- -1
    }
    alp <- atan(abs(ysection)/abs(xsection))
    if(!is.na(alp)){
      x_new <- xmult*cos(alp)*r + x1
      y_new <- ymult*sin(alp)*r + y1
    }else{
      x_new <- xmult*1*r + x1
      y_new <- ymult*1*r + y1
    }
    out <- round(cbind("x1"=x_new, "y1" = y_new),8)
    return(out)
  }
  
  # angle of the arrows for text rotation
  angle <- function(x0,y0,x1=0,y1 = 0){
    alpha <- atan((y0-y1)/(x0-x1))
    alpha <- alpha/(2*pi)  *360
    return(alpha)
  }
  
  # How well is voi described?
  data <- data.frame(X); colnames(data) <- colnames(X)
  colnames(data) <- gsub("[`]","",colnames(data))
  formula <- paste0("m_voi <- lm(data = data,`", gsub("[`]","",voi) ,"`~.-1 )")
  eval(parse(text = formula))
  
  beta <- round(summary(m_voi)$coef[,"Estimate"],2)
  wald <- round(summary(m_voi)$coef[,"t value"],2)
  dd$label[dd$label==gsub("[`]","",voi)] <- paste0(dd$label[dd$label==gsub("[`]","",voi)]," \n R2: ", round(summary(m_voi)$r.squared,3) )
  
  label <- paste0("b: ", beta, "\n t: ", wald)
  label <- cbind("from" = gsub("\\`","", names(coef(m_voi)))  , label)
  label <- data.frame(label)
  label <- cbind(label, "wald" = wald)
  label <- merge(dd_arr[dd_arr$to==gsub("[`]","",voi),], label , by = "from") 
  label$xloc <- (label$x0+label$x1)/2# x-location of text
  label$yloc <- (label$y0+label$y1)/2# y-location of text
  label$srt <- angle(x0 = label$x0,y0 = label$y0,x1=0,y1 = 0) # angle of text
  label$srt[is.na(label$srt)] <- 0
  label[,c("x1","y1")] <- t(apply(cbind(label$x0, label$y0,radius_circle ),1, function(x){trimmer(x[1],x[2],r=x[3])}))
  
  # linewidth of egde to voi (lwdtovoi  = maximal lwd)
  label$lwd <- abs(label$wald)
  label$lwd[label$lwd > lwdtovoi ] <- lwdtovoi
  
  if(subR2){
    # R2 from the explanatory variables:
    label$R2 <- NA
    for(i in 1:length(label$from)){
      vari <- label$from[i]
      formula <- paste0("m <- lm(data = data,`", vari ,"`~.-1-`",gsub("[`]","",voi),"` )")
      eval(parse(text = formula))
      label$R2[i] <- summary(m)$r.squared
      dd$label[dd$label==vari] <- paste0(dd$label[dd$label==vari]," \n R2: ", round(summary(m)$r.squared,3) )
    }
  }
  
  #sorting order with respect to wald statistics
  label <- label[match( rownames(dd)[-1], label$from ),]# order so that circle is ordered!
  
  label$from <- label$from[order(abs(label$wald), decreasing = TRUE)]
  label$lwd <- label$lwd[order(abs(label$wald), decreasing = TRUE)]
  label$label <- label$label[order(abs(label$wald), decreasing = TRUE)]
  label$wald <- label$wald[order(abs(label$wald), decreasing = TRUE)]
  
  s <- c(1,match( label$from,  rownames(dd) ))
  rownames(dd) <- rownames(dd)[s]
  dd$label <- dd$label[s]
  
  # plotting
  par(pty="s",cex = 1,mar = mar)
  plot(x = dd[,1], y = dd[,2], col = NA,
       ylim = c(-1-radius_circle, 1+radius_circle), xlim = c(-1-radius_circle, 1+radius_circle),
       type='n',axes=FALSE, ylab = "",xlab ="", xaxt = "n", yaxt ="n",
       main = main,cex.main = cex.main
  )
  for(i in 1:nrow(label)){
    with(label, arrows(x0 = x0[i], y0 = y0[i],x1 = x1[i], y1 = y1[i],
                       col = col.edge.line, code = 2, length=0.3, angle = 5, lwd = lwd[i] ) )
    with(label, text(x = xloc[i], y = yloc[i], labels = label[i], col = col.edge.text,
                     srt = srt[i], cex = cex.tovoi) )
  }
  symbols(x = dd[,"x"], y = dd[,"y"], circles = rep(radius_circle, nrow(dd)),
          add=T, inches=F, fg=c("black"), bg = dd$color )
  text(x = dd[,"x"], y = dd[,"y"], labels = dd$label, cex = cex.node )
  
}
