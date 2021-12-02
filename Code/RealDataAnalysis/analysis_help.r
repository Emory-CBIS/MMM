# ------------------------
# This script includes useful function for real data analysis
# -------------------------

# function to plot probability of latent structural connectivity
PlotPiAlpha <- function(pi.a, G, M, n=2, figure = TRUE){
  # Args:
  #   pi.a: estimated latent structural connectivity matrix with latent structural state as 1
  #   G: number of groups
  #   n: integer indicating the number of decimal places
  # Return: 
  #   heatmap figure:  (G*M) * M
  #   combined structural connectivity matrix
  mat <- NULL
  for(g in 1:G){
    mat <- rbind(mat, pi.a[g,,]) 
  }
    group <- rep(c(1:G),rep(M,G)) #split variable Group
    mycol <- colorRamp2( c(0, 1), c("white", "red")) #color bar
    ht <- Heatmap(mat, 
            rect_gp = gpar(col = "white"), 
            cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(round(mat[i, j],n), x = x, y = y)
            },
            col = mycol,
            split = paste0("Group", group), 
            column_title = "Anatomical Prior Probablity", 
            name = "value", # legend value
            cluster_rows = FALSE, 
            cluster_columns = FALSE)
    if (figure == TRUE) {
      draw(ht)
    }
  return(mat)
}

# function to plot probability of latent functional connectivity
PlotPiF <- function(p_1, p0, p1, G, M, n=3, combine=TRUE, figure=TRUE){
  # Args:
  #   p_1: estimated latent functional connectivity matrix with latent functional state as -1
  #   p0: estimated latent functional connectivity matrix with latent functional state as 0
  #   p1: estimated latent functional connectivity matrix with latent functional state as 1
  #   G: number of group
  #   M: number of module
  #   n: integer indicating the number of decimal places
  #   combine: combine all figures or not
  #   figure: plot figure or not
  # Returns:
  #   heatmap figure:  (G*M) * M
  #   combined functional connectivity matrix
  D <- dim(p1)
  mat <- matrix(0, nrow=D[1]*D[2], ncol=D[3]*D[4]*3)
  for(g in 1:G){
    for(j in 1:2){
      mat[((g-1)*D[2]+1):(g*D[2]),((j-1)*3*D[2]+1):(j*3*D[2])] <- cbind(p_1[g,,,j], p0[g,,,j], p1[g,,,j])
    }
  }
  group <- rep(c(1:G),rep(M,G)) #split variable Group
  mycol <- colorRamp2( c(0, 1), c("white", "red")) #color bar
  df = data.frame(F = c(rep("-1", M), rep("0", M), rep("1", M)))
  ha1 = HeatmapAnnotation(df = df, col = list(F = c("-1" =  "blue", "0" = "grey", "1" = "red"))) # heatmap annotaion
  ha2 = HeatmapAnnotation(df = df, col = list(F = c("-1" =  "blue", "0" = "grey", "1" = "red")), show_legend = FALSE) # heatmap annotaio
  ht1 <- Heatmap(mat[,1:(M*3)], 
          rect_gp = gpar(col = "white"), 
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(round(mat[,1:(M*3)][i, j],n), x = x, y = y)
          },
          col = mycol,
          split = paste0("Group", group), 
          column_title = "Functional Prior Probablity \n Anatomical = 0", 
          name = "value", # legend value
          top_annotation = ha1,
          cluster_rows = FALSE, 
          cluster_columns = FALSE)
  ht2 <- Heatmap(mat[,(M*3+1):(M*6)], 
                 rect_gp = gpar(col = "white"), 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(round(mat[,(M*3+1):(M*6)][i, j],n), x = x, y = y)
                 },
                 col = mycol,
                 split = paste0("G", group), 
                 column_title = "Functional Prior Probablity \n Anatomical = 1",
                 show_heatmap_legend = FALSE,
                 name = "value", # legend value
                 top_annotation = ha2,
                 cluster_rows = FALSE, 
                 cluster_columns = FALSE)
  if (combine == TRUE) {
    if (figure == TRUE){
      draw(ht1 + ht2)
    }
    return(mat)
  }
  else {
    if (figure == TRUE) {
      draw(ht1)
    }
    return(mat)
  }
}

# function to plot probability of latent structural connectivity separately for each group
PlotPiAlphaSplit<-function(pi.a, G, M, n=2, sim = FALSE, order = FALSE, figure = TRUE, s = c(1,2,3,5,6,4,7,8,9)){
  # Args:
  #   pi.a: estimated latent structural connectivity matrix with latent structural state as 1
  #   G: group number
  #   M: number of module
  #   n: integer indicating the number of decimal places
  #   sim: simulation data or real data 
  #   order: reorder DMN module or not (only used when sim = FALSE)
  #   figure: plot figure or not
  #   s: new order of moduels (only used when order = TRUE)
  # Return: 
  #   heatmap figure:  M * M
  #   structural connectivity matrix
  mat <- NULL
  mat <- rbind(mat, pi.a[G,,]) 
  
  mycol <- colorRamp2( c(0, 1), c("white", "red")) #color bar
  if (sim == FALSE) {
    if (figure == TRUE& order == FALSE) {
      ht <- Heatmap(mat, 
                    rect_gp = gpar(col = "grey"), 
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(round(mat[i, j],n), x = x, y = y, gp=gpar(col="black", fontsize=20, fontface="bold"))
                    },
                    col = mycol,
                    name = "value", # legend value
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE)
      draw(ht, show_heatmap_legend = FALSE)
      return(mat)
    }
    
    if (figure == TRUE& order == TRUE) {
      mat <- mat[s,s]
      ht <- Heatmap(mat, 
                    rect_gp = gpar(col = "grey"), 
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(round(mat[i, j],n), x = x, y = y, gp=gpar(col="black", fontsize=20, fontface="bold"))
                    },
                    col = mycol,
                    name = "value", # legend value
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE)
      draw(ht, show_heatmap_legend = FALSE)
      return(mat)
    }
    
    if (figure == FALSE& order == FALSE) {
      return(mat)
    }
    
    if (figure == FALSE& order == TRUE) {
      mat <- mat[s,s]
      return(mat)
    }
  } else {
    if (fig == TRUE) {
      ht <- Heatmap(mat, 
                    rect_gp = gpar(col = "grey"), 
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(round(mat[i, j], n), x = x, y = y, gp=gpar(col="black", fontsize=20, fontface="bold"))
                    },
                    col = mycol,
                    name = "value", # legend value
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE)
      draw(ht, show_heatmap_legend = FALSE)
      return(mat)
    }
    if (fig == FALSE) {
      return(mat)
    }
  }
}

# function to plot probability of latent functional connectivity separately for each group
PlotPiFSplit<-function(p_1, p0, p1, G, M, n=3,  sim = FALSE, order = TRUE, figure = TRUE, combine = FALSE, fig = 1, s = c(1,2,3,5,6,4,7,8,9)){
  # Args:
  #  p_1: estimated latent functional connectivity matrix with latent functional state as -1
  #  p0: estimated latent functional connectivity matrix with latent functional state as 0
  #  p1: estimated latent functional connectivity matrix with latent functional state as 1
  #  G: group number
  #  M: number of module
  #  n: integer indicating the number of decimal places
  #  sim: simulation data or real data 
  #  order: reorder DMN module or not (only used when sim = FALSE)
  #  figure: plot figure or not
  #  combine: draw combined figure or not
  #  fig: if combine is FALSE, which one to draw: 1 -> structural connectivity 0; 2 -> structural connectivity 1
  #  s: new order of moduels (only used when order = TRUE)
  # Return: 
  #   heatmap figure:  M * (3*M)
  #   functional connectivity matrix
  D <- dim(p1)
  mat <- matrix(0, nrow=D[2], ncol=D[3]*D[4]*3)
  
  g <- 1  
  for(j in 1:2){
    mat[((g-1)*D[2]+1):(g*D[2]),((j-1)*3*D[2]+1):(j*3*D[2])] <- cbind(p_1[G,,,j], p0[G,,,j], p1[G,,,j])
  }

  mycol <- colorRamp2( c(0, 1), c("white", "red")) #color bar
  df = data.frame(F = c(rep("-1", M), rep("0", M), rep("1", M)))
  ha = HeatmapAnnotation(df = df, col = list(F = c("-1" =  "blue", "0" = "grey", "1" = "red"))) # heatmap annotaion
  if (sim == FALSE) {
    if (order == FALSE) {
      ht1 <- Heatmap(mat[,1:(M*3)], 
                     rect_gp = gpar(col = "grey"), 
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(round(mat[,1:(M*3)][i, j],n), x = x, y = y,gp=gpar(col="black", fontsize=15, fontface="bold"))
                     },
                     col = mycol,
                     name = "value", # legend value
                     show_heatmap_legend = FALSE,
                     top_annotation = ha,
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE)
      ht2 <- Heatmap(mat[,(M*3+1):(M*6)], 
                     rect_gp = gpar(col = "grey"), 
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(signif(mat[,(M*3+1):(M*6)][i, j],2), x = x, y = y,gp=gpar(col="black", fontsize=15, fontface="bold"))
                     },
                     col = mycol,
                     show_heatmap_legend = FALSE,
                     name = "value", # legend value
                     top_annotation = ha,
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE)
    }
    else {
      #s <- order(o)
      mat11 <- mat[,1:M][s,s]
      mat12 <- mat[,(M+1):(2*M)][s,s]
      mat13 <- mat[,(2*M+1):(3*M)][s,s]
      mat1 <- cbind(mat11, mat12, mat13)
      ht1 <- Heatmap(mat1, 
                     rect_gp = gpar(col = "grey"), 
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(round(mat1[i, j], n), x = x, y = y,gp=gpar(col="black", fontsize=15, fontface="bold"))
                     },
                     col = mycol,
                     name = "value", # legend value
                     show_heatmap_legend = FALSE,
                     top_annotation = ha,
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE)
      
      mat21 <- mat[,(3*M+1):(4*M)][s,s]
      mat22 <- mat[,(4*M+1):(5*M)][s,s]
      mat23 <- mat[,(5*M+1):(6*M)][s,s]
      mat2 <- cbind(mat21, mat22, mat23)
      ht2 <- Heatmap(mat2, 
                     rect_gp = gpar(col = "grey"), 
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(round(mat2[i, j], n), x = x, y = y,gp=gpar(col="black", fontsize=15, fontface="bold"))
                     },
                     col = mycol,
                     show_heatmap_legend = FALSE,
                     name = "value", # legend value
                     top_annotation = ha,
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE)
      mat <- cbind(mat1, mat2)
    }
  } else {
    ht1 <- Heatmap(mat[,1:(M*3)], 
                   rect_gp = gpar(col = "grey"), 
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(round(mat[,1:(M*3)][i, j], n), x = x, y = y,gp=gpar(col="black", fontsize=15, fontface="bold"))
                   },
                   col = mycol,
                   name = "value", # legend value
                   show_heatmap_legend = FALSE,
                   top_annotation = ha,
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE)
    ht2 <- Heatmap(mat[,(M*3+1):(M*6)], 
                   rect_gp = gpar(col = "grey"), 
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(round(mat[,(M*3+1):(M*6)][i, j], n), x = x, y = y,gp=gpar(col="black", fontsize=15, fontface="bold"))
                   },
                   col = mycol,
                   show_heatmap_legend = FALSE,
                   name = "value", # legend value
                   top_annotation = ha,
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE)
  }
  if (combine == TRUE) {
    if (figure == TRUE){
      draw(ht1 + ht2)
    }
    return(mat)
  }
  else {
    if (fig == 1 ){
      if (figure == TRUE){
        draw(ht1)
      }
      return(mat[,1:(M*3)])
    }
    if (fig == 2){
      if (figure == TRUE){
        draw(ht2)
      }
      return(mat[,(M*3+1):(M*6)])
    }
  }
}

# function to create heatmap for module-wise estimates
Heatmap9by9.signif <- function(mat, col, legend=FALSE, size=20, n=3, signif, pmat){
  # Args:
  #  mat: value matrix 
  #  col: color bar for heatmap
  #  legend: plot legend or not
  #  size: font size
  #  n: integer indicating the number of decimal places  
  #  signif: significance level matrix
  #  pmat: p-value matrix
  # Return: 
  #   figure
  #   matrix
  pmat2 <- pmat
  pmat2[signif==0] <- 1
  pmat2 <- 1 - pmat2
  pmat2[upper.tri(pmat2)] <- NA
  ht <- Heatmap(pmat2, 
                na_col = "white", 
                rect_gp = gpar(col = "white"), 
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(signif[i, j] == 1){
                    if(mat[i,j] >=0){
                      if(i >=j){
                        grid.text(round(mat[i, j], n), x = x, y = y, gp=gpar(col="red", fontsize=size, fontface="bold"))
                      }
                    } else {
                      if(i >= j){
                        grid.text(round(mat[i, j], n), x = x, y = y, gp=gpar(col="royalblue1", fontsize=size, fontface="bold"))
                      }
                    }
                  }
                  if(signif[i, j] == 0){
                    if(i >=j){
                      grid.text(round(mat[i, j], n), x = x, y = y, gp=gpar(col="black", fontsize=size, fontface="bold"))
                    }
                  }
                },
                col = col,
                name = "value",
                cluster_rows = FALSE, 
                cluster_columns = FALSE)
  if(legend==FALSE) {
    draw(ht, show_heatmap_legend = FALSE)
  }
  else {
    draw(ht)
  }
}

# function to create circle graph of brain connectivity
circlegraph <- function(mat.beta, mat.p, filename = "../Results/Circle Graph/Anatomical Diff.pdf", alpha  = 0.05, factor = 100, color.limit = 0.2){
  # Args:
  #  mat.beta: value matrix
  #  mat.p: pvalue matrix
  #  filename: output filename
  #  alpha: statistical significance level 
  #  factor: factor parameter for line width
  #  color.limit: color parameter
  #  p0: estimated latent functional connectivity matrix with latent functional state as 0
  # Return: 
  #   figure
  label = c(0, 0, 0, 0, 0, 1, 1, 1, 1)
  cols = c(rep("turquoise",sum(label==0)), rep("khaki",sum(label==1)))
  s <- c(4, 2, 1, 3, 5, 8, 6, 7, 9)
  node.names = c( "SM", "OP Vis", "Med Vis", "Lat Vis", "Aud", "FPL", "DMN", "EC", "FPR")
  node.names = factor(node.names,levels=unique(node.names))
  
  # preprocess of mat
  p = nrow(mat.beta)
  mat.p <- round(mat.p, 5)
  mat.p[mat.p == 0] <- 0.0001
  mat.beta[mat.p > alpha] = 0
  mat.p[mat.p > alpha] = 0
  mat.p <- mat.p[s, s]
  mat.beta <- mat.beta[s, s]
  
  #circos.clear()
  pdf(file = filename, width = 6, height = 6)
  colormap = colorRamp2(c(-color.limit,0,color.limit),c("blue","white","red"))
  circos.par("track.height" = 0.2, start.degree= 10)
  circos.initialize(factors = node.names,xlim=c(0,1))
  circos.track(bg.col=cols,ylim=c(0,1),
               panel.fun=function(x,y){
                 circos.text(0.5,CELL_META$ylim[1]+uy(5,'mm'),CELL_META$sector.index,
                             facing='inside',niceFacing = T,
                             cex=1)
               })
  for(i in 2:p){
    for(j in 1:(i-1)){
      if(mat.p[i,j]!=0){
        circos.link(node.names[j], 0.5, node.names[i], 0.5, col = colormap(mat.beta[i,j]), lwd= factor * log10(mat.p[i,j]))
      }
    }
  }
  dev.off()
}
