library(circlize)
circlegraph <- function(mat.beta, mat.p, filename = "../Results/Circle Graph/Anatomical Diff.pdf", alpha  = 0.05, factor = 100, color.limit = 0.2){
  # defaulf setup
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
