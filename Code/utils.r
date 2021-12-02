# ------------------------
# This file includes two functions are useful in real data analysis
# ------------------------


# ------------------------
# Transform mat to vec
# ------------------------
mat2vec <- function(Mat)
{
  return(Mat[lower.tri(Mat)])
}

# ------------------------
# Transform mat to vec
# ------------------------
NodeMod <- function(N,Q,M_st){
  # Args:
  #   N: total number of connections
  #   Q: number of nodes
  #   M_st: a vector of module index
  Node_Mod = data.frame(Connection=c(1:N))
  a<-NULL
  b<-NULL
  for(i in 1:Q){
    a = c(a,rep(i,Q-i));
    b = c(b, c((i+1):Q));  
  }
  b = b[1:N]
  Node_Mod$Node1 = a;
  Node_Mod$Node2 = b;
  rm(a);rm(b)
  for(i in 1:N){
    node1<-Node_Mod$Node1[i]
    node2<-Node_Mod$Node2[i]
    Node_Mod$Mod1[i]<-M_st[node1]
    Node_Mod$Mod2[i]<-M_st[node2]
  }
  return(Node_Mod)
}