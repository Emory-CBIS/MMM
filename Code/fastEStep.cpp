// C++ implementation of raddlog, used to deal with the numerical issue of logSum
#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double Raddlog_C(double x, double y) {
  double a = x;
  double b = y;
  double result=0;  
  bool idx1 = a > (b + 200);
  if(idx1){
    result = a;
  }
  bool idx2 = b > (a + 200);
  if(idx2){
    result = b;
  }
  bool idx0= !(idx1 || idx2);
  if(idx0){
    result = a + log1p(exp(b - a));
  }
  return result;
}

// [[Rcpp::export]]
// C++ implementation of Estep and Mstep
List Update_group_module_log(NumericVector II, NumericMatrix III,NumericVector I_G, arma::cube pi_a, arma::cube pi_a0, List pi_f1, List pi_f0, List pi_f_1, 
                             NumericMatrix Node_Mod, NumericMatrix D, NumericMatrix B, NumericVector Rho, NumericVector X2, NumericVector Xi2,NumericMatrix Gamma, 
                             NumericMatrix Mu, NumericMatrix Sigma2){
  NumericVector dim = B.attr("dim");
  int G = II.size() - 1; 
  arma::cube P_hat(G, dim[1], 6); // 2 is the group number 
  arma::cube p_hat(G, dim[1], 6);
  arma::cube p_hat1(G, dim[1], 6);
  P_hat.fill(0);
  p_hat.fill(0);
  p_hat1.fill(0);
  double b=0;
  double a=0;
  for(int n=0; n<dim[1]; n++){
    int s = Node_Mod(n,3)-1;
    int t = Node_Mod(n,4)-1;
    for(int j=0; j<2; j++){
      for(int k=0; k<3; k++){
        int l = 3*j + k;
        for(int g=0; g < G; g++){
          arma::cube pi_f1_c = as<arma::cube>(pi_f1[g]);
          arma::cube pi_f0_c = as<arma::cube>(pi_f0[g]);
          arma::cube pi_f_1_c = as<arma::cube>(pi_f_1[g]);
          a = (pi_a(g,s,t)) * (1-j) + (pi_a0(g,s,t)) * j + (pi_f1_c(s,t,j))*(k-1)*(k+1)*k/6 + pi_f0_c(s,t,j)*k*(k-2)*(-1) + (pi_f_1_c(s,t,j))*(k+1)*(k-1)*(k-2)/2;
          for(int i= II[g]; i<II[g+1]; i++){
            if(D(i,n)==0){
              b = std::log(Rho[j])+std::log(R::dnorm(B(i,n),Mu(j,k),sqrt(Sigma2(j,k)),0));
              P_hat(g,n,l) =  P_hat(g,n,l) + b;
            }
            else{
              b = std::log(1-Rho[j]) + std::log(R::dnorm(B(i,n),Mu(j,k),sqrt(Sigma2(j,k)),0)) + std::log(R::dnorm(D(i,n),X2(j,0),sqrt(Xi2(j,0)),0)*Gamma(j,0) + R::dnorm(D(i,n),X2(j,1),sqrt(Xi2(j,1)),0)*Gamma(j,1) + R::dnorm(D(i,n),X2(j,2),sqrt(Xi2(j,2)),0)*Gamma(j,2));
              P_hat(g,n,l) = P_hat(g,n,l) + b;
            }
          }
          P_hat(g,n,l) = P_hat(g,n,l) + a;
        }
      }
    }
  }

  
  // probability calculated by Raddlog 
  for(int g=0; g<G; g++){
    for(int n=0; n<dim[1]; n++){
      for(int l=0; l<6; l++){
        double logt=P_hat(g,n,0);
        for(int z=1; z<6; z++){
          logt = Raddlog_C(logt,P_hat(g,n,z));
        }
        p_hat(g,n,l) = exp(P_hat(g,n,l) - logt);
      }
    }
  }
  
  arma::cube w(dim[0],dim[1],6);
  for(int n=0; n<dim[1]; n++){
    for(int g=0; g < G; g++){
      for(int i= II[g]; i<II[g+1]; i++){
        for(int j=0; j<2; j++){
          for(int l=0; l<3; l++){
            int o = 3*j + l;
            if(j==0){
              w(i,n,o) = (p_hat(g,n,0) + p_hat(g,n,1) + p_hat(g,n,2))* Gamma(j,l)*R::dnorm(D(i,n),X2(j,l),sqrt(Xi2(j,l)),0)/(R::dnorm(D(i,n),X2(j,0),sqrt(Xi2(j,0)),0)*Gamma(j,0)+ R::dnorm(D(i,n),X2(j,1),sqrt(Xi2(j,1)),0)*Gamma(j,1) + R::dnorm(D(i,n),X2(j,2),sqrt(Xi2(j,2)),0)*Gamma(j,2));
            }
            if(j==1){
              w(i,n,o) = (p_hat(g,n,3) + p_hat(g,n,4) + p_hat(g,n,5))* Gamma(j,l)*R::dnorm(D(i,n),X2(j,l),sqrt(Xi2(j,l)),0)/(R::dnorm(D(i,n),X2(j,0),sqrt(Xi2(j,0)),0)*Gamma(j,0)+ R::dnorm(D(i,n),X2(j,1),sqrt(Xi2(j,1)),0)*Gamma(j,1) + R::dnorm(D(i,n),X2(j,2),sqrt(Xi2(j,2)),0)*Gamma(j,2));
            }
          }
        }
      }
    }
  }
  
  NumericVector Rho_new(2);
  NumericVector nume(2);
  NumericVector domi(2);
  for(int j=0; j<2; j++){
    for(int n=0; n<dim[1]; n++){
      for(int g=0; g<G; g++){
        if(j==0){
          nume[j] = nume[j]+(p_hat(g,n,0)+p_hat(g,n,1)+p_hat(g,n,2))*III(g,n);
          domi[j] = domi[j]+(p_hat(g,n,0)+p_hat(g,n,1)+p_hat(g,n,2))*I_G[g];
        }
        if(j==1){
          nume[j] = nume[j]+(p_hat(g,n,3)+p_hat(g,n,4)+p_hat(g,n,5))*III(g,n);
          domi[j] = domi[j]+(p_hat(g,n,3)+p_hat(g,n,4)+p_hat(g,n,5))*I_G[g];
        }
      }
    }
    Rho_new[j] =nume[j]/domi[j];
  }
  
  return List::create(Named("P.hat")=P_hat,Named("p.hat")=p_hat,Named("w")=w,Named("Rho.new")=Rho_new);
}


