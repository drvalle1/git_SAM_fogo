#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
int whichLessDVPresence(double value, NumericVector prob) {
  int res=-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples cs's from a categorical distribution
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector cs(prob.nrow());

  for(int i=0; i<prob.nrow();i++){
    cs[i]=whichLessDVPresence(randu[i],prob(i,_));
  }
  return cs;
}

//' This function converts v's into theta's
// [[Rcpp::export]]
NumericVector convertSBtoNormal(NumericVector v) {
  NumericVector res(v.size());
  
  res[0]=v[0];
  double prod=1-v[0];
  for(int j=1; j<v.size();j++){
    res(j)=v[j]*prod;    
    prod=prod*(1-v[j]);
  }
  
  return (res);
}