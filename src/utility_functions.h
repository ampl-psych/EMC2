#ifndef utility_h
#define utility_h

#include <Rcpp.h>
using namespace Rcpp;

LogicalVector contains(CharacterVector sv, std::string txt) {
  LogicalVector res(sv.size());
  for (int i = 0; i < sv.size(); i ++) {
    res[i] = (sv[i] == txt);
  }
  return res;
}

NumericVector pnorm_multiple(NumericVector x){
  NumericVector out(x.size());
  for(int i = 0; i < x.size(); i++){
    out[i] = R::pnorm(x[i], 0, 1, TRUE, FALSE);
  }
  return out;
}

LogicalVector contains_multiple(CharacterVector sv, CharacterVector inputs) {
  LogicalVector res(sv.size());
  for (int i = 0; i < sv.size(); i ++) {
    int k = 0;
    for (int j = 0; j < inputs.size(); j++){
      if (sv[i] == inputs[j]){
        k++;
      }
    }
    res[i] = k > 0;
  }
  return res;
}

NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition) {
  int n=X.nrow(), k=X.ncol();
  NumericMatrix out(sum(condition),k);
  for (int i = 0, j = 0; i < n; i++) {
    if(condition[i]) {
      out(j,_) = X(i,_);
      j = j+1;
    }
  }
  return(out);
}

#endif


