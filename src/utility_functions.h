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

#endif


