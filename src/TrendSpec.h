// // TrendPhase marker
// enum class TrendPhase { Premap, Pretransform, Posttransform };
//
// inline TrendPhase parse_phase(const std::string& ph) {
//   if (ph == "premap")        return TrendPhase::Premap;
//   if (ph == "pretransform")  return TrendPhase::Pretransform;
//   if (ph == "posttransform") return TrendPhase::Posttransform;
//   Rcpp::stop("Unknown trend phase: '%s'", ph.c_str());
// }
//
// // One trend operation, as parsed from one element of R's `trend` list
// struct TrendOp {
//   TrendPhase phase;
//   std::string kernel;                 // <- now a string
//   std::string target_param;
//   Rcpp::CharacterVector input_params;
//   Rcpp::List spec;                    // full R spec if you need extra fields
//
//   // Construct directly from the R list element
//   TrendOp(const Rcpp::List& cur) {
//     using Rcpp::as;
//     phase        = parse_phase(as<std::string>(cur["phase"]));
//     kernel       = as<std::string>(cur["kernel"]);   // e.g. "deltaRule", "rw", ...
//     target_param = as<std::string>(cur["target"]);
//     input_params = cur["params"];                    // CharacterVector of names
//     spec         = cur;
//   }
// };
//
//
// struct TrendSpec {
//   std::vector<TrendOp> premap_ops;
//   std::vector<TrendOp> pretransform_ops;
//   std::vector<TrendOp> posttransform_ops;
//
//   bool has_premap() const       { return !premap_ops.empty(); }
//   bool has_pretransform() const { return !pretransform_ops.empty(); }
//   bool has_posttransform() const{ return !posttransform_ops.empty(); }
//
//   static TrendSpec from_list(const Rcpp::List& trend) {
//     TrendSpec ts;
//     for (int i = 0; i < trend.size(); ++i) {
//       Rcpp::List cur = trend[i];
//       TrendOp op(cur);  // uses the ctor above
//
//       switch (op.phase) {
//       case TrendPhase::Premap:        ts.premap_ops.push_back(op); break;
//       case TrendPhase::Pretransform:  ts.pretransform_ops.push_back(op); break;
//       case TrendPhase::Posttransform: ts.posttransform_ops.push_back(op); break;
//       }
//     }
//     return ts;
//   }
// };
