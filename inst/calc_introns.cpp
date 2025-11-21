#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame calc_all_introns(IntegerVector start, IntegerVector end,
                           CharacterVector seqnames,
                           CharacterVector strand,
                           CharacterVector transcript_id) {
  int n = start.size();
  std::vector<int> gs, ge;
  std::vector<std::string> gseq, gstrand, gtx;

  for(int i = 0; i < n - 1; i++) {
    // 如果当前 exon 和下一个 exon 属于同一个 transcript
    if(transcript_id[i] == transcript_id[i+1]) {
      int s = end[i] + 1;
      int e = start[i+1] - 1;
      if(s <= e) {
        gs.push_back(s);
        ge.push_back(e);
        gseq.push_back(as<std::string>(seqnames[i]));
        gstrand.push_back(as<std::string>(strand[i]));
        gtx.push_back(as<std::string>(transcript_id[i]));
      }
    }
  }

  return DataFrame::create(
    _["seqnames"] = gseq,
    _["start"] = gs,
    _["end"] = ge,
    _["strand"] = gstrand,
    _["transcript_id"] = gtx
  );
}
