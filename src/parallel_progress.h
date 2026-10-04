// Thread-safe progress reporting and interrupt handling for OpenMP loops.
//
// Worker threads only increment an atomic counter. The R API (progress bar
// output, user-interrupt check) is touched exclusively by the master thread
// (thread 0 of the team, i.e. the R main thread), and at a throttled rate.
// RcppProgress' own increment() is not used inside parallel regions because
// its master/worker updates of the counter race with each other.

#ifndef CGEI_PARALLEL_PROGRESS_H
#define CGEI_PARALLEL_PROGRESS_H

#include <Rcpp.h>

#include <atomic>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "boxfilter.h"  // cgei::omp_thread_num()
#include "eta_progress_bar.h"

namespace cgei {

// Number of threads to use for `requested` cores (1 without OpenMP).
inline int resolve_threads(const int requested) {
#ifdef _OPENMP
  return requested < 1 ? 1 : requested;
#else
  (void)requested;
  return 1;
#endif
}

class ParallelProgress {
 public:
  ParallelProgress(const unsigned long total, const bool display)
      : total_(total), display_(display), step_(total / 200 > 0 ? total / 200 : 1),
        bar_(), pb_(total, display, bar_) {}

  // Called by any thread after finishing n work items.
  void tick(const unsigned long n = 1) {
    const unsigned long done = done_.fetch_add(n, std::memory_order_relaxed) + n;
    if (omp_thread_num() == 0) poll(done);
  }

  // True once the user interrupted the computation.
  bool aborted() const { return aborted_.load(std::memory_order_relaxed); }

  // To be called by the master thread outside of parallel regions.
  void finish() {
    if (!aborted() && Progress::check_abort()) aborted_.store(true);
    if (display_ && !aborted() && total_ > 0) pb_.update(total_);
    if (aborted()) throw Rcpp::internal::InterruptedException();
  }

 private:
  void poll(const unsigned long done) {
    if (++polls_ % 16 == 0 && Progress::check_abort()) aborted_.store(true);
    if (display_ && done - shown_ >= step_) {
      pb_.update(done);
      shown_ = done;
    }
  }

  const unsigned long total_;
  const bool display_;
  const unsigned long step_;
  std::atomic<unsigned long> done_{0};
  std::atomic<bool> aborted_{false};
  unsigned long polls_ = 0;  // master thread only
  unsigned long shown_ = 0;  // master thread only
  ETAProgressBar bar_;
  Progress pb_;
};

}  // namespace cgei

#endif
