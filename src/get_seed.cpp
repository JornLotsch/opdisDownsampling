#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>

using namespace Rcpp;

// Function to set seed and get RNG state
SEXP getRngStateForSeed(int seed) {
    // Set seed using R's RNG
    Function setSeed("set.seed");
    setSeed(seed);

    // Get .Random.seed from global environment
    Environment globalEnv = Environment::global_env();
    if (!globalEnv.exists(".Random.seed")) {
        return R_NilValue;
    }

    // Get the SEXP directly and duplicate it instead of using clone()
    SEXP randomSeed = globalEnv[".Random.seed"];
    return Rf_duplicate(randomSeed);
}

// Fast comparison of integer vectors (RNG states)
bool compareRngStates(SEXP state1, SEXP state2) {
    if (TYPEOF(state1) != INTSXP || TYPEOF(state2) != INTSXP) {
        return false;
    }

    IntegerVector v1(state1);
    IntegerVector v2(state2);

    if (v1.size() != v2.size()) {
        return false;
    }

    for (int i = 0; i < v1.size(); i++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }

    return true;
}

// Main search function - exported to R
// [[Rcpp::export]]
IntegerVector findMatchingSeedsCpp(IntegerVector candidates,
                                  SEXP targetState,
                                  bool verbose = false,
                                  int progressEvery = 10000) {

    std::vector<int> matches;
    int nCandidates = candidates.size();

    if (verbose) {
        Rcout << "Searching " << nCandidates << " candidate seeds..." << std::endl;
    }

    for (int i = 0; i < nCandidates; i++) {
        // Progress reporting
        if (verbose && (i + 1) % progressEvery == 0) {
            Rcout << "Processed " << (i + 1) << "/" << nCandidates << " candidates" << std::endl;
        }

        // Check for user interruption every 1000 iterations
        if (i % 1000 == 0) {
            checkUserInterrupt();
        }

        int seed = candidates[i];
        SEXP currentState = getRngStateForSeed(seed);

        if (currentState != R_NilValue) {
            if (compareRngStates(currentState, targetState)) {
                matches.push_back(seed);

                if (verbose) {
                    Rcout << "Found matching seed: " << seed << std::endl;
                }

                // Return immediately after first match
                break;
            }
        }
    }

    return wrap(matches);
}

// Batch processing version
// [[Rcpp::export]]
IntegerVector findMatchingSeedsBatchCpp(IntegerVector candidates,
                                       SEXP targetState,
                                       int batchSize = 1000,
                                       bool verbose = false) {

    int nCandidates = candidates.size();
    int nBatches = (nCandidates + batchSize - 1) / batchSize;

    if (verbose) {
        Rcout << "Processing " << nCandidates << " candidates in " << nBatches << " batches" << std::endl;
    }

    for (int batch = 0; batch < nBatches; batch++) {
        int startIdx = batch * batchSize;
        int endIdx = std::min(startIdx + batchSize, nCandidates);

        IntegerVector batchCandidates(endIdx - startIdx);
        for (int i = 0; i < batchCandidates.size(); i++) {
            batchCandidates[i] = candidates[startIdx + i];
        }

        if (verbose) {
            Rcout << "Processing batch " << (batch + 1) << "/" << nBatches << std::endl;
        }

        IntegerVector batchResult = findMatchingSeedsCpp(batchCandidates, targetState, false);

        if (batchResult.size() > 0) {
            return batchResult;
        }

        // Force garbage collection between batches
        R_gc();
    }

    return IntegerVector();
}