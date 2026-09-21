#ifndef _FASERDATADIR_HH_
#define _FASERDATADIR_HH_ 1

#include <cstdlib>
#include <stdexcept>
#include <string>

#include <TSystem.h>

// Resolves FASER's consolidated simulation/reconstruction data tree,
// $FASERDATA (exported by setup.sh/common_setup.sh - see those scripts),
// instead of the hardcoded relative "output/"/"input/" paths and
// inter-directory symlinks (FASERG4/output, Batch/input, EvDisplay/input,
// data/batch/input) this replaces. Every executable that reads or writes
// FASERG4/batchreco data should go through this instead of building its
// own relative or hand-typed path, so there's exactly one place that
// knows where the data actually lives - and so a checkout can point
// FASERDATA somewhere else entirely (a scratch area, an EOS path, a
// per-target location) just by setting the environment variable, with no
// code change here.
namespace FASER {

// Returns "$FASERDATA/subdir" (or plain "$FASERDATA" if subdir is empty),
// with no trailing slash, creating it - and any missing parent
// directories - if it doesn't already exist (like `mkdir -p`, via ROOT's
// gSystem so the behaviour is identical on macOS/Linux/lxplus).
//
// Throws std::runtime_error, rather than silently falling back to a
// relative path, if FASERDATA isn't set at all (source setup.sh first)
// or the directory can't be created (e.g. a permissions problem) -
// exactly the class of silent-wrong-location bug this replaces.
inline std::string GetDataDir(const std::string& subdir = "") {
    const char* base = std::getenv("FASERDATA");
    if (!base || *base == '\0') {
        throw std::runtime_error(
            "FASERDATA is not set. Source setup.sh (or common_setup.sh) "
            "before running this executable, or export FASERDATA yourself "
            "to point at FASER's consolidated data directory.");
    }

    std::string path = base;
    if (!subdir.empty()) {
        path += "/" + subdir;
    }

    // gSystem->AccessPathName(...) follows POSIX access() convention: it
    // returns kFALSE (0) if the path exists and is accessible, so "true"
    // here means "does NOT exist yet".
    if (gSystem->AccessPathName(path.c_str())) {
        // Second argument requests recursive creation of missing parents,
        // like `mkdir -p`. A nonzero return can also just mean "another
        // process created it first" under a race, so only treat this as a
        // real failure if the directory still doesn't exist afterwards.
        if (gSystem->mkdir(path.c_str(), true) != 0 && gSystem->AccessPathName(path.c_str())) {
            throw std::runtime_error("FASER::GetDataDir: could not create directory: " + path);
        }
    }

    return path;
}

} // namespace FASER

#endif
