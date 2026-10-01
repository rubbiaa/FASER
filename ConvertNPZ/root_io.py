"""Use the ROOT dictionary and data directories provided by FASER."""

import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


def data_directory():
    """Match common_setup.sh and the FASER reconstruction runner."""
    return Path(os.environ.get("FASERDATA", str(REPO_ROOT / "data"))).resolve()


def dictionary_path():
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    default = REPO_ROOT / "CoreUtils" / f"libCoreUtilsDict{suffix}"
    return Path(os.environ.get("FASER_NPZ_LIBRARY", str(default))).resolve()


def load_dictionary(path=None):
    """Load the normal FASER build's library; return the PyROOT module."""
    path = Path(path).resolve() if path is not None else dictionary_path()
    if not path.is_file():
        raise FileNotFoundError(
            f"FASER dictionary not found: {path}\n"
            f"Build FASER first: cmake --build {REPO_ROOT / 'build'} --target CoreUtilsDict -j 8")
    try:
        import ROOT
    except ImportError as exc:
        raise RuntimeError(
            "PyROOT is unavailable in this Python environment. Source FASER's setup "
            "and use the Python version that your ROOT installation was built for.") from exc
    if ROOT.gSystem.Load(str(path)) < 0:
        raise RuntimeError(f"Could not load {path}; check the FASER ROOT/GenFit environment")
    return ROOT
