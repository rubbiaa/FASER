#!/usr/bin/env python3
"""
Authors: Dr. Saul Alonso-Monsalve, Fabio Cufino
Convert neutrino and single-particle reconstruction to the original NPZ formats.

Both modes share hit decoding, CSR matching, segmentation, kinematics and
NPZ writing. Neutrinos retain the 41-field format, including aligned MDT truth
and extended muon diagnostics. Particles retain the 37-field format, including
the primary momentum. The default cut is 20 non-ghost PS voxels in both modes.

Build the dictionary against the enclosing FASER checkout and run
preflight_check_reco.py to check the legacy TMuTrack streamer layout.
"""

import os
import sys

import glob
import numpy as np
import tqdm
import argparse
from root_io import data_directory, load_dictionary

# =============================================================================
# CONFIGURATION AND SETUP
# =============================================================================

ROOT = None
TFile = None
tporeco_event = None
tcal_event = None


def initialize_root(geometry=None):
    """Load the matching dictionary and, for the particle reader, its geometry."""
    global ROOT, TFile, tporeco_event, tcal_event
    ROOT = load_dictionary()
    TFile = ROOT.TFile
    # Keep schema/streamer errors visible while suppressing ROOT info chatter.
    ROOT.gErrorIgnoreLevel = ROOT.kError
    if geometry is not None:
        if not os.path.isfile(geometry):
            raise FileNotFoundError(f"Geometry file not found: {geometry}")
        if not ROOT.TGeoManager.Import(geometry):
            raise RuntimeError(f"Could not load geometry: {geometry}")
    tcal_event = ROOT.TcalEvent()
    tporeco_event = ROOT.TPORecoEvent()

# =============================================================================
# CONSTANTS
# =============================================================================

# true_hits columns: track_id, parent_id, primary_id, pdg, x, y, z, module,
# energy, is_primary, is_secondary, is_tau_decay, is_charm_decay
TRACK_ID = 0
PARENT_ID = 1
PRIMARY_ID = 2
PDG = 3
ENERGY = 8

MUONIC_PDGS = [-13, 13]
ELECTROMAGNETIC_PDGS = [-11, 11, -15, 15, 22]

MAXMUTRACKS = 10
N_LAYERS_Z = 20

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_discrete_coord_from_id(ID, N_LAYERS_Z=20):
    """
    Decode a channel ID into discrete voxel coordinates (x, y, z).
    Works without geometry by unpacking the encoded ID.
    Returns a flattened z index combining layer and intra-layer depth.
    """
    hittype = ID // 100_000_000_000

    if hittype == 0:
        ix = ID % 1000
        iy = (ID // 1000) % 1000
        iz = (ID // 1_000_000) % 1000
        ilayer = ID // 1_000_000_000
        z = iz + ilayer * N_LAYERS_Z
        return ix, iy, z

    elif hittype == 1:
        ix = ID % 10000
        iy = (ID // 10000) % 10000
        ilayer = (ID // 100_000_000) % 100
        icopy = (ID // 10_000_000_000) % 10
        z = ilayer * 2 + icopy
        return ix, iy, z
    else:
        raise ValueError(f"Unknown hittype {hittype}")


def get_module_from_id(ID):
    """
    Derive the PS module index directly from the encoded channel ID.
    Replaces TcalEvent.getChannelModulefromID() which is unavailable in v9.
    """
    hittype = ID // 100_000_000_000
    if hittype == 0:
        return int(ID // 1_000_000_000)
    elif hittype == 1:
        return int((ID // 100_000_000) % 100)
    else:
        raise ValueError(f"Unknown hittype {hittype} for ID {ID}")


def load_caldata_v9(cal_dir, run_number, event_id, cal_path=None):
    """
    Load ECAL, AHCAL, true-track, and MDT-track data for one FASERCALDATA v9 event.

    Bypasses TcalEvent.Load_event() and binds only the branches needed for the NPZ.
    mdttracks provides the truth-validation block. Keep cal_file open while
    consuming the returned vectors; ROOT owns the pointed objects backing
    tracks_vec/mdttracks_vec.

    Returns:
        (cal_file, rearcal_vec, rearhcal_vec, tracks_vec, mdttracks_vec)
        Caller must call cal_file.Close() when done.
        Returns (None, empty, empty, empty, empty) if the file is missing.
    """
    if cal_path is None:
        cal_path = os.path.join(cal_dir, f"FASERG4-Tcalevent_{run_number}_{event_id}.root")
    if not os.path.exists(cal_path):
        empty_dep = ROOT.std.vector("TcalEvent::REARCALDEPOSIT")()
        empty_trk = ROOT.std.vector("DigitizedTrack*")()
        empty_mdt = ROOT.std.vector("MDTTrack*")()
        return None, empty_dep, empty_dep, empty_trk, empty_mdt

    cal_file = ROOT.TFile(cal_path, "read")
    tree = cal_file.Get("calEvent")

    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("rearcal*",  1)
    tree.SetBranchStatus("rearhcal*", 1)
    tree.SetBranchStatus("tracks*",   1)
    tree.SetBranchStatus("mdttracks*", 1)

    rearcal_vec  = ROOT.std.vector("TcalEvent::REARCALDEPOSIT")()
    rearhcal_vec = ROOT.std.vector("TcalEvent::REARCALDEPOSIT")()
    tracks_vec   = ROOT.std.vector("DigitizedTrack*")()
    mdttracks_vec = ROOT.std.vector("MDTTrack*")()

    tree.SetBranchAddress("rearcal",  rearcal_vec)
    tree.SetBranchAddress("rearhcal", rearhcal_vec)
    tree.SetBranchAddress("tracks",   tracks_vec)
    tree.SetBranchAddress("mdttracks", mdttracks_vec)
    tree.GetEntry(0)

    return cal_file, rearcal_vec, rearhcal_vec, tracks_vec, mdttracks_vec


def get_true_hits(tracks_vec, po_event, is_tau: bool, is_charmed: bool, skip_zero_id=True):
    """
    Build the original float32 (N, 13) true-hit array from DigitizedTrack objects.
    tracks_vec comes from direct CAL branches or TcalEvent.getfTracks(). The
    channel helpers decode IDs without using geometry. skip_zero_id preserves
    the different ID-0 behavior of the two production readers.
    """
    rows, hit_ids_list = [], []
    particles  = po_event.POs
    tau_decay  = po_event.taudecay
    charm_decay = po_event.charmdecay

    po_geant_ids = {trk.geanttrackID for trk in particles}

    tau_decay_geant_id = -1
    if is_tau:
        tau_parent_ids = [trk.m_trackid_in_particle[0] for trk in tau_decay if trk.nparent == 1]
        if tau_parent_ids:
            tau_parent_id = tau_parent_ids[0]
            tau_decay_geant_id = next(
                (trk.geanttrackID for trk in particles
                 if getattr(trk, "nparent", None) == 1
                 and trk.m_trackid_in_particle[0] == tau_parent_id), -1)

    charm_decay_geant_id = -1
    if is_charmed:
        charm_parent_ids = [trk.m_trackid_in_particle[0] for trk in charm_decay if trk.nparent == 1]
        if charm_parent_ids:
            charm_parent_id = charm_parent_ids[0]
            charm_decay_geant_id = next(
                (trk.geanttrackID for trk in particles if trk.m_track_id == charm_parent_id), -1)

    for trk in tracks_vec:
        track_id  = trk.ftrackID
        parent_id = trk.fparentID
        primary_id = trk.fprimaryID
        pdg        = trk.fPDG

        for hid, energy in zip(trk.fhitIDs, trk.fEnergyDeposits):
            # Preserve the neutrino reader's ID-0 cut. The original particle
            # reader retained it, so that mode passes skip_zero_id=False.
            if (skip_zero_id and hid == 0) or tcal_event.getChannelTypefromID(hid) != 0 or energy == 0:
                continue

            x, y, z = get_discrete_coord_from_id(hid, N_LAYERS_Z=N_LAYERS_Z)
            module   = tcal_event.getChannelModulefromID(hid)

            is_primary_flag     = (track_id == primary_id) and (parent_id == 0)
            is_secondary_flag   = parent_id in po_geant_ids
            is_tau_decay_flag   = parent_id == tau_decay_geant_id
            is_charm_decay_flag = parent_id == charm_decay_geant_id

            rows.append([
                track_id, parent_id, primary_id, pdg,
                x, y, z, module, energy,
                is_primary_flag, is_secondary_flag, is_tau_decay_flag, is_charm_decay_flag
            ])
            hit_ids_list.append(hid)

    if not rows:
        return None, np.empty((0,), dtype=np.int64)

    hits    = np.asarray(rows, dtype=np.float32)
    hit_ids = np.asarray(hit_ids_list, dtype=np.int64)
    return hits, hit_ids


def _build_true_index_by_id(true_ids):
    """Build an index mapping hit IDs to their positions in the true_hits array."""
    index = {}
    for i, hid in enumerate(true_ids):
        index.setdefault(int(hid), []).append(i)
    for k, v in index.items():
        index[k] = np.asarray(v, dtype=np.int32)
    return index


def get_reco_hits_and_csr_map(fPORecoEvent, true_hits, true_ids):
    """
    Processes reconstructed hits and builds CSR-style mapping from reco -> true.
    The module is decoded from the channel ID; geometry is not needed.

    Returns:
        reco_hits: (N_reco, 6) float32 array [x,y,z,module,RawEnergy,ghost_flag]
            - x, y, z are discrete voxel coordinates
        true_index: 1D int32 array of concatenated true indices
        indptr: 1D int32 array of length N_reco+1 (CSR row pointer)
        ghost_mask: 1D bool array (True for ghosts; real unmatched voxels stay False)
        link_weight: 1D float32 array, same length as true_index
                     (fraction of reco energy attributed to each true hit)
    """
    num_voxels = len(fPORecoEvent.PSvoxelmap)
    reco_hits = np.zeros((num_voxels, 6), dtype=np.float32)
    ghost_mask = np.zeros(num_voxels, dtype=bool)

    no_true = (true_hits is None) or (true_hits.shape[0] == 0)
    id_index = _build_true_index_by_id(true_ids) if not no_true else {}

    matched_lists, weight_lists = [], []

    for i, (voxel_id, psvoxel_3d) in enumerate(fPORecoEvent.PSvoxelmap):
        x, y, z = get_discrete_coord_from_id(voxel_id, N_LAYERS_Z=N_LAYERS_Z)
        module = get_module_from_id(voxel_id)

        reco_hits[i, 0:3] = [x, y, z]
        reco_hits[i, 3] = module
        reco_hits[i, 4] = psvoxel_3d.RawEnergy
        reco_hits[i, 5] = float(psvoxel_3d.ghost)

        if psvoxel_3d.ghost == 0:
            if no_true:
                reco_hits[i, 5] = 2.0
                matched_lists.append(np.empty(0, dtype=np.int32))
                weight_lists.append(np.empty(0, dtype=np.float32))
            else:
                matches = id_index.get(int(voxel_id))
                if matches is None or matches.size == 0:
                    reco_hits[i, 5] = 2.0
                    matched_lists.append(np.empty(0, dtype=np.int32))
                    weight_lists.append(np.empty(0, dtype=np.float32))
                else:
                    matched_lists.append(matches)
                    energies = true_hits[matches, ENERGY].astype(np.float32, copy=False)
                    total = energies.sum()
                    weights = energies / total if total > 0 else np.full_like(energies, 1.0 / len(energies))
                    weight_lists.append(weights)
        else:
            ghost_mask[i] = True
            matched_lists.append(np.empty(0, dtype=np.int32))
            weight_lists.append(np.empty(0, dtype=np.float32))

    counts = np.fromiter((m.size for m in matched_lists), count=num_voxels, dtype=np.int32)
    indptr = np.empty(num_voxels + 1, dtype=np.int32)
    indptr[0] = 0
    np.cumsum(counts, out=indptr[1:])
    total = int(indptr[-1])
    true_index = np.empty(total, dtype=np.int32)
    link_weight = np.empty(total, dtype=np.float32)

    offset = 0
    for m, w in zip(matched_lists, weight_lists):
        n = m.size
        if n:
            true_index[offset:offset+n] = m
            link_weight[offset:offset+n] = w
        offset += n

    return reco_hits, true_index, indptr, ghost_mask, link_weight


def process_labels_csr(true_index, indptr, ghost_mask, true_hits, out_lepton_pdg, is_cc, link_weight=None):
    """
    Computes seg_labels with optional weighted contributions.

    Returns:
        seg_labels: (num_hits, 4) float32
            [:,0] -> Ghost label (1 if ghost/unmatched, 0 otherwise)
            [:,1] -> (muonic + electromagnetic) E_dep, minus primary-lepton E
            [:,2] -> hadronic E_dep, minus primary-lepton E
            [:,3] -> primary-lepton E_dep
    """
    num_hits = indptr.size - 1
    seg_labels = np.zeros((num_hits, 4), dtype=np.float32)
    no_true = (true_hits is None) or (true_hits.shape[0] == 0)

    for i in range(num_hits):
        if ghost_mask[i] or no_true:
            seg_labels[i] = [1.0, 0.0, 0.0, 0.0]
            continue

        sl = slice(indptr[i], indptr[i+1])
        if sl.start == sl.stop:
            seg_labels[i] = [1.0, 0.0, 0.0, 0.0]
            continue

        matched = true_hits[true_index[sl]]
        pdgs = matched[:, PDG].astype(np.int32, copy=False)
        energies = matched[:, ENERGY].astype(np.float32, copy=False)
        if link_weight is not None:
            energies = energies * link_weight[sl]

        mu_mask = np.isin(pdgs, MUONIC_PDGS)
        em_mask = np.isin(pdgs, ELECTROMAGNETIC_PDGS)
        muem_mask = mu_mask | em_mask
        had_mask = ~muem_mask

        primary_mask = (
            (matched[:, TRACK_ID].astype(np.int32) == matched[:, PRIMARY_ID].astype(np.int32)) &
            (matched[:, PARENT_ID].astype(np.int32) == 0) &
            np.isin(pdgs, out_lepton_pdg)
        ) if is_cc else np.zeros_like(pdgs, dtype=bool)

        muem_sum = energies[muem_mask].sum() - energies[muem_mask & primary_mask].sum()
        had_sum = energies[had_mask].sum() - energies[had_mask & primary_mask].sum()
        prim_sum = energies[primary_mask].sum()

        seg_labels[i] = [0.0, muem_sum, had_sum, prim_sum]

    return seg_labels


def get_muon_spectrometer_truth(fMuTracks, mdttracks):
    """Build the MDT truth block (muspec_true) aligned to RECO tracks by ftrackID.

    Returns an (11, ntracks) array with the v8 row order the ML pipeline
    expects: charge, npoints, px, py, pz, p, chi2, ndof, pval, fpErr, fipErr.
    Momenta are unsmeared CAL truth; chi2=0, ndof=nhits, pval=1 and
    fpErr = fipErr = 0 are placeholders.
    """
    ntracks = min(len(fMuTracks), MAXMUTRACKS)
    true_info = np.zeros((11, ntracks), dtype=np.float32)
    if ntracks == 0:
        return true_info

    mdt_by_id = {m.ftrackID: m for m in mdttracks}

    for i, trk in enumerate(fMuTracks[:ntracks]):
        m = mdt_by_id.get(trk.ftrackID)
        if m is None or m.mom.size() == 0:
            continue
        pdg = m.fPDG
        charge = -1.0 if pdg == 13 else (1.0 if pdg == -13 else 0.0)
        nhits = m.pos.size()
        m0 = m.mom[0]
        px, py, pz = m0.x() / 1000.0, m0.y() / 1000.0, m0.z() / 1000.0  # MeV -> GeV
        p = float(np.sqrt(px * px + py * py + pz * pz))
        if p <= 0:
            continue

        true_info[:, i] = [charge, nhits, px, py, pz, p, 0.0, nhits, 1.0, 0.0, 0.0]

    return true_info


# Diagnostic members that exist only in the locally-extended TMuTrack build.
_MUSPEC_DIAG_FIELDS = (
    "fanalytic_ok", "fkalman_ok", "ffit_status", "fnStations",
    "fQOverPAnalytic", "fQOverPAnalyticErr", "fchi2Analytic", "fnDoFAnalytic",
)


def get_muon_spectrometer(fMuTracks, extended=True):
    """Extract muon spectrometer track information (same layout as v8).

    Reads the RECO TPORecoEvent.fMuTracks (ReconstructMDT + GenFit Kalman fit).
    muspec_info rows: charge, npoints, px, py, pz, p, chi2, ndof, pval, fpErr, fipErr.

    Also returns muspec_extra with rows fit_ok, p_analytic, qoverp, qoverp_err
    and row 4 charge_mode. Failed fits (fit_ok=0) have npoints=0,
    ndof=0 and placeholder fit values but a valid p_analytic.

    charge_mode records how ReconstructMDT decided the sign: 0 ambiguous /
    both hypotheses failed, 1 both converged with a clear chi2/NDF winner,
    2 only mu- converged, 3 only mu+ converged, 4 slope-change tiebreaker,
    5 rescue path. Row 4 is appended, so the historical rows 0-3 keep their
    positions for existing readers.

    trackid is returned separately.
    """
    ntracks = min(len(fMuTracks), MAXMUTRACKS)
    rows = np.array([
        [t.fcharge for t in fMuTracks[:ntracks]],
        [t.fpos.size() for t in fMuTracks[:ntracks]],
        [t.fpx for t in fMuTracks[:ntracks]],
        [t.fpy for t in fMuTracks[:ntracks]],
        [t.fpz for t in fMuTracks[:ntracks]],
        [t.fp for t in fMuTracks[:ntracks]],
        [t.fchi2 for t in fMuTracks[:ntracks]],
        [t.fnDoF for t in fMuTracks[:ntracks]],
        [t.fpval for t in fMuTracks[:ntracks]],
        [t.fpErr for t in fMuTracks[:ntracks]],
        [t.fipErr for t in fMuTracks[:ntracks]],
    ], dtype=np.float32)

    if not extended:
        return ntracks, rows

    extra = np.array([
        [float(t.ffit_ok) for t in fMuTracks[:ntracks]],
        [t.fpAnalytic for t in fMuTracks[:ntracks]],
        [t.fQOverP for t in fMuTracks[:ntracks]],
        [t.fQOverPErr for t in fMuTracks[:ntracks]],
        [float(t.fchargeMode) for t in fMuTracks[:ntracks]],
    ], dtype=np.float32)

    # NaN-filled against upstream TMuTrack, populated against the extended build.
    diagnostics = np.array([
        [float(getattr(t, field, np.nan)) for t in fMuTracks[:ntracks]]
        for field in _MUSPEC_DIAG_FIELDS
    ], dtype=np.float32).reshape(len(_MUSPEC_DIAG_FIELDS), ntracks)

    trackid = np.array(
        [[t.ftrackID for t in fMuTracks[:ntracks]]], dtype=np.int32
    ).reshape(1, ntracks)

    return ntracks, rows, extra, diagnostics, trackid


def build_pair_array(objs, dtype=None):
    """Build a structured NumPy array with fields: g4_id, track_id, parent_id, pdg."""
    if dtype is None:
        dtype = np.dtype([
            ('g4_id', np.int32),
            ('track_id', np.int32),
            ('parent_id', np.int32),
            ('pdg', np.int32),
        ])

    def row(o):
        parent_id = o.m_trackid_in_particle[0] if o.nparent == 1 else -1
        return (o.geanttrackID, o.m_track_id, parent_id, o.m_pdg_id)

    return np.fromiter((row(o) for o in objs), dtype=dtype)


# HCAL
def getChannelXYZRearHCal(moduleID):
    """Extract XYZ coordinates from rear AHCal module ID."""
    x = moduleID % 1000
    y = (moduleID // 1000) % 1000
    z = (moduleID // 1000000) % 1000
    return x, y, z


# ECAL
def getChannelXYZRearCal(moduleID):
    """Extract XYZ coordinates from rear ECAL module ID."""
    x = moduleID % 1000
    y = (moduleID // 1000) % 1000
    z = (moduleID // 1000000) % 1000
    return x, y, z


def th2d_to_dense_numpy(hist):
    """Convert a ROOT TH2 histogram to a dense array with shape (y_bins, x_bins)."""
    if not hist:
        return np.empty((0, 0), dtype=np.float32)

    n_bins_x = hist.GetNbinsX()
    n_bins_y = hist.GetNbinsY()
    values = np.zeros((n_bins_y, n_bins_x), dtype=np.float32)
    for ix in range(1, n_bins_x + 1):
        for iy in range(1, n_bins_y + 1):
            values[iy - 1, ix - 1] = hist.GetBinContent(ix, iy)
    return values


def extract_ps_views(tporeco_event, geom_detector):
    """Extract XZ, YZ, and the physical XY PS views from TPORecoEvent."""
    n_modules = int(geom_detector.NRep)
    voxel = float(geom_detector.fScintillatorVoxelSize)
    nz_view = int(float(geom_detector.fSandwichLength) / voxel) if voxel > 0 else N_LAYERS_Z

    xz_view = th2d_to_dense_numpy(tporeco_event.Get2DViewXPS())
    yz_view = th2d_to_dense_numpy(tporeco_event.Get2DViewYPS())

    xy_views = []
    z_views = tporeco_event.zviewPS
    for layer in range(n_modules):
        hist = z_views[layer] if layer < len(z_views) else None
        xy_views.append(th2d_to_dense_numpy(hist))

    return {
        "view_xz": xz_view,
        "view_yz": yz_view,
        "view_xy": np.stack(xy_views, axis=0) if xy_views else np.empty((0, 0, 0), dtype=np.float32),
        "view_n_modules": np.asarray(n_modules, dtype=np.int32),
        "view_zviewps_size": np.asarray(len(z_views), dtype=np.int32),
        "view_nz_bins_per_module": np.asarray(nz_view, dtype=np.int32),
        "view_nz_scint_bins_per_module": np.asarray(N_LAYERS_Z, dtype=np.int32),
    }


def divide_list_into_chunks(input_list, num_chunks=1):
    chunk_size, remainder = divmod(len(input_list), num_chunks)
    chunks, start = [], 0
    for i in range(num_chunks):
        end = start + chunk_size + (1 if i < remainder else 0)
        chunks.append(input_list[start:end])
        start = end
    return chunks


# =============================================================================
# MAIN PROCESSING FUNCTION
# =============================================================================

def generate_events(number, chunks, disable, max_events=None, output_dir=None, reco_paths=None,
                    version=None, base_path=None, include_views=False, min_non_ghost=20,
                    skip_existing=False, tcal_dir=None, mode="neutrino", geometry=None):
    """Process RECO files with shared calculations and the selected NPZ format."""
    if reco_paths is None:
        raise ValueError("reco_paths must be provided")
    if mode not in ("neutrino", "single-particle"):
        raise ValueError(f"Unknown mode: {mode}")
    single_particle = mode == "single-particle"
    if single_particle and geometry is None:
        geometry = os.environ.get(
            "FASER_NPZ_GEOMETRY", str(data_directory() / "GDML" / "FASERCAL_V10.gdml"))
    initialize_root(geometry=geometry)

    _version = version or "v9.0_6000"
    if not _version.startswith("v"):
        _version = "v" + _version

    root_base_path = base_path or str(data_directory())
    out_dir = output_dir if output_dir is not None else str(data_directory() / "npz")
    os.makedirs(out_dir, exist_ok=True)

    # Support both the original chunk directories and current flat faserG4
    # output, including interaction suffixes such as _nueCC.root.
    cal_root = tcal_dir or os.path.join(
        root_base_path, f"FASERCALDATA_{_version}" if base_path or version else "faserG4")
    cal_paths_by_event = {}
    cal_paths = glob.glob(os.path.join(cal_root, "FASERG4-Tcalevent_*_*.root"))
    cal_paths += glob.glob(os.path.join(cal_root, "chunk_*", "FASERG4-Tcalevent_*_*.root"))
    for cal_path in sorted(cal_paths):
        try:
            parts = os.path.basename(cal_path).removesuffix(".root").split("_")
            cal_key = (int(parts[1]), int(parts[2]))
        except (IndexError, ValueError):
            continue
        previous = cal_paths_by_event.get(cal_key)
        if previous is not None and not os.path.samefile(previous, cal_path):
            raise ValueError(f"Multiple CAL files for run/event {cal_key}: {previous}, {cal_path}")
        cal_paths_by_event[cal_key] = cal_path

    if number == 0:
        print("--- File Discovery Check ---")
        print(f"Mode: {mode}")
        print(f"Total RECO files found: {len(reco_paths)}")
        print(f"CAL folder: {cal_root}")
        print(f"Output folder: {out_dir}")
        print(f"CAL files indexed by run/event ID: {len(cal_paths_by_event)}")
        print(f"Minimum non-ghost PS voxels: {min_non_ghost}")
        print("--- End Check ---\n")

    all_chunks = divide_list_into_chunks(reco_paths, num_chunks=chunks)
    if number >= len(all_chunks):
        print(f"Error: Chunk number {number} is out of range for {len(all_chunks)} chunks.")
        return

    chunk = all_chunks[number]
    if not chunk:
        print(f"Chunk {number} is empty. Skipping.")
        return

    t = tqdm.tqdm(enumerate(chunk), total=len(chunk), disable=disable)
    events_saved_in_chunk = 0

    for i, reco_file_path in t:
        reco_file = TFile(reco_file_path, "read")
        if not reco_file or reco_file.IsZombie():
            raise OSError(f"Could not open RECO file: {reco_file_path}")
        reco_tree = reco_file.Get("RecoEvent")
        if not reco_tree:
            empty = reco_file.GetListOfKeys().GetSize() == 0
            # BatchReco creates these diagnostics before processing events,
            # and creates RecoEvent only after the first accepted event.
            energy_hist = reco_file.Get("h_fullevent_Evis")
            diagnostics = [reco_file.Get(name) for name in ("ParticleGun", "MuonSpectrometer")]
            empty = empty or bool(
                energy_hist and energy_hist.InheritsFrom("TH1") and energy_hist.GetEntries() == 0
                and all(tree and tree.InheritsFrom("TTree") and tree.GetEntries() == 0
                        for tree in diagnostics))
            reco_file.Close()
            if empty:
                print(f"Skipping empty reconstruction file: {reco_file_path}")
                continue
            raise ValueError(f"No RecoEvent tree in {reco_file_path}")
        if not single_particle:
            reco_tree.SetBranchStatus("*", 0)
            for branch_pattern in ("fTPOEvent*", "PSvoxelmap*", "fMuTracks*"):
                reco_tree.SetBranchStatus(branch_pattern, 1)
            if include_views:
                # geom_detector is only needed to bin the PS views
                for branch_pattern in ("xviewPS*", "yviewPS*", "zviewPS*", "geom_detector*"):
                    reco_tree.SetBranchStatus(branch_pattern, 1)
        reco_tree.SetBranchAddress("TPORecoEvent", tporeco_event)
        total_entries = reco_tree.GetEntries()

        for entry_idx in range(total_entries):
            reco_tree.GetEntry(entry_idx)

            po_event = tporeco_event.GetPOEvent()

            run_number, event_id = po_event.run_number, po_event.event_id

            # Resume: an npz for this event already exists, so skip it before
            # opening the CAL file or running the fit -- that is where the time
            # goes. Lets an interrupted pass pick up instead of redoing work.
            if skip_existing and os.path.exists(
                    f'{out_dir}/run_{run_number}_event_{event_id}.npz'):
                continue

            cal_path = cal_paths_by_event.get((int(run_number), int(event_id)))
            if cal_path is None:
                raise FileNotFoundError(
                    f"No CAL file indexed for run {run_number}, event {event_id} under {cal_root}")
            cal_dir = os.path.dirname(cal_path)
            is_cc = bool(po_event.isCC)
            is_es = bool(po_event.isES())
            is_tau = bool(po_event.istau)
            is_charmed = bool(po_event.isCharmed())

            po = build_pair_array(po_event.POs)
            tau_decay = build_pair_array(po_event.taudecay)
            charm_decay = build_pair_array(po_event.charmdecay)

            primary_vertex = np.array([po_event.prim_vx.x(), po_event.prim_vx.y(), po_event.prim_vx.z()])
            if single_particle:
                if len(po_event.POs) > 0:
                    primary = po_event.POs[0]
                    primary_momentum = np.array([primary.m_px, primary.m_py, primary.m_pz])
                else:
                    primary_momentum = np.zeros(3)
            e_vis = po_event.Evis
            sp_momentum = np.array([po_event.spx, po_event.spy, po_event.spz])
            vis_sp_momentum = np.array([po_event.vis_spx, po_event.vis_spy, po_event.vis_spz])
            jet_momentum = np.array([po_event.jetpx, po_event.jetpy, po_event.jetpz])
            pt_miss = po_event.ptmiss

            in_neutrino, out_lepton = po_event.in_neutrino, po_event.out_lepton
            in_neutrino_pdg = in_neutrino.m_pdg_id
            in_neutrino_momentum = np.array([in_neutrino.m_px, in_neutrino.m_py, in_neutrino.m_pz])
            in_neutrino_energy = in_neutrino.m_energy
            out_lepton_pdg = out_lepton.m_pdg_id
            out_lepton_momentum = np.array([out_lepton.m_px, out_lepton.m_py, out_lepton.m_pz])
            out_lepton_energy = out_lepton.m_energy

            tau_vis_momentum = np.array([po_event.tauvis_px, po_event.tauvis_py, po_event.tauvis_pz])
            tau_decay_mode = int(po_event.tau_decaymode)
            tau_decay_length = float(po_event.tauDecaylength())
            tau_kink_angle = float(po_event.tauKinkAngle())

            # Extract PS views (only when --include-views is set; disabled by default)
            view_data = extract_ps_views(tporeco_event, tporeco_event.geom_detector) if include_views else {}

            # Keep the original loading methods: the particle reader also
            # refreshes po_event from CAL before extracting truth-hit labels.
            cal_file = None
            if single_particle:
                if tcal_event.Load_event(cal_dir + os.sep, run_number, event_id, 0, po_event) != 0:
                    raise RuntimeError(f"Could not load particle CAL data for run {run_number}, event {event_id}")
                rearcal, rearhcal = tcal_event.rearCalDeposit, tcal_event.rearHCalDeposit
                tracks_vec = tcal_event.getfTracks()
            else:
                cal_file, rearcal, rearhcal, tracks_vec, mdttracks = load_caldata_v9(
                    cal_dir, run_number, event_id, cal_path=cal_path)

            ecal_hits = np.zeros((rearcal.size(), 4), dtype=np.float32)
            for h_idx, x in enumerate(rearcal):
                ecal_hits[h_idx, :3] = getChannelXYZRearCal(x.moduleID)
                ecal_hits[h_idx, 3] = x.energyDeposit

            ahcal_hits = np.zeros((rearhcal.size(), 4), dtype=np.float32)
            for h_idx, x in enumerate(rearhcal):
                ahcal_hits[h_idx, :3] = getChannelXYZRearHCal(x.moduleID)
                ahcal_hits[h_idx, 3] = x.energyDeposit

            true_hits, true_ids = get_true_hits(
                tracks_vec, po_event, is_tau, is_charmed, skip_zero_id=not single_particle)
            # Model-facing muon spectrometer: corrected direct fMuTracks in the
            # established v8 11-row layout. Keep aligned CAL truth separately.
            if single_particle:
                muspec_ntracks, muspec_info = get_muon_spectrometer(tporeco_event.fMuTracks, extended=False)
            else:
                (muspec_ntracks, muspec_reco_info, muspec_reco_extra,
                 muspec_reco_diagnostics,
                 muspec_reco_trackid) = get_muon_spectrometer(tporeco_event.fMuTracks)
                muspec_true = get_muon_spectrometer_truth(tporeco_event.fMuTracks, mdttracks)
                muspec_info = muspec_reco_info

            if cal_file:
                cal_file.Close()

            # --- Reco hits and CSR mapping ---
            reco_hits, true_index, indptr, ghost_mask, link_weight = get_reco_hits_and_csr_map(
                tporeco_event, true_hits, true_ids
            )

            n_non_ghost = int((~ghost_mask).sum())
            if min_non_ghost > 0 and n_non_ghost < min_non_ghost:
                continue

            seg_labels = process_labels_csr(
                true_index, indptr, ghost_mask, true_hits,
                out_lepton_pdg, is_cc, link_weight=link_weight
            )

            event_data = {
                'run_number': run_number,
                'event_id': event_id,
                'is_cc': is_cc,
                'is_es': is_es,
                'is_tau': is_tau,
                'is_charmed': is_charmed,
                'po': po,
                'tau_decay': tau_decay,
                'charm_decay': charm_decay,
                'e_vis': e_vis,
                'sp_momentum': sp_momentum,
                'vis_sp_momentum': vis_sp_momentum,
                'jet_momentum': jet_momentum,
                'pt_miss': pt_miss,
                'primary_vertex': primary_vertex,
                'true_hits': true_hits,
                'reco_hits': reco_hits,
                'true_index': true_index,
                'indptr': indptr,
                'ghost_mask': ghost_mask,
                'link_weight': link_weight,
                'in_neutrino_pdg': in_neutrino_pdg,
                'in_neutrino_momentum': in_neutrino_momentum,
                'in_neutrino_energy': in_neutrino_energy,
                'out_lepton_pdg': out_lepton_pdg,
                'out_lepton_momentum': out_lepton_momentum,
                'out_lepton_energy': out_lepton_energy,
                'tau_vis_momentum': tau_vis_momentum,
                'tau_decay_mode': tau_decay_mode,
                'tau_decay_length': tau_decay_length,
                'tau_kink_angle': tau_kink_angle,
                'ecal_hits': ecal_hits,
                'ahcal_hits': ahcal_hits,
                'seg_labels': seg_labels,
                'muspec_ntracks': muspec_ntracks,
                'muspec_info': muspec_info,
            }
            if single_particle:
                event_data['primary_momentum'] = primary_momentum
            else:
                if true_hits is None:
                    event_data['true_hits'] = np.empty((0, 13), dtype=np.float32)
                event_data.update({
                    'muspec_true': muspec_true,
                    'muspec_reco_info': muspec_reco_info,
                    'muspec_reco_extra': muspec_reco_extra,
                    'muspec_reco_diagnostics': muspec_reco_diagnostics,
                    'muspec_reco_trackid': muspec_reco_trackid,
                })
            event_data.update(view_data)

            output_filename = f'{out_dir}/run_{run_number}_event_{event_id}.npz'
            np.savez_compressed(output_filename, **event_data)

            events_saved_in_chunk += 1
            t.set_description(f"File: {os.path.basename(reco_file_path)} | Saved (Run:{run_number}, Event:{event_id})")

            if max_events is not None and events_saved_in_chunk >= max_events:
                reco_file.Close()
                return events_saved_in_chunk

        reco_file.Close()

    print(f"\nSuccessfully processed chunk {number}/{chunks}.")
    print(f"Saved {events_saved_in_chunk} NPZ files in {out_dir}")
    return events_saved_in_chunk


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Convert FASER reconstruction to neutrino or single-particle NPZ files.")
    parser.add_argument("--mode", choices=("neutrino", "single-particle"), required=True,
                        help="Preserve the corresponding production NPZ format")
    parser.add_argument('--number', type=int, default=0, help="Chunk number to process (0-based; default: 0)")
    parser.add_argument('--chunks', type=int, default=1, help="Total number of chunks (default: 1)")
    parser.add_argument("--disable", action="store_true", default=False, help="Disable progress bar")
    parser.add_argument("--max-events", type=int, default=None, help="Stop after saving this many events")
    parser.add_argument("--output-dir", type=str, default=None, help="Output directory (default: $FASERDATA/npz)")
    parser.add_argument("--base-path", type=str, default=None, help="Legacy directory containing FASERCALDATA_* and FASERCALRECODATA_*")
    parser.add_argument("--reco-dir", type=str, default=None,
                        help="RECO directory (default: $FASERDATA/batch)")
    parser.add_argument("--reco-file", type=str, default=None, help="Convert one RECO file instead of a directory")
    parser.add_argument("--tcal-dir", type=str, default=None,
                        help="CAL directory, flat or with chunk_* subdirectories (default: $FASERDATA/faserG4)")
    parser.add_argument("--version", type=str, default=None, help="Use the legacy dataset layout (e.g. v10.0_10000)")
    parser.add_argument("--min-non-ghost", type=int, default=20, help="Minimum number of non-ghost PS voxels required to save an event; use 0 to disable")
    parser.add_argument("--skip-existing", action="store_true", default=False,
                        help="Skip events whose npz already exists (resume an interrupted pass)")
    parser.add_argument("--include-views", action="store_true", help="Store xviewPS, yviewPS, and physical zviewPS XY views")
    parser.add_argument("--energy", type=float, default=None,
                        help="Single-particle mode: append an energy subdirectory (e.g. 50GeV) to --output-dir")
    parser.add_argument("--geometry", type=str, default=None,
                        help="Single-particle GDML (default: $FASERDATA/GDML/FASERCAL_V10.gdml; FASER_NPZ_GEOMETRY overrides it)")
    args = parser.parse_args()

    if not (0 <= args.number < args.chunks):
        parser.error(f"number must be in [0, {args.chunks-1}]")
    if args.reco_file and args.reco_dir:
        parser.error("choose --reco-file or --reco-dir")
    if args.max_events is not None and args.max_events < 1:
        parser.error("--max-events must be positive")
    if args.mode != "single-particle" and (args.energy is not None or args.geometry is not None):
        parser.error("--energy and --geometry apply to single-particle mode")

    version_val = args.version or ("v9.0_6000" if args.base_path else None)
    if version_val and not version_val.startswith("v"):
        version_val = f"v{version_val}"

    reco_dir = os.path.abspath(args.reco_dir) if args.reco_dir else os.path.join(
        args.base_path or str(data_directory()),
        f"FASERCALRECODATA_{version_val}" if version_val else "batch")
    # Production chunk directories also contain one-file-per-event input
    # symlinks. Only BatchReco outputs are RECO files.
    _reco_paths = [os.path.abspath(args.reco_file)] if args.reco_file else sorted(
        glob.glob(os.path.join(reco_dir, "Batch-TPORecevent_*.root")) +
        glob.glob(os.path.join(reco_dir, "chunk_*", "Batch-TPORecevent_*.root")))

    if not _reco_paths:
        print(f"Error: No RECO ROOT files found in {reco_dir}")
        sys.exit(1)

    out_dir = args.output_dir or str(data_directory() / "npz")
    if args.energy is not None:
        out_dir = os.path.join(out_dir, f"{args.energy:.0f}GeV")

    saved = generate_events(args.number, args.chunks, args.disable, max_events=args.max_events,
                    output_dir=out_dir, reco_paths=_reco_paths, version=version_val,
                    base_path=args.base_path, include_views=args.include_views,
                    min_non_ghost=args.min_non_ghost, skip_existing=args.skip_existing,
                    tcal_dir=args.tcal_dir, mode=args.mode, geometry=args.geometry)
    return 1 if args.mode == "single-particle" and not saved and not args.skip_existing else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as exc:
        sys.exit(f"error: {exc}")
