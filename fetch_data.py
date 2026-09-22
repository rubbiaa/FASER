#!/usr/bin/env python3
"""
fetch_data.py -- download named input-data files from a public CERNBox
share into $FASERDATA, skipping any that are already present.

Why: input samples like FASERMC-PO-Run10000-0_53954_3DCAL.root are too
big to commit to git (they used to live in FASERG4/ before being moved to
$FASERDATA/GENIE/ -- see docs/REGRESSION_TESTS.md and the run_faserps.py
docstring) and are *.root-gitignored, so a fresh `git clone` of this repo
has none of them. This script is the other half of that: a small, growable
manifest (REMOTE_FILES below) of "here's a public URL, here's where it
goes locally, here's its known sha256" -- add an entry here whenever a new
large input sample needs distributing, rather than inventing a new
one-off download step each time. See the module docstring's "Growing this
manifest" note below for what that looks like in practice.

The download source is a CERNBox public link (Share -> "Can view", no
password/login), not a personal EOS token: unlike a token, a public link
is *meant* to be shared and committed -- see the commit that introduced
this file for why an EOS token was the wrong tool here. CERNBox's own
docs (https://cernbox.docs.cern.ch/web/sharing/public-share/) give the
direct-download URL pattern this script builds:
    https://cernbox.cern.ch/remote.php/dav/public-files/<share token>/<path>
authenticated via HTTP Basic auth with the share token as the username
and an empty password (the standard ownCloud/Nextcloud WebDAV public-share
convention CERNBox is built on).

STATUS: the URL construction and sha256-verification logic have NOT been
exercised against a live download -- cernbox.cern.ch isn't reachable from
either of the sandboxes this was written in (both hit
"blocked-by-allowlist" from their egress proxy; confirmed against
github.com working fine from the same shells, so it's a proxy allowlist
gap, not a general outage). Please run this for real once
(`python3 fetch_data.py --force`) and let me know if it doesn't work --
most likely failure mode if the URL pattern above is subtly wrong is an
HTTP error or an HTML login/error page instead of the ROOT file, which
the sha256 check below should catch rather than silently writing garbage
into data/GENIE/.

Usage:
    python3 fetch_data.py                # fetch anything missing, skip the rest
    python3 fetch_data.py --force         # re-download everything regardless
    python3 fetch_data.py --list          # list manifest entries and exit
    python3 fetch_data.py --dry-run       # show what would be fetched, don't fetch
"""
import argparse
import base64
import hashlib
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

CERNBOX_DAV_BASE = "https://cernbox.cern.ch/remote.php/dav/public-files"

# One entry per distributable input file. "share_token" is the CERNBox
# public-link token (safe to commit -- see module docstring), "remote_path"
# is the file's path *inside* that share, "local_relpath" is where it goes
# under $FASERDATA, and "sha256"/"size_bytes" are the known-good values
# (from the copy already in this repo's history) used to verify a download
# actually succeeded instead of silently accepting a truncated file or an
# HTML error page.
#
# Growing this manifest: add a new dict here for the next big input file
# (same CERNBox "GENIE_V10" share, a different one, doesn't matter) --
# nothing else in this script needs to change. That's the point of having
# this be a small data table instead of one-off code per file.
REMOTE_FILES = [
    {
        "share_token": "xlFcHS50ZsH2Vmb",
        "remote_path": "FASERMC-PO-Run10000-0_53954_3DCAL.root",
        "local_relpath": "GENIE/FASERMC-PO-Run10000-0_53954_3DCAL.root",
        "size_bytes": 4598075,
        "sha256": "5b0a67941fd267371a877005a50b15e3ac373346acdd82d4a03683d6ec4923ab",
    },
]


def faserdata_dir():
    """Same $FASERDATA-or-REPO_ROOT/data fallback run_faserps.py's
    _default_genie_input_file() and run_batchreco.py's
    _default_geometry_file() use, so this resolves consistently with them
    whether or not setup.sh has been sourced yet."""
    return Path(os.environ.get("FASERDATA", str(REPO_ROOT / "data")))


def sha256_of(path: Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def download_one(entry: dict, *, force: bool, dry_run: bool) -> bool:
    """Returns True if the file is present and verified (whether freshly
    downloaded, already-present-and-valid, or dry-run-would-fetch)."""
    dest = faserdata_dir() / entry["local_relpath"]
    label = entry["local_relpath"]

    if dest.is_file() and not force:
        actual = dest.stat().st_size
        if actual == entry.get("size_bytes"):
            print(f"[fetch_data] {label}: already present ({actual} bytes) -- skipping")
            return True
        print(f"[fetch_data] {label}: present but wrong size "
              f"(have {actual}, expected {entry.get('size_bytes')}) -- re-fetching")

    if dry_run:
        print(f"[fetch_data] {label}: would fetch (dry run)")
        return True

    url = f"{CERNBOX_DAV_BASE}/{entry['share_token']}/{entry['remote_path']}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_dest = dest.with_suffix(dest.suffix + ".part")

    print(f"[fetch_data] {label}: fetching from {url}")
    auth = base64.b64encode(f"{entry['share_token']}:".encode()).decode()
    request = urllib.request.Request(url, headers={"Authorization": f"Basic {auth}"})
    try:
        with urllib.request.urlopen(request, timeout=120) as response, open(tmp_dest, "wb") as out:
            content_type = response.headers.get("Content-Type", "")
            if "html" in content_type.lower() or "text" in content_type.lower():
                # A public CERNBox link that's misconfigured, expired, or
                # password-protected typically serves an HTML page here
                # instead of the file -- catch that up front rather than
                # writing an "HTML page with a .root extension" and letting
                # ROOT produce a confusing error much later.
                out.close()
                tmp_dest.unlink(missing_ok=True)
                sys.exit(
                    f"error: {label}: server returned Content-Type '{content_type}', "
                    f"expected a binary file -- the public link may be expired, "
                    f"password-protected, or the remote_path in fetch_data.py's "
                    f"REMOTE_FILES is wrong. Response came from: {url}"
                )
            total = 0
            while chunk := response.read(1 << 20):
                out.write(chunk)
                total += len(chunk)
    except urllib.error.HTTPError as e:
        tmp_dest.unlink(missing_ok=True)
        sys.exit(f"error: {label}: HTTP {e.code} fetching {url} -- {e.reason}")
    except urllib.error.URLError as e:
        tmp_dest.unlink(missing_ok=True)
        sys.exit(
            f"error: {label}: could not reach {url} ({e.reason}). "
            f"Check network access to cernbox.cern.ch, or fetch the file manually "
            f"and place it at {dest}."
        )

    expected_sha = entry.get("sha256")
    if expected_sha:
        actual_sha = sha256_of(tmp_dest)
        if actual_sha != expected_sha:
            bad_dest = dest.with_suffix(dest.suffix + ".bad")
            tmp_dest.replace(bad_dest)
            sys.exit(
                f"error: {label}: sha256 mismatch after download "
                f"(expected {expected_sha}, got {actual_sha}). "
                f"Kept the bad download at {bad_dest} for inspection instead of "
                f"overwriting {dest}. This usually means the public link's "
                f"content changed, or the download was corrupted/truncated -- "
                f"not necessarily a bug in this script."
            )

    tmp_dest.replace(dest)
    print(f"[fetch_data] {label}: downloaded and verified ({total} bytes) -> {dest}")
    return True


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force", action="store_true", help="Re-download even if already present.")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be fetched, don't fetch.")
    parser.add_argument("--list", action="store_true", help="List manifest entries and exit.")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.list:
        for entry in REMOTE_FILES:
            print(f"{entry['local_relpath']}  ({entry.get('size_bytes', '?')} bytes)")
        return 0

    ok = True
    for entry in REMOTE_FILES:
        if not download_one(entry, force=args.force, dry_run=args.dry_run):
            ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
