#!/usr/bin/env python3
"""Local functional test for fetch_data.py's download/verify logic, using
a throwaway local HTTP server instead of the real (network-blocked from
here) cernbox.cern.ch. Exercises: skip-if-present, successful download +
sha256 verify, sha256-mismatch handling, and HTML-error-page rejection.
Does NOT test the real CERNBox URL/auth -- that needs a real network path,
see fetch_data.py's own docstring.
"""
import hashlib
import http.server
import importlib.util
import os
import shutil
import sys
import tempfile
import threading
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
spec = importlib.util.spec_from_file_location("fetch_data", REPO_ROOT / "fetch_data.py")
fetch_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fetch_data)


class Handler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        # Route by the share token embedded in the path, so different
        # test cases can get different canned responses.
        if "/BADTOKEN/" in self.path:
            body = b"<html><body>Not authorized</body></html>"
            self.send_response(200)
            self.send_header("Content-Type", "text/html")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        elif "/MISMATCHTOKEN/" in self.path:
            body = b"this is not the file you expected" * 100
            self.send_response(200)
            self.send_header("Content-Type", "application/octet-stream")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        elif "/GOODTOKEN/" in self.path:
            body = GOOD_CONTENT
            self.send_response(200)
            self.send_header("Content-Type", "application/octet-stream")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        else:
            self.send_response(404)
            self.end_headers()

    def log_message(self, fmt, *args):
        pass  # quiet


GOOD_CONTENT = os.urandom(50_000)
GOOD_SHA = hashlib.sha256(GOOD_CONTENT).hexdigest()

server = http.server.HTTPServer(("127.0.0.1", 0), Handler)
port = server.server_port
thread = threading.Thread(target=server.serve_forever, daemon=True)
thread.start()

fetch_data.CERNBOX_DAV_BASE = f"http://127.0.0.1:{port}"

results = []


def check(name, condition):
    results.append((name, condition))
    print(("PASS" if condition else "FAIL") + f": {name}")


with tempfile.TemporaryDirectory() as td:
    os.environ["FASERDATA"] = td
    fasdata = Path(td)

    # --- Test 1: successful download + sha256 verify ---
    entry_good = {
        "share_token": "GOODTOKEN", "remote_path": "file.root",
        "local_relpath": "GENIE/good.root",
        "size_bytes": len(GOOD_CONTENT), "sha256": GOOD_SHA,
    }
    ok = fetch_data.download_one(entry_good, force=False, dry_run=False)
    dest = fasdata / "GENIE" / "good.root"
    check("good download succeeds", ok and dest.is_file())
    check("good download content matches", dest.read_bytes() == GOOD_CONTENT)
    check("no leftover .part file", not (fasdata / "GENIE" / "good.root.part").exists())

    # --- Test 2: skip-if-present-and-correct-size (no second network call needed) ---
    mtime_before = dest.stat().st_mtime
    ok2 = fetch_data.download_one(entry_good, force=False, dry_run=False)
    check("second call skips (mtime unchanged)", ok2 and dest.stat().st_mtime == mtime_before)

    # --- Test 3: --force re-downloads even if present ---
    ok3 = fetch_data.download_one(entry_good, force=True, dry_run=False)
    check("force re-downloads", ok3 and dest.is_file())

    # --- Test 4: sha256 mismatch is caught, bad file quarantined, exits nonzero ---
    entry_mismatch = {
        "share_token": "MISMATCHTOKEN", "remote_path": "file.root",
        "local_relpath": "GENIE/mismatch.root",
        "size_bytes": 999999, "sha256": "0" * 64,
    }
    pid = os.fork()
    if pid == 0:
        # child: this should sys.exit(...) with an error
        devnull = open(os.devnull, "w")
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            fetch_data.download_one(entry_mismatch, force=False, dry_run=False)
            os._exit(0)  # should not reach here
        except SystemExit as e:
            os._exit(1 if e.code else 0)
    else:
        _, status = os.waitpid(pid, 0)
        exited_nonzero = os.WIFEXITED(status) and os.WEXITSTATUS(status) == 1
        check("sha256 mismatch exits with error", exited_nonzero)
        bad_file = fasdata / "GENIE" / "mismatch.root.bad"
        check("bad download quarantined as .bad, not overwriting dest", bad_file.is_file())
        check("dest never created for mismatched file", not (fasdata / "GENIE" / "mismatch.root").exists())

    # --- Test 5: HTML content-type response is rejected, not saved as the file ---
    entry_bad = {
        "share_token": "BADTOKEN", "remote_path": "file.root",
        "local_relpath": "GENIE/bad.root",
        "size_bytes": 123, "sha256": "1" * 64,
    }
    pid = os.fork()
    if pid == 0:
        devnull = open(os.devnull, "w")
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            fetch_data.download_one(entry_bad, force=False, dry_run=False)
            os._exit(0)
        except SystemExit as e:
            os._exit(1 if e.code else 0)
    else:
        _, status = os.waitpid(pid, 0)
        exited_nonzero = os.WIFEXITED(status) and os.WEXITSTATUS(status) == 1
        check("HTML response rejected", exited_nonzero)
        check("no file written for HTML response", not (fasdata / "GENIE" / "bad.root").exists())
        check("no .part leftover for HTML response", not (fasdata / "GENIE" / "bad.root.part").exists())

server.shutdown()

print()
n_fail = sum(1 for _, ok in results if not ok)
print(f"{len(results) - n_fail}/{len(results)} passed")
sys.exit(1 if n_fail else 0)
