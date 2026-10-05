#!/usr/bin/env python3
"""Build a HTCondor queue list of (chunk, start_evt, nevt) triples.

Splits a total number of events into fixed-size chunks and writes one
"<chunk> <start> <nevt>" line per chunk, suitable for:

    queue chunk, start, nevt from jobs.list

Example:
    ./make_jobs_list.py --total-events 350432 --chunk-size 5000 --out jobs_r6003.list
"""
import argparse
import sys


def build_rows(total_events: int, chunk_size: int):
    rows = []
    start = 0
    chunk = 0
    while start < total_events:
        nevt = min(chunk_size, total_events - start)
        rows.append((chunk, start, nevt))
        start += chunk_size
        chunk += 1
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--total-events", type=int, required=True, help="Total number of events to process")
    p.add_argument("--chunk-size", type=int, default=5000,
                   help="Events per chunk (default: 5000; pick 10000-50000, aim for ~10-20 chunks)")
    p.add_argument("--out", default="jobs.list", help="Output list file (default: jobs.list)")
    args = p.parse_args()

    if args.total_events <= 0:
        sys.exit("ERROR: --total-events must be a positive integer")
    if args.chunk_size <= 0:
        sys.exit("ERROR: --chunk-size must be a positive integer")

    rows = build_rows(args.total_events, args.chunk_size)

    with open(args.out, "w") as f:
        for chunk_id, start_evt, nevt in rows:
            f.write(f"{chunk_id} {start_evt} {nevt}\n")

    n_chunks = len(rows)
    print(f"Wrote {args.out} with {n_chunks} chunk(s) "
          f"(total_events={args.total_events}, chunk_size={args.chunk_size})")
    if n_chunks < 5 or n_chunks > 30:
        print(f"NOTE: {n_chunks} chunks is outside the usual 10-20 sweet spot; "
              f"consider adjusting --chunk-size.", file=sys.stderr)


if __name__ == "__main__":
    main()
