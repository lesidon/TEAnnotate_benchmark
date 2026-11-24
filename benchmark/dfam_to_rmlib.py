#!/usr/bin/env python3
"""Turn a Dfam FASTA ('>ACCESSION name') into a RepeatMasker-ready library
('>name#Type/SubType') by attaching Dfam classifications.

Dfam's FASTA API (format=fasta) omits the RepeatMasker class tag, so a library fed
straight to `RepeatMasker -lib` is annotated as 'Unknown' and lib-test.py cannot
split it into TE categories (LTR/LINE/SINE/…). This joins each family to its Dfam
`repeat_type_name` / `repeat_subtype_name` — fetched from the summary API (paginated),
or read from a cached --summary-json — and rewrites the headers.

Usage:
  # fetch classifications from the API for clade 9606 (human) + descendants:
  python dfam_to_rmlib.py --fasta dfam_human.raw.fasta --out human.curated.lib.fa --clade 9606
  # or use a cached summary JSON (offline / reproducible):
  python dfam_to_rmlib.py --fasta dfam_human.raw.fasta --out human.curated.lib.fa \
      --summary-json human_summary.json
"""
import argparse
import base64
import gzip
import json
import sys
import time
import urllib.parse
import urllib.request

API = "https://www.dfam.org/api/families"


def _get_json(url, retries=4):
    for n in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except Exception as e:                       # noqa: BLE001 (network is best-effort)
            if n == retries - 1:
                raise
            sys.stderr.write(f"[dfam_to_rmlib] retry {n + 1}/{retries} ({e})\n")
            time.sleep(5)


def fetch_summary(clade, descendants=True, page=1000):
    """Return {accession_base: (type, subtype, name)} for a clade via the API."""
    families, start = {}, 0
    while True:
        params = {"clade": clade, "format": "summary", "limit": page, "start": start}
        if descendants:
            params["clade_relatives"] = "descendants"
        data = _get_json(f"{API}?{urllib.parse.urlencode(params)}")
        results = data.get("results", []) if isinstance(data, dict) else data
        if not results:
            break
        _ingest(results, families)
        got = len(results)
        start += got
        total = data.get("total_count") if isinstance(data, dict) else None
        if got < page or (total is not None and start >= total):
            break
        time.sleep(0.5)
    return families


def _ingest(results, families):
    for fam in results:
        acc = (fam.get("accession") or "").split(".")[0]
        if acc:
            families[acc] = (fam.get("repeat_type_name") or "",
                             fam.get("repeat_subtype_name") or "",
                             fam.get("name") or "")


def rm_class(rtype, subtype):
    rtype = (rtype or "").strip()
    subtype = (subtype or "").strip()
    if not rtype:
        return "Unknown"
    return f"{rtype}/{subtype}" if subtype else rtype


def parse_header(line):
    """'>DF000000053.4 AluYa5 …' -> (accession_base, name)."""
    parts = line[1:].strip().split()
    acc = parts[0] if parts else ""
    name = parts[1] if len(parts) > 1 else acc
    return acc.split(".")[0], name


def load_fasta_text(path):
    """Return FASTA text from a raw FASTA, a .gz, or Dfam's download=true envelope.

    The Dfam API with download=true returns JSON: {"attachment":..., "body":"<b64>"}
    where <b64> is base64(gzip(fasta)). Also transparently handles a gzip file, and
    rejects an HTML error page (e.g. a GitHub blob URL saved by mistake)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    head = raw.lstrip()[:1]
    if head == b"{":                                  # Dfam JSON envelope
        obj = json.loads(raw.decode("utf-8", "replace"))
        body = obj.get("body", "")
        if not body:
            raise ValueError(f"{path}: Dfam JSON envelope has no 'body'")
        blob = base64.b64decode(body)
        if blob[:2] == b"\x1f\x8b":
            blob = gzip.decompress(blob)
        return blob.decode("utf-8", "replace")
    if raw[:2] == b"\x1f\x8b":                         # bare gzip
        return gzip.decompress(raw).decode("utf-8", "replace")
    if head == b"<":                                  # HTML page, not a library
        raise ValueError(f"{path} looks like an HTML page, not a Dfam FASTA")
    return raw.decode("utf-8", "replace")


def rewrite(fasta_in, out_path, classmap):
    n = unknown = 0
    text = load_fasta_text(fasta_in)
    with open(out_path, "w", encoding="utf-8") as out:
        for line in text.splitlines():
            if line.startswith(">"):
                n += 1
                acc, name = parse_header(line)
                rtype, subtype, _ = classmap.get(acc, ("", "", ""))
                cls = rm_class(rtype, subtype)
                if cls == "Unknown":
                    unknown += 1
                # RepeatMasker library entry: sanitise whitespace in the id.
                out.write(f">{name.replace(' ', '_')}#{cls}\n")
            else:
                out.write(line + "\n")
    return n, unknown


def load_summary_json(path):
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)
    results = raw.get("results", raw) if isinstance(raw, dict) else raw
    families = {}
    _ingest(results, families)
    return families


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="Dfam FASTA to reclassify")
    ap.add_argument("--out", help="RepeatMasker-ready library FASTA (required unless --check)")
    ap.add_argument("--clade", default="9606", help="NCBI taxon id (default 9606 = human)")
    ap.add_argument("--no-descendants", dest="descendants", action="store_false",
                    help="do not include clade descendants (default: include)")
    ap.add_argument("--summary-json", default=None,
                    help="use a cached Dfam summary JSON instead of calling the API")
    ap.add_argument("--check", action="store_true",
                    help="only validate that --fasta fully decodes to FASTA (exit 0/1); no "
                         "API call, no --out. Used to detect a truncated Dfam download.")
    ap.set_defaults(descendants=True)
    args = ap.parse_args(argv)

    if args.check:
        # Decode-only: a truncated download raises here (e.g. base64 'Incorrect padding'),
        # so a non-zero exit tells the caller to re-download rather than reuse a bad file.
        try:
            text = load_fasta_text(args.fasta)
        except Exception as e:                       # noqa: BLE001
            sys.stderr.write(f"[dfam_to_rmlib] --check failed: {e}\n")
            return 1
        return 0 if ">" in text else 1
    if not args.out:
        ap.error("--out is required unless --check is given")

    classmap = (load_summary_json(args.summary_json) if args.summary_json
                else fetch_summary(args.clade, args.descendants))
    if not classmap:
        sys.stderr.write("[dfam_to_rmlib] WARNING: no classifications obtained; "
                         "headers will be '#Unknown' (total-level scoring only)\n")
    n, unknown = rewrite(args.fasta, args.out, classmap)
    sys.stderr.write(f"[dfam_to_rmlib] {n} families -> {args.out}; "
                     f"{len(classmap)} classifications, {unknown} left Unknown\n")
    if n == 0:
        sys.stderr.write("[dfam_to_rmlib] ERROR: no FASTA records parsed — the input did "
                         "not decode to FASTA (JSON envelope / gzip / HTML?)\n")
        return 1
    if classmap and unknown == n:
        # Library is still written (usable at total level); flag the join problem.
        sys.stderr.write("[dfam_to_rmlib] WARNING: every family is Unknown despite "
                         "classifications — accession join likely failed\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
