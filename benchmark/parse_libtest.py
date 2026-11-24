#!/usr/bin/env python3
"""Parse an EDTA lib-test.py report into a tidy per-category metrics TSV.

lib-test.py writes ``<test>.repeatmasker.lib.report`` containing, per TE category,
a summary line of the form (tab-separated):

    <test>.<category>.lib.report<TAB>sens spec accu prec FDR F1 TP TN FP FN

Categories: ltr, nonltr, line, sine, tir, mite, helitron, others, any.

Output columns:
    kingdom species method scoring_type category
    sensitivity specificity accuracy precision FDR F1 TP TN FP FN
"""
import argparse
import re
import sys

CATEGORIES = ("ltr", "nonltr", "line", "sine", "tir", "mite", "helitron", "others", "any")
# field[0] looks like "<something>.<category>.lib.report"
_CAT_RE = re.compile(r"\.(%s)\.lib\.report\s*$" % "|".join(CATEGORIES), re.IGNORECASE)
_METRIC_COLS = ["sensitivity", "specificity", "accuracy", "precision", "FDR", "F1",
                "TP", "TN", "FP", "FN"]


def _num(tok):
    tok = tok.strip()
    if tok in ("", "NA", "nan", "NaN", "-"):
        return ""
    try:
        f = float(tok)
        return str(int(f)) if f.is_integer() else repr(f)
    except ValueError:
        return ""


def parse_report(path):
    """Yield (category, [10 metric strings]) for each summary line found."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if ".lib.report" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            m = _CAT_RE.search(fields[0])
            if not m:
                continue
            cat = m.group(1).lower()
            vals = [_num(x) for x in fields[1:11]]
            vals += [""] * (10 - len(vals))       # pad if TP..FN absent
            yield cat, vals[:10]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, help="<test>.repeatmasker.lib.report")
    ap.add_argument("--species", required=True)
    ap.add_argument("--kingdom", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--out", required=True, help="output tidy metrics TSV")
    args = ap.parse_args(argv)

    rows = list(parse_report(args.report))
    if not rows:
        sys.stderr.write(
            f"[parse_libtest] WARNING: no category summary lines in {args.report}\n")

    header = ["kingdom", "species", "method", "scoring_type", "category"] + _METRIC_COLS
    with open(args.out, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for cat, vals in rows:
            out.write("\t".join([args.kingdom, args.species, args.method, "gold", cat]
                                + vals) + "\n")
    sys.stderr.write(f"[parse_libtest] {args.method}/{args.species}: "
                     f"{len(rows)} categories -> {args.out}\n")


if __name__ == "__main__":
    main()
