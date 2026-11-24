"""Colour-vision check of the figure palette (Python port of the dataviz skill's validate_palette.js
CVD and normal-vision Delta E checks: OKLab x100, Machado et al. 2009 severity-1.0 protan/deutan).

    python palette_check.py "#2a78d6,#eb6834,..."   (all pairs; prints the worst pairs)
"""
import itertools, math, sys

MACHADO = {
    'protan': [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
    'deutan': [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
}


def lin(h):
    h = h.lstrip('#')
    s = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    return [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in s]


def oklab(rgb):
    r, g, b = rgb
    l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)


def sim(h, kind):
    r = lin(h)
    if kind is None:
        return r
    M = MACHADO[kind]
    return [min(1, max(0, sum(M[i][k] * r[k] for k in range(3)))) for i in range(3)]


def de(a, b, kind=None):
    return 100 * math.dist(oklab(sim(a, kind)), oklab(sim(b, kind)))


def contrast(a, b='#ffffff'):
    def lum(h):
        r, g, bb = lin(h)
        return 0.2126 * r + 0.7152 * g + 0.0722 * bb
    hi, lo = sorted([lum(a), lum(b)], reverse=True)
    return (hi + 0.05) / (lo + 0.05)


if __name__ == '__main__':
    pal = [c.strip() for c in sys.argv[1].split(',')]
    worst_cvd = min((min(de(a, b, 'protan'), de(a, b, 'deutan')), a, b) for a, b in itertools.combinations(pal, 2))
    worst_n = min((de(a, b), a, b) for a, b in itertools.combinations(pal, 2))
    print(f'worst all-pairs CVD dE {worst_cvd[0]:.1f} ({worst_cvd[1]} vs {worst_cvd[2]}); target >= 8, floor 6')
    print(f'worst all-pairs normal dE {worst_n[0]:.1f} ({worst_n[1]} vs {worst_n[2]}); floor 15')
    for c in pal:
        print(f'  {c} contrast vs white {contrast(c):.2f}')
