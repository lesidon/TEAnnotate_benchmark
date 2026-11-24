"""Shared style and data loading for the manuscript figures (matplotlib, print/light only).

Colours: a winter palette (midnight navy, glacier blue, aurora violet, slate, frost silver) on cool
blue-grey neutrals. All five series colours pass the all-pairs colour-vision checks of
palette_check.py (worst CVD dE 11.5, glacier vs violet; worst normal-vision dE 16.1, midnight vs
slate). Every series also carries a distinct marker or hatch, so identity never rests on colour alone.
"""
import csv, os
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FIGDIR = os.path.dirname(HERE)

INK, INK2, MUTED = '#0d1b2a', '#465263', '#788495'   # text: 17.4, 7.9 and 3.8:1 on white
GRID, AXIS = '#e2e8ef', '#b8c3cf'
MIDNIGHT, GLACIER, AURORA = '#123666', '#4ba5d2', '#8b62d4'
SLATE, FROST = '#555f6f', '#bcc4cd'

# method -> (label, colour, marker, filled, size)
METHODS = {
    'teannotate':    ('TEAnnotate', MIDNIGHT, 'o', True, 34),
    'teannotate_wo': ('TEAnnotate, no TEtrimmer', MIDNIGHT, 'o', False, 34),
    'earlgrey':      ('Earl Grey', GLACIER, 's', True, 24),
    'hite':          ('HiTE', AURORA, '^', True, 30),
    'edta':          ('EDTA', SLATE, 'D', True, 20),
    'heliano':       ('HELIANO', FROST, 'v', True, 24),
}
TOOL_STYLE = {  # discovery tools in Fig. 4, same identity as Fig. 2
    'earlgrey': METHODS['earlgrey'], 'hite': METHODS['hite'], 'edta': METHODS['edta'], 'heliano': METHODS['heliano'],
}

# development genomes first, then held-out
ORDER = ['arabidopsis_thaliana', 'oryza_sativa', 'drosophila_melanogaster', 'caenorhabditis_elegans',
         'zymoseptoria_tritici', 'septoria_linicola',
         'chlamydomonas_reinhardtii', 'zea_mays', 'caenorhabditis_briggsae', 'gallus_gallus', 'taeniopygia_guttata',
         'mus_musculus', 'danio_rerio', 'blumeria_hordei', 'ectocarpus_sp7']
PENDING = set()


def italic(display):
    """Matplotlib mathtext label with the binomial in italics ('Ectocarpus sp. 7' keeps 'sp. 7' upright)."""
    if display.startswith('Ectocarpus'):
        return r'$\it{Ectocarpus}$ sp. 7'
    return r'$\it{' + display.replace(' ', r'\ ') + '}$'


def setup():
    mpl.rcParams.update({
        'font.family': 'Arial', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7,
        'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5,
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
        'axes.edgecolor': AXIS, 'axes.linewidth': 0.6, 'axes.labelcolor': INK, 'text.color': INK,
        'xtick.color': INK2, 'ytick.color': INK2, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.5, 'axes.axisbelow': True,
        'axes.spines.top': False, 'axes.spines.right': False,
        'legend.frameon': False, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 400,
    })


def read(name):
    with open(os.path.join(DATA, name), newline='') as fh:
        return list(csv.DictReader(fh))


def num(x):
    return float(x) if x not in (None, '') else None


def save(fig, stem):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(FIGDIR, f'{stem}.{ext}'), bbox_inches='tight', pad_inches=0.02)
    print('wrote', stem)
