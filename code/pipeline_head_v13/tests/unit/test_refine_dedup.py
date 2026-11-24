"""bin/refine_dedup.py: a user-library sequence stays its CD-HIT cluster's
representative (review item 8); chimera pieces get unused IDs and a lineage row
(round-2 S03, S01)."""
import argparse
import subprocess
import sys

from conftest import BIN, headers, load_script

A, C, G = 'A' * 240, 'C' * 244, 'G' * 300


def setup(tmp_path):
    (tmp_path / 'combined.fa').write_text(
        f'>rnd-1_family-1#LTR/Gypsy\n{A}\n>other#DNA/hAT\n{G}\n'
        f'>rnd-1_family-1_user#LTR/Copia\n{C}\n')
    (tmp_path / 'dedup.fa').write_text(f'>rnd-1_family-1#LTR/Gypsy\n{A}\n>other#DNA/hAT\n{G}\n')
    # CD-HIT kept the de novo sequence; the user's is a member of its cluster.
    (tmp_path / 'clstr').write_text(
        '>Cluster 0\n'
        '0\t240nt, >rnd-1_family-1#LTR/Gypsy... *\n'
        '1\t244nt, >rnd-1_family-1_user#LTR/Copia... at 1:240:1:244/+/99.00%\n'
        '>Cluster 1\n'
        '0\t300nt, >other#DNA/hAT... *\n')
    (tmp_path / 'lineage.tsv').write_text(
        'id\torigin\tparent\tsource_id\n'
        'rnd-1_family-1\tdenovo\t-\trnd-1_family-1\n'
        'other\tdenovo\t-\tother\n'
        'rnd-1_family-1_user\tuser\t-\trnd-1_family-1\n')


def refine(run_script, tmp_path, *extra):
    run_script('refine_dedup', '--combined', tmp_path / 'combined.fa', '--dedup', tmp_path / 'dedup.fa',
               '--clstr', tmp_path / 'clstr', '--out', tmp_path / 'out.fa',
               '--report', tmp_path / 'rep.tsv', '--no-split', *extra)
    return headers(tmp_path / 'out.fa')


def test_user_sequence_becomes_representative(run_script, tmp_path):
    setup(tmp_path)
    out = refine(run_script, tmp_path, '--lineage', tmp_path / 'lineage.tsv',
                 '--out-lineage', tmp_path / 'out_lineage.tsv')
    assert out == ['rnd-1_family-1_user#LTR/Copia', 'other#DNA/hAT']
    assert 'user_representative' in (tmp_path / 'rep.tsv').read_text()
    assert 'rnd-1_family-1_user\tuser\t-\trnd-1_family-1' in (tmp_path / 'out_lineage.tsv').read_text()


def test_without_lineage_cd_hit_choice_stands(run_script, tmp_path):
    setup(tmp_path)
    assert refine(run_script, tmp_path) == ['rnd-1_family-1#LTR/Gypsy', 'other#DNA/hAT']


def test_duplicate_ids_stop_the_run(tmp_path):
    setup(tmp_path)
    (tmp_path / 'combined.fa').write_text(f'>x#LTR/Gypsy\n{A}\n>x#DNA/hAT\n{G}\n')
    r = subprocess.run([sys.executable, str(BIN / 'refine_dedup.py'), '--combined', str(tmp_path / 'combined.fa'),
                        '--dedup', str(tmp_path / 'dedup.fa'), '--clstr', str(tmp_path / 'clstr'),
                        '--out', str(tmp_path / 'out.fa'), '--report', str(tmp_path / 'rep.tsv'), '--no-split'],
                       capture_output=True, text=True)
    assert r.returncode != 0 and 'duplicate IDs' in r.stderr


def test_lineage_reader_keeps_a_sequence_called_id(tmp_path):
    """Round-3 T02: only the first line is the header."""
    rd = load_script('refine_dedup')
    p = tmp_path / 'lin.tsv'
    p.write_text('id\torigin\tparent\tsource_id\nid\tuser\t-\tid\nx\tdenovo\tid\tid\n')
    lin = rd.read_lineage(p)
    assert lin['id']['origin'] == 'user' and lin['x']['parent'] == 'id'


def test_chimera_pieces_get_unused_ids_and_a_parent():
    rd = load_script('refine_dedup')
    host, other, donor = 'T' * 1500, 'C' * 900, 'G' * 600
    lib = [('family#LTR/Gypsy', host), ('family_part1#LTR/Gypsy', other), ('donor#LINE/L1', donor)]
    # A 600 bp LINE segment (bp 501-1100) inside the 1500 bp LTR: the LTR is the host.
    rd.blast_pairs = lambda q, s, threads, ident: [('s0', 's2', 95.0, 600, (501, 1100), (1, 600))]
    args = argparse.Namespace(threads=1, min_identity=80.0, min_segment=500, max_foreign_frac=0.5,
                              min_fragment=80, unlabel_conflicts=False)
    lineage = {'family': {'origin': 'user', 'parent': '-', 'source_id': 'fam0'}}
    taken = {'family', 'family_part1', 'donor'}
    out = rd.chimera_step(lib, args, [], rd.defaultdict(int), taken, lineage)
    names = [n for n, _ in out]
    # family_part1 already exists: the pieces become _part2 and _part3, not a second _part1
    assert names == ['family_part2#LTR/Gypsy', 'family_part3#LTR/Gypsy',
                     'family_part1#LTR/Gypsy', 'donor#LINE/L1']
    assert lineage['family_part2'] == {'origin': 'user', 'parent': 'family', 'source_id': 'fam0'}
    assert 'family_part1' not in lineage          # the pre-existing sequence gains no parent
