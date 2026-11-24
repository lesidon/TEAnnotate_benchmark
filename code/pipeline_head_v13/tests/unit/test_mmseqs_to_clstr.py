"""bin/mmseqs_to_clstr.py: CD-HIT's greedy rule on an all-vs-all MMseqs2 hit table, on
both strands, written in CD-HIT's .clstr format."""
from conftest import headers

LONG, MID, SHORT, OTHER = 'A' * 1000, 'C' * 900, 'G' * 500, 'T' * 950
# MMseqs2 names a sequence by its whole first header word, label included
ID = {'short': 'short#DNA', 'long': 'long#LTR/Gypsy', 'mid': 'mid#LTR/Gypsy', 'other': 'other#LINE'}


def setup(tmp_path, hits):
    (tmp_path / 'lib.fa').write_text(
        f'>short#DNA\n{SHORT}\n>long#LTR/Gypsy\n{LONG}\n>mid#LTR/Gypsy\n{MID}\n>other#LINE\n{OTHER}\n')
    # query target fident qstart qend qlen tstart tend tlen
    (tmp_path / 'hits.tsv').write_text(
        ''.join('\t'.join(map(str, (ID[h[0]], ID[h[1]]) + h[2:])) + '\n' for h in hits))


def run(run_script, tmp_path, *extra):
    run_script('mmseqs_to_clstr', '--fasta', tmp_path / 'lib.fa', '--hits', tmp_path / 'hits.tsv',
               '--out', tmp_path / 'dedup.fa', *extra)
    return headers(tmp_path / 'dedup.fa'), (tmp_path / 'dedup.fa.clstr').read_text()


def test_minus_strand_member_joins_with_its_coordinates(run_script, tmp_path):
    # mid on the minus strand of long, reported only as long -> mid
    setup(tmp_path, [('long', 'mid', 0.9, 51, 950, 1000, 900, 1, 900)])
    reps, clstr = run(run_script, tmp_path)
    assert reps == ['short#DNA', 'long#LTR/Gypsy', 'other#LINE']
    assert '900nt, >mid#LTR/Gypsy... at 1:900:51:950/-/90.00%' in clstr


def test_best_earlier_representative_wins(run_script, tmp_path):
    # short hits long (85%) and other (95%); other is a representative before it, so it wins
    setup(tmp_path, [('short', 'long', 0.85, 1, 500, 500, 1, 500, 1000),
                     ('short', 'other', 0.95, 1, 500, 500, 101, 600, 950)])
    reps, clstr = run(run_script, tmp_path)
    assert 'short#DNA' not in reps
    block = clstr.split('>Cluster')[3]
    assert '>other#LINE... *' in block and '>short#DNA... at 1:500:101:600/+/95.00%' in block


def test_member_of_a_member_stays_a_representative_if_unlinked(run_script, tmp_path):
    # mid joins long; short aligns only to mid (not a representative) so it stands alone
    setup(tmp_path, [('mid', 'long', 0.9, 1, 900, 900, 1, 900, 1000),
                     ('short', 'mid', 0.99, 1, 500, 500, 1, 500, 900)])
    reps, _ = run(run_script, tmp_path)
    assert reps == ['short#DNA', 'long#LTR/Gypsy', 'other#LINE']


def test_identity_and_coverage_thresholds(run_script, tmp_path):
    setup(tmp_path, [('mid', 'long', 0.79, 1, 900, 900, 1, 900, 1000),   # identity too low
                     ('short', 'long', 0.9, 1, 350, 500, 1, 350, 1000)])  # covers 70% of short
    reps, _ = run(run_script, tmp_path)
    assert len(reps) == 4
    reps, _ = run(run_script, tmp_path, '--identity', '0.75', '--coverage', '0.7')
    assert reps == ['long#LTR/Gypsy', 'other#LINE']
