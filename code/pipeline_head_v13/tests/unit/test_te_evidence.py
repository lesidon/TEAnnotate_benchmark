"""bin/te_evidence.py: genomic copies counted per locus (review R10), trust and tool
from the lineage only (R11, round-2 S01), and a failed domain search (R12)."""
import pytest

SEQ = 'ACGT' * 250   # 1000 bp


def hit(q, s, a, b):
    return f'{q}\t{s}\t98.0\t1000\t1\t1000\t1000\t{a}\t{b}\n'


@pytest.fixture
def evidence(run_script, tmp_path):
    (tmp_path / 'lib.fa').write_text(
        f'>dupHSP#Unknown\n{SEQ}\n'               # two identical HSP rows, one locus
        f'>palin#DNA\n{SEQ}\n'                    # forward + reverse HSP at one locus
        f'>twoLoci#Unknown\n{SEQ}\n'              # two distinct loci on one chromosome
        f'>rnd-1_family-1#Unknown\n{SEQ}\n'       # de novo, same ID as a user family
        f'>rnd-1_family-1_user#LTR/Gypsy\n{SEQ}\n'   # the user family after CD_HIT_EST's rename
        f'>foo_part1#DNA/hAT\n{SEQ}\n'            # a de novo family that happens to be called foo_part1
        f'>bar_part2#DNA/hAT\n{SEQ}\n'            # a piece actually cut from the user family bar
        f'>rnd-5_family-3_part1#LTR\n{SEQ}\n')    # a piece cut from an Earl Grey consensus
    (tmp_path / 'raw.fa').write_text(
        '>dupHSP#Unknown\nA\n>palin#DNA\nA\n>twoLoci#Unknown\nA\n>rnd-1_family-1#Unknown\nA\n'
        '>foo_part1#DNA/hAT\nA\n>rnd-5_family-3#LTR\nA\n')
    (tmp_path / 'prov.tsv').write_text('foo_part1\thite\nrnd-5_family-3\tearlgrey\n')
    (tmp_path / 'lineage.tsv').write_text(
        'id\torigin\tparent\tsource_id\n'
        'rnd-1_family-1\tdenovo\t-\trnd-1_family-1\n'
        'rnd-1_family-1_user\tuser\t-\trnd-1_family-1\n'
        'foo_part1\tdenovo\t-\tfoo_part1\n'
        'foo\tuser\t-\tfoo\n'
        'bar\tuser\t-\tbar\n'
        'bar_part2\tuser\tbar\tbar\n'
        'rnd-5_family-3\tdenovo\t-\trnd-5_family-3\n'
        'rnd-5_family-3_part1\tdenovo\trnd-5_family-3\trnd-5_family-3\n')
    (tmp_path / 'cls.tsv').write_text('')
    (tmp_path / 'support.tsv').write_text('')
    rows = (hit('dupHSP', 'chr1', 100, 1099) + hit('dupHSP', 'chr1', 100, 1099) +
            hit('palin', 'chr2', 500, 1499) + hit('palin', 'chr2', 1499, 500) +
            hit('twoLoci', 'chr1', 100, 1099) + hit('twoLoci', 'chr1', 5000, 5999))
    (tmp_path / 'genome9.tsv').write_text(rows)
    (tmp_path / 'genome7.tsv').write_text(
        ''.join('\t'.join(l.split('\t')[:7]) + '\n' for l in rows.splitlines()))

    def run(genome='genome9.tsv', *extra):
        run_script('te_evidence', '--library', tmp_path / 'lib.fa', '--raw', tmp_path / 'raw.fa',
                   '--provenance', tmp_path / 'prov.tsv', '--lineage', tmp_path / 'lineage.tsv',
                   '--tesorter-cls', tmp_path / 'cls.tsv',
                   '--support-hits', tmp_path / 'support.tsv', '--genome-hits', tmp_path / genome,
                   '--out-library', tmp_path / 'keep.fa', '--out-nonte', tmp_path / 'nonte.fa',
                   '--out-removed', tmp_path / 'rm.fa', '--evidence', tmp_path / 'ev.tsv', *extra)
        lines = [l.rstrip('\n').split('\t') for l in open(tmp_path / 'ev.tsv')]
        return {r[0]: dict(zip(lines[0], r)) for r in lines[1:]}
    return run


def test_overlapping_hsps_are_one_locus(evidence):
    ev = evidence()
    assert ev['dupHSP']['full_copies'] == '1' and ev['dupHSP']['action'] == 'remove'
    assert ev['palin']['full_copies'] == '1'
    assert ev['twoLoci']['full_copies'] == '2' and ev['twoLoci']['tier'] == '3'


def test_user_family_is_trusted_after_rename(evidence):
    ev = evidence()
    assert ev['rnd-1_family-1_user']['tier'] == 'T'
    assert ev['rnd-1_family-1']['tier'] != 'T'


def test_a_part_suffix_alone_gives_no_trust(evidence):
    ev = evidence()
    assert ev['foo_part1']['tier'] != 'T'          # unrelated to the user's foo
    assert ev['foo_part1']['tool'] == 'hite'
    assert ev['bar_part2']['tier'] == 'T'          # recorded as cut from the user's bar


def test_a_piece_takes_its_parents_tool(evidence):
    assert evidence()['rnd-5_family-3_part1']['tool'] == 'earlgrey'


def test_old_seven_column_table_counts_rows(evidence):
    assert evidence('genome7.tsv')['dupHSP']['full_copies'] == '2'


def test_homology_is_not_a_domain(run_script, tmp_path):
    """Round-3 T01: TEsorter's second pass (Domains 'none') is sequence support (tier 2)
    that may classify an Unknown, never a domain (tier 1) or a label override."""
    (tmp_path / 'lib.fa').write_text(
        f'>homUnk#Unknown\n{SEQ}\n>homConflict#DNA/hAT\n{SEQ}\n>domConflict#DNA/hAT\n{SEQ}\n'
        f'>id#DNA/hAT\n{SEQ}\n')
    (tmp_path / 'cls.tsv').write_text(
        '#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\n'
        'homUnk#Unknown\tLTR\tGypsy\tunknown\tnone\t?\tnone\n'
        'homConflict#DNA/hAT\tLTR\tGypsy\tunknown\tnone\t?\tnone\n'
        'domConflict#DNA/hAT\tLTR\tGypsy\tTekay\tyes\t+\tGAG|Ty3-GAG RT|Ty3-RT\n'
        'id#DNA/hAT\tLTR\tGypsy\tTekay\tyes\t+\tRT|Ty3-RT\n')
    (tmp_path / 'lineage.tsv').write_text('id\torigin\tparent\tsource_id\nid\tuser\t-\tid\n')
    for f in ('raw.fa', 'support.tsv', 'genome.tsv'):
        (tmp_path / f).write_text('')
    run_script('te_evidence', '--library', tmp_path / 'lib.fa', '--raw', tmp_path / 'raw.fa',
               '--lineage', tmp_path / 'lineage.tsv', '--tesorter-cls', tmp_path / 'cls.tsv',
               '--support-hits', tmp_path / 'support.tsv', '--genome-hits', tmp_path / 'genome.tsv',
               '--relabel-unknown', '--relabel-conflicts',
               '--out-library', tmp_path / 'keep.fa', '--out-nonte', tmp_path / 'nonte.fa',
               '--out-removed', tmp_path / 'rm.fa', '--evidence', tmp_path / 'ev.tsv')
    lines = [l.rstrip('\n').split('\t') for l in open(tmp_path / 'ev.tsv')]
    ev = {r[0]: dict(zip(lines[0], r)) for r in lines[1:]}
    assert ev['homUnk']['tier'] == '2' and ev['homUnk']['label_out'] == 'LTR/Gypsy'
    assert ev['homConflict']['tier'] == '2' and ev['homConflict']['label_out'] == 'DNA/hAT'
    assert ev['homConflict']['homology'] == 'LTR/Gypsy' and ev['homConflict']['domain'] == ''
    assert ev['domConflict']['tier'] == '1' and ev['domConflict']['label_out'] == 'LTR/Gypsy'
    # round-3 T02: a user record called "id" is trusted, and keeps its label
    assert ev['id']['tier'] == 'T' and ev['id']['label_out'] == 'DNA/hAT'


def test_failed_domain_search_removes_nothing(evidence, tmp_path):
    ev = evidence('genome9.tsv', '--domains-unassessed')
    assert ev['dupHSP']['action'] == 'keep_unassessed'
    assert (tmp_path / 'rm.fa').read_text() == ''
