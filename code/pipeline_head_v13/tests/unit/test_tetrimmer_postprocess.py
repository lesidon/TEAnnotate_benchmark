"""bin/tetrimmer_postprocess.py: the user's records leave TEtrimmer as given (round-3
T03), inputs kept unassessed by the gate are re-injected (T05), and the output lineage
names each consensus's input. Unknown-labelled consensi skip the label re-check, so
no blastn is needed."""
from conftest import headers

S = 'ACGT' * 100


def setup(tmp_path):
    (tmp_path / 'input.fa').write_text(
        f'>u1#LTR/Gypsy\n{S}\n>d1#Unknown\n{S}\n>d2#Unknown\n{S}\n>d3#Unknown\n{S}\n>d4#Unknown\n{S}\n')
    t = tmp_path / 'tt'
    t.mkdir()
    # TEtrimmer rebuilt u1 (the user's record) and d1; d2-d4 have no consensus
    (t / 'TEtrimmer_consensus.fasta').write_text(
        f'>u1_cons#Unknown\tTEtrimmer\tReliable\n{S}\n>d1_cons#Unknown\tTEtrimmer\tReliable\n{S}\n')
    (t / 'Sequence_name_mapping.txt').write_text(
        'original_name\tmodified_name\nu1#LTR/Gypsy\tseq1\nd1#Unknown\tseq2\n')
    (t / 'summary.txt').write_text('input_name,output_name\nseq1,u1_cons\nseq2,d1_cons\n')
    (tmp_path / 'evidence.tsv').write_text(
        'name\ttier\taction\n'
        'u1\tT\tkeep\nd1\t3\tkeep\nd2\t4\tkeep_unassessed\nd3\t3\tkeep\nd4\t2\tkeep\n')
    (tmp_path / 'lineage.tsv').write_text(
        'id\torigin\tparent\tsource_id\n'
        'u1\tuser\t-\tu1_original\nd1\tdenovo\t-\td1\nd2\tdenovo\t-\td2\n'
        'd3\tdenovo\t-\td3\nd4\tdenovo\t-\td4\n')


def test_user_records_unassessed_inputs_and_lineage(run_script, tmp_path):
    setup(tmp_path)
    run_script('tetrimmer_postprocess', '--tetrimmer-dir', tmp_path / 'tt', '--input', tmp_path / 'input.fa',
               '--evidence', tmp_path / 'evidence.tsv', '--reinject-tiers', '1,2,T',
               '--lineage', tmp_path / 'lineage.tsv', '--out', tmp_path / 'out.fa',
               '--out-user', tmp_path / 'user.fa', '--out-lineage', tmp_path / 'out_lineage.tsv',
               '--report', tmp_path / 'rep.tsv')
    # d1's consensus; d2 (gate kept it unassessed) and d4 (tier 2) re-injected; d3 (tier 3) not;
    # the consensus rebuilt from the user's u1 dropped
    assert headers(tmp_path / 'out.fa') == ['d1_cons#Unknown', 'd2#Unknown', 'd4#Unknown']
    # the user's record, exactly as given
    assert (tmp_path / 'user.fa').read_text() == f'>u1#LTR/Gypsy\n{S[:60]}\n' + '\n'.join(
        S[i:i + 60] for i in range(60, len(S), 60)) + '\n'
    rows = [l.rstrip('\n').split('\t') for l in open(tmp_path / 'out_lineage.tsv')][1:]
    lin = {r[0]: r[1:] for r in rows}
    assert lin['d1_cons'] == ['denovo', 'd1', 'd1']
    assert lin['d2'] == ['denovo', '-', 'd2']
    assert lin['u1'] == ['user', '-', 'u1_original']
    rep = (tmp_path / 'rep.tsv').read_text()
    assert 'drop_user_rebuilt' in rep and 'kept unassessed' in rep


def test_without_lineage_user_records_are_reinjected_as_before(run_script, tmp_path):
    setup(tmp_path)
    run_script('tetrimmer_postprocess', '--tetrimmer-dir', tmp_path / 'tt', '--input', tmp_path / 'input.fa',
               '--evidence', tmp_path / 'evidence.tsv', '--reinject-tiers', '1,2,T',
               '--out', tmp_path / 'out.fa', '--report', tmp_path / 'rep.tsv')
    assert headers(tmp_path / 'out.fa') == ['u1_cons#Unknown', 'd1_cons#Unknown', 'd2#Unknown', 'd4#Unknown']
