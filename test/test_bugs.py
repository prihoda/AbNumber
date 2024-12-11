import pytest
from abnumber import Chain, ChainParseError, Position
from abnumber.germlines import HUMAN_IMGT_IG_V, HUMAN_IMGT_IG_J


def test_imgt_61A():
    assert Position.from_string('61A', 'H', 'imgt') < Position.from_string('61', 'H', 'imgt')
    seq = 'EVQLVESGGGLVQPGGSLRLSCAASGIILDYYPIGWFRQAPGKEREGVAFITNSDDSTIYTNYADSVKGRFTISRDKNSLYLQMNSLRAEDTAVYYCSSKASFLIGKDDQGIDAGEYDYWGQGTMVTVSS'
    chain = Chain(seq, 'imgt')
    assert chain.seq == seq

def test_imgt_33A():
    assert Position.from_string('33A', 'H', 'imgt') < Position.from_string('33', 'H', 'imgt')
    seq = 'EVQLVESGGGLVQPGGSLRLSCAASGIILELVISDYYPIGWFRQAPGKEREGVAFITNSDDSTIYTNYADSVKGRFTISRDKNSLYLQMNSLRAEDTAVYYCSSKASFLIGKDDQGIDAGEYDYWGQGTMVTVSS'
    chain = Chain(seq, 'imgt')
    assert chain.seq == seq


def test_kabat_insertion():
    # QVQLVESGGGLVQAGGSLSLSCAYSDSGRAFATVVMAWFRQPPGKDRDFVAGIRRSTNTYYADSVKGRFTISRDNAKNTVYLHINLKPEDTAVYYCAAXXXXXXXXXXXXXXAWGQGTQVTVSS
    #                                                                                     x

    seq = 'QVQLVESGGGLVQAGGSLSLSCAYSDSGRAFATVVMAWFRQPPGKDRDFVAGIRRSTNTYYADSVKGRFTISRDNAKNTVYLHINLKPEDTAVYYCAAXXXXXXXXXXXXXXAWGQGTQVTVSS'
    chain = Chain(seq, 'kabat')
    assert chain.seq == seq

    seq = 'QVQLVESGGGLVQAGGSLSLSCAYSDSGRAFATVVMAWFRQPPGKDRDFVAGIRRSTNTYYADSVKGRFTISRDNAKNTVYLHINALKPEDTAVYYCAAXXXXXXXXXXXXXXAWGQGTQVTVSS'
    chain = Chain(seq, 'kabat')
    assert chain.seq == seq


def test_kappa_imgt_21():
    from anarci import run_hmmer
    sequence = 'DVVMTQSPLSLPVTLGQPASISCRSSQSLVYSDGNTYLNWFQQRPGQSPRRLIYKVSNRDSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCMQGTFGQGTKVEIK'
    alignment = run_hmmer([('id', sequence)], hmm_database="ALL", hmmerpath="", ncpu=None, bit_score_threshold=80)[0]
    print(len(alignment))
    print('---')
    for pos in alignment:
        print(pos)
    print('---')
    for pos, aa in zip(alignment, sequence):
        print(pos, aa)


def test_light_chain_IMGT_position_21():
    # Check bug from ANARCI 2021.02.04
    # When numbering full Kappa chains, position IMGT 21 contains a gap
    # When numbering V gene only, position IMGT 21 contains an amino acid as expected
    # Test against this by making sure that same numbering is assigned when numbering V gene and VJ genes concatenated
    # https://github.com/oxpig/ANARCI/issues/17
    for germline in HUMAN_IMGT_IG_V['K']['aligned_sequences']:
        v_seq = HUMAN_IMGT_IG_V['K']['aligned_sequences'][germline].replace('-', '')
        first_j_gene = list(HUMAN_IMGT_IG_J['K']['aligned_sequences'].keys())[0]
        j_seq = HUMAN_IMGT_IG_J['K']['aligned_sequences'][first_j_gene].replace('-', '')
        vj_seq = v_seq + j_seq
        try:
            v_chain = Chain(v_seq, 'imgt')
            vj_chain = Chain(vj_seq, 'imgt')
        except Exception as e:
            print(e)
            continue
        v_positions = [str(p) for p in v_chain.positions]
        vj_positions = [str(p) for p in vj_chain.positions]

        len_limit = len(v_seq) - 20
        assert ','.join(v_positions[:len_limit]) == ','.join(vj_positions[:len_limit])
