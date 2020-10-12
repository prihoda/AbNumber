from abnumber import Chain


def test_chain_from_str():
    var = 'ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK'
    tail = 'RTV'
    chain = Chain.from_str(var + tail)
    assert chain.seq == var

    expected_format = '''
    ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK
                          ^^^^^^                 ^^^                                    ^^^^^^^^^          
    '''.strip()

    assert chain.format() == expected_format

    assert chain.cdr1_seq == 'CRASQG'
    assert chain.cdr2_seq == 'LLI'
    assert chain.cdr3_seq == 'IYYCQQFNS'
