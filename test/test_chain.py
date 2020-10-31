from abnumber import Chain


def test_light_chain_from_str():
    var = 'ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK'
    tail = 'RTV'
    chain = Chain(var + tail, scheme='imgt')
    assert chain.seq == var
    assert chain.tail == tail

    expected_format = '''
ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK
                          ^^^^^^                 ^^^                                    ^^^^^^^^^          
    '''

    assert chain.format().strip() == expected_format.strip()

    assert chain.cdr1_seq == 'QGISSA'
    assert chain.cdr2_seq == 'DAS'
    assert chain.cdr3_seq == 'QQFNSYPLT'


def test_heavy_chain_from_str():
    var = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
    tail = 'AKTTAPSVYPLA'
    chain = Chain(var + tail, scheme='imgt')
    assert chain.seq == var
    assert chain.tail == tail

    expected_format = '''
QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                         ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^           
'''

    assert chain.format().strip() == expected_format.strip()

    assert chain.cdr1_seq == 'GYTFTRYT'
    assert chain.cdr2_seq == 'INPSRGYT'
    assert chain.cdr3_seq == 'ARYYDDHYCLDY'
