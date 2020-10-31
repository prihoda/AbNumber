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


def test_heavy_chain_slice():
    seq = 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS'
    chain = Chain(seq, scheme='imgt')

    assert str(chain[:'5']) == 'EVQLQ'
    assert str(chain['H111':]) == 'DDERGHYCLDYWGQGTTLTVSS'
    assert str(chain['3':'2']) == ''
    assert str(chain['1000':'2000']) == ''

    expected_slice = ', '.join([f'{p}={aa}' for p, aa in chain['111':'112']])
    assert expected_slice == 'H111=D, H111A=D, H111B=E, H112B=R, H112A=G, H112=H'


def test_heavy_alignment_slice():
    seq = 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS'
    chain = Chain(seq, scheme='imgt')
    seq2 = 'QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYDDYLDRWGQGTTLTVSSAKTTAP'
    chain2 = Chain(seq2, scheme='imgt')
    alignment = chain.align(chain2)

    expected_format = '''
EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPS-RGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS
+|||.||||||.||||+|||||||||||.||||||||||||||||+||||||||.|.||||||||||||||||||||||||||.+||||||||||||||||..........||.|||||||||||
QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYD-------DYLDRWGQGTTLTVSS
                         ^^^^^^^^                 ^^^^^^^^^                                      ^^^^^^^^^^^^^^^^^
'''
    assert expected_format.strip() == alignment.format().strip()

    assert alignment['111B'] == ('E', '-')
    assert alignment.raw[105] == ('E', '-')
