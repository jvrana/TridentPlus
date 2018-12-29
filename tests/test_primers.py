import pytest
import tridentplus
import os
from tridentplus.utils import saved_browser
import jdna


@pytest.fixture(scope='session')
def gfp():
    return 'atggccgatgatgaagttgccgccctcgctgcagccccggtagaaaaaatgagtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttctgttatggtgttcaatgcttttcaagatacccagatcatatgaaacggcatgactttttcaagagtgccatgcccgaaggttatgtacaggaaagaactatatttttcaaagatgacgggaactacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaaggtattgattttaaagaagatggaaacattcttggacacaaattggaatacaactataactcacacaatgtatacatcatggcagacaaacaaaagaatggaatcaaagttaacttcaaaattagacacaacattgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttcgaaagatcccaacgaaaagagagaccacatggtccttcttgagtttgtaacagctgctgggattacacatggcatggatgaactatacaaatag'


@pytest.fixture(scope='function')
def fake_primer():
    class FakePrimer():
        def __init__(self, seq, name='fake'):
            self.name = name
            self.seq = seq

    return FakePrimer


@pytest.fixture(scope='session')
def primers(session, datadir):
    filepath = os.path.join(datadir, 'browser.pkl')

    with saved_browser.SavedBrowser(filepath, session) as browser:
        if browser.model_cache:
            primer_type = browser.cached_where('SampleType', {'name': 'Primer'})[0]
            primers = browser.cached_where('Sample', {'sample_type_id': primer_type.id})
        else:
            primers = tridentplus.primers.load_primers(browser)

    return primers


def test_find_init_bindings(fake_primer, gfp):

    p1 = fake_primer(gfp[10:50].lower(), name='p1')
    p2 = fake_primer(gfp[100:120].upper(), name='p2')
    p3 = fake_primer(jdna.alphabet.AmbiguousDNA.rc(p2.seq), name='p2_rc')

    template = gfp[:200] + p1.seq + gfp[200:]

    primers = [p1, p2, p3]

    bindings = tridentplus.primers.find_initial_bindings(primers, template, 10)
    for b in bindings:
        print(b.primer.name)
    assert len(bindings) == 3


def test_primer_df(fake_primer, gfp):

    p1 = fake_primer(gfp[10:50].lower(), name='p1')
    p2 = fake_primer(gfp[100:120].upper(), name='p2')
    p3 = fake_primer(jdna.alphabet.AmbiguousDNA.rc(p2.seq), name='p2_rc')

    template = gfp[:200] + p1.seq + gfp[200:]

    primers = [p1, p2, p3]

    df = tridentplus.primers.primer_binding_df(primers, template)
    assert len(df) == 4
