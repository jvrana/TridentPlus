import re

import jdna
import pandas as pd
import primer3
from pydent.browser import Browser
from pydent import AqSession
from collections import namedtuple
from collections import OrderedDict
from warnings import warn
from tridentplus.utils.saved_browser import SavedBrowser
import os
from tqdm import tqdm
from itertools import product
from primer3plus import Primer3Design

reverse_complement = jdna.alphabet.AmbiguousDNA.rc
named_primer = namedtuple('Primer', ['sequence', 'name'])
here = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.join(here, 'data')


def _valid_dna_sequence(seq):
    char = set(seq)
    valid_characters = set(jdna.alphabet.UnambiguousDNA.characters())
    if char.difference(valid_characters):
        return False
    return True


def valid_primers(primers, verbose=False):
    valid_primers = []
    for p in primers:
        if _valid_dna_sequence(p.seq):
            valid_primers.append(p)
        else:
            if verbose:
                warn("Primer {} (id={}) does not have a valid sequence '{}'".format(p.name, p.id, p.seq))
    return valid_primers


def load_primers(browser: Browser, verbose=False) -> list:
    """
    Find and load all primers from Aquarium, grabbing the anneal and overhang sequences, concatenating them,
    and saving to the `seq` attribute for each primer

    :param browser: pydent Browser instance
    :type browser: pydent.browser.Browser
    :return: list of primers with the parsed sequence in the 'seq' field
    :rtype: list
    """
    print("Retrieving primers from Aquarium...")
    primers = browser.where({}, sample_type="Primer")
    print("Retrieving properties for {} primers".format(len(primers)))
    browser.recursive_retrieve(primers, {'sample_type': 'field_types', 'field_values': {}})

    for p in primers:
        props = p.properties
        anneal = props['Anneal Sequence']
        overhang = props['Overhang Sequence']
        if anneal is None:
            anneal = ''
        if overhang is None:
            overhang = ''
        seq = overhang + anneal
        seq = re.sub('\s+', '', seq)
        p.seq = seq

    # filter primers

    return valid_primers(primers, verbose)


def cache_primers(session, filepath=None, overwrite=False) -> list:
    if filepath is None:
        filepath = os.path.join(data_dir, 'primers.pkl')
    with SavedBrowser(filepath, session) as browser:
        if browser.model_cache and not overwrite:
            primer_type = browser.cached_where('SampleType', {'name': 'Primer'})[0]
            primers = browser.cached_where('Sample', {'sample_type_id': primer_type.id})
        else:
            primers = load_primers(browser)
    return valid_primers(primers)


Binding = namedtuple('Binding', ['matchseq', 'primerseq', 'primer', 'direction'])
Match = namedtuple('Match', ['anneal', 'overhang', 'start', 'end'])


def find_initial_bindings(primers: list, template: str, min_bases: int) -> list:
    """
    Quickly finds initial bindings to primers

    :param primers: list of loaded primers (with sequence in the 'seq' attribute) or a tuple of the (sequence, name)
    :type primers: list
    :param template: template sequence
    :type template: basestring
    :param min_bases: minimum number of bases to find matches
    :type min_bases: int
    :return: list of matches at namedtuples (matching_sequence, full primer sequence, primer model, 1 or -1 for direction)
    :rtype: list
    """

    jseq = jdna.Sequence(template)
    template = str(jseq)
    template_rc = str(jseq.copy().rc())
    assert isinstance(primers, list)
    matches = []

    for primer in primers:
        if isinstance(primer, tuple):
            primer = named_primer(*primer)
        try:
            test_seq = primer.seq[min_bases:]
            if len(test_seq) < min_bases:
                continue
            fwd_matches = re.findall(test_seq, template, re.IGNORECASE)
            rev_matches = re.findall(test_seq, template_rc, re.IGNORECASE)
            for _m in fwd_matches:
                matches.append(Binding(test_seq, primer.seq, primer, 1))
            for _m in rev_matches:
                matches.append(Binding(test_seq, primer.seq, primer, -1))
        except Exception as e:
            raise e
    return list(set(matches))


def _extend_match(matchseq: str, primer: str, template: str) -> list:
    """
    Extends a primer binding match

    :param matchseq: initial matching sequence
    :type matchseq: basestring
    :param primer: complete primer sequence
    :type primer: basestring
    :param template: template
    :type template: basestring
    :return: list of matches [(annealing, overhang, start of annealing, and end of annealing), ...]
    :rtype: list
    """
    matches = []
    for match in re.finditer(matchseq.upper(), template.upper()):
        e = match.end()
        i = 1
        __i = i
        tmp_seq = template[e - i:e]
        tmp_matched = primer[-i:]

        __seq = tmp_seq
        __matched = tmp_matched
        while (i <= len(primer) and tmp_seq.upper() == tmp_matched.upper()):
            __i = i
            __seq = tmp_seq
            __matched = tmp_matched
            i += 1
            tmp_seq = template[e - i:e]
            tmp_matched = primer[-i:]
        annealing = primer[-__i:]
        overhang = primer[:-__i]
        matches.append(Match(annealing, overhang, e - __i, e))
    return matches


def primer_bindings(primers: list, template: str, min_bases=10) -> pd.DataFrame:
    """
    Generate a primer binding dataframe from a list of primers

    :param primers: list of loaded primers (with sequence in the 'seq' attribute)
    :type primers: list
    :param template: template sequence
    :type template: basestring
    :return: data frame of primer binding sites
    :rtype: pandas.DataFrame
    """

    bindings = find_initial_bindings(primers, template, min_bases)

    rows = []

    for binding in bindings:
        t = template
        if binding.direction == -1:
            t = reverse_complement(t)
        matches = _extend_match(binding.matchseq, binding.primerseq, t)
        for match in matches:
            if binding.direction == 1:
                abs_start = match.start
                abs_end = match.end
            else:
                abs_start = len(template) - match.start
                abs_end = len(template) - match.end
            row = OrderedDict()
            row['name'] = binding.primer.name
            row['sequence'] = binding.primerseq
            row['direction'] = binding.direction
            row['overhang'] = match.overhang
            row['annealing'] = match.anneal
            row['start'] = match.start
            row['end'] = match.end
            row['abs_start'] = abs_start
            row['abs_end'] = abs_end
            row['Tm'] = round(primer3.calcTm(match.anneal[-60:].upper(), dv_conc=15), 2)
            row['match'] = binding.matchseq
            rows.append(row)
    return rows


def primer_binding_df(primers: list, template: str, min_bases=10) -> list:
    """
    Generate a primer binding dataframe from a list of primers

    :param primers: list of loaded primers (with sequence in the 'seq' attribute)
    :type primers: list
    :param template: template sequence
    :type template: basestring
    :return: list of OrderedDict
    :rtype: list
    """

    rows = primer_bindings(primers, template, min_bases=min_bases)
    if rows:
        return pd.DataFrame(rows, columns=rows[0].keys())


def pick_primers(session: AqSession, template: str, size_range: tuple, target=None, primer_list=None, max_anneal_len=36, ):
    if primer_list is None:
        primer_list = cache_primers(session)
    bindings = primer_bindings(primer_list, str(template))
    bindings_by_annealing = {}
    for b in bindings:
        bindings_by_annealing.setdefault(b['annealing'][-max_anneal_len:].upper(), []).append(b)
    fwd = [b['annealing'][-max_anneal_len:] for b in bindings if b['direction'] == 1]
    rev = [b['annealing'][-max_anneal_len:] for b in bindings if b['direction'] == -1]

    print("Found {} forward bindings and {} reverse bindings".format(len(fwd), len(rev)))

    pairs = list(product(fwd, rev))
    designer = Primer3Design()
    all_results = []
    reasons = []
    for f, r in tqdm(pairs, ):
        results = designer.check_pcr_primers(str(template), f, r, size_range=size_range, target=target,
                                             max_iterations=15)
        if results[0]:
            all_results += list(results[0].values())
        else:
            reasons.append(results[1])
    all_results = sorted(all_results, key=lambda x: x['PAIR']['PENALTY'])

    # add original binding information
    for r in all_results:
        left_bindings = bindings_by_annealing[r['LEFT']['SEQUENCE'].upper()]
        right_bindings = bindings_by_annealing[r['RIGHT']['SEQUENCE'].upper()]
        r['LEFT']['bindings'] = left_bindings
        r['RIGHT']['name'] = right_bindings
    return all_results