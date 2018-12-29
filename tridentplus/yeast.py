import os
from itertools import product

import jdna
import pyblast
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intermine.webservice import Service

here = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.join(here, 'data')

service = Service("https://yeastmine.yeastgenome.org:443/yeastmine/service")


def _chromosome_list():
    """Return a list of chromosomes"""
    query = service.new_query("Chromosome")
    query.add_view(
        "primaryIdentifier",
        "featureType", "length"
    )
    return query


def _chromosome_sequence(identifier=None):
    """Return the sequence of a specific chromosome"""
    query = service.new_query("Chromosome")
    query.add_view(
        'primaryIdentifier', 'length', 'featureType', 'sequence.residues'
    )
    if identifier:
        query.add_constraint("primaryIdentifier", "=", identifier, code="A")
    return query


default_genome_path = os.path.join(data_dir, 'SCer.fasta')


def yeast_genome(name=None, overwrite=False):
    """
    Return the full sequence of the yeast genome from intermine. Save it to 'name'. Filepath at 'name' already exists,
    load that file. If 'overwrite', always pull data and overwrite file even if it exists.
    """
    if name is None:
        name = default_genome_path
    filepath = os.path.abspath(name)
    if overwrite or not os.path.isfile(filepath):
        rows = _chromosome_sequence().rows()
        seqrecords = []
        for row in rows:
            seq = SeqRecord(Seq(row['sequence.residues']), id=row['primaryIdentifier'])
            seqrecords.append(seq)
        SeqIO.write(seqrecords, name, format='fasta')
        return seqrecords, name
    else:
        return list(SeqIO.parse(name, format='fasta')), name


def yeast_genome_json():
    """Return json object representing the yeast genome"""
    seqs, filepath = yeast_genome()
    return [{
        'name': seq.id,
        'bases': str(seq.seq),
        'circular': False
    } for seq in seqs]


def chromosome(iden):
    for seq in yeast_genome()[0]:
        if seq.id == iden:
            return seq


def align(queries: list, subjects: list):
    query_json = [{
        'bases': q,
        'name': 'query_{}'.format(i),
        'circular': False
    } for i, q in enumerate(queries)]

    subject_json = [{
        'bases': q,
        'name': 'subject_{}'.format(i),
        'circular': False
    } for i, q in enumerate(subjects)]

    aligner = pyblast.JSONBlast(query_json=query_json, subject_json=subject_json, span_origin=False)
    aligner.quick_blastn()
    return aligner


def align_to_genome(sequences: list):
    """Preform a blast alignment to the yeast genome """
    assert isinstance(sequences, list)

    query_json = []

    for seq in sequences:
        if isinstance(seq, str):
            query_json.append({
                'bases': str(seq),
                'name': 'query',
                'circular': False
            })
        elif isinstance(seq, jdna.Sequence):
            query_json.append({
                'bases': str(seq),
                'name': seq.name,
                'circular': False
            })
        else:
            raise TypeError("Type {} not recognized".format(type(seq)))

    aligner = pyblast.JSONBlast(query_json=query_json, subject_json=yeast_genome_json(), span_origin=False)
    aligner.quick_blastn()
    return aligner


def alignments_to_integration_sites(alignments):
    # group alignments by chromosome
    group_by_subject_name = {}
    for align in alignments:
        group_by_subject_name.setdefault(align['subject']['name'], []).append(align)

    sites = []

    for subject_name, aligns in group_by_subject_name.items():

        # pair up alignments
        pairs = list(product(aligns, aligns))

        # filter out pairs that use different queries
        pairs = [(p1, p2) for p1, p2 in pairs if p1['query']['sequence_id'] == p2['query']['sequence_id']]

        # filter out pairs that have same start and end
        pairs = [(p1, p2) for p1, p2 in pairs if p1['query']['start'] != p2['query']['start'] and
                 p1['query']['end'] != p2['query']['end']]

        for a1, a2 in pairs:
            s1 = a1['subject']
            s2 = a2['subject']
            q1 = a1['query']
            q2 = a2['query']

            print('{} {} {} -> {} {} {}'.format(s1['start'], s1['end'], s1['strand'], s2['start'], s2['end'],
                                                s2['strand']))
            print('{} {} {} -> {} {} {}'.format(q1['start'], q1['end'], q1['strand'], q2['start'], q2['end'],
                                                q2['strand']))

            if q1['sequence_id'] != q2['sequence_id']:
                print('REASON: diff queries')
                continue
            if q1['start'] == q2['start']:
                print('REASON: same qstart')
                continue
            if s1['strand'] != s2['strand']:
                print('REASON: diff strands')
                continue
            if q1['end'] > q2['start']:
                print("REASON: qend > qstart")
                continue

            direction = None
            s1_max = max([s1['start'], s1['end']])
            s1_min = min([s1['start'], s1['end']])
            s2_max = max([s2['start'], s2['end']])
            s2_min = min([s2['start'], s2['end']])

            if s1['strand'] == 'plus' and s1_max <= s2_min:
                direction = s2_min - s1_max
            elif s1['strand'] == 'minus' and s2_max <= s1_min:
                direction = s2_max - s1_min
            if direction is not None:
                sites.append((direction, a1, a2))
    return sites


def _aligner_to_integrations(aligner, integrant=None, template=None, min_homology=100, flanking_bps=500,
                             max_distance=5000):
    alignments = aligner.results.alignments
    sites = alignments_to_integration_sites(alignments)
    valid_sites = []
    invalid_sites = []
    for site in sites:
        if site[1]['meta']['alignment_length'] < min_homology:
            invalid_sites.append({
                'REASON': 'homology 1 less than min homology ({} < {})'.format(site[1]['alignment_length'],
                                                                               min_homology),
                'SITE': site
            })
            continue
        if site[2]['meta']['alignment_length'] < min_homology:
            invalid_sites.append({
                'REASON': 'homology 2 less than min homology ({} < {})'.format(site[2]['alignment_length'],
                                                                               min_homology),
                'SITE': site
            })
            continue
        if abs(site[0]) > max_distance:
            invalid_sites.append({
                'REASON': 'max distance between homologies is too great ({} > {})'.format(abs(site[0]), max_distance),
                'SITE': site
            })
            continue
        valid_sites.append(site)

    integrations = []
    for site in valid_sites:
        if integrant is None:
            integrant = aligner.seq_dict[site[1]['query']['sequence_id']]['bases']
        if template is None:
            template = aligner.seq_dict[site[1]['subject']['sequence_id']]['bases']
        query_range = (site[1]['query']['end'], site[2]['query']['start'])
        integrant_seq = integrant[query_range[0]:query_range[1]]
        direction = 1
        if site[0] > 0:
            left = (site[1]['subject']['start'], site[1]['subject']['end'])
            right = (site[2]['subject']['start'], site[2]['subject']['end'])
        else:
            integrant_seq.rc()
            direction = -1
            left = (site[2]['subject']['end'], site[2]['subject']['start'])
            right = (site[1]['subject']['end'], site[1]['subject']['start'])
        flank_left = (left[0] - flanking_bps, left[0])
        flank_right = (right[0], right[0] + flanking_bps)

        flank_left_seq = template[flank_left[0]:flank_left[1]]
        left_seq = template[left[0]:left[1]]
        right_seq = template[right[0]:right[1]]
        flank_right_seq = template[flank_right[0]:flank_right[1]]

        integration_info = {
            'subject_name': site[1]['subject']['name'],
            'query_name': site[1]['query']['name'],
            'subject': template,
            'query': integrant,
            'subject_pos': '{}:{}-{}'.format(site[1]['subject']['name'], left[0], right[1]),
            'query_pos': '{}:{}-{}'.format(site[1]['query']['name'], query_range[0], query_range[1]),
            'left_homology_pos': '{}:{}-{}'.format(site[1]['subject']['name'], left[0], left[1]),
            'right_homology_pos': '{}:{}-{}'.format(site[1]['subject']['name'], right[0], right[1]),
            'left_homology_range': left,
            'right_homology_range': right,
            'flank_left_range': flank_left,
            'flank_right_range': flank_right,
            'query_range': query_range,
            'sequence': {
                'flank_left': flank_left_seq,
                'left_homology': left_seq,
                'integrant': integrant_seq,
                'right_homology': right_seq,
                'flank_right': flank_right_seq
            },
            'direction': direction,
            'alignment': (site[1], site[2])
        }

        integrations.append(integration_info)
    return integrations


def genomic_integration(integrant: jdna.Sequence, min_homology=100, flanking_bps=500, max_distance=5000):
    alignments = align_to_genome([str(integrant)])
    return _aligner_to_integrations(alignments, integrant=integrant, min_homology=min_homology,
                                    flanking_bps=flanking_bps, max_distance=max_distance)


def homologous_recombination(integrant: jdna.Sequence, template: jdna.Sequence, min_homology=100, flanking_bps=500,
                             max_distance=5000):
    alignments = align([str(integrant)], [str(template)])
    return _aligner_to_integrations(alignments, integrant=integrant, template=template, min_homology=min_homology,
                                    flanking_bps=flanking_bps, max_distance=max_distance)
