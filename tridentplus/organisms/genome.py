"""
Genome services
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tridentplus.primers import pick_pcr_primers
from pydent import AqSession
from itertools import product

import jdna
import pyblast


class BlastUtils(object):
    """Blast utilities"""

    @classmethod
    def align(cls, queries: list, subjects: list):
        """
        Performs an alignment between queries and subjects. Subject sequences will be collated into
        a single FASTA file and local blast installation will build a blastdb in a temporary
        directory and perform a query.

        :param queries: list of query sequences (as strings)
        :type queries: list
        :param subjects: list of subject sequences (as strings) to build a blastdb out of
        :type subjects:
        :return: aligner instance, with results available via `aligner.results` or `aligner.results.alignments`. Sequence
                    dictionary is available via `aligner.seq_dict`. Raw results available via `aligner.raw_results`.
        :rtype: pyblast.JSONBlast
        """
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

    @classmethod
    def alignments_to_integration_sites(cls, alignments):
        """
        Parse alignments (from :class:`pyblast.JSONBlast` instance) into integration sites.

        :param alignments: alignment list
        :type alignments: list
        :return: list of sites compose of tuple of (integration_direction, alignment1, alignment2)
        :rtype: list
        """
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

                if q1['sequence_id'] != q2['sequence_id']:
                    continue
                if q1['start'] == q2['start']:
                    continue
                if s1['strand'] != s2['strand']:
                    continue
                if q1['end'] > q2['start']:
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

    @classmethod
    def aligner_to_integrations(cls, aligner, integrant=None, template=None, min_homology=100, flanking_bps=500,
                                 max_distance=5000):
        """
        Parse a **run** :class:`pyblast.JSONBlast` instance into valid integration sites. Returns the following info
        for each integration:

        .. code-block::

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
                'flanking_bps': flanking_bps,
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


        :param aligner: aligner
        :type aligner: pyblast.JSONBlast
        :param integrant: integrant sequence
        :type integrant: jdna.Sequence
        :param template: template sequence
        :type template: jdna.Sequence
        :param min_homology: minimum number of homologous bps to consider a valid recombination event (default: 100)
        :type min_homology: int
        :param flanking_bps: the number of outside flanking bps to return (default: 500)
        :type flanking_bps: int
        :param max_distance: maxmum distance (bps) between homologous regions to consider a valid recombination event (default: 10000)
        :type max_distance: int
        :return: list of integration events
        :rtype: list
        """

        alignments = aligner.results.alignments
        sites = cls.alignments_to_integration_sites(alignments)
        valid_sites = []
        invalid_sites = []
        for site in sites:
            if site[1]['meta']['alignment_length'] < min_homology:
                invalid_sites.append({
                    'REASON': 'homology 1 less than min homology ({} < {})'.format(site[1]['meta']['alignment_length'],
                                                                                   min_homology),
                    'SITE': site
                })
                continue
            if site[2]['meta']['alignment_length'] < min_homology:
                invalid_sites.append({
                    'REASON': 'homology 2 less than min homology ({} < {})'.format(site[2]['meta']['alignment_length'],
                                                                                   min_homology),
                    'SITE': site
                })
                continue
            if abs(site[0]) > max_distance:
                invalid_sites.append({
                    'REASON': 'max distance between homologies is too great ({} > {})'.format(abs(site[0]),
                                                                                              max_distance),
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
                'flanking_bps': flanking_bps,
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

    @staticmethod
    def integration_sequence(integration_results: dict):
        """Parse an integration result into a :class:`jdna.Sequence` instance."""

        flank_left = jdna.Sequence(str(integration_results['sequence']['flank_left']))
        flank_left.annotate(None, None, integration_results['subject_name'])

        left = jdna.Sequence(str(integration_results['sequence']['left_homology']))
        left.annotate(None, None, integration_results['left_homology_pos'])

        flank_right = jdna.Sequence(str(integration_results['sequence']['flank_right']))
        flank_right.annotate(None, None, integration_results['subject_name'])

        right = jdna.Sequence(str(integration_results['sequence']['right_homology']))
        right.annotate(None, None, integration_results['right_homology_pos'])

        integration_seq = flank_left + left + integration_results['sequence']['integrant'] + right + flank_right

        return integration_seq

    # TODO: allow a provided list of primers to select from.
    # TODO: allow a design of left or right primers
    @classmethod
    def integration_qc_primers(cls, integration_results, primer_list, size_range, session: AqSession, pick_left=False,
                               pick_right=False):

        """Parse an integration result into a :class:`jdna.Sequence` instance and select PCR primers
        from the provided Aquarium session.

        :param integration_results: integration result
        :type integration_results: dict
        :param size_range: size range of pcr products
        :type size_range: tuple
        :param session: session
        :type session: AqSession
        :param pick_left: whether to pick left primers
        :type pick_left: bool
        :param pick_right: whether to pick right primers
        :type pick_right: bool
        :return: (seq, left_qc_primers, right_qc_primers)
        :rtype: tuple
        """

        integration_seq = cls.integration_sequence(integration_results)
        flanking_bps = integration_results['flanking_bps']

        left_qc_primers = []
        right_qc_primers = []
        if pick_left:
            left_homology_len = len(integration_results['sequence']['left_homology'])
            left_target = (flanking_bps, left_homology_len)
            left_qc_primers = pick_pcr_primers(primer_list, integration_seq[:flanking_bps * 2], size_range,
                                               target=left_target)
        if pick_right:
            _integration_seq = integration_seq[-flanking_bps * 2:]
            right_homology_len = len(integration_results['sequence']['right_homology'])
            right_target = (len(_integration_seq) - flanking_bps - right_homology_len, right_homology_len)
            right_qc_primers = pick_pcr_primers(primer_list, _integration_seq, size_range,
                                                target=right_target)
        return integration_seq, left_qc_primers, right_qc_primers


    @classmethod
    def homologous_recombination(cls, integrant: jdna.Sequence, template: jdna.Sequence, min_homology=100, flanking_bps=500,
                                 max_distance=5000):
        """
        Parse a **run** :class:`pyblast.JSONBlast` instance into valid integration sites. Returns the following info
        for each integration:

        .. code-block::

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
                'flanking_bps': flanking_bps,
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


        :param integrant: integrant sequence
        :type integrant: jdna.Sequence
        :param template: template sequence
        :type template: jdna.Sequence
        :param min_homology: minimum number of homologous bps to consider a valid recombination event (default: 100)
        :type min_homology: int
        :param flanking_bps: the number of outside flanking bps to return (default: 500)
        :type flanking_bps: int
        :param max_distance: maxmum distance (bps) between homologous regions to consider a valid recombination event (default: 10000)
        :type max_distance: int
        :return: list of integration events
        :rtype: list
        """

        alignments = cls.align([str(integrant)], [str(template)])
        results = cls.aligner_to_integrations(alignments, integrant=integrant, template=template,
                                           min_homology=min_homology,
                                           flanking_bps=flanking_bps, max_distance=max_distance)
        return results

class Genome(BlastUtils):
    """Genome utilities"""

    genome_name = 'Unknown'
    service_url = "interminewebservice"
    service = None
    default_genome_path = None

    @classmethod
    def chromosomes(cls):
        """Return a list of chromosomes"""
        query = cls.service.new_query("Chromosome")
        query.add_view(
            "primaryIdentifier",
            "featureType", "length"
        )
        return query

    @classmethod
    def _chromosome_sequence(cls, identifier=None):
        """
        Return the sequence of a specific chromosome

        :param identifier: chromosome identifier
        :type identifier: basestring
        :return:
        :rtype:
        """
        query = cls.service.new_query("Chromosome")
        query.add_view(
            'primaryIdentifier', 'length', 'featureType', 'sequence.residues'
        )
        if identifier:
            query.add_constraint("primaryIdentifier", "=", identifier, code="A")
        return query

    @classmethod
    def genome(cls, filepath=None, overwrite=False) -> list:
        """
        Return genome.

        :param filepath: filepath to save the genome
        :type filepath: basestring
        :param overwrite: if True, overwrite existing file
        :type overwrite: bool
        :return: list of SeqRecords
        :rtype: list
        """
        if filepath is None:
            filepath = cls.default_genome_path
        filepath = os.path.abspath(filepath)
        if overwrite or not os.path.isfile(filepath):
            rows = cls._chromosome_sequence().rows()
            seqrecords = []
            for row in rows:
                seq = SeqRecord(Seq(row['sequence.residues']), id=row['primaryIdentifier'])
                seqrecords.append(seq)
            SeqIO.write(seqrecords, filepath, format='fasta')
            return seqrecords
        else:
            return list(SeqIO.parse(filepath, format='fasta'))

    @classmethod
    def genome_json(cls):
        """Return json object representing the genome"""
        seqs = cls.genome()
        return [{
            'name': seq.id,
            'bases': str(seq.seq),
            'circular': False
        } for seq in seqs]

    @classmethod
    def chromosome(cls, iden):
        """Return a chromosome"""
        for seq in cls.genome()[0]:
            if seq.id == iden:
                return seq

    @classmethod
    def genome_aligner(cls, sequences: list):
        """
        Return a blast aligner to the genome

        :param sequences: list of sequences to align to genome
        :type sequences: list
        :return: aligner
        :rtype: pyblast.JSONBlast
        """
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

        aligner = pyblast.JSONBlast(query_json=query_json, subject_json=cls.genome_json(), span_origin=False)
        return aligner

    @classmethod
    def align_to_genome(cls, sequences: list):
        """
        Align sequences to genome using blast.

        :param sequences: list of sequences
        :type sequences: list
        :return: aligner
        :rtype: pyblast.JSONBlast
        """
        aligner = cls.genome_aligner(sequences)
        aligner.quick_blastn()
        return aligner

    @classmethod
    def genomic_integration(cls, integrant: jdna.Sequence, min_homology=100, flanking_bps=500, max_distance=5000):
        """
        Simulate a genomic integration event.

        :param integrant: integrant sequence
        :type integrant: jdna.Sequence
        :param min_homology: minimum homology to consider
        :type min_homology: int
        :param flanking_bps:
        :type flanking_bps:
        :param max_distance:
        :type max_distance:
        :return:
        :rtype:
        """
        alignments = cls.align_to_genome([str(integrant)])
        results = cls.aligner_to_integrations(alignments, integrant=integrant, min_homology=min_homology,
                                           flanking_bps=flanking_bps, max_distance=max_distance)
        return results
