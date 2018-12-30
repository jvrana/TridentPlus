"""
Organisms
"""

import os
from tridentplus.organisms.genome import Genome
from intermine.webservice import Service
import logging


def reset_logger():
    root = logging.getLogger()
    list(map(root.removeHandler, root.handlers[:]))
    list(map(root.removeFilter, root.filters[:]))

# necessary since importing intermine results in global changes in logging out of our control
reset_logger()

here = os.path.abspath(os.path.dirname(__file__))
genome_dir = os.path.join(here, '..', 'data', 'genomes')


class Yeast(Genome):
    """Yeast genome services"""

    genome_name = 'Saccharomyces_cerevisiae'
    service_url = "https://yeastmine.yeastgenome.org:443/yeastmine/service"
    service = Service(service_url)
    default_genome_path = os.path.join(genome_dir, genome_name)
