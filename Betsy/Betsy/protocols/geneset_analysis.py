#geneset_analysis.py

from Betsy.protocol_utils import Parameter

PRETTY_NAME="Find enriched gene sets."

SCORE_PATHWAY = 'Score Pathway parameters'
CATEGORIES = [SCORE_PATHWAY]
#input predicates
INPUTS = [
    'gse_id',
    'gse_id_and_platform',
    'cel_files',
    'geneset_file']

#output predicates
OUTPUTS = 'make_geneset_report'

#convert predicates to prolog facts
predicate2arguments = {
    'gse_id': ([], '[]'),
    'gse_id_and_platform': ([], '[]'),
    'cel_files': (['version', 'unknown_version'], '[]'),
    'geneset_file':([],'[]')}

#parameter objects
PARAMETERS=[Parameter('automatch', pretty_name='Auto Match',choices=['yes_automatch','no_automatch'],
                      category=SCORE_PATHWAY),
            Parameter('geneset', pretty_name='Gene Set',type='string',category=SCORE_PATHWAY),
            Parameter('allgenes', pretty_name='All genes',choices=['yes_allgenes','no_allgenes'],category=SCORE_PATHWAY)]


    
