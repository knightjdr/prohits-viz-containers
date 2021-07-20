from saint_domain_enrich.parse_args import parse_args
from saint_domain_enrich.enrich import enrich

'''
Usage:

python3 saint_domain_enrich.py \
-d domains.json \
-f 0.01 \
-g genemap.json \
-i refseqp \
-s saint.txt \
-t 0

output: domain-enrichment-saint.xlsx
'''

args = parse_args()
enrich(args)
