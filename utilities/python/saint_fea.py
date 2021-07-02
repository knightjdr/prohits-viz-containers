from saint_fea.parse_args import parse_args
from saint_fea.enrich import enrich

'''
Usage:

python3 saint_fea.py \
-f 0.01 \
-s saint.txt \
-t 0

output: enrichment-saint.xlsx
'''

args = parse_args()
enrich(args)
