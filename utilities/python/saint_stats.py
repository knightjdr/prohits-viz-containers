from saint_stats.parse_args import parse_args
from saint_stats.summarize import summarize

'''
Usage:

python3 saint_stats.py \
-f 0.01 \
-s saint.txt

output: saint-statistics.txt
'''

args = parse_args()
summarize(args)
