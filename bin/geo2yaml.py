import sys
import argparse
import re

import numpy as np
import pandas as pd
import pysradb


# def unescaped_str(arg_str):
#     """
#     Borrowed from https://stackoverflow.com/questions/34145686/handling-argparse-escaped-character-as-option
#     """
#     return codecs.decode(str(arg_str), 'unicode-escape')

class SmartFormatter(argparse.HelpFormatter):
    '''
         Custom Help Formatter used to split help text when '\n' was 
         inserted in it.
    '''

    def _split_lines(self, text, width):
        r = []
        for t in text.splitlines(): r.extend(argparse.HelpFormatter._split_lines(self, t, width))
        return r

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Create the input section for distiller\'s project.yml from GEO/ENA/SRA accessions.',
        # formatter_class=argparse.RawDescriptionHelpFormatter
        formatter_class=SmartFormatter
        )
    parser.add_argument(
        'accessions', 
        metavar='N', 
        type=str, 
        nargs='+',
        help='GEO/SRA/ENA accession with a Hi-C project. Multiple values are allowed.')
    parser.add_argument(
        '--title_sub', 
        nargs=2, 
        action='append',
#         type=unescaped_str,
        default = [],
        help='A list of regular expression substitutions to clean up the experiment titles. '
        'Multiple sequential substitutions are allowed. ' 
        'Each substitution must be provided using a separate flag --title_sub followed by '
        'a pair of regular expressions pat repl, separated by a space, '
        'where pat is the matching pattern and repl is the replacement string. '
        'Internally, these expressions are then provided to pandas.Series.str.replace() or re.sub(). '
        'The default substitutions (1) replace spaces with underscores and (2) remove characters not matching '
        'A–Z a–z 0–9 ._- (a.k.a. the POSIX portable file name character set):'
        '\n'
        '--title_sub \'\\s\' \'_\' --title_sub \'[^\\w_.-]\' \'\''
    )
    parser.add_argument(
        '--group_sub', 
        nargs=2, 
        action='append',
#         type=unescaped_str,
        default = [],
        help='A list of regular expression substitutions to convert experiment titles into groups. '
        'The usage is same as above. The default substitution removes patterns like _R1/_R2/_rep1/-R1/R1 '
        'at the end of the experiment title:'
        '\n'
        '--group_sub \'[_-](R|rep)[\\d+]$\' \'\''
    )

    parser.add_argument(
        '--filter_pre', 
        nargs=1, 
        action='append',
        default = [],
        type=str,
        help='A regular expression to filter datasets by their *unedited* name. '
        'If multiple filters are provided, select datasets that satisfy at least one of the filters. '
        '--filter \'[Hh][Ii]-?[Cc]\''
    )

    parser.add_argument(
        '--filter_post', 
        action='append',
        default = [],
        type=str,
        help='A regular expression to filter datasets by their *edited* name. '
        'If multiple filters are provided, select datasets that satisfy at least one of the filters. '
        '--filter \'[Hh][Ii]-?[Cc]\''
    )
    
    return parser.parse_args(args)

def to_downloadable(queries):
    out_queries = []
    for q in queries:
        if q.startswith('GSE'):
            out_queries += list(
                pysradb.SRAweb()
                .gse_to_srp(q)
                .study_accession
            )
        else: 
            out_queries.append(q)
    return out_queries

DEFAULT_TITLE_SUB = [
        ('\s', '_'),
        ('[^\w_.-]', '') # the first character cannot be a hyphen!!
    ]
DEFAULT_GROUP_SUB =  [
    ('[_-](R|rep)_?[\d+]$', '')
]

TAB_CHAR = '    '

args = parse_args(sys.argv[1:])

db = pysradb.SRAweb()

queries = to_downloadable(args.accessions)
 
srr_table = pd.concat([    
    db.sra_metadata(q)
    for q in queries
])

srr_table = srr_table[['experiment_title', 'run_accession']]

srr_table['experiment_title'] = (
    srr_table['experiment_title']
    .str.split(';')
    .str.get(0)
    .str.split(':')
    .str.get(1)
    .str.strip()
)

if args.filter_pre:
    mask = np.logical_or.reduce([
        srr_table['experiment_title'].str.contains(fltr, regex=True) 
        for fltr in list(args.filter_pre)])
    srr_table = srr_table[mask]

for re_sub in (args.title_sub if args.title_sub else DEFAULT_TITLE_SUB):
    srr_table['experiment_title'] = (
            srr_table.experiment_title
            .str.replace(re_sub[0], re_sub[1], regex=True)
        )

if args.filter_post:
    mask = np.logical_or.reduce([
        srr_table['experiment_title'].str.contains(fltr, regex=True) 
        for fltr in list(args.filter_post)])
    srr_table = srr_table[mask]


srr_table=srr_table.sort_values(['experiment_title','run_accession'])

srr_table['lane'] = (
    'lane'
    + (srr_table.groupby('experiment_title').cumcount()+1)
        .astype('str')
)

group = srr_table.experiment_title
for sub in (args.group_sub if args.group_sub else DEFAULT_GROUP_SUB):
    group = group.str.replace(sub[0], sub[1])
srr_table['group'] = group

# Keeping this code in case YAML structures will become useful:

# out_raw_reads_paths = {}
# for title, grouped in srr_table.groupby('experiment_title'):
#     out_raw_reads_paths[title] = {
#         row.lane:f'- sra:{row.run_accession}'
#         for _,row in grouped.iterrows()
#     }

# out_library_groups = {}
# for group, grouped in srr_table.groupby('group'):
#     experiment_titles = list(grouped.experiment_title.unique())
#     if len(experiment_titles) > 1:
#         out_library_groups[group] = experiment_titles


out_raw_reads_paths = [f'{TAB_CHAR}raw_reads_paths:']
for title, grouped in srr_table.groupby('experiment_title'):
    out_raw_reads_paths.append(f'{TAB_CHAR}{TAB_CHAR}{title}:')
    for _, row in grouped.iterrows():
        out_raw_reads_paths.append(f'{TAB_CHAR}{TAB_CHAR}{TAB_CHAR}{row.lane}:')
        out_raw_reads_paths.append(f'{TAB_CHAR}{TAB_CHAR}{TAB_CHAR}{TAB_CHAR}- sra:{row.run_accession}')

out_library_groups = [f'{TAB_CHAR}library_groups:']
for group, grouped in srr_table.groupby('group'):
    experiment_titles = grouped.experiment_title.unique()
    if len(experiment_titles) > 1:
        out_library_groups.append(f'{TAB_CHAR}{TAB_CHAR}{group}:')
        out_library_groups += [f'{TAB_CHAR}{TAB_CHAR}{TAB_CHAR}- {title}' 
                               for title in experiment_titles]

out = '\n'.join(['input:']+out_raw_reads_paths+out_library_groups) 

print(out)

