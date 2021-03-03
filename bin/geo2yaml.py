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
        formatter_class=lambda prog: SmartFormatter(prog,
            indent_increment=2,
            max_help_position=8)
        )

    parser.add_argument(
        'accessions', 
        metavar='N', 
        type=str, 
        nargs='+',
        help='GEO/SRA/ENA accession with a Hi-C project. Multiple values are allowed.')
    
    parser.add_argument(
        '--print-srr-table',
        action='store_true',
        help='If provided, print the table of SRR accessions instead of a distiller yaml config.'
    )

    parser.add_argument(
        '--library-eval',
        type=str,
        action='append', 
        default = [],
        help='A list of functions to generate library IDs from an SRR table. '
        'Internally, each supplied string is fed into pandas.DataFrame.eval(srr_df, s, engine="python"). '
        'The result is stored in a column d_lib, which can be reused after the first function. '
        'The default values are: \n'
        '--library-eval \'experiment_title.str.extract(": ([^;]*);")\' \\\n'
        '--library-eval \'d_lib.str.replace("\s", "_", regex=True)\' \\\n'
        '--library-eval \'d_lib.str.replace("[^\w_.-]", "", regex=True)\' \\\n',
    )

    parser.add_argument(
        '--group-eval', 
        type=str,
        action='append',
        default = [],
        help='A list of functions to generate library IDs from an SRR table. '
        'Internally, each supplied string is fed into pandas.DataFrame.eval(srr_df, s, engine="python"). '
        'The result is stored in a column d_group, which can be reused after the first function. '
        'The default values are: \n'
        '--group-eval \'d_lib.str.replace("[_-](R|rep)_?[\d+]$", "", regex=True)\''
    )

    parser.add_argument(
        '--filter-srr', 
        action='append',
        default = [],
        type=str,
        help='A list of functions to filter the SRR table. '
        'Internally, each supplied string is fed into pandas.DataFrame.query(srr_df, s, engine="python"). '
        'By default, no filtering is applied. Example: \n'
        '--filter-srr d_lib.str.contains("[Hh][Ii]-?[Cc]", regex=True)'
        
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

DEFAULT_LIBRARY_EVAL = [
    'experiment_title.str.extract(": ([^;]*);")',
    'd_lib.str.replace("\s", "_", regex=True)',
    'd_lib.str.replace("[^\w_.-]", "", regex=True)',
]

DEFAULT_GROUP_EVAL = [
    'd_lib.str.replace("[_-](R|rep)_?[\d+]$", "", regex=True)',
]

DEFAULT_SRR_FILTER = []

TAB_CHAR = '    '

args = parse_args(sys.argv[1:])

db = pysradb.SRAweb()

queries = to_downloadable(args.accessions)
 
srr_df = pd.concat([    
    db.sra_metadata(q, detailed=True)
    for q in queries
])

for library_eval in (args.library_eval if args.library_eval else DEFAULT_LIBRARY_EVAL):
    srr_df['d_lib'] = srr_df.eval(library_eval, engine='python')

for group_eval in (args.group_eval if args.group_eval else DEFAULT_GROUP_EVAL):
    srr_df['d_group'] = srr_df.eval(group_eval, engine='python')

for srr_filter in (args.filter_srr if args.filter_srr else DEFAULT_SRR_FILTER):
    srr_df = srr_df.query(srr_filter, engine='python')

srr_df = srr_df.sort_values(['d_lib', 'run_accession'])

srr_df['lane'] = (
    'lane'
    + (srr_df.groupby('d_lib').cumcount()+1).astype('str')
)


# Keeping this code in case YAML structures will become useful:

# out_raw_reads_paths = {}
# for title, grouped in srr_df.groupby('d_lib'):
#     out_raw_reads_paths[title] = {
#         row.lane:f'- sra:{row.run_accession}'
#         for _,row in grouped.iterrows()
#     }

# out_library_groups = {}
# for group, grouped in srr_df.groupby('group'):
#     distiller_libraries = list(grouped.d_lib.unique())
#     if len(distiller_libraries) > 1:
#         out_library_groups[group] = distiller_libraries


if args.print_srr_table:
    srr_df.to_csv(sys.stdout, sep='\t', index=False)
else:
    out_raw_reads_paths = [f'{TAB_CHAR}raw_reads_paths:']
    for title, grouped in srr_df.groupby('d_lib'):
        out_raw_reads_paths.append(f'{TAB_CHAR}{TAB_CHAR}{title}:')
        for _, row in grouped.iterrows():
            out_raw_reads_paths.append(f'{TAB_CHAR}{TAB_CHAR}{TAB_CHAR}{row.lane}:')
            out_raw_reads_paths.append(f'{TAB_CHAR}{TAB_CHAR}{TAB_CHAR}{TAB_CHAR}- sra:{row.run_accession}')

    out_library_groups = [f'{TAB_CHAR}library_groups:']
    for group, grouped in srr_df.groupby('d_group'):
        distiller_libraries = grouped['d_lib'].unique()
        if len(distiller_libraries) > 1:
            out_library_groups.append(f'{TAB_CHAR}{TAB_CHAR}{group}:')
            out_library_groups += [f'{TAB_CHAR}{TAB_CHAR}{TAB_CHAR}- {title}' 
                                for title in distiller_libraries]

    out = '\n'.join(['input:']+out_raw_reads_paths+out_library_groups) 

    print('# generated by geo2yaml.py:')
    print('# python ' + ' '.join(sys.argv))
    print(out)

