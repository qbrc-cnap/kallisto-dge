import argparse 
import pandas as pd

DESEQ_RESULTS = 'deseq_results'
GENE_COUNTS = 'gene_counts'
TRANSCRIPT_TO_GENE_MAP = 't2g'
HUMAN_TAG = 'human_tag'
MOUSE_TAG = 'mouse_tag'

HUMAN_PREFIX = 'ENST'

def parse_input():
    '''
    Parses the commandline input, returns a dict
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', required=True, dest=DESEQ_RESULTS)
    parser.add_argument('-c', required=True, dest=GENE_COUNTS)
    parser.add_argument('-t', required=True, dest=TRANSCRIPT_TO_GENE_MAP)
    parser.add_argument('-ht', required=True, dest=HUMAN_TAG)
    parser.add_argument('-mt', required=True, dest=MOUSE_TAG)

    args = parser.parse_args()
    return vars(args)


if __name__ == '__main__':
    
    args = parse_input()

    # the transcript to gene map has three columns: 'ensg', 'enst', and 'name'
    t2g_file = args[TRANSCRIPT_TO_GENE_MAP]
    df = pd.read_csv(t2g_file, sep='\t')

    # split the table into human and mouse components
    human_rows = df['enst'].apply(lambda x: x.startswith(HUMAN_PREFIX))
    human_tx = df.loc[human_rows]
    mouse_tx = df.loc[~human_rows]

    # get the human and mouse gene symbols
    all_human_genes = set(human_tx['name'])
    all_mouse_genes = set(mouse_tx['name'])

    # open the gene-level count table and split that based on whether human or mouse gene
    counts = pd.read_csv(args[GENE_COUNTS], sep='\t', index_col=0)
    human_genes = all_human_genes.intersection(counts.index)
    mouse_genes = all_mouse_genes.intersection(counts.index)
    human_counts = counts.loc[human_genes]
    mouse_counts = counts.loc[mouse_genes]
    contrast, suffix = args[GENE_COUNTS].split('.')
    human_counts_file = '%s.%s.%s' % (contrast, args[HUMAN_TAG], suffix)
    human_counts.to_csv(human_counts_file, sep='\t', index_label='gene')
    mouse_counts_file = '%s.%s.%s' % (contrast, args[MOUSE_TAG], suffix)
    mouse_counts.to_csv(mouse_counts_file, sep='\t', index_label='gene')
    del counts
    del human_counts
    del mouse_counts

    # open the gene-level deseq results table 
    # and split that based on whether human or mouse gene
    deseq_results = pd.read_csv(args[DESEQ_RESULTS], sep='\t', index_col=0)
    human_genes = all_human_genes.intersection(deseq_results.index)
    mouse_genes = all_mouse_genes.intersection(deseq_results.index)
    human_table = deseq_results.loc[human_genes]
    mouse_table = deseq_results.loc[mouse_genes]
    contrast, suffix = args[DESEQ_RESULTS].split('.')
    human_file = '%s.%s.%s' % (contrast, args[HUMAN_TAG], suffix)
    human_table.to_csv(human_file, sep='\t', index_label='gene')
    mouse_file = '%s.%s.%s' % (contrast, args[MOUSE_TAG], suffix)
    mouse_table.to_csv(mouse_file, sep='\t', index_label='gene')


