import click
import pandas as pd
from bgreference import refseq #hg19, hg38
from collections import defaultdict
from tqdm import tqdm

tqdm.pandas()

def create_snv_class(df):
    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    x = df['TRIPLET']

    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])

    return out


def return_ordered_matrix(df):
    samples_dict = defaultdict(dict)
    order = snvs_order()

    for sample, data in df.groupby(by='SAMPLE'):
        dic_count = data['VARIANT_CLASS'].value_counts().to_dict()
        for i in order:
            samples_dict[sample][i] = dic_count.get(i, 0)

    matrix = pd.DataFrame.from_dict(samples_dict)
    matrix = matrix.loc[order]

    return matrix


# generate order of SBS
def snvs_order():
    order = []
    first = ['A', 'C', 'G', 'T']
    pyr = ['C', 'T']
    for p in pyr:
        for mut in first:
            if mut != p:
                for f in first:
                    for f2 in first:
                        comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                        order.append(comb)
    return order


def get_triplet(df, genomeref):
    try:
        triplet = refseq(genomeref, df['CHROM'], df['POS'] - 1, 3)

    except:
        triplet = 'NOT_IN_CHROM/genome reference not properly installed?'

    return triplet


def create_snvs_matrix(file, genome, output):

    print('Loading file...')
    df = pd.read_csv(file, sep='\t')

    print('Selecting SNVs...')
    # select whether we have SNVs or others
    df['len_alt'] = df['ALT'].str.len()

    # number of characters in ref
    df['len_ref'] = df['REF'].str.len()

    # select snvs
    snv_df = df[(df['len_alt'] == 1) & (df['len_ref'] == 1) & (df['ALT'] != '-') & (df['REF'] != '-')]

    print('Getting contexts...')

    snv_df['TRIPLET'] = snv_df.progress_apply(get_triplet, args=(genome,), axis=1)

    #if genome == 'hg19':
        #snv_df['TRIPLET'] = snv_df.progress_apply(get_triplet, args=(genome,), axis=1)

    #elif genome == 'hg38':
    #    snv_df['TRIPLET']  = snv_df.progress_apply(get_triplet, args=(hg38,), axis=1)

    # remove those variants that do not fall in normal chromosomes
    snv_df = snv_df[snv_df['TRIPLET'] != 'NOT_IN_CHROM']

    snv_df['REF_GENOME'] = snv_df.TRIPLET.map(lambda x: x[1])
    bad_ref = snv_df[snv_df['REF'] != snv_df['REF_GENOME']]
    wrong_mapped = len(bad_ref)

    if wrong_mapped > 0:
        print('WARNING! you have {} variants that do not mapp with the genome reference'.format(wrong_mapped))

    snv_df['VARIANT_CLASS'] = snv_df.apply(create_snv_class, axis=1)

    matrix = return_ordered_matrix(snv_df[['VARIANT_CLASS', 'SAMPLE']])

    #Add index labels (trinucleotides):

    channels = ['ACAA', 'ACCA', 'ACGA', 'ACTA', 'ACAG', 'ACCG', 'ACGG', 'ACTG', 'ACAT', 'ACCT', 'ACGT', 'ACTT', 'ATAA', 'ATCA', 'ATGA', 'ATTA', 'ATAC', 'ATCC', 
                'ATGC', 'ATTC', 'ATAG', 'ATCG', 'ATGG', 'ATTG', 'CCAA', 'CCCA', 'CCGA', 'CCTA', 'CCAG', 'CCCG', 'CCGG', 'CCTG', 'CCAT', 'CCCT', 'CCGT', 'CCTT', 
                'CTAA', 'CTCA', 'CTGA', 'CTTA', 'CTAC', 'CTCC', 'CTGC', 'CTTC', 'CTAG', 'CTCG', 'CTGG', 'CTTG', 'GCAA', 'GCCA', 'GCGA', 'GCTA', 'GCAG', 'GCCG', 
                'GCGG', 'GCTG', 'GCAT', 'GCCT', 'GCGT', 'GCTT', 'GTAA', 'GTCA', 'GTGA', 'GTTA', 'GTAC', 'GTCC', 'GTGC', 'GTTC', 'GTAG', 'GTCG', 'GTGG', 'GTTG', 
                'TCAA', 'TCCA', 'TCGA', 'TCTA', 'TCAG', 'TCCG', 'TCGG', 'TCTG', 'TCAT', 'TCCT', 'TCGT', 'TCTT', 'TTAA', 'TTCA', 'TTGA', 'TTTA', 'TTAC', 'TTCC', 
                'TTGC', 'TTTC', 'TTAG', 'TTCG', 'TTGG', 'TTTG']


    matrix['channels'] = channels
    matrix = matrix.set_index('channels')
    matrix.index.name = None

    #Transpose matrix
    #matrix = matrix.T


    matrix.to_csv(output, index = True, header = True, sep ='\t')


@click.command()
@click.option('--input_file',
                type=click.Path(exists=True),
                help="Input data",
                required=True)
@click.option('--genome_reference',
                #type=click.Choice(['hg19', 'hg38']),
                help="Genome reference",
                required = True)
@click.option('--output_file',
                type=click.Path(),
                help="output file",
                required=True)
def run(input_file, genome_reference, output_file):
    create_snvs_matrix(file = input_file, genome = genome_reference, output = output_file)

if __name__ == '__main__':
    run()
