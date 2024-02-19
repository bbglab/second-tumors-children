import click
import pandas as pd
from aux_functions import read_vcf
from tqdm import tqdm
from multiprocessing import Pool
from bgreference import hg38
tqdm.pandas()

def most_damaging(string):
    """
    in case of more than one consequence type per row that is coma-separated it returns the most damaging or severe
    """
    if "," in string:
        string_list = string.split(",")
        for dam in DAMAGING_CONSEQ:
            if dam in string_list:
                return dam
                break
    else:
        return string


def get_cannonical(group):
    """
    For each variant, that is unique combination of #CHROM, POS, REF, ALT, it keeps CANONICAL transcript annotations.
    If more than consequence type is annotated, it keeps the most damaging or severe either if there is more than one
    entry (row) per variant or either is coma-separated in the same row
    :param group: groupby object (~df slice) of each variant and all the reported annotations by VEP92
    :return:a one row data frame for each variant
    """
    name, dff = group
    if len(dff) > 1:
        if 'YES' in list(dff['CANONICAL'].unique()):  # get canonical in case of more than one annotation
            dff = dff[dff['CANONICAL'] == 'YES']
            if len(dff) > 1:
                conseqs = dff['Consequence'].tolist()  # get the most severe
                conseqs = [most_damaging(x) for x in
                           conseqs]  # in case of two consequence type in the same row get the most damaging
                for cons in DAMAGING_CONSEQ:
                    if cons in conseqs:
                        dff = dff[dff['Consequence'].str.contains(cons)]
                        if len(dff) > 1:
                            dff = dff.drop_duplicates(subset=['#Uploaded_variation'], keep='first')
                            break
                        else:
                            break
            else:# get the most severe
                conseqs = dff['Consequence'].tolist()
                conseqs = [most_damaging(x) for x in
                           conseqs]  # in case of two consequence type in the same row get the most damaging
                for cons in DAMAGING_CONSEQ:
                    if cons in conseqs:
                        dff = dff[dff['Consequence'].str.contains(cons)]
                        if len(dff) > 1:
                            dff = dff.drop_duplicates(subset=['#Uploaded_variation'], keep='first')
        else:  # get the most severe
            conseqs = dff['Consequence'].tolist()
            conseqs = [most_damaging(x) for x in
                       conseqs]  # in case of two consequence type in the same row get the most damaging
            for cons in DAMAGING_CONSEQ:
                if cons in conseqs:
                    dff = dff[dff['Consequence'].str.contains(cons)]
                    if len(dff) > 1:
                        dff = dff.drop_duplicates(subset=['#Uploaded_variation'], keep='first')

        return dff
    else:
        conseqs = dff['Consequence'].tolist()
        conseqs = [most_damaging(x) for x in
                   conseqs]  # in case of two consequence type in the same row get the most damaging
        for cons in DAMAGING_CONSEQ:
            if cons in conseqs:
                dff = dff[dff['Consequence'].str.contains(cons)]
                if len(dff) > 1:
                    dff = dff.drop_duplicates(subset=['#Uploaded_variation'], keep='first')
        return dff


def process_af(rw):
    """
    The column gnomADg_AF contains the allele frequency population information but it has mixed data types.
    :param rw: row of the dataframe (one variant)
    :return: a unique value of gnomADg_AF for each variant
    """
    if (rw['gnomADg_AF'] == '-') or (rw['gnomADg_AF'] == '.'):
        rw['gnomADg_AF'] = 0
    elif "," in str(rw['gnomADg_AF']):
        lista = rw['gnomADg_AF'].split(",")
        try:
            lista = [float(x) for x in lista]
        except ValueError:
            lista = list(set(lista))
            if "." in lista:
                lista.remove(".")
            if "-" in lista:
                lista.remove("-")
            if len(lista) != 0:
                lista = [float(x) for x in lista]
            else:
                lista = [0]
        mini = min(lista)
        rw['gnomADg_AF'] = float(mini)
    else:
        do_nothing="do_nothing"
    return rw


def separate(rw):
    """
    Separate information from id variant column #Uploaded_variation
    :param rw: row of the dataframe (one variant)
    :return: row with recovered information in new columns
    """

    if "_" in rw['#Uploaded_variation']:
        rw['#CHROM'] = rw['#Uploaded_variation'].split("_")[0]
        rw['POS'] = rw['#Uploaded_variation'].split("_")[1]
        rw['Change'] = rw['#Uploaded_variation'].split("_")[2]
        rw['REF'] = rw['Change'].split("/")[0]
        rw['ALT'] = rw['Change'].split("/")[1]
    else:
        rw['#CHROM'] = str(rw['Location'].split(":")[0])
        rw['POS'] = rw['Location'].split(":")[1]
        if "-" in rw['POS']:
            rw['POS'] = int(rw['POS'].split("-")[0])
        else:
            rw['POS'] = int(rw['POS'])
        rw['REF'] = hg19(str(rw['#CHROM']), rw['POS'], 1)
        rw['ALT'] = rw['Allele']
    return rw


DAMAGING_CONSEQ= ['transcript_ablation',
'splice_acceptor_variant',
'splice_donor_variant',
'stop_gained',
'frameshift_variant',
'stop_lost',
'start_lost',
'transcript_amplification',
'inframe_insertion',
'inframe_deletion',
'missense_variant',
'protein_altering_variant',
'splice_region_variant',
'incomplete_terminal_codon_variant',
'start_retained_variant',
'stop_retained_variant',
'synonymous_variant',
'coding_sequence_variant',
'mature_miRNA_variant',
'5_prime_UTR_variant',
'3_prime_UTR_variant',
'non_coding_transcript_exon_variant',
'intron_variant',
'NMD_transcript_variant',
'non_coding_transcript_variant',
'upstream_gene_variant',
'downstream_gene_variant',
'TFBS_ablation',
'TFBS_amplification',
'TF_binding_site_variant',
'regulatory_region_ablation',
'regulatory_region_amplification',
'feature_elongation',
'regulatory_region_variant',
'feature_truncation',
'intergenic_variant'
]

@click.command()
@click.option('--input',
              '-i',
              required = True,
              help="File that contains the tab delimited maf annotated by vep101")
@click.option('--output',
              '-o',
              required = True,
              help="File name for the output. E.g. 'AQ5174_chr1_vep_filt.tsv'")
@click.option('--cores',
              '-c',
              required = True,
              help="Number of cores to process table")


def cli(input, output, cores):
    """
    It takes the TAB table format of the output of VEP end filters by the canonical and damaging consequence type
    """
        
    # read vep annotated file
    df_ann = read_vcf(input)
    
    # get one annotation per variant
    df_ann_filt = pd.DataFrame()

    with Pool(int(cores)) as pool:
        for df_res in tqdm(pool.imap(get_cannonical, df_ann.groupby("#Uploaded_variation")),
                           total=len(df_ann.groupby("#Uploaded_variation"))):
            df_ann_filt = df_ann_filt.append(df_res, ignore_index=True)

    # prepare tab file for merge

    df_ann_filt = df_ann_filt.apply(lambda x: separate(x), axis=1)

    # process ambiguous values and turn them to 0 and keep the minimum AF in case of two values
    df_ann_filt = df_ann_filt.apply(lambda x: process_af(x), axis=1)
    df_ann_filt['gnomADg_AF'] = df_ann_filt['gnomADg_AF'].apply(lambda x: 0 if x == "." or x == "-" else float(x))

    print(len(df_ann_filt))

    # write results
    df_ann_filt.to_csv(output, sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    cli()