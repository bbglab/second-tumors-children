import pandas as pd
import os
import click
from multiprocessing import Pool
from aux_functions import read_vcf, fix_indels
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
from tqdm import tqdm
tqdm.pandas()


def process_haplocaller(row,nor_col,allele):
    '''
    :param row: each row is a variant of the VCF
    :param tum_col: tumor sample id
    :param nor_col: normal sample id
    :return: row with more columns with info from VCF
    '''
    row['n_ref_reads'] = int(row[nor_col].split(':')[1:2][0].split(',')[0])
    if allele == 0:
        row['n_alt_reads'] = int(row[nor_col].split(':')[1:3][0].split(',')[1])
        row['n_AF'] = float(row['INFO'].split(';')[1].split('=')[1])
    elif allele == 1:
        row['n_alt_reads'] = int(row[nor_col].split(':')[1:3][0].split(',')[1])
        row['n_AF'] = float(row['INFO'].split(';')[1].split('=')[1].split(',')[0])
        row['ALT'] = str(row['ALT'].split(',')[0])
    elif allele == 2:
        row['n_alt_reads'] = int(row[nor_col].split(':')[1:3][0].split(',')[2])
        row['n_AF'] = float(row['INFO'].split(';')[1].split('=')[1].split(',')[1])
        row['ALT'] = str(row['ALT'].split(',')[1])
    else:
        print('Specify which allele must be annotated: 1 or 2. If only 1 allele, indicate it as 0.')
 
    row['DP_normal'] = int(row[nor_col].split(':')[1:3][1].split(',')[0])
    row[['ALT']] = row[['ALT']].astype(str)
    
    # extract info genotype
    row['GT_normal'] = row[nor_col].split(':')[0]
    
    # infer mutation type
    if len(row['ALT']) != len(row['REF']):
        row['mut_type'] = 'indel'
    elif (len(row['ALT']) == len(row['REF'])) & (len(row['REF']) > 1):
        row['mut_type'] = 'mnv'
    else:
        row['mut_type'] = 'snv'
        
    return row

def process_chr_table (input_data):
    
    chrom = input_data[0]
    output_dir = input_data[1]
    normal_sample_id = input_data[2]

    print(chrom)
    df = pd.read_csv(output_dir + 'tmp/' + chrom + '.tsv.gz',sep='\t')

    df0 = df[~df['ALT'].str.contains(',')]
    df1 = df[df['ALT'].str.contains(',')]
    df2 = df[df['ALT'].str.contains(',')]

    df0 = df0.progress_apply(lambda x: process_haplocaller(x,normal_sample_id,allele=0), axis=1)
    df1 = df1.progress_apply(lambda x: process_haplocaller(x,normal_sample_id,allele=1), axis=1)
    df2 = df2.progress_apply(lambda x: process_haplocaller(x,normal_sample_id,allele=2), axis=1)

    df = pd.concat([df0,df1,df2],ignore_index=True)

    # CHANGE NORMAL/TUMOR COLUMNS NAME
    df = df.rename(columns={normal_sample_id:'NORMAL'})

    # CHANGE INDEL FORMAT
    df['ref_original'] = df['REF']
    df['alt_original'] = df['ALT']
    df['pos_original'] = df['POS']  
    df = df.progress_apply(lambda rw: fix_indels(rw),axis=1)
    df.drop(columns=['ref_original','alt_original','pos_original'],inplace=True)

    # ENSURE USEFUL COLUMNS      
    df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL' ,
         'n_ref_reads','n_alt_reads','n_AF','DP_normal','GT_normal','mut_type']]

    # ENSURE DATA TYPE
    df[['POS','n_ref_reads','n_alt_reads']] = df[['POS','n_ref_reads','n_alt_reads']].astype(int)     

    df.to_csv(output_dir+normal_sample_id+"_"+chrom+'.maf.gz', sep='\t', index=False, compression='gzip')

    #REMOVE TMP FILE
    os.remove(output_dir + 'tmp/' + chrom + '.tsv.gz')
    
    
@click.command()

@click.option('--input_vcf',
              '-i',
              required = True,
              help="path to vcf file from caller: mutect or sage")
@click.option('--output_dir',
              '-o',
              required = True,
              help="path to output folder")
@click.option('--normal_sample_id',
              '-n_id',
              required=True,
              help="Normal sample id. Should be the same as the column in VCF")
@click.option('--cpus',
              '-c',
              required=True,
              default=1,
              type=int,
              help="Number of cpus to be used")

def cli(input_vcf, output_dir, normal_sample_id, cpus):

    """
    Process VCF files from HaplotypeCaller and make filtered MAF, by chromosome
    """
 
    chroms = list(range(1,23))
    chroms = ['chr'+str(chrom) for chrom in chroms]
    sex_chroms = ['chrX','chrY']
    chroms = chroms + sex_chroms
    
    # CHECK/CREATE tmp/ FOLDER
    if not os.path.exists(output_dir + 'tmp'):
        os.makedirs(output_dir + 'tmp')
    
    # SPLIT FILE INTO ONE FILE PER CHROMOSOME
    if len(os.listdir(output_dir + 'tmp/')) != 24:
        
        # READ VCF and FILTER PASS variants
        df = read_vcf(input_vcf)

        df = df[df['FILTER'].str.contains('PASS')]
        
        print('save chr tables in tmp/ folder')
        for chrom in chroms:
            df_chr = df[df['#CHROM']==chrom]
            df_chr.to_csv(output_dir + 'tmp/' + chrom + '.tsv.gz',sep='\t',index=None,compression='gzip')
           
    # PREPARE INPUT DATA TO PARALELLISE 
    input_data = [[chrom,output_dir,normal_sample_id] for chrom in chroms]
    
    with Pool(cpus) as p:
        p.map(process_chr_table, input_data)
    
    #REMOVE TMP FOLDER (should be empty)
    os.rmdir(output_dir + 'tmp/')
    
if __name__ == "__main__":
    cli()
