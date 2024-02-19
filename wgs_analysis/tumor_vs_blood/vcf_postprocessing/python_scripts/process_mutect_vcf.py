import pandas as pd
import click
from aux_functions import read_vcf, fix_indels
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
from tqdm import tqdm
tqdm.pandas()


def process_mutect2(row, tum_col, nor_col):
    '''
    :param row: each row is a variant of the VCF
    :param tum_col: tumor sample id
    :param nor_col: normal sample id
    :return: row with more columns with info from VCF
    '''
    # extract info of reads
    row['t_AF'] = row[tum_col].split(':')[2]
    row['n_AF'] = row[nor_col].split(':')[2]
    row['t_ref_reads'] = int(row[tum_col].split(':')[1:2][0].split(',')[0])
    row['t_alt_reads'] = int(row[tum_col].split(':')[1:2][0].split(',')[1])
    row['DP_tumor'] = int(row['t_ref_reads']+row['t_alt_reads'])
    row['n_ref_reads'] = int(row[nor_col].split(':')[1:2][0].split(',')[0])
    row['n_alt_reads'] = int(row[nor_col].split(':')[1:2][0].split(',')[1])
    row['DP_normal'] = int(row[nor_col].split(':')[1:2][0].split(',')[0])
    row[['ALT']] = row[['ALT']].astype(str)
    # extract info genotype
    row['GT_normal'] = row[nor_col].split(':')[0]
    row['GT_tumor'] = row[tum_col].split(':')[0]

    # infer mutation type
    if len(row['ALT']) != len(row['REF']):
        row['mut_type'] = 'indel'
    elif (len(row['ALT']) == len(row['REF'])) & (len(row['REF']) > 1):
        row['mut_type'] = 'mnv'
    else:
        row['mut_type'] = 'snv'

    return row

def process_mutect2_tumor_only(row, tum_col):
    '''
    :param row: each row is a variant of the VCF
    :param tum_col: tumor sample id
    :return: row with more columns with info from VCF
    '''
    # extract info of reads
    row['AF'] = row[tum_col].split(':')[2]
    row['ref_reads'] = int(row[tum_col].split(':')[1:2][0].split(',')[0])
    row['alt_reads'] = int(row[tum_col].split(':')[1:2][0].split(',')[1])
    row['DP'] = int(row['ref_reads']+row['alt_reads'])
    row[['ALT']] = row[['ALT']].astype(str)
    # extract info genotype
    row['GT'] = row[tum_col].split(':')[0]

    # infer mutation type
    if len(row['ALT']) != len(row['REF']):
        row['mut_type'] = 'indel'
    elif (len(row['ALT']) == len(row['REF'])) & (len(row['REF']) > 1):
        row['mut_type'] = 'mnv'
    else:
        row['mut_type'] = 'snv'

    return row


@click.command()

@click.option('--input_vcf',
              '-i',
              required = True,
              help="path to vcf file from caller: mutect or sage")
@click.option('--output_dir',
              '-o',
              required = True,
              help="path to output folder")
@click.option('--tumor_sample_id',
              '-t_id',
              required=True,
              help="Tumor sample id. Should be the same as the column in VCF"
              )
@click.option('--other_tumor_sample_id',
              '-ot_id',
              required=False,
              default=None,
              help="Tumor sample id from the other tumor, from which mutations will be rescued."
              "Should be the same as the column in VCF"
              "If not provided it will just select PASS variants from the tumor file"
              )
@click.option('--panel_of_normals',
              '-pon',
              is_flag = True,
              required=False,
              default=False,
              help="If true, include variants present in the panel of normals from mutect"
              "This is set up for the analysis with Mutect2 tumor only mode"
              "If not provided it will just select PASS variants from the tumor file"
              )
@click.option('--normal_sample_id',
              '-n_id',
              required=False,
              default=None,
              help="Normal sample id. Should be the same as the column in VCF"
              "If not provided, will run as Mutect2 tumor only mode"
              )
@click.option('--germline/--somatic',
              '-g/-s',
              required=False,
              default=False,
              help="Whether the vcf file is for germline alterations."
              )

def cli(input_vcf, output_dir, tumor_sample_id, other_tumor_sample_id, panel_of_normals, normal_sample_id,germline):

    """
    Process VCF files from Mutect2 or SAGE and make filtered MAF
    """
       
    # READ VCF and FILTER variants
    df = read_vcf(input_vcf)

    if other_tumor_sample_id is not None: # (include rescued mutations from the other tumor from the same patient)
    	input_vcf2 = input_vcf.replace(tumor_sample_id,other_tumor_sample_id)
    	
    	df_o = read_vcf(input_vcf2)
    	df_o = df_o[['#CHROM','POS','REF','ALT']][df_o['FILTER']=='PASS']
    
    	df = pd.merge(df,df_o,how='left',indicator='in_other_tumor')
    	df['in_other_tumor'] = df['in_other_tumor'].replace('left_only',False)
    	df['in_other_tumor'] = df['in_other_tumor'].replace('both',True)
    	df = df[(df['FILTER'].str.contains('PASS'))|(df['in_other_tumor']==True)]

    elif panel_of_normals == True: # include variants filtered by the panel of normals from mutect only tumor
      df = df[(df['FILTER']=='PASS')|(df['FILTER']=='panel_of_normals')]

    else: # include only Filter == PASS
      df = df[df['FILTER']=='PASS']

    if normal_sample_id is not None:
      df = df.progress_apply(lambda x: process_mutect2(x, tumor_sample_id, normal_sample_id), axis=1)

      # CHANGE NORMAL/TUMOR COLUMNS NAME
      df = df.rename(columns={normal_sample_id:'NORMAL',tumor_sample_id:'TUMOR'})


      # ENSURE DATA TYPE
      df[['n_alt_reads', 'n_ref_reads',
              't_alt_reads', 't_ref_reads', 'POS']] = df[['n_alt_reads', 'n_ref_reads', 't_alt_reads',
                                                          't_ref_reads', 'POS']].astype(int)

      # ENSURE USEFUL COLUMNS
      df = df[
          ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','NORMAL' ,'TUMOR' ,
           't_AF', 'n_AF', 'DP_tumor', 't_alt_reads', 't_ref_reads', 'DP_normal', 'n_alt_reads', 'n_ref_reads', 'mut_type',
           'GT_normal','GT_tumor']]        
      
      # WRITE RESULTS
      if germline == False:
          output_file = output_dir + tumor_sample_id + '_vs_' + normal_sample_id + '_process.maf.gz'
          df.to_csv(output_file, sep='\t', index=False,compression='gzip')
          
      elif germline == True:
          # SAVE TABLES DIVIDED BY CHROMOSOME
          chrom_groups = df.groupby("#CHROM")

          for chr_ in chrom_groups.groups:
              df_chr = chrom_groups.get_group(chr_)            
              df_chr.to_csv(output_dir+normal_sample_id+"_"+chr_+'.maf.gz', sep='\t', index=False,compression='gzip')

    else:

      df = df.progress_apply(lambda x: process_mutect2_tumor_only(x, tumor_sample_id), axis=1)

      # CHANGE NORMAL/TUMOR COLUMNS NAME
      df = df.rename(columns={tumor_sample_id:'TUMOR'})

      # ENSURE DATA TYPE
      df[['alt_reads', 'ref_reads', 'POS']] = df[['alt_reads','ref_reads', 'POS']].astype(int)

      # ENSURE USEFUL COLUMNS
      df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','TUMOR',
            'AF', 'DP', 'alt_reads', 'ref_reads', 'mut_type','GT']]        
      
      # WRITE RESULTS
      output_file = output_dir + tumor_sample_id + '_process.maf.gz'
      df.to_csv(output_file, sep='\t', index=False,compression='gzip')
              
if __name__ == "__main__":
    cli()
