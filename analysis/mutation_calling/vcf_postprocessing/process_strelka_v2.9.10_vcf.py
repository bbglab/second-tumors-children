import click
import os
from aux_functions import read_vcf, snvs_processing, indels_processing, calculate_AF
from tqdm import tqdm
import pandas as pd
tqdm.pandas()


@click.command()
# @click.option('--patient',
#               '-pt',
#               required = True,
#               help="Patient id: pt1 - pt6")

@click.option('--tumor_sample_id',
              '-t_id',
              required = True,
              help="Tumor sample id from CNAG")
@click.option('--other_tumor_sample_id',
              '-ot_id',
              required=False,
              default=None,
              help="Tumor sample id from the other tumor, from which mutations will be rescued."
              "Should be the same as the column in VCF"
              )
@click.option('--normal_sample_id',
              '-n_id',
              required = True,
              help="Normal sample id from CNAG")
@click.option('--input_vcf_snvs',
              '-i',
              required = True,
              help="Normal sample id from CNAG")
@click.option('--output_dir',
              '-o',
              required = True,
              help="Normal sample id from CNAG")

def cli(input_vcf_snvs, output_dir, tumor_sample_id, other_tumor_sample_id, normal_sample_id):
    """
    It processes the results of Strelka and extracts information of the calls to return processed VCF in the form of
    TSV file
    """
       
    # read strelka results
    df_snvs = read_vcf(input_vcf_snvs)
    input_vcf_indels = input_vcf_snvs.replace('snvs','indels')
    df_indels = read_vcf(input_vcf_indels)

    if other_tumor_sample_id is not None:
      # read strelka results from the other tumor
      input_vcf_snvs_o = input_vcf_snvs.replace(tumor_sample_id,other_tumor_sample_id)
      input_vcf_indels_o = input_vcf_snvs_o.replace('snvs','indels')
    
      #Rescue muts from other tumor
      #snvs
      df_snvs_o = read_vcf(input_vcf_snvs_o)
      df_snvs_o = df_snvs_o[['#CHROM','POS','REF','ALT']][df_snvs_o['FILTER'].str.contains('PASS')]
    
      df_snvs = pd.merge(df_snvs,df_snvs_o,how='left',indicator='in_other_tumor')
      df_snvs['in_other_tumor'] = df_snvs['in_other_tumor'].replace('left_only',False)
      df_snvs['in_other_tumor'] = df_snvs['in_other_tumor'].replace('both',True)
      df_snvs = df_snvs[(df_snvs['FILTER'].str.contains('PASS'))|(df_snvs['in_other_tumor']==True)]    
      #indels
      df_indels_o = read_vcf(input_vcf_indels_o)
      df_indels_o = df_indels_o[['#CHROM','POS','REF','ALT']][df_indels_o['FILTER'].str.contains('PASS')]
      
      df_indels = pd.merge(df_indels,df_indels_o,how='left',indicator='in_other_tumor')
      df_indels['in_other_tumor'] = df_indels['in_other_tumor'].replace('left_only',False)
      df_indels['in_other_tumor'] = df_indels['in_other_tumor'].replace('both',True)
      df_indels = df_indels[(df_indels['FILTER'].str.contains('PASS'))|(df_indels['in_other_tumor']==True)]
    else:
      df_snvs = df_snvs[df_snvs['FILTER'].str.contains('PASS')]
      df_indels = df_indels[df_indels['FILTER'].str.contains('PASS')]

    # filter out weird chromosomes
    chroms = list(range(1,23))
    chroms = ['chr'+str(chrom) for chrom in chroms]
    chroms = chroms + ['chrX','chrY']

    df_snvs['#CHROM'].astype(str)
    df_snvs = df_snvs[df_snvs["#CHROM"].isin(chroms)]
    
    df_indels['#CHROM'].astype(str)
    df_indels = df_indels[df_indels["#CHROM"].isin(chroms)]

    # mut type processing
    df_snvs_proc = df_snvs.progress_apply(lambda x: snvs_processing(x), axis=1)
    df_indels_proc = df_indels.apply(lambda x: indels_processing(x), axis=1)

    # write results
    df_snvs_proc[['n_alt_reads', 'n_ref_reads',
        't_alt_reads', 't_ref_reads', 'POS']] = df_snvs_proc[['n_alt_reads', 'n_ref_reads', 't_alt_reads',
                                                    't_ref_reads', 'POS']].astype(int)
    df_snvs_proc['t_AF'] = df_snvs_proc.apply(lambda row: calculate_AF(row['t_alt_reads'],row['t_ref_reads']),axis=1)
    df_snvs_proc['n_AF'] = df_snvs_proc.apply(lambda row: calculate_AF(row['n_alt_reads'],row['n_ref_reads']),axis=1)
    
    df_indels_proc[['n_alt_reads', 'n_ref_reads',
        't_alt_reads', 't_ref_reads', 'POS']] = df_indels_proc[['n_alt_reads', 'n_ref_reads', 't_alt_reads',
                                                    't_ref_reads', 'POS']].astype(int)
    df_indels_proc['t_AF'] = df_indels_proc.apply(lambda row: calculate_AF(row['t_alt_reads'],row['t_ref_reads']),axis=1)
    df_indels_proc['n_AF'] = df_indels_proc.apply(lambda row: calculate_AF(row['n_alt_reads'],row['n_ref_reads']),axis=1)
    
    df_snvs_proc = df_snvs_proc[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR', 
                                 't_AF','n_AF','DP_tumor','t_alt_reads', 't_ref_reads', 'DP_normal', 'n_alt_reads',
                                 'n_ref_reads', 'mut_type','GT_normal','GT_tumor']]
    df_indels_proc = df_indels_proc[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL',
                                     'TUMOR', 't_AF','n_AF','DP_tumor','t_alt_reads', 't_ref_reads', 'DP_normal',
                                     'n_alt_reads', 'n_ref_reads', 'mut_type','GT_normal','GT_tumor']]
    
    df_proc = pd.concat([df_snvs_proc,df_indels_proc],ignore_index=True)

    sample = tumor_sample_id + '_vs_' + normal_sample_id
    df_proc.to_csv(os.path.join(output_dir, sample + "_process.maf.gz"), sep='\t', index=False,compression='gzip')

if __name__ == '__main__':
    cli()