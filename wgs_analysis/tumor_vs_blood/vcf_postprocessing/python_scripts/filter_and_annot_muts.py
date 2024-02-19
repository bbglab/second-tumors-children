import pandas as pd
import os
import click
from tqdm import tqdm
from functions_filterings_snvs_after_vep import cnv_snvs_intersect, calculate_VAF, calculate_ccf
tqdm.pandas()

most_damaging = ['transcript_ablation',
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
'stop_retained_variant']

variants_dict = {'truncating':['transcript_ablation','splice_acceptor_variant','splice_donor_variant',
                'stop_gained','frameshift_variant','stop_lost','start_lost','transcript_amplification'],
                'miss_inframe':['inframe_insertion','inframe_deletion','missense_variant'],
                'other':['protein_altering_variant','splice_region_variant','incomplete_terminal_codon_variant',
                'start_retained_variant','stop_retained_variant']}

def damaging_col(x):
    x = str(x)
    x_list = x.split(',')
    y = list(set(x_list) & set(most_damaging))
    if y != []:       
        return True
    else:
        return False
    
def aa_change_col(row):
    pos = str(row['Protein_position'])
    aa = str(row['Amino_acids'])
    if '/' in aa:
        aa_ref = aa.split('/')[0]
        aa_alt = aa.split('/')[1]
        aa_change = aa_ref + pos + aa_alt
    elif aa == '-':
        aa_change = aa
    else:
        aa_change = aa + pos
    return aa_change

def real_af (row):
    alt = row['n_alt_reads']
    ref = row['n_ref_reads']
    if alt == 0:
        af = 0
    else:
        af = round(alt/(alt+ref),3)
    return af

def variant_type(x,variants_dict):
    if isinstance(x, str):
        if ',' in x:
            x = x.split(',')
            for c in x:
                if c in variants_dict['truncating']:
                    return 'truncating'
                elif c in variants_dict['miss_inframe']:
                    return 'miss_inframe'
                elif c in variants_dict['other']:
                    return 'other'
                else:
                    return None
        else:
            if x in variants_dict['truncating']:
                return 'truncating'
            elif x in variants_dict['miss_inframe']:
                return 'miss_inframe'
            elif x in variants_dict['other']:
                return 'other'
            else:
                return None

    else:
        return None
    
@click.command()

@click.option('--input_dir',
              '-i',
              required = True,
              help="path to folder from process_vep output")
@click.option('--output_dir',
              '-o',
              required = True,
              help="path to output folder")
@click.option('--tumor_sample_id',
              '-t_id',
              required=False,
              help="Tumor sample id. Should be the same as the column in VCF."
              "Only needed for 'intersect' caller (somatic)")
@click.option('--normal_sample_id',
              '-n_id',
              required=True,
              help="Normal sample id. Should be the same as the column in VCF")
@click.option('--caller',
              '-c',
              required=True,
              help="Name of the variant caller used: 'sage' or 'hc' (haplotype caller) for germline,"
              "'intersect' for somatic")
@click.option('--gnomad_threshold',
              '-gn',
              required=True,
              default=1,
              type=float,              
              help="Gnomad AF threshold to be applied. If none, specify '1'")
@click.option('--cnv_file',
              '-cnv',
              required=False,              
              help="Path to cnv file from purple or ascat")   
@click.option('--purity_file',
              '-pur',
              required=False,              
              help="Path to purity file from purple")
@click.option('--ccf_threshold_clonals',
              '-ccf',
              required=False,
              default=None,              
              help="CCF threshold to select clonals")   


def prepare_table_annotations(input_dir,output_dir,tumor_sample_id,normal_sample_id,caller,gnomad_threshold,
                              cnv_file,purity_file,ccf_threshold_clonals):
    
    print('Loading table')
    
    files = os.listdir(input_dir)
    
    if caller == 'hc':
        pt1_df = pd.DataFrame() 
        for file in files:
            df1 = pd.read_csv(input_dir+file,sep='\t')
            pt1_df = pd.concat([pt1_df,df1],ignore_index=True)
    elif caller == 'sage' or caller == 'intersect':
        file = files[0]
        pt1_df = pd.read_csv(input_dir+file,sep='\t')
    else:
        print('Specify the caller: sage, hc or intersect')
    
    print('Select gnomADg_AF<'+str(gnomad_threshold))    
    pt1_df = pt1_df[pt1_df['gnomADg_AF']<gnomad_threshold]
    
    print('Annotate damaging')
    pt1_df['Damaging'] = pt1_df['Consequence'].progress_apply(lambda x: damaging_col(x))
    pt1_df['mut'] = pt1_df['#CHROM']+':'+pt1_df['POS'].astype(str)+':'+pt1_df['REF']+'>'+pt1_df['ALT']
    pt1_df.drop_duplicates(subset='mut',inplace=True,keep='first')
    
    print('Annotate aa_change')
    pt1_df['aa_change'] = pt1_df.progress_apply(lambda row: aa_change_col(row),axis=1)
    
    print('Annotate real AF')
    pt1_df['n_AF_real'] = pt1_df.progress_apply(lambda row: real_af(row),axis=1)
    
    print('Annotate intogen drivers')
    drivers_df = pd.read_csv('/workspace/datasets/intogen/runs/v2023/20230224_release2023/run/intogen_analysis/unique_drivers.tsv',sep='\t')
    drivers = drivers_df['SYMBOL'].tolist()
    pt1_df['intogen'] = pt1_df['SYMBOL'].progress_apply(lambda x: True if x in drivers else False)
    
    print('Annotate germline genes from Rahman et al Nature 2014')
    germ_old_df = pd.read_excel('/workspace/projects/sjd_pediatric_tumors/data/germline_muts_data/41586_2014_BFnature12981_MOESM96_ESM.xlsx')
    germ_old = germ_old_df['Gene  Symbol'].tolist()
    pt1_df['germline'] = pt1_df['SYMBOL'].progress_apply(lambda x: True if x in germ_old else False)
    
    print('Annotate germline genes from MSKCC, Fiala et al Nature 2021')
    germ_new_df = pd.read_excel('/workspace/projects/sjd_pediatric_tumors/data/germline_muts_data/43018_2021_172_MOESM2_ESM.xlsx',sheet_name='Supp Table 3',skiprows=1)
    germ_new = germ_new_df['Gene'].tolist()
    pt1_df['germline_mskcc'] = pt1_df['SYMBOL'].progress_apply(lambda x: True if x in germ_new else False)
    
    print('Annotate germline genes from Akhavanfard et al NatCom 2020')
    germ_akh_df = pd.read_excel('/workspace/projects/sjd_pediatric_tumors/data/germline_muts_data/41467_2020_16067_MOESM5_ESM.xlsx',skiprows=2)
    germ_akh = germ_akh_df['Gene  Symbol'].tolist()
    pt1_df['germline_akh'] = pt1_df['SYMBOL'].apply(lambda x: True if x in germ_akh else False)
    
    print('Annotate cancer gene role\n')
    role_df = pd.read_csv('/workspace/projects/sjd_pediatric_tumors/data/consensus_role_cancer_genes.tsv',sep='\t')
    role_df = role_df[['Gene.Symbol','consensus_role']]
    role_df = role_df.rename(columns={'Gene.Symbol':'SYMBOL','consensus_role':'role'})
    pt1_df = pd.merge(pt1_df,role_df,how='left')
    
    print('Annotate variant_type')
    pt1_df['variant_type'] = pt1_df['Consequence'].progress_apply(lambda x: variant_type(x,variants_dict))
    
    print(caller)
    if caller == 'sage' or caller =='hc':
        pt1_df.to_csv(output_dir+normal_sample_id+'_filt.maf.gz',sep='\t',compression='gzip',index=None)
    elif caller == 'intersect':
        #Calculate CCF
        #Get cnv and purity values from purple
        cnv_df = pd.read_csv(cnv_file,sep='\t')
        purity = pd.read_csv(purity_file,sep='\t')
        purity = purity['purity'][0]
        
        #Get df with cn and snv intersect
        if 'purple' in cnv_file:
            caller_cn = 'purple'
        else:
            caller_cn = 'ascat'
        pt1_df = cnv_snvs_intersect(cnv_df,pt1_df,caller_cn)

        #Calculate AF and CCF
        pt1_df['t_AF'] = pt1_df.apply(lambda row: calculate_VAF(row['t_alt_reads'],row['t_ref_reads']),axis=1)
        pt1_df['n_AF'] = pt1_df.apply(lambda row: calculate_VAF(row['n_alt_reads'],row['n_ref_reads']),axis=1)
        pt1_df['t_CCF'] = pt1_df.apply(lambda row: calculate_ccf(row['t_AF'],row['CN'],purity),axis=1)
        pt1_df['n_CCF'] = pt1_df.apply(lambda row: calculate_ccf(row['n_AF'],row['CN'],purity),axis=1)
        pt1_df['clonal'] = pt1_df['t_CCF'].apply(lambda x: True if x > float(ccf_threshold_clonals) else False)
        
        pt1_df.to_csv(output_dir+tumor_sample_id+'_vs_'+normal_sample_id+'_filt.maf.gz',sep='\t',compression='gzip',index=None)

        
if __name__ == "__main__":
    prepare_table_annotations()