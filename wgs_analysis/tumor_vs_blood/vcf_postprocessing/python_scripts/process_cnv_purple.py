import pandas as pd
import click


def cna_annotation(cn):
    if cn >= 2.2:
        cna = 'amp'
    elif cn <= 1.8:
        cna = 'del'
    else:
        cna = '-'
    return cna

def gene_info_annotations(pt1_df):
    
    print('Annotate intogen drivers')
    drivers_df = pd.read_csv('/workspace/datasets/intogen/runs/20210108/intogen_20210108/unique_drivers.tsv',sep='\t')
    drivers = drivers_df['SYMBOL'].tolist()
    pt1_df['intogen'] = pt1_df['SYMBOL'].apply(lambda x: True if x in drivers else False)
    
    print('Annotate germline genes')
    germ_old_df = pd.read_excel('/workspace/projects/sjd_pediatric_tumors/data/germline_muts_data/41586_2014_BFnature12981_MOESM96_ESM.xlsx')
    germ_old = germ_old_df['Gene  Symbol'].tolist()
    pt1_df['germline'] = pt1_df['SYMBOL'].apply(lambda x: True if x in germ_old else False)
    
    print('Annotate germline genes from MSKCC')
    germ_new_df = pd.read_excel('/workspace/projects/sjd_pediatric_tumors/data/germline_muts_data/43018_2021_172_MOESM2_ESM.xlsx',sheet_name='Supp Table 3',skiprows=1)
    germ_new = germ_new_df['Gene'].tolist()
    pt1_df['germline_mskcc'] = pt1_df['SYMBOL'].apply(lambda x: True if x in germ_new else False)
    
    print('Annotate cancer gene role\n')
    role_df = pd.read_csv('/workspace/projects/sjd_pediatric_tumors/data/consensus_role_cancer_genes.tsv',sep='\t')
    role_df = role_df[['Gene.Symbol','consensus_role']]
    role_df = role_df.rename(columns={'Gene.Symbol':'SYMBOL','consensus_role':'role'})
    pt1_df = pd.merge(pt1_df,role_df,how='left')

    return pt1_df


@click.command()

@click.option('--input_tsv',
              '-i',
              required = True,
              help="path to tsv file from purple cn, with gene annotation")
@click.option('--output_dir',
              '-o',
              required = True,
              help="path to output folder")
@click.option('--tumor_sample_id',
              '-t_id',
              required=True,
              help="Tumor sample id")
@click.option('--normal_sample_id',
              '-n_id',
              required=True,
              help="Normal sample id")


def cnv_table (input_tsv,output_dir,tumor_sample_id,normal_sample_id):
    df1 = pd.read_csv(input_tsv,sep='\t')
    df1['CNA'] = df1['maxCopyNumber'].apply(lambda x: cna_annotation(x))
    df1 = df1.rename(columns={'chromosome':'CHROM','gene':'SYMBOL','maxCopyNumber':'CN','transcriptId':'Feature',
                              'minMinorAlleleCopyNumber':'CN_min_allele'})
    
    #Gene info annotations (role, intogen, germline)
    df1 = gene_info_annotations(df1)
    
    #Annotate chrom and cytoband:
    
    df1['cytoband'] = df1['CHROM']+':'+df1['chromosomeBand']
    
    #Get useful columns
    df1 = df1[['SYMBOL','germline', 'germline_mskcc', 'intogen', 'role',
            'CNA', 'CN','CN_min_allele','cytoband']]
    
    df1.to_csv(output_dir+tumor_sample_id+'_vs_'+normal_sample_id+'.maf.gz',sep='\t',compression='gzip',index=None)

if __name__ == "__main__":
    cnv_table()