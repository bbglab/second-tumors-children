import pandas as pd
import click
import allel
import pybedtools

def get_positions (row):
    alt_1 = row['ALT_1']
    if isinstance(alt_1,str):     
        if '[' in alt_1:
            chrom = alt_1.split('[')[1].split(':')[0]
            pos = alt_1.split('[')[1].split(':')[1]
        elif ']' in alt_1:
            chrom = alt_1.split(']')[1].split(':')[0]
            pos = alt_1.split(']')[1].split(':')[1]
        else:
            chrom = row['CHROM']
            pos = row['POS']
    else:
        chrom = '-'
        pos = '-'
    row['chrom_alt'] = chrom
    row['pos_alt'] = pos
    return row

def get_fusions(row,positions):
    chrom_1 = str(row['chrom_alt'])
    pos_1 = str(row['pos_alt'])
    symbol_1 = str(row['Hugo_Symbol'])
    
    if symbol_1=='nan':
        symbol_1 = '-'
        
    for pos in positions:
        chrom_2 = str(pos[0])
        pos_2 = str(pos[1])
        symbol_2 = str(pos[2])
        
        if symbol_2 == 'nan':
            symbol_2 = '-'
        if chrom_2 == chrom_1 and pos_2 == pos_1:
            if symbol_1 == symbol_2 and symbol_1 != '-':
                pos_ref = row['POS']
                if int(pos_ref) == int(pos_1):
                    fusion = symbol_1 + '-ins'
                else:
                    fusion = symbol_1 + '-del'

            else:
                fusion = symbol_1 + '/' + symbol_2
            break
        else:
            fusion = '-'
            
    return fusion


file = '/workspace/projects/sjd_pediatric_tumors/data/hg38.chrom.sizes'
with open(file) as f:
    rows = (line.split('\n')[0] for line in f)
    rows = (line.split('\t') for line in rows)
    chrom_sizes = {row[0]:row[1] for row in rows}


def get_distance(row):    
    chrom1 = row['CHROM']
    chrom2 = row['chrom_alt']
    if chrom1!='nan' and chrom2!='-':     
        if chrom1 == chrom2:
            pos1 = row['POS']
            pos2 = row['pos_alt']
            dist = int(pos2) - int(pos1)
            chrom_size = int(chrom_sizes[chrom1])
            dist_rel = round(dist / chrom_size, 2)
        else:
            dist = '-'
            dist_rel = '-'
    else:
        dist = '-'
        dist_rel = '-'
        
    row['distance'] = dist
    row['distance_rel'] = dist_rel
    return row

def annot_sv(row):
    chr1 = row['CHROM']
    chr2 = row['chrom_alt']
    fusion = row['fusion']
    if 'del' in fusion:
        sv_type = 'del'
    elif 'ins' in fusion:
        sv_type = 'ins'
    elif chr1 == chr2:
        sv_type = 'inv'
    elif '-' not in fusion:
        sv_type = 'fusion'
    else:
        sv_type = 'other'
    return sv_type
    
def annot_chrom (row):
    chr1 = row['CHROM']
    chr2 = row['chrom_alt']
    if chr1 == chr2:
        chrom = chr1
    else:
        chrom = chr1 + '/' + chr2
        
    return chrom


def gene_info_annotations(pt1_df):
    
    #Annotate mut column
    print(pt1_df.dtypes)
    pt1_df['mut'] = pt1_df['CHROM'].astype(str)+':'+pt1_df['POS'].astype(str)+':'+pt1_df['REF']+'>'+pt1_df['ALT']
    
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
    
    print('Annotate cancer gene census translocations')
    cgc_transl_df = pd.read_csv('/workspace/projects/sjd_pediatric_tumors/data/Census_transFri_Apr_9_10_22_06_2021.tsv',sep='\t')
    cgc_transl = cgc_transl_df['Gene Symbol'].tolist()
    pt1_df['cgc_transl'] = pt1_df['SYMBOL'].apply(lambda x: True if x in cgc_transl else False)

    return pt1_df

def annotate_sv_genes (sv_df,pos_df):

    #Preapre cnv_df
    sv_bt_df = sv_df.copy(deep=True)
    sv_bt_df['POS2'] = sv_df['POS']+1
    sv_bt_df = sv_bt_df[['CHROM','POS','POS2','ID','REF','ALT_1','ALT_2','ALT_3','QUAL','FILTER_PASS']]

    #Prepare snvs_df
    pos_bt_df = pos_df.copy(deep=True)
    pos_bt_df['chr'] = 'chr' + pos_bt_df['chr']
    pos_bt_df = pos_bt_df[['chr', 'gene_start', 'gene_end','Hugo_Symbol', 'Feature']]

    #pybedtools cnv + snvs
    sv_bt = pybedtools.BedTool.from_dataframe(sv_bt_df)
    pos_bt = pybedtools.BedTool.from_dataframe(pos_bt_df)

    sv_and_pos = pos_bt.intersect(sv_bt)

    #Pass intersection to dataframe
    sv_pos_df = sv_and_pos.to_dataframe()
    sv_pos_df = sv_pos_df.rename(columns={'chrom':'CHROM','start':'POS','name':'Hugo_Symbol','score':'Feature'})
    sv_pos_df = sv_pos_df[['CHROM','POS','Hugo_Symbol','Feature']]
    sv_final_df = pd.merge(sv_df,sv_pos_df,how='left')
    
    #Annotate fusions
    sv_final_df = sv_final_df.apply(lambda row: get_positions(row),axis=1)
    positions = sv_final_df[['CHROM','POS','Hugo_Symbol']].values.tolist()
    sv_final_df['fusion'] = sv_final_df.apply(lambda row: get_fusions(row,positions),axis=1)
    
    sv_final_df = sv_final_df.apply(lambda row: get_distance(row),axis=1)
    
    #Fix dels
    dels = sv_final_df['fusion'][sv_final_df['fusion'].str.contains('-del')].drop_duplicates().tolist()
    del_df = pd.DataFrame()
    for d in dels:
        d_df = sv_final_df[sv_final_df['fusion']==d]
        d_pos = d_df['POS'].tolist()
        d_pos = [str(p) for p in d_pos]
        d_pos_alt = d_df['pos_alt'].tolist()
        for i in range(0,len(d_pos)):
            if i == len(d_pos)-1:
                break
            if (d_pos[i] == d_pos_alt[i+1]) and (d_pos[i+1] == d_pos_alt[i]):
                row = d_df[i:i+1]
                del_df = pd.concat([del_df,row],ignore_index=True)
    sv_final_df = sv_final_df[~sv_final_df['fusion'].str.contains('-del')]
    sv_final_df = pd.concat([sv_final_df,del_df],ignore_index=True)
    
    #Annotate chrom info:
    sv_final_df['chr/chr'] = sv_final_df.apply(lambda row: annot_chrom(row),axis=1)
    
    #Order by chrom and pos
    chroms = list(range(1,22))
    chroms = ['chr'+str(c) for c in chroms]
    sex_chroms = ['chrX','chrY']
    chroms = chroms + sex_chroms
    t = pd.CategoricalDtype(categories=chroms, ordered=True)
    sv_final_df['CHROM']=pd.Series(sv_final_df['CHROM'], dtype=t)
    sv_final_df.sort_values(by=['CHROM','POS'],inplace=True)
    
    #Rename some columns
    sv_final_df = sv_final_df.rename(columns={'ALT_1':'ALT','Hugo_Symbol':'SYMBOL'})
    
    #Annotate sv info:
    sv_final_df['sv_type'] = sv_final_df.apply(lambda row: annot_sv(row),axis=1)
    
    #Gene info annotations (role, intogen, germline)
    sv_final_df = gene_info_annotations(sv_final_df)
    
    #Return usefull info
    sv_final_df = sv_final_df[['SYMBOL', 'germline', 'germline_mskcc', 'intogen', 'role','mut',
           'fusion','cgc_transl','chr/chr','sv_type', 'distance','distance_rel']]
  
    return sv_final_df

@click.command()

@click.option('--input_vcf',
              '-i',
              required = True,
              help="path to vcf file from caller: gridds")
@click.option('--output_dir',
              '-o',
              required = True,
              help="path to output folder")
@click.option('--genomic_positions_file',
              '-gp',
              required=True,
              help="path to file with the annotation of the genomic positions of each gene")
@click.option('--canonical_transcripts_file',
              '-ct',
              required=True,
              help="path to file with the annotation of the canonical transcripts with ensembl code")
@click.option('--tumor_sample_id',
              '-t_id',
              required=False,
              default = None,
              help="Tumor sample id")
@click.option('--normal_sample_id',
              '-n_id',
              required=True,
              help="Normal sample id")
@click.option('--is_germline',
              required=False,
              is_flag = True,
              help="GRIPPS germline analysis")

def cli(input_vcf, output_dir, genomic_positions_file, canonical_transcripts_file, tumor_sample_id, normal_sample_id,is_germline):
    
    #Prepare gene info annotations
    pos_df = pd.read_csv(genomic_positions_file,sep='\t')
    pos_df = pos_df.rename(columns={'Transcript stable ID':'Feature','Chromosome/scaffold name':'chr','Gene start (bp)':'gene_start','Gene end (bp)':'gene_end'})
    pos_df = pos_df[['Feature','chr','gene_start','gene_end']].drop_duplicates()
    can_df = pd.read_csv(canonical_transcripts_file,header=None,sep='\t')
    can_df = can_df.rename(columns={0:'Gene_id',1:'Feature',2:'Hugo_Symbol'})
    can_df = can_df.drop(columns='Gene_id')
    pos_df = pd.merge(pos_df,can_df,how='left')
    pos_df = pos_df[['Hugo_Symbol','Feature','chr','gene_start','gene_end']]
    
    #Read vcf file and filter PASS
    sv_df = allel.vcf_to_dataframe(input_vcf)
    sv_df = sv_df[sv_df['FILTER_PASS']==True]
    
    #Run annotation function
    t1_df = annotate_sv_genes(sv_df,pos_df)
    
    #Save results
    if is_germline == True:
    	t1_df.to_csv(output_dir+normal_sample_id+'.maf.gz',compression='gzip',sep='\t',index=None)
    else:
    	t1_df.to_csv(output_dir+tumor_sample_id+'_vs_'+normal_sample_id+'.maf.gz',compression='gzip',sep='\t',index=None)
    

if __name__ == "__main__":
    cli()