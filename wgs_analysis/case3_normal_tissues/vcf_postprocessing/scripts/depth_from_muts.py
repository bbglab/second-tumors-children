import pandas as pd
import gzip
from tqdm import tqdm
import numpy as np
import json
import click


@click.command()
@click.option('--input_somatic_muts',
                help='path to file containing the positions of the mutations. Has columns CHROM and POS')
@click.option('--input_folder', 
                help='path to the folder containing the depth, output from samtools.')
@click.option('--output_file', 
                help='path to the folder output: it will generate a tsv file with depth by position fo mutations')
                
def main(input_somatic_muts,input_folder, output_file):
    '''Reads a file with mutations and gets the depth counts
    from samtools file
     input: path to mutations file, path to the depth file (output from samtools)
     output: json file with depth per mutation position'''

    mut_df = pd.read_csv(input_somatic_muts,sep='\t')

    chroms = list(range(1,23))
    chroms = ['chr'+str(chrom) for chrom in chroms]
    chroms.append('chrX')
    chroms.append('chrY')

    all_df = pd.DataFrame()
    
    for chrom in tqdm(chroms):  
        
        positions = mut_df['POS'][mut_df['CHROM']==chrom].tolist()
        input_file = input_folder + chrom + '_depth.txt.gz'  
        
        df = pd.read_csv(input_file,sep='\t',header=None)
        df = df[df[1].isin(positions)]
        all_df = pd.concat([all_df,df],ignore_index=True)

    all_df.to_csv(output_file+'.gz',sep='\t',index=None,header=None,compression='gzip')


if __name__ == '__main__':
    main()
