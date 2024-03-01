import click
import json
import pysam
import pandas as pd
from tqdm import tqdm
tqdm.pandas()

def extract_indels_info(row,suffix,sam_file):
    chrom = row['CHROM']
    pos1 = int(row['POS'])-2
    isdel = 0
    isins = 0
    no_indel = 0
    for pileupcolumn in sam_file.pileup(chrom, pos1, pos1+1):
        if pileupcolumn.pos == pos1:
            
            for pileupread in pileupcolumn.pileups:
                if pileupread.indel != 0:
                    if pileupread.indel < 0:
                        isdel = isdel +1
                    elif pileupread.indel > 0:
                        isins = isins +1
                else:
                    no_indel = no_indel +1
    row['ref_reads'+suffix] = no_indel
    if isins == 0:
        alt_reads = isdel
        indel_type = 'del'
        depth = no_indel + isdel
    else:
        alt_reads = isins
        indel_type = 'ins'
        depth = no_indel + isins
        
    if depth == 0:
        AF = 0
    else:
        AF = round((alt_reads / depth),4) 
        
    row['alt_reads'+suffix] = alt_reads
    row['indel_type'] = indel_type
    row['DP'+suffix] = depth       
    row['AF'+suffix] = float(AF)
    return row

@click.command()
@click.option('--muts_file',
               help='path to file containing the mutations, with these columns: ["CHROM", "POS", "REF", "ALT","tissue"]')
@click.option('--samples_json',
               help='path to json file with the sample dictionary, e.g. "AQ5174":"blood"')
@click.option('--suffixes_json',
               help='path to json file with the suffixes dictionary (for table columns), "AQ5174":"_b"')
@click.option('--base_path',
               help='path to the output from HMF pipeline')
@click.option('--sarek_path',
               help='path to the output from sarek pipeline',
               default=None)
@click.option('--output_file',
               help='path to the output file')


def main (muts_file,samples_json,suffixes_json,base_path,sarek_path,output_file):
    '''recount reads from indel mutations directly from bams using pysam.
    input: table with mutations, sample dictionary and column suffix dictionary, path to HMF or sarek pipeline
    output: table with added columns, alt_reads, ref_reads, DP, AF, tissue'''
    
    muts_df = pd.read_csv(muts_file,sep='\t')
    indel_df = muts_df[muts_df['mut_type']=='indel']

    with(open(samples_json,'r')) as f:
        samples_dict = json.load(f)
    with(open(suffixes_json,'r')) as f:
        suffixes_dict = json.load(f)

    for sample in tqdm(samples_dict.keys()):
        tissue = samples_dict[sample]
        if tissue in ['tumor1','blood']:
            path_to_bam = base_path+sample+'/aligner/'+sample+'.bam'
        elif tissue == 'tumor2':
            base_path = base_path.replace('-t1','-t2')
            path_to_bam = base_path+sample+'/aligner/'+sample+'.bam'
        else:
            path_to_bam = sarek_path+'preprocessing/recalibrated/'+sample+'/'+sample+'.recal.bam'

        sam_file = pysam.AlignmentFile(path_to_bam, "rb")
        suffix = suffixes_dict[sample]
        indel_df = indel_df.progress_apply(lambda row: extract_indels_info(row,suffix,sam_file),axis=1)
        indel_df[tissue+'_old'] = indel_df[tissue]
        indel_df[tissue] = indel_df['alt_reads'+suffix].apply(lambda x: True if x != 0 else False)
        
    indel_df.to_csv(output_file,sep='\t',compression='gzip',index=None)

if __name__ == '__main__':
    main()
