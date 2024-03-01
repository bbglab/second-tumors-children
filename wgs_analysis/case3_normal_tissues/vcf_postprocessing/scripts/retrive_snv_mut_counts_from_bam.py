import click
import json
import pysam
import pandas as pd
from tqdm import tqdm

tqdm.pandas()

def extract_variant_info(row, sample, suffix, sam):
    
    chrom = row['CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']
    
    #extract read counts from sam file
    mut_counts = sam.count_coverage(chrom,int(pos)-1,int(pos), quality_threshold = 0)
    mut_dict = {}
    mut_dict['A'] = mut_counts[0][0]
    mut_dict['C'] = mut_counts[1][0]
    mut_dict['G'] = mut_counts[2][0]
    mut_dict['T'] = mut_counts[3][0]
    ref_reads = mut_dict[ref]
    alt_reads = mut_dict[alt]
    depth = mut_counts[0][0] + mut_counts[1][0] + mut_counts[2][0] + mut_counts[3][0]
    if depth == 0:
    	AF = 0
    else:
    	AF = round((alt_reads / depth),4)

    # extract info of reads
    row['AF'+suffix] = float(AF)
    row['ref_reads'+suffix] = int(ref_reads)
    row['alt_reads'+suffix] = int(alt_reads)
    row['DP'+suffix] = int(depth)

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
    '''recount reads from snv mutations directly from bams using pysam.
    input: table with mutations, sample dictionary and column suffix dictionary, path to HMF or sarek pipeline
    output: table with added columns, alt_reads, ref_reads, DP, AF, tissue'''

    muts_df = pd.read_csv(muts_file,sep='\t')
    snv_df = muts_df[muts_df['mut_type']=='snv']

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
        
        #extract variant info
        suffix = suffixes_dict[sample]
        snv_df = snv_df.progress_apply(lambda row: extract_variant_info(row, sample, suffix, sam_file),axis=1)
        snv_df[tissue+'_old'] = snv_df[tissue]
        snv_df[tissue] = snv_df['alt_reads'+suffix].apply(lambda x: True if x != 0 else False)

    snv_df.to_csv(output_file,sep='\t',compression='gzip',index=None)

if __name__ == '__main__':
    main()
