import click
import json
import pysam
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool

tqdm.pandas()

def extract_snvs_info(row, suffix, sam_file, tissue):   

    chrom = row['CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']
    
    #extract read counts from sam file
    mut_counts = sam_file.count_coverage(chrom,int(pos)-1,int(pos), quality_threshold = 0)
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

    return row[['CHROM','POS','REF','ALT','AF'+suffix, 'ref_reads'+suffix, 'alt_reads'+suffix,'DP'+suffix, tissue ]]

def extract_indels_info(row,suffix,sam_file,tissue):
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
    return row[['CHROM','POS','REF','ALT','AF'+suffix, 'ref_reads'+suffix, 'alt_reads'+suffix,'DP'+suffix, tissue ]]

#create "intermediate" functions, because in the paralelized function we cannot put multiple variables in the function
def snvs_info (group):
    name, chr_df = group

    sam_file = pysam.AlignmentFile(path_to_bam, "rb")

    chr_df = chr_df.apply(lambda row: extract_snvs_info (row, suffix, sam_file, tissue),axis=1)
    return chr_df

def indel_info (group):
    name, chr_df = group
    
    sam_file = pysam.AlignmentFile(path_to_bam, "rb")

    chr_df = chr_df.apply(lambda row: extract_indels_info(row, suffix, sam_file, tissue),axis=1)
    return chr_df

@click.command()
@click.option('--muts_file',
    help='path to file containing the mutations, with these columns: ["CHROM", "POS", "REF", "ALT","tissue"]')
@click.option('--samples_json',
    help='path to json file with the sample dictionary, e.g. "AQ5174":"blood"')
@click.option('--suffixes_json',
    help='path to json file with the suffixes dictionary (for table columns), "AQ5174":"_b"')
@click.option('--platinum_path',
    help='path to the output from HMF pipeline',
    default=None)
@click.option('--sarek_path',
    help='path to the output from sarek pipeline',
    default=None)
@click.option('--bams_json',
    help='path to json file with path to bam for each sample',
    default=None)
@click.option('--variant_type',
    help='variant type to be annotated: snv or indel')
@click.option('--output_file',
    help='path to the output file')
@click.option('--cpus',
              '-c',
              required=True,
              default=1,
              type=int,
              help="Number of cpus to be used")

def main (muts_file,samples_json,suffixes_json,platinum_path,sarek_path,bams_json,output_file,variant_type,cpus):
    '''recount reads from snv mutations directly from bams using pysam.
    input: table with mutations, sample dictionary and column suffix dictionary, path to HMF or sarek pipeline
    output: table with added columns, alt_reads, ref_reads, DP, AF, tissue'''

    muts_df = pd.read_csv(muts_file,sep='\t')
    variant_df = muts_df[muts_df['mut_type']==variant_type].drop_duplicates(subset=['CHROM','POS','REF','ALT'],keep='first')

    with(open(samples_json,'r')) as f:
        samples_dict = json.load(f)
    with(open(suffixes_json,'r')) as f:
        suffixes_dict = json.load(f)

    final_df = pd.DataFrame(columns={'CHROM','POS','REF','ALT'})

    for sample in tqdm(samples_dict.keys()):
        
        #define the path to the bam file
        global path_to_bam
        global tissue

        tissue = samples_dict[sample]

        if bams_json != None:
            with(open(bams_json,'r')) as f:
                bams_dict = json.load(f)           
            path_to_bam = bams_dict[sample]
        else:
            if tissue in ['tumor1','blood']:
                path_to_bam = platinum_path+sample+'/aligner/'+sample+'.bam'
            elif tissue == 'tumor2':
                platinum_path = platinum_path.replace('-t1','-t2')
                path_to_bam = platinum_path+sample+'/aligner/'+sample+'.bam'
            else:
                path_to_bam = sarek_path+'preprocessing/recalibrated/'+sample+'/'+sample+'.recal.bam'
        
        #define the suffix to use in the columns per sample
        global suffix
        suffix = suffixes_dict[sample]
        
        #define the function per variant type
        if variant_type == 'snv':
            f = snvs_info
        elif variant_type == 'indel':
            f = indel_info    
        else:
            print('Please specify the variant type: snv or indel')
        
        #run the function in a paralelized mode: group the df per chromosome, run function in paralell and then concat
        sample_df = pd.DataFrame()
        with Pool(int(cpus)) as pool:

            for chr_df in tqdm(pool.imap(f, variant_df.groupby("CHROM")),
                               total=len(variant_df.groupby("CHROM"))):
                sample_df = pd.concat([sample_df,chr_df], ignore_index=True)

        final_df = pd.merge(final_df,sample_df, how='outer')

        print(len(final_df))


    for sample in tqdm(samples_dict.keys()):

        tissue = samples_dict[sample]

        suffix = suffixes_dict[sample]

        final_df[tissue+'_old'] = final_df[tissue]
        final_df[tissue] = final_df['alt_reads'+suffix].apply(lambda x: True if x != 0 else False)


    final_df.to_csv(output_file,sep='\t',compression='gzip',index=None)

if __name__ == '__main__':
    main()
