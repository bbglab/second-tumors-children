import click
import pandas as pd
import json

sorter_context = ['ACAA', 'ACCA', 'ACGA', 'ACTA', 'CCAA', 'CCCA', 'CCGA', 'CCTA', 'GCAA', 'GCCA', 'GCGA', 'GCTA', 'TCAA', 'TCCA',
	'TCGA', 'TCTA', 'ACAG', 'ACCG', 'ACGG', 'ACTG', 'CCAG', 'CCCG', 'CCGG', 'CCTG', 'GCAG', 'GCCG', 'GCGG', 'GCTG', 'TCAG', 'TCCG', 'TCGG',
	'TCTG', 'ACAT', 'ACCT', 'ACGT', 'ACTT', 'CCAT', 'CCCT', 'CCGT', 'CCTT', 'GCAT', 'GCCT', 'GCGT', 'GCTT', 'TCAT', 'TCCT', 'TCGT', 'TCTT',
	'ATAA', 'ATCA', 'ATGA', 'ATTA', 'CTAA', 'CTCA', 'CTGA', 'CTTA', 'GTAA', 'GTCA', 'GTGA', 'GTTA', 'TTAA', 'TTCA', 'TTGA', 'TTTA', 'ATAC',
	'ATCC', 'ATGC', 'ATTC', 'CTAC', 'CTCC', 'CTGC', 'CTTC', 'GTAC', 'GTCC', 'GTGC', 'GTTC', 'TTAC', 'TTCC', 'TTGC', 'TTTC', 'ATAG', 'ATCG',
	'ATGG', 'ATTG', 'CTAG', 'CTCG', 'CTGG', 'CTTG', 'GTAG', 'GTCG', 'GTGG', 'GTTG', 'TTAG', 'TTCG', 'TTGG', 'TTTG']

@click.command()
@click.option('--pcawg_sigs_file',
                type=click.Path(exists=True),
                help="PCAWG signatures file",
                required=False,
                default=None)
@click.option('--cosmic_sigs_file',
                type=click.Path(exists=True),
                help="COSMIC signatures file",
                required=False,
                default=None)
@click.option('--out_file',
                type=click.Path(),
                help="output file",
                required=True)
@click.option('--sigs',
				type=click.STRING,
				help='Signatures selected, sepparated by comma: e.g. SBS1,SBS5,SBS18',
				required=True)
@click.option('--adjust_contexts',
				is_flag=True,
				default=False,
				help='Specify if the signatures have to be adjusted to a subset of genomic regions (e.g. panel)',
				required=False)
@click.option('--context_panel',
				type=click.Path(exists=True),
				default=None,
				help='Contexts count file from panel, contexts as rows, no column',
				required=False)
@click.option('--context_genome',
				type=click.Path(exists=True),
				default=None,
				help='Contexts count file from genome, contexts as rows, no column',
				required=False)

def run(pcawg_sigs_file,cosmic_sigs_file, out_file, sigs, adjust_contexts, context_genome, context_panel):

	if pcawg_sigs_file != None:
		#Read table with pcawg signatures
		sigs_df = pd.read_csv(pcawg_sigs_file,sep='\t')
		sigs_df = sigs_df.T

		#Reformat index
		channels = sigs_df.index.tolist()
		channels = [c[0]+c[2]+c[6]+c[4] for c in channels]

		sigs_df['index'] = channels
		sigs_df = sigs_df.set_index('index')
		sigs_df = sigs_df.reindex(sorter_context)

	elif cosmic_sigs_file != None:
		#Read table with signatures from cosmic
		sigs_df = pd.read_csv(cosmic_sigs_file,sep='\t')
		sigs_df['index'] = sigs_df['Type'].apply(lambda c: c[0]+c[2]+c[6]+c[4])
		sigs_df.drop('Type',inplace=True,axis=1)
		sigs_df.set_index('index',inplace=True)
		sigs_df = sigs_df.reindex(sorter_context)
	else:
		print('Please specify the signature file: from pcawg or cosmic')

	#Select columns with corresponding signatures
	sigs = sigs.split(',')
	selected_sigs_df = sigs_df[sigs]

	#Transform signatures proportional to the context in the subset of genomic regions (e.g. panel)
	if adjust_contexts == True:

	    with open (context_genome) as f:
	        context_genome_dict = json.load(f)
	    with open (context_panel) as f:
	        context_panel_dict = json.load(f)

	    cols = list(selected_sigs_df)
	    
	    selected_sigs_df['channel'] = selected_sigs_df.index
	    
	    def correct_context(row,col,context_genome_dict,context_panel_dict):
	        #gcontext = context_genome_dict[x.index]
	        channel = row['channel']
	        context_genome_count = context_genome_dict[channel]
	        context_panel_count = context_panel_dict[channel]
	        value = row[col]
	        new_value = value*context_panel_count/context_genome_count
	        return new_value
	    
	    for col in cols:
	        selected_sigs_df[col] = selected_sigs_df.apply(lambda row: correct_context(row,col,context_genome_dict,context_panel_dict),axis=1)
	    
	    selected_sigs_df.drop(columns='channel',inplace=True)

	#Fix float number; all channels from each signature must be == 1

	col_dict = selected_sigs_df.sum().to_dict()

	def fix_proportion_num (x,col):
		total_sum = col_dict[col]
		new_col = x/total_sum
		return new_col

	cols = list(selected_sigs_df)
	
	for col in cols:
		selected_sigs_df[col] = selected_sigs_df[col].apply(lambda x:fix_proportion_num(x,col))

	#Save table
	selected_sigs_df.to_csv(out_file,sep='\t',header=True,index=True)

if __name__ == '__main__':
	run()
