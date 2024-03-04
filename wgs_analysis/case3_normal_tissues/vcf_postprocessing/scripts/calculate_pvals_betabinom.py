import click
import json
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors

import scipy.stats as stats
import numpy as np
import os
from functools import partial

tqdm.pandas()

def betabinom_pval(row, alt_blood_cols, DP_blood_cols, suffixes_dict, tolerance):

	alternate_allele_count = row[alt_blood_cols].to_numpy()
	alternate_allele_count = tf.convert_to_tensor(alternate_allele_count, dtype=tf.float32)

	depth_count = row[DP_blood_cols].to_numpy()
	depth_count = tf.convert_to_tensor(depth_count, dtype=tf.float32)

	# define the variables for the optimization problem

	alpha = tfp.util.TransformedVariable(.5, tfp.bijectors.Log(), name='alpha')

	def get_model(total_counts, alpha_var):
		return tfd.BetaBinomial(total_counts, alpha_var, tolerance - alpha_var)

	# define model

	model = get_model(depth_count, alpha)

	# define convergence criterion    

	convergence_criterion = tfp.optimizer.convergence_criteria.LossNotDecreasing(rtol=0.25, window_size=5, min_num_steps=10, name=None)

	# train model (MLE)

	num_steps = 100
	res = tfp.math.minimize(
		loss_fn=lambda: -tf.reduce_sum(model.log_prob(alternate_allele_count)),
		num_steps=num_steps,
		optimizer=tf.optimizers.Adam(learning_rate=1e-1), 
		trainable_variables=model.trainable_variables,
		convergence_criterion=convergence_criterion
	)

	# probability distribution of alternate allele counts given the obs_depth
	alpha_hat = model.trainable_variables[0].numpy()
	beta_hat = tf.cast(tolerance - alpha_hat, dtype=tf.float32)
    
	for sample in suffixes_dict.keys():

		suffix = suffixes_dict[sample]
		alt_read_col = 'alt_reads'+ suffix
		DP_col = 'DP'+suffix

		obs_alt = row[alt_read_col]
		obs_depth = row[DP_col]

		mle_betabinom = tfd.BetaBinomial(obs_depth, alpha_hat, beta_hat)

		## p-value: Prob(X >= obs_alt) ##
		# sum of discrete probs higher or equal than obs_alt

		onesided_pval = sum([mle_betabinom.prob(x) for x in range(int(obs_alt), int(obs_depth) + 1)])

		row['pval_sa'+suffix] = onesided_pval.numpy()

	return row

def process_chunk(chunk_filename):
	chunk = pd.read_csv(chunk_filename,sep='\t')
	# Do some processing on the chunk
	processed_chunk = chunk.apply(lambda row: betabinom_pval(row, alt_blood_cols, DP_blood_cols, suffixes_dict, tol),axis=1)

	# Save the processed chunk to a temporary file
	processed_chunk.to_csv(chunk_filename, index=False,sep='\t')


@click.command()
@click.option('--muts_file',
	help='path to file containing the mutations, with these columns: ["CHROM", "POS", "REF", "ALT"], all the alt_reads and DP cols from all tissues'
	'and all the alt_reads and DP cols from the cohort (blood samples')
@click.option('--suffixes_json',
	help='path to json file with the suffixes dictionary (for table columns), "AQ5174":"_b"')
@click.option('--suffixes_blood_json',
	help='path to json file with the suffixes dictionary from the cohort/blood samples (for table columns), "AQ5174":"_b1"')
@click.option('--tolerance',
	help='tolerance level to perform the betabinomial test. Higher the tolerance, more samples are significant, and less artifacts are flagged',
	default=30,
	type=int)
@click.option('--output_file',
	help='path to the output file')
@click.option('--chunk_size',
	required=False,
	default=1000,
	type=int,
	help='size of the chunk of rows to paralelize')

def main (muts_file,suffixes_json,suffixes_blood_json,tolerance,output_file,chunk_size):
	'''performs betabinomial test.
	input: table with mutations info, column suffix dictionary of both tissues and blood cohort, tolerance, cpus
	output: table with added columns: pvals_sa_<suffix>'''

	muts_df = pd.read_csv(muts_file,sep='\t')

	global suffixes_dict

	with(open(suffixes_json,'r')) as f:
		suffixes_dict = json.load(f)
	with(open(suffixes_blood_json,'r')) as f:
		suffixes_blood_dict = json.load(f)

	global alt_blood_cols
	global DP_blood_cols

	alt_blood_cols = ['alt_reads'+suffixes_blood_dict[sample] for sample in suffixes_blood_dict.keys() if sample not in ['AQ5174']]

	DP_blood_cols = ['DP'+suffixes_blood_dict[sample] for sample in suffixes_blood_dict.keys() if sample not in ['AQ5174']]

	global tol

	tol = tolerance

	final_df = pd.DataFrame()

	#run the function in a paralelized mode: group the df per chromosome, run function in paralell and then concat


	def parallelize_dataframe(temp_dir, df, chunk_size, num_processes=mp.cpu_count()):
		# Create temporary directory for storing processed chunks
		os.makedirs(temp_dir, exist_ok=True)

		# Split the dataframe into chunks based on number of processes
		print('Saving temporary files...')

		for i in tqdm(range(0, len(df), chunk_size)):
			chunk = df.iloc[i:i + chunk_size]
			chunk_filename = os.path.join(temp_dir, f"{chunk.index[0]}_{chunk.index[-1]}.tsv")
			chunk.to_csv(chunk_filename, index=False,sep='\t')

		chunk_filenames = [temp_dir +"/"+ filename for filename in os.listdir(temp_dir)]


		print('Multiprocess betabinomial per file....')

		# Create pool of processes to process each chunk in parallel
		with mp.Pool(num_processes) as pool, tqdm(total=len(chunk_filenames)) as progress_bar:
			for _ in pool.imap_unordered(process_chunk, chunk_filenames):
				progress_bar.update(1)

		print('Combine all tables...')

		# Combine the processed chunks into a single dataframe
		processed_chunks = [pd.read_csv(os.path.join(temp_dir, filename),sep='\t') for filename in os.listdir(temp_dir)]
		processed_df = pd.concat(processed_chunks, ignore_index=True)

		print('Remove all temporary files...')
		# Delete the temporary directory and files
		for filename in os.listdir(temp_dir):
			os.remove(os.path.join(temp_dir, filename))
		os.rmdir(temp_dir)

		return processed_df

	temp_dir = muts_file.split('/')[-1].split('.')[0]
	
	final_df = parallelize_dataframe(temp_dir=temp_dir,df=muts_df,chunk_size=chunk_size)


	final_df.to_csv(output_file,sep='\t',compression='gzip',index=None)

if __name__ == '__main__':
	main()
