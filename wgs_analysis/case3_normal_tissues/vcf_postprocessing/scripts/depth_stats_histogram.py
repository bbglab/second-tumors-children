import gzip
from tqdm import tqdm
import numpy as np
from scipy.stats import t
import json
import matplotlib.pyplot as plt
import click


@click.command()
@click.option('--input_file', 
                help='path to the file depth, output from samtools. Depth is the 3rd column')
@click.option('--output_file_stats', 
                help='path to the file output: json file with stats info')
@click.option('--output_file_hist', 
                help='path to the file output: png file with histogram')
                
def main(input_file, output_file_stats,output_file_hist):
    '''Reads a file with depth counts and retrieves mean, sd
     and confidence intervals at .95, .9, .8, .5
     input: path to the depth file, output from samtools
     output: json file with stats'''
    counts = []
    with gzip.open(input_file,'r') as file:
        for line in tqdm(file):
            line = str(line)
            count = line.split('\\t')[2].split('\\n')[0]
            if int(count) < 400:
                counts.append(int(count))
        file.close()
        
        mean = np.mean (counts)
        sd = np.std (counts)
        
        dof = len(counts)-1 
        confidence = 0.95
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        left_ci_95 = mean-sd*t_crit/np.sqrt(len(count))
        right_ci_95 = mean+sd*t_crit/np.sqrt(len(count))
        confidence = 0.90
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        left_ci_90 = mean-sd*t_crit/np.sqrt(len(count))
        right_ci_90 = mean+sd*t_crit/np.sqrt(len(count))
        confidence = 0.80
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        left_ci_80 = mean-sd*t_crit/np.sqrt(len(count))
        right_ci_80 = mean+sd*t_crit/np.sqrt(len(count))   
        confidence = 0.50
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        left_ci_50 = mean-sd*t_crit/np.sqrt(len(count))
        right_ci_50 = mean+sd*t_crit/np.sqrt(len(count)) 

        #save stats in a json file
        values_dict = {'mean':mean,
    		'sd': sd,
    		'left_ci_95':left_ci_95,
    		'right_ci_95':right_ci_95,
    		'left_ci_90':left_ci_90,
    		'right_ci_90':right_ci_90,
    		'left_ci_80':left_ci_80,
    		'right_ci_80':right_ci_80,
    		'left_ci_50':left_ci_50,
    		'right_ci_50':right_ci_50}

        with open(output_file_stats,'w') as f:
    	    json.dump(values_dict,f)
        f.close()
        
        #historgam
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 4)
        sample = input_file.split('/')[8]
        chrom = input_file.split('/')[9].split('_')[0]
        title = sample + ' ' + chrom
        ax.hist(counts, bins=100, range=(0,400),color='lightblue')
        plt.title(title,size=20,pad=20)
        plt.xlim(left=-50,right=400)
        #plt.axvline(60)
        #plt.axvline(200)
        plt.axvline(mean,color='black')

        plt.axvline(left_ci_95,color='red')
        plt.axvline(right_ci_95,color='red')
        plt.text(x=.9,y=1.05,s='CI=0.95',transform = ax.transAxes,color='red')

        plt.axvline(mean-2*sd,color='blue')
        plt.axvline(mean+2*sd,color='blue')
        plt.text(x=.9,y=1.10,s='2SD',transform = ax.transAxes,color='blue')

        
        plt.savefig(output_file_hist, bbox_inches="tight")


if __name__ == '__main__':
    main()
