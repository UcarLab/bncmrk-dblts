import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd 

#def scr_pipe(file_name,gene_name,exp_rate=0.1):
def scr_pipe(file_name,exp_rate=0.1):
    counts_matrix = scipy.io.mmread(file_name).T.tocsc()
    if "/" in file_name:
        home_dir="/".join(file_name.split("/")[0:len(file_name.split("/"))-5])
        pr_name=file_name.split("/")[-2]
        file_name=file_name.split("/")[-1]
    else :
        home_dir="/".join(os.getcwd().split("/")[0:len(os.getcwd().split("/"))-3])
        pr_name=os.getcwd().split("/")[-1]
 
    name=file_name.split(".raw")[0]
    out_dir=home_dir+"/Scrublet/output/"+pr_name+"/"
    #genes = np.array(scr.load_genes(home_dir+data_dir+pr_name + 'genes.tsv', delimiter='\t', column=1))
    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=exp_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                       min_cells=3, 
                                                       min_gene_variability_pctl=85, 
                                                       n_prin_comps=30)
    print(out_dir+name+".Scr.doublet.scores."+str(exp_rate)+"exp_rate.csv")
    pd.DataFrame(doublet_scores).to_csv(out_dir+name+".Scr.doublet.scores."+str(exp_rate)+"exp_rate.csv")