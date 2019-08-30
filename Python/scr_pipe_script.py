import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import glob

home_dir="/projects/ucar-lab/danaco/bncmrk-dblts"

os.chdir(home_dir+"/Python")

from scr_pipe import *

pr_names=list(["CZI.PBMC","PBMC.8.HTO","Four.Cell.12.HTO"])

for i in pr_names:
    os.chdir(home_dir+"/Scrublet/input/"+i)
    fl=glob.glob("*raw.counts.for.Scrublet.mtx")
    for j in fl:
        scr_pipe(file_name=j)
        
        

    