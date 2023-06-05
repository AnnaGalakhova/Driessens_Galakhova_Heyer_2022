# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 09:43:56 2022

proces large matrices of raw counts into log10CPM values 
and write to a new matrix 

@author: stand  & annagalakhova
"""

import pandas as pd
import numpy as np
import dask.dataframe as dd
from pathlib import Path
import os 

filepath = Path(r'C/Users/annagalakhova/PhD INF CNCR VU/DATA/transcriptomics/RNAseq/VISp_mouse')

df = pd.read_csv(r"/Users/annagalakhova/PhD INF CNCR VU/DATA/transcriptomics/RNAseq/VISp_mouse/mouse_VISp_2018-06-14_intron-matrix.csv", chunksize=2000)
filename = "/Users/annagalakhova/PhD INF CNCR VU/DATA/transcriptomics/RNAseq/VISp_mouse/mouse_VISp_2018-06-14_intron-matrix_CPMlog10.csv"
print('iterator made')
c=0
header = True
for chunk in df:
    print(c)
    libsize = np.sum(chunk.iloc[:,1:].values, axis=1)
    logcpm = np.log10(((chunk.iloc[:,1:].T / libsize).T *1_000_000) + 1)
    logcpm = logcpm.astype('float32') #reduce memory from float64
    logcpm.insert(0, 'sample_name', chunk.iloc[:,0])
    if c == 0:
        logcpm.to_csv(filename, index = False)
        c+=1
    else:
        logcpm.to_csv(filename, mode='a', header=False, index = False)
        c+=1

