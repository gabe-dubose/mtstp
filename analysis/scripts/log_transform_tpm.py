#!/usr/bin/env python3

import pandas as pd
import numpy as np

#load data
tpm_counts = pd.read_csv('/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/dpl_tpm_counts_kallisto.csv', index_col=0)

#add pseudo count
tpm_counts = tpm_counts + 1

log_tpm_counts = np.log(tpm_counts)

log_tpm_counts.to_csv('/home/gabe/Desktop/mtstp/data/intermediate_data/count_tables/dpl_log_tpm_counts_kallisto.csv')