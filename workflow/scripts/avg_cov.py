import pandas as pd

data = pd.read_table(snakemake.input[0], header=None)
avg_cov = round(data.iloc[:, 2] / (data.iloc[:, 1] / 1000), 3) #truncated to 3 decimals
data.loc[:,4] =  avg_cov
data.to_csv(snakemake.output[0], sep="\t", header=None)