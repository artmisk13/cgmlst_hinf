import pandas as pd
import numpy as np
import statsmodels.api as sm

#input file (infile): the dataframe from correlation-test.R script (i.e. compare_data2) as csv
#outfile: user-defined output file (.txt)

infile = sys.argv[1]
outfile = sys.argv[2]

df1 = pd.read_csv(infile)

x = np.array(df1['df_allele.value']).reshape((-1, 1))
y = np.array(df1['df_nt.log10value'])

x1 = sm.add_constant(x) 

#regression model based on ordinary least squares (OLS): statsmodels.regression.linear_model.OLS
model1 = sm.OLS(y, x1)
results = model1.fit()

with open(outfile, 'w') as f:
  f.write(results.summary())
