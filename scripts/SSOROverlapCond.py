import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

results = pd.read_csv("../results/results.csv")
results = results.dropna()


# noinspection SqlNoDataSourceInspection,SqlDialectInspection
queried_data = sqldf('''
SELECT p_value, MIN(cond_number) as CN, overlap
FROM results
WHERE mesh = 'mesh-cube-20' AND preconditioner = 'ssor'  
GROUP BY p_value, overlap
ORDER BY p_value
'''
                     )



plt.figure()
for overlap in queried_data['overlap'].unique():
    subset = queried_data[queried_data["overlap"] == overlap]
    plt.plot(subset['p_value'], subset['CN'], label=f'overlap = {overlap}')

plt.xlabel('p')
plt.ylabel('condition number')
plt.yscale('log')
plt.title('SSOR p_value vs condition number with different overlap values')
plt.grid()
plt.legend()
plt.savefig('../results/plot/SSOROverlapCond')
plt.close()
