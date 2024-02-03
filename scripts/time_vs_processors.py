import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

results = pd.read_csv("../results/results.csv")
results = results.dropna()

# noinspection SqlNoDataSourceInspection,SqlDialectInspection
queried_data = sqldf('''
SELECT mesh, n, preconditioner, MIN(elapsed_time) as MinTime
FROM results
GROUP BY preconditioner, n
ORDER BY preconditioner
'''
                     )

plt.figure()
for preconditioner in queried_data['preconditioner'].unique():
    subset = queried_data[queried_data["preconditioner"] == preconditioner]
    plt.plot(subset['n'], subset['MinTime'], label=preconditioner)

plt.xlabel('n')
plt.ylabel('time [ms]')
plt.yscale('log')
plt.title('Execution time in ms vs number of processors')
plt.grid()
plt.legend()
plt.savefig('../results/plot/time_vs_processors')
plt.close()
