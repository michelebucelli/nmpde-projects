import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd


results = pd.read_csv("../results/results.csv")
results = results.dropna()


queried_data = sqldf('''
    SELECT mesh, p_value, preconditioner, MIN(iterations) as MinIterations
    FROM results
    GROUP BY preconditioner, p_value
    ORDER BY preconditioner
''')


plt.figure()
for preconditioner in queried_data['preconditioner'].unique():
    subset = queried_data[queried_data["preconditioner"] == preconditioner]
    plt.plot(subset['p_value'], subset['MinIterations'], label=preconditioner)

plt.xlabel('p')
plt.ylabel('Number of Iterations')
plt.yscale('log')
plt.title('Number of Iterations vs Sphere Exponent for Different Preconditioners')
plt.grid()
plt.legend()
plt.savefig('../results/plot/iterations_vs_p')
plt.close()
