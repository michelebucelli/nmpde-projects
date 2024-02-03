import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

# Read the results from the CSV file
results = pd.read_csv("../results/results.csv")

# Use pandasql to group by preconditioner and r_value, finding the minimum number of iterations
queried_data = sqldf('''
    SELECT mesh, r, preconditioner, MIN(iterations) as MinIterations
    FROM results
    GROUP BY preconditioner, r
    ORDER BY preconditioner
''')

# Plotting
plt.figure()
for preconditioner in queried_data['preconditioner'].unique():
    subset = queried_data[queried_data["preconditioner"] == preconditioner]
    plt.plot(subset['r'], subset['MinIterations'], label=preconditioner)

plt.xlabel('r')
plt.ylabel('Number of Iterations')
plt.yscale('log')  # Adjust the scale if needed
plt.title('Number of Iterations vs r Value for Different Preconditioners')
plt.grid()
plt.legend()
plt.savefig('../results/plot/iterations_vs_r')
plt.close()
