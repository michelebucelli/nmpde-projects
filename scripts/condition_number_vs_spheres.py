import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

# Read the results from the CSV file
results = pd.read_csv("../results/results.csv")

# Use pandasql to group by preconditioner and r_value, finding the minimum number of iterations
queried_data = sqldf('''
    SELECT mesh, spheres, preconditioner, MIN(cond_number) as CN
    FROM results
    WHERE mesh = 'mesh-cube-40'
    GROUP BY preconditioner, spheres
    ORDER BY preconditioner
'''
                     )

# Plotting
plt.figure()
for preconditioner in queried_data['preconditioner'].unique():
    subset = queried_data[queried_data["preconditioner"] == preconditioner]
    plt.plot(subset['spheres']**3, subset['CN'], label=preconditioner)

plt.xlabel('Number of spheres')
plt.ylabel('Condition number')
plt.yscale('log')  # Adjust the scale if needed
plt.title('Number of sphere vs condition number')
plt.grid()
plt.legend()
plt.savefig('../results/plot/cond_vs_spheres')
plt.close()
