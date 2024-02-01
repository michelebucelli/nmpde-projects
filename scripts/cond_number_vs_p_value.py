import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd


results = pd.read_csv("../results/results.csv")
results = results.dropna()


queried_data = sqldf('''
    SELECT mesh, p_value, preconditioner, MIN(cond_number) as CN
    FROM results
    GROUP BY preconditioner, p_value
    ORDER BY preconditioner
''')



plt.figure()
for preconditioner in queried_data['preconditioner'].unique():
    subset = queried_data[queried_data["preconditioner"] == preconditioner]
    plt.plot(subset['p_value'], subset['CN'], label=preconditioner)

plt.xlabel('p')
plt.ylabel('Condition number')
plt.yscale('log')
plt.title('Condition number vs p')
plt.grid()
plt.legend()
plt.yscale('log')
plt.savefig('../results/plot/cond_numb_vs_p')
plt.close()


