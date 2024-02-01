import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

results = pd.read_csv("../results/results.csv")
results = results.dropna()


# noinspection SqlNoDataSourceInspection,SqlDialectInspection
queried_data = sqldf('''
SELECT p_value, MIN(cond_number) as MinTime, sweeps
FROM results
WHERE mesh = 'mesh-cube-20' AND preconditioner = 'amg'  AND p_value < 8 
GROUP BY p_value, sweeps
ORDER BY p_value
'''
                     )

print(queried_data)

plt.figure()
for sweep in queried_data['sweeps'].unique():
    subset = queried_data[queried_data["sweeps"] == sweep]
    plt.plot(subset['p_value'], subset['MinTime'], label=f'sweep = {sweep}')

plt.xlabel('p')
plt.ylabel('condition number')
plt.yscale('log')
plt.title('AMG p_value vs condition number with different sweep value')
plt.grid()
plt.legend()
plt.savefig('../results/plot/AMGsweepCond7')
plt.close()

######################Plotting for P=8####################################

# noinspection SqlNoDataSourceInspection,SqlDialectInspection
queried_data = sqldf('''
SELECT p_value, MIN(cond_number) as MinTime, sweeps
FROM results
WHERE mesh = 'mesh-cube-20' AND preconditioner = 'amg' 
GROUP BY p_value, sweeps
ORDER BY p_value
'''
                     )


plt.figure()
for sweep in queried_data['sweeps'].unique():
    subset = queried_data[queried_data["sweeps"] == sweep]
    plt.plot(subset['p_value'], subset['MinTime'], label=f'sweep = {sweep}')

plt.xlabel('p')
plt.ylabel('condition number')
plt.yscale('log')
plt.title('AMG p_value vs condition number with different sweep value')
plt.grid()
plt.legend()
plt.savefig('../results/plot/AMGsweepCond8')
plt.close()
