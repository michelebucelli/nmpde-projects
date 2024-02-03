import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

results = pd.read_csv("../results/results.csv")
results = results.dropna()

# noinspection SqlNoDataSourceInspection,SqlDialectInspection
queried_data = sqldf('''
SELECT p_value, MIN(elapsed_time) as MinTime, sweeps
FROM results
WHERE mesh = 'mesh-cube-20' AND preconditioner = 'amg' AND p_value < 8 
GROUP BY p_value, sweeps
ORDER BY p_value
'''
                     )


plt.figure()
for sweep in queried_data['sweeps'].unique():
    subset = queried_data[queried_data["sweeps"] == sweep]
    plt.plot(subset['p_value'], subset['MinTime'], label=f'sweep = {sweep}')

plt.xlabel('p')
plt.ylabel('time [ms]')
plt.title('Execution time in ms vs sphere exponent for different preconditioners')
plt.grid()
plt.legend()
plt.savefig('../results/plot/AMGsweepTime7')
plt.close()

######################P_value free############################

# noinspection SqlNoDataSourceInspection,SqlDialectInspection
queried_data = sqldf('''
SELECT p_value, MIN(elapsed_time) as MinTime, sweeps
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
plt.ylabel('time [ms]')
plt.title('Execution time in ms vs sphere exponent for different preconditioners')
plt.grid()
plt.legend()
plt.savefig('../results/plot/AMGsweepTime8')
plt.close()