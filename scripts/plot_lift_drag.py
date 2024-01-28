
import matplotlib.pyplot as plt
import numpy as np

# Open the file at ../results/lift_drag.csv
with open('../results/lift_drag.csv', 'r') as file:
    # Read the file
    data = file.read()

# Split the data into lines
lines = data.split('\n')

# Remove the first line (header)
lines.pop(0)

# Initialize the lists
time_step = []
lift_coefficient = []
drag_coefficient = []
reynolds_number = []

# Iterate over the lines
for line in lines[:-1]:
    # Split the line into values
    values = line.split(',')

    # Append the values to the lists
    time_step.append(float(values[0]))
    lift_coefficient.append(float(values[1]))
    drag_coefficient.append(float(values[2]))
    reynolds_number.append(float(values[3]))

# Convert the lists to numpy arrays
time_step = np.array(time_step)
lift_coefficient = np.array(lift_coefficient)
drag_coefficient = np.array(drag_coefficient)
reynolds_number = np.array(reynolds_number)

# Plot the lift coefficient
plt.figure(1)
plt.plot(time_step, lift_coefficient, 'k-')
plt.xlabel('Time step')
plt.ylabel('Lift coefficient')
plt.grid(True)

# Plot the drag coefficient
plt.figure(2)
plt.plot(time_step, drag_coefficient, 'k-')
plt.xlabel('Time step')
plt.ylabel('Drag coefficient')
plt.grid(True)

# Plot the Reynolds number
plt.figure(3)
plt.plot(time_step, reynolds_number, 'k-')
plt.xlabel('Time step')
plt.ylabel('Reynolds number')
plt.grid(True)

# Show the plots
plt.show()
