
import sys

import matplotlib.pyplot as plt
import numpy as np

time_step_start = 0
time_step_end = 0

# Read the optional two integer arguments
# These are the number of time steps to cut off the beginning and end of the
# simulation
if len(sys.argv) == 3:
    # Set the time steps
    time_step_start = int(sys.argv[1])
    time_step_end = int(sys.argv[2])

# Open the file at ../results/lift_drag.csv
with open('../results/lift_drag.csv', 'r') as file:
    # Read the file
    data = file.read()

# Split the data into lines
lines = data.split('\n')

# Check the header to see if weak forces are included, then remove it
header = lines.pop(0)
if "weak" in header:
    plot_twice = True
else:
    plot_twice = False

# Initialize the lists
time_step = []
lift_coefficient = []
drag_coefficient = []
weak_lift_coefficient = []
weak_drag_coefficient = []
reynolds_number = []

# Iterate over the lines
for line in lines[time_step_start:-1-time_step_end]:
    # Split the line into values
    values = line.split(',')

    # Append the values to the lists
    time_step.append(float(values[0]))
    lift_coefficient.append(float(values[1]))
    drag_coefficient.append(float(values[2]))
    if plot_twice:
        weak_lift_coefficient.append(float(values[3]))
        weak_drag_coefficient.append(float(values[4]))
        reynolds_number.append(float(values[5]))
    else:
        reynolds_number.append(float(values[3]))

# Convert the lists to numpy arrays
time_step = np.array(time_step)
lift_coefficient = np.array(lift_coefficient)
drag_coefficient = np.array(drag_coefficient)
weak_lift_coefficient = np.array(weak_lift_coefficient)
weak_drag_coefficient = np.array(weak_drag_coefficient)
reynolds_number = np.array(reynolds_number)

# Plot the lift coefficient
plt.figure(1)
plt.plot(time_step, lift_coefficient, 'k-')
if plot_twice:
    plt.plot(time_step, weak_lift_coefficient, 'r-')
plt.xlabel('Time step')
plt.ylabel('Lift coefficient')
plt.grid(True)

# Plot the drag coefficient
plt.figure(2)
plt.plot(time_step, drag_coefficient, 'k-')
if plot_twice:
    plt.plot(time_step, weak_drag_coefficient, 'r-')
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
