"""
Created on Wed Sep 25 07:17:42 2024

@author: Morufdeen ATILOLA

"""
# Import the required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import time
from build_Y_bus import build_Y_bus
from build_Jacobian  import build_Jacobian
from calculate_system_loss import calculate_system_loss


# Define the function for the backward_forward method for distribution network analysis. 
def nrpi_method(load_data, line_data, slack_bus_voltage=1.05, tolerance=1e-6, max_iter = 100, runs = 10):
    # arguments includes the constant slack bus voltage, convergence value, and max. iterations, all of which can be changed when the function is called
    prompt = (
        "Ensure the argument data i.e., load and line data are in this format:\n"
        "load data: column1 - bus index, column2 - real power (P), and column3 - reactive power (Q)\n"
        "line data: column1 - sending bus index, column2 - receiving bus index, column3 - resistance (R), and reactance (X)\n"
        "\n"
        "The default values are slack_bus_voltage = 1.5, tolerance = 1e-6, and max_iter = 100\n"
        "If you decide to change any of the constant argument values, such as slack_bus_voltage, tolerance, or max_iter, just pass the new value in their respective position\n"
        "If you will jump any of the order, kindly put the constant value of the one you are jumping in its position and input the new value of the argument you are changing\n"
        "For example: [v, iteration] = bfs_method(load_data, line_data, 1.5, 1e-4), this is because, I want to keep the slack bus the same but want to change the tolerance\n"
        "I did not put the value of max_iter because, it is the last one and I want to keep it the same.\n"
        )
    
    print(prompt)
    
    user_input = input("Press Enter to continue or type 'quit' to exit: \n");
    
    if user_input.lower() == 'quit':
        print("Exiting the function as per user request probably because the data does not conform with the requirement")
        iteration = []  # Return an empty array
        v = []
        return v, iteration
    
    # System base values
    V_base = 12.66  # Nominal voltage in kV
    S_base = 1  # Base power in MVA
    Z_base = (V_base ** 2) / S_base  # Base impedance in ohms

    num_buses = load_data.shape[0]  # Number of buses
    P_load = load_data[:, 1] / 1000  # Real power in p.u.
    Q_load = load_data[:, 2] / 1000  # Reactive power in p.u.

    # Preallocate arrays to improve performance
    max_voltage_errors = np.zeros(max_iter)
    cumulative_iter_times = np.zeros(max_iter)
    computation_times = np.zeros(runs)

    # Begin NR method for multiple runs
    for run in range(runs):
        # Initialize bus voltages
        v = np.ones(num_buses, dtype=complex)  # Initial guess
        v[0] = slack_bus_voltage  # Slack bus voltage

        # Construct the Y_bus matrix using the helper function script
        Y_bus = build_Y_bus(line_data, num_buses, Z_base)

        # Measure the computational time
        start_time = time.time()

        # Start NR Iterations
        for iteration in range(max_iter):
            # Start timing for the iteration
            iter_start_time = time.time()

            # Initialize mismatch vectors
            P_mismatch = np.zeros(num_buses)
            Q_mismatch = np.zeros(num_buses)

            # Calculate the current injection and power mismatches
            for i in range(1, num_buses):  # Skip the slack bus
                V_i = abs(v[i])
                delta_i = np.angle(v[i])

                P_cal = 0
                Q_cal = 0
                for j in range(num_buses):
                    V_j = abs(v[j])
                    delta_j = np.angle(v[j])
                    Y_ij = abs(Y_bus[i, j])
                    theta_ij = np.angle(Y_bus[i, j])

                    P_cal += V_i * V_j * Y_ij * np.cos(delta_i - delta_j - theta_ij)
                    Q_cal += V_i * V_j * Y_ij * np.sin(delta_i - delta_j - theta_ij)

                # Power mismatch
                P_mismatch[i] = P_load[i] - P_cal
                Q_mismatch[i] = Q_load[i] - Q_cal

            # Build the Jacobian matrix using the helper function script
            Jacobian = build_Jacobian(v, Y_bus, num_buses)

            # Solve for corrections
            mismatch = np.concatenate([P_mismatch[1:], Q_mismatch[1:]])
            corrections = np.linalg.solve(Jacobian, mismatch)

            # Update the voltages
            delta_mismatch = corrections[:num_buses-1]
            V_mismatch = corrections[num_buses-1:]

            V_new = abs(v[1:]) + V_mismatch
            delta_new = np.angle(v[1:]) + delta_mismatch

            v[1:] = V_new * np.exp(1j * delta_new)  # Polar to rectangular

            # End timing the iteration
            iter_time = time.time() - iter_start_time

            # Check for max_difference
            max_diff = np.max(np.abs(corrections))

            # Store max voltage difference and cumulative iteration time
            max_voltage_errors[iteration] = max_diff

            if iteration == 0:
                cumulative_iter_times[iteration] = iter_time
            else:
                cumulative_iter_times[iteration] = cumulative_iter_times[iteration - 1] + iter_time

            # Print max voltage difference at an iteration and its respective
            # computation time till that iteration.
            print(f'Iteration {iteration + 1}: max voltage difference = {max_diff:.10f}')
            print(f'Time taken till {iteration + 1} iteration: {iter_time:.4f} seconds')

            if max_diff <= tolerance:
                break  # Break the loop if convergence is achieved

        # End computational time measurement
        computation_times[run] = time.time() - start_time

    # Prepare results and print them
    results = [] # Initialize an empty list to store the results.
    print(f"Converged in {iteration+1} iterations.") # Print the number of iterations it took to converge.
    print("Bus Voltages:") # Print a header for the bus voltages.

    for i, voltage in enumerate(v): # Iterate over each bus voltage.
        magnitude = abs(voltage) # Calculate the magnitude of the voltage.
        angle = np.angle(voltage, deg=True)  # Calculate the angle of the voltage in degrees.
        rectangular_form = f"{voltage.real:.4f} + {voltage.imag:.4f}j"
        polar_form = f"{magnitude:.4f} ∠ {angle:.2f}°"
        print(f"Bus {i+1}: Rectangular form: {rectangular_form} p.u., Polar form: {polar_form}") # Print the voltage in both forms.
        results.append((i+1, voltage, polar_form)) # Append the bus number, original voltage, and polar form to the results list.
        
    # Compute system loss, substation power, and other results
    # Calculate system loss
    total_active_loss, total_reactive_loss = calculate_system_loss(num_buses, line_data, v, Z_base)

    print(f'Total active power loss: {total_active_loss:.4f} p.u. = {total_active_loss * 1000:.4f} kW')
    print(f'Total reactive power loss: {total_reactive_loss:.4f} p.u. = {total_reactive_loss * 1000:.4f} kVAR')

    # Calculate substation power
    substation_active_power = np.sum(P_load) + total_active_loss
    substation_reactive_power = np.sum(Q_load) + total_reactive_loss
    
    print(f"Substation active power: {substation_active_power:.4f} p.u. = {substation_active_power:.4f} MW")
    print(f"Substation reactive power: {substation_reactive_power:.4f} p.u. = {substation_reactive_power:.4f} MVAR")
    
    # Find minimum and maximum voltages and their corresponding bus indices
    min_voltage = np.min(np.abs(v))
    min_index = np.argmin(np.abs(v))  # Index of the bus with minimum voltage
    
    max_voltage = np.max(np.abs(v))
    max_index = np.argmax(np.abs(v))  # Index of the bus with maximum voltage
    
    # Print the minimum and maximum voltages along with the bus indices
    print(f'Minimum voltage = {min_voltage:.4f} pu at bus {min_index + 1}')  # Adding 1 for 1-based indexing
    print(f'Maximum voltage = {max_voltage:.4f} pu at bus {max_index + 1}')  # Adding 1 for 1-based indexing
    
    # Create a formatted string representation of the computation times so as to return it as a formated list of strings in 4dp.
    formatted_times = [f"{time:.4f}" for time in computation_times]
    
    # Print the formatted list of computation times
    print(f"Computation times for each run: {formatted_times}")
    
    # Calculate and print the average computation time
    average_time = sum(computation_times) / len(computation_times)
    print(f"Average computation time for the model after {runs} runs: {average_time:.4f} seconds")

    # Plot the voltage magnitudes and angle
    voltage_magnitudes = np.abs(v)
    voltage_angle = np.angle(v)
    bus_indices = np.arange(1, num_buses + 1)

    # Create a figure with 2 subplots
    plt.figure(figsize=(12, 8))
    
    # First subplot for Voltage Magnitude
    plt.subplot(1, 2, 1)  # 1 row, 2 columns, 1st subplot
    plt.plot(bus_indices, voltage_magnitudes, marker='o', linestyle='-', color='b', label='Voltage Magnitude')
    plt.xlabel('Bus Index')
    plt.ylabel('Voltage Magnitude (p.u.)')
    plt.title('Bus Voltage Magnitudes')
    plt.grid(True)
    plt.legend()
    
    # Second subplot for Voltage Angle
    plt.subplot(1, 2, 2)  # 1 row, 2 columns, 2nd subplot
    plt.plot(bus_indices, voltage_angle, marker='o', linestyle='-', color='r', label='Voltage Angle')
    plt.xlabel('Bus Index')
    plt.ylabel('Voltage Angle (radian)')
    plt.title('Bus Voltage Angles')
    plt.grid(True)
    plt.legend()
    
    # Show the combined plots
    plt.tight_layout()
    plt.show()
    
    # Plot the maximum voltage error vs. computation time till iteration
    plt.figure(figsize=(10, 6))
    plt.plot(cumulative_iter_times, max_voltage_errors, marker='o', linestyle='-', color='b')
    plt.xlabel('Computation Time (seconds)')
    plt.ylabel('Maximum Error')
    plt.title('Maximum Error vs. Computation Time')
    plt.grid(True)  # Add grid if desired
    plt.show()
    
    # Write results to a text file with UTF-8 encoding
    with open(f'bus_voltage_{len(load_data)}_bus.txt', 'w', encoding='utf-8') as txt_file:
        txt_file.write(f"Converged in {iteration} iterations.\n")
        txt_file.write("Bus Voltages:\n")
        txt_file.write("Bus Number\tRectangular Form\tPolar Form\n")
        for result in results:
            txt_file.write(f"{result[0]}\t{result[1].real:.4f} + {result[1].imag:.4f}j\t{result[2]}\n")
    
    # Write results to a CSV file with UTF-8 encoding
    with open(f'bus_voltage_{len(load_data)}_bus.csv', 'w', newline='', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["Bus Number", "Rectangular Form", "Polar Form"])
        for result in results:
            csv_writer.writerow([result[0], f"{result[1].real:.4f} + {result[1].imag:.4f}j", result[2]])

    return v, iteration+1


