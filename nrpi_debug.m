% Define the function for the Newton Raphson Current Injection Method for distribution network analysis.
function [v, iteration] = nrpi_debug(load_data, line_data, slack_bus_voltage, tolerance, max_iter)
    % Arguments include the constant slack bus voltage, convergence value, and max. iterations,
    % all of which can be changed when the function is called.

    prompt = {'Ensure the argument data i.e, load and line data are in this format:'
              'load data: column1 - bus index, column2 - real power (P), and column3 - reactive power (Q)'
              'line data: column1 - sending bus index, column2 - receiving bus index, column3 - resistance (R), and reactance (X)'
              '\n'
              'The default values are slack_bus_voltage = 1.5, tolerance = 1e-6, and max_iter = 100'
              'If you decide to change any of the constant argument values, such as slack_bus_voltage, tolerance, or max_iter, just pass the new value in their respective positions.'
              'If you skip any argument, maintain the default values in the skipped position.'
              'For example: [v, iteration] = nr_power_injection_method(load_data, line_data, 1.5, 1e-4)'
              'I did not put the value of max_iter because it is the last one and I want to keep it the same.'
              };
    disp(prompt);
    user_input = input('Press Enter to continue or type ''quit'' to exit: ', 's');

    % Check if the user wants to quit
    if strcmp(user_input, 'quit')
        disp('Exiting the function as per user request.');
        iteration = [];
        v = [];
        return;
    end

    % Set default values if arguments are not provided
    if nargin < 3, slack_bus_voltage = 1; end
    if nargin < 4, tolerance = 1e-6; end
    if nargin < 5, max_iter = 100; end
    
    % System base values
    V_base = 12.66; % Nominal voltage in kV
    S_base = 1; % Base power in MVA
    Z_base = (V_base^2) / S_base; % Base impedance in ohms

    num_buses = size(load_data, 1); % Number of buses
    P_load = load_data(:, 2) / 1000; % Real power in p.u.
    Q_load = load_data(:, 3) / 1000; % Reactive power in p.u.
      
    % Preallocate arrays to improve performance
    max_voltage_errors = zeros(1, max_iter);
    cumulative_iter_times = zeros(1, max_iter);
   
    % Initialize bus voltages
    v = ones(num_buses, 1) + 1j * zeros(num_buses, 1); % Initial guess
    v(1) = slack_bus_voltage; % Slack bus voltage

    % Construct the Y_bus matrix using the helper function script
    Y_bus = build_Y_bus(line_data, num_buses, Z_base);
        
    % Start NR Iterations
    for iteration = 1:max_iter % Iterate up to max_iter times to perform the BFS algorithm.
            % Start timing for the iteration
            iter_start_time = tic;
            
            % Initialize mismatch vectors
            P_mismatch = zeros(num_buses, 1);
            Q_mismatch = zeros(num_buses, 1);

            % Calculate the current injection and power mismatches
            for i = 2:num_buses % Skip the slack bus
                V_i = abs(v(i));
                delta_i = angle(v(i));

                P_cal = 0;
                Q_cal = 0;
                for j = 1:num_buses
                    V_j = abs(v(j));
                    delta_j = angle(v(j));
                    Y_ij = abs(Y_bus(i, j));
                    theta_ij = angle(Y_bus(i, j));

                    P_cal = P_cal + V_i * V_j * Y_ij * cos(delta_i - delta_j - theta_ij);
                    Q_cal = Q_cal + V_i * V_j * Y_ij * sin(delta_i - delta_j - theta_ij);
                    
                end
                
                % Power mismatch
                P_mismatch(i) = - P_load(i) - P_cal;
                Q_mismatch(i) = - Q_load(i) - Q_cal;
                             
            end

            % Build the Jacobian matrix using the help function script
            Jacobian = build_Jacobian(v, Y_bus, num_buses);
            % Jacobian_inv = inv(Jacobian);
            
            % J1 = Jacobian_inv(1:num_buses-1, 1:num_buses-1);
            % J2 = Jacobian_inv(1:num_buses-1, num_buses:end);
            % J3 = Jacobian_inv(num_buses:end, 1:num_buses-1);
            % J4 = Jacobian_inv(num_buses:end, num_buses:end);
            
            % Solve for corrections
            mismatch = [P_mismatch(2:end); Q_mismatch(2:end)];
            corrections = Jacobian \ mismatch;
            
            % delta_mismatch = J1*P_mismatch(2:end) + J2*Q_mismatch(2:end);
            % V_mismatch = J3*P_mismatch(2:end) + J4*Q_mismatch(2:end);

            % Update the voltages
            delta_mismatch = corrections(1:num_buses-1);
            V_mismatch = corrections(num_buses:end);

            V_new = abs(v(2:end)) + V_mismatch;
            delta_new = angle(v(2:end)) + delta_mismatch;

            v(2:end) = V_new .* exp(1j * delta_new); % Polar to rectangular

            % End timing the iteration
            iter_time = toc(iter_start_time);
            
            % Check for max_difference
            max_diff = max(abs(V_mismatch));
            
            % Store max voltage difference and cumulative iteration time
            max_voltage_errors(iteration) = max_diff;  % Assign value directly
        
            if iteration == 1
                cumulative_iter_times(iteration) = iter_time;
            else
                cumulative_iter_times(iteration) = cumulative_iter_times(iteration - 1) + iter_time;
            end

            % Print max voltage difference at an iteration and its respective
            % computation time till that iteration.
            fprintf('Iteration %d: max voltage difference = %.10f\n', iteration, max_diff);
            fprintf('Time taken till %d iteration: %.4f seconds\n', iteration, iter_time);
            
            if max(abs(V_mismatch)) <= tolerance % Check if the maximum voltage difference between iterations is less than the tolerance.
                break; % Break the loop if convergence is achieved.
            end
     end
        
        
      % Prepare results and print them
    fprintf('Converged in %d iterations.\n', iteration); % Print the number of iterations it took to converge.

    disp("Jacobian");
    disp(Jacobian);
    disp("inv(Jacobian)");
    disp(inv(Jacobian));
    disp("mismatch");
    disp(mismatch);
    disp("corrections");
    disp(corrections);
    disp(real(v));