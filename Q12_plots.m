
close all
clear all
figure;
line_colors = {'b', 'g', 'r', 'c'};
load("Stochastic_Simulation_data.mat")
for sim_num = 1:4
    % Construct the filename
    filename = sprintf('tau_leap_j%d.mat', sim_num);
    tau_list = [1,0.5,0.25,0.125];
    % Load data from file
    load(filename);  % Assuming the file contains variables 'time' and 'data'
    hold on;
    error = abs(mean_A_tau_leap - mean(mean_A));
     plot(time_vec_tau_leap, error, 'k-', 'LineWidth', 2,'Color',line_colors{sim_num}, 'DisplayName', sprintf('tau mean absolute error = %f', tau_list(sim_num)));
    xlabel('Time');
    ylabel('Data Values');
    title('Absolute Mean Error comparison for different values of tau');
    legend('show');
end

hold off;