function main()
    % Define parameters
    k1nu=0.1; % k1/v 
    k_minus1=1; % k_-1
    k2=0.1; %k_2
    k_minus2nu=1; %(sec^{-1})

    % Set initial conditions
    V0 = 0;
    M0 = 1;

    % Set time span for the simulation
    tspan = [0 100]; % Adjust the time span as needed

    % Solve the system of differential equations using ode15s
    [t, Y] = ode15s(@(t,X)RHS_Norm(t, X, k1nu, k2, k_minus1, k_minus2nu),tspan,[V0,M0]);

    % Extract the solutions
    secondMoment = Y(:, 1);
    M = Y(:, 2);

    % Plot the results
    figure;
    subplot(2, 1, 1);
    plot(t, secondMoment, 'b', 'LineWidth', 2);
    xlabel('Time');
    ylabel('<n^2>');
    title('Solution of d<n^2>/dt');


    subplot(2, 1, 2);
    plot(t, M, 'r', 'LineWidth', 2);
    xlabel('Time');
    ylabel('M(t)');
    title('Solution of dM/dt');

end

function dydt = RHS_Norm(t, X, k1mu, k2, k_minus1, k_minus2mu)
    % X(1) corresponds to <n^2>, and X(2) corresponds to M
    dydt(1) = (1)*((-2*k1mu)*(3.*X(1).*X(2)-2.*(X(2))^3) + (3 * k1mu + 2 * k_minus1 - 2 * k2) .* X(1) + (k2 - k1mu + k_minus1 + 2 * k_minus2mu) * X(2) + k_minus2mu);
    dydt(2) = (k1mu - k2 + k_minus1) * X(2) - k1mu * X(1) + k_minus2mu;
dydt=dydt';
end
