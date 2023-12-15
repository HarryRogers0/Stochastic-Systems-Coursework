% Code modified for use from Stochastic methods course.
% Original Author Kit Yates.

clear all
close all

%Define the rate constant of degradation
k1mu=1*10^(-1); % k1/v
kminus1=1; % k_-1
k2=0.1; %k_2
kminus2mu=1; %(sec^{-1}


% Define the number of repeats we will do
M=1000;

% Define the final time we will simulate to
T_final=100;

%Define the recording time interval
rec_step=0.1;

%Define a vector of the time_points
time_vec=[0:rec_step:T_final];

%Calculate the number of recording steps required
num_rec_steps=T_final/rec_step;

%Define the vector which will reocrd the number of particles  at each time
%point
rec_vector_A=zeros(num_rec_steps+1,M);

%Define the initial number of particles
A_init=1;

%Define the initial number of molecules forthe recording vector
rec_vector_A(1,:)=A_init*ones(1,M);

%Define a vector which will hold the propensity functions
a=zeros(4,1);

%Run through a for loop for each of the repeats
for i=1:M
    
    %Define the initial time to be zero
    t=0;
    
    %Initialise the times to help recording
    t_after=0;
    
    %   initialise the number of particles for this repeat
    A=A_init;
    
    %Run through a for loop for each of the time-steps
    %Enter the while loop
    while t<T_final
        
        %Define the propensity functions
        a(1)=A*(A-1)*k1mu;
        a(2)=A*kminus1;
        a(3)=A*k2;
        a(4)=kminus2mu;
        
        %Calculate the cumulative sum of a
        cumsuma=cumsum(a);
        
        %Calculate the sum of the propensity functions
        a0=cumsuma(end);
        
        %Determine the time for the next reaction
        tau=(1/a0)*log(1/rand);
        
        %Update the time
        t=t+tau;
        
        %Write the current value of A to A_old
        A_old=A;
        
        %multiply the random number by the sum of the propensities and
        %store so we can reuse
        ra0=rand*a0;
        
        %if we have not run out of molecules
        if tau<inf;
            %Determine which reaction happened
            if ra0<cumsuma(1)
                %Implement the production reaction
                A=A_old-1;
            elseif ra0<cumsuma(2)
                %Implement the degradation reaction
                A=A_old+1;
            elseif ra0<cumsuma(3)
                % Implement a production reaction
                A=A_old-1;
            else
                % Implement the other production reaction
                A=A_old+1;
            end
        end
        %Calculate the times for recording
        t_before=t_after;
        t_after=t;
        
        %Calculate the indices of the time step before and the current time
        %step in terms of recording
        ind_before=ceil((t_before+eps)/rec_step);
        ind_after=min(floor(t_after/rec_step),num_rec_steps);
        
        %Find out how many time-steps to write to
        steps_to_write=ind_after-ind_before+1;
        
        if steps_to_write>0 && steps_to_write~=Inf
            rec_vector_A(ind_before+1:ind_after+1,i)=(A_old)*ones(steps_to_write,1);
        end
        
    end
    
end

%Solve the ODEs for the deterministic system
%new = rhs(1,1,k1mu,kminus1,k2,kminus2mu);
[t_deterministic,X_deterministic]=ode15s(@(t,X)(k1mu - k2 + kminus1)*X(1)-k1mu*X(1)^2 + kminus2mu,[0,T_final],[1]);

%Write A and B from X.
A_deterministic=X_deterministic(:,1)';

%Calculate the mean over all 1000 realisations
mean_A=mean(rec_vector_A,2)';
%Calculate the variance over all 1000 realisations
var_A = var(rec_vector_A')';

%Plot the first 30 trajectories for A and B
f1=figure;
f2 = figure;
f3 = figure;

%     Change to figure f1
figure(f1)
plot(time_vec,rec_vector_A(:,1:30));
figure(f2)
plot(time_vec,mean_A);
figure(f3)
plot(time_vec,var_A);

save('Stochastic_Simulation_data','rec_vector_A','mean_A','time_vec')
%     Change to figure f1
figure(f1)
hold on
%Also plot on the mean over 1000 realisations
[hA1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the analytically determined mean
[hA2]=plot(t_deterministic,A_deterministic,'k--','linewidth',5); %Comment line for Q4

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of A particles')

legend([hA1,hA2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs'); % Comment out for Q4
%legend([hA1],['Mean of M=',num2str(M),' simulations']); %uncomment for Q4

hold off

%     Change to figure f1
figure(f2)
hold on
% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('mean number of A particles')

%     Change to figure f2
figure(f2)

%Save as a .fig as well


hold off


% Calculating the stationary distribution histogram (Question 11)


%Define the time after which we consider the process to be stationary
t_stat=10;

%Disregard the first t_stat seconds worth of data
rec_vector_A_cut=rec_vector_A(t_stat/rec_step:end,:);

%Set the maximum value of X
max_data_A=max(rec_vector_A_cut(:));

%Decide the range into which the data will fall
X_A=[0:1:max_data_A];

%extract histogram data from each of the repeats after 40 seconds when the
%process should have roughly reached steady state
[N,X_A]=hist(rec_vector_A_cut(:),X_A);

%Normalised the data to make it a PDF
N=N/sum(N);

% Plot the histogram of the stationary frequencies
[h1]=bar(X_A,N);

%correct the axes
axis([0 max_data_A+0.5 0 ceil(max(N)*100)/100])


%Calculate the simulated mean and variance
mu=mean(rec_vector_A_cut(:));
sim_var=var(rec_vector_A_cut(:));
sigma = sqrt(sim_var);

% Calculate the PDF normal distribution using simulated mean and S.D
pdf_values = exp(-(X_A - mu).^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi));

hold on
%Also plot the analytically determined mean
[h2]=plot(X_A,pdf_values,'r','linewidth',5);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('number of particles')
ylabel('stationary distribution')

legend([h1,h2],['Long-run SSA'],'PGF analysis');



       