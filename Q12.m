% Code modified for use from Stochastic methods course.
% Original Author Kit Yates.

clear all
close all

k1mu=1*10^(-1); % k1/v
kminus1=1; % k_-1
k2=0.1; %k_2
kminus2mu=1; %(sec^{-1}


% Define the final time we will simulate to
T_final=128;

%Define the number of different taus to test
num_taus=4;

%What value should we add to increase the lower index from 1
index_increment= 6;

%index the power (of 2) of the number of steps to take
tau_power_indices=1+index_increment:num_taus+index_increment;

%Define the list of taus to use
tau_list=T_final*(1/2).^tau_power_indices;

%Define a vector to hold the time for each tau
tau_time=zeros(num_taus,1);

%Define a vector to hold the number of paths for each tau
tau_paths=zeros(num_taus,1);

%Define a vector which will hold the number of molecules of A at the end of
%the simulation
time_final=zeros(num_taus,1);

%Determine th enumber of reactions
num_reactions=4;

%create the vector which tels us the number of firings of each reaction
num_firings=zeros(num_reactions,1);
    
%Enter the for loop which allows us to run the tau leaping simulation for
%each different value of tau
for j=tau_power_indices
    j
    %Extract the releavent value of tau from the list
    tau=tau_list(j-index_increment);
    
    %Define a vector of the time_points
    time_vec_tau_leap=0:tau:T_final;
    
    %Calculate the number of recording steps required
    num_rec_steps=2^j;
    
    %Define the vector which will reocrd the number of particlesat each time
    %point
    rec_vector_A_tau_leap=zeros(num_rec_steps+1,1);
    
    %Define the initial number of particles
    A_init=1;
    
    %Define the tolerance (how small the variance has to get)
    var_tol=0.01;
    
    %Calculate sample variance
    sample_variance=var(rec_vector_A_tau_leap)/1;
    
    %Initialise the indexing variable i
    i=0;
    
    %(re-)start the timer
    tic
    

    %Run through a for loop for each of the repeats until the sample mean
    %is less than out tolerance (providing we have a t leat two readings)
    while  (sample_variance)>var_tol || i<2
        
        i=i+1;
        
        if i>1
            %Add another row to rec_vec
            rec_vector_A_tau_leap=[rec_vector_A_tau_leap,zeros(num_rec_steps+1,1)];
        end
        
        %Define the initial number of molecules for the recording vector
        rec_vector_A_tau_leap(1,i)=A_init;
        
        %Define the initial time to be zero
        t=0;
        
        %   initialise the number of particles for this repeat
        A=A_init;
        
        %Run through a for loop for each of the time-steps
        %Enter the while loop
        for k=1:2^j
            
            %Generate the number of firings of each reaction
            num_firings(1)=poissrnd(A*(A-1)*k1mu*tau);
            num_firings(2)=poissrnd(A*kminus1*tau);
            num_firings(3)=poissrnd(A*k2*tau);
            num_firings(4)=poissrnd(kminus2mu*tau);
                       
            %Update the species numbers by adding or sutracting poisson distributed
            %random variates with the correct mean.
            A=A-num_firings(1)-num_firings(3)+num_firings(2)+num_firings(4);

            % Deal with edge case where A may be less than zero.
            if A < 0
                A = 0;
            end
            %Record the current value of the species numbers
            rec_vector_A_tau_leap(k+1,i)=A;
            
        end
        
        %Calculate the sample variance of A
        sample_variance=var(rec_vector_A_tau_leap(end,:))/i;
        
    end
    %output the time taken
    tau_time(j-index_increment)=toc;
    %output the number of paths used
    tau_paths(j-index_increment)=i;
    
    %Calculate the mean over all realisations
    mean_A_tau_leap=mean(rec_vector_A_tau_leap,2);
    
    %Extract the final value of A
    mean_A_final(j-index_increment)=mean_A_tau_leap(end);
    %Save the relevant data for this value of tau
    save(['tau_leap_j',num2str(j-index_increment)],'mean_A_tau_leap','rec_vector_A_tau_leap','num_taus','tau','index_increment','time_vec_tau_leap')
end

%Number of paths are as follows:
disp(['tau=',num2str(tau_list(1)),' tau=',num2str(tau_list(2)),' tau=',num2str(tau_list(3)),' tau=',num2str(tau_list(4))])
disp(['paths=',num2str(tau_paths(1)),' paths=',num2str(tau_paths(2)),' paths=',num2str(tau_paths(3)),' paths=',num2str(tau_paths(4))])

%Final numbers are as follows
disp(['tau=',num2str(tau_list(1)),' tau=',num2str(tau_list(2)),' tau=',num2str(tau_list(3)),' tau=',num2str(tau_list(4))])
disp(['mean_A=',num2str(mean_A_final(1)),' mean_A=',num2str(mean_A_final(2)),' mean_A=',num2str(mean_A_final(3)),' mean_A=',num2str(mean_A_final(4))])