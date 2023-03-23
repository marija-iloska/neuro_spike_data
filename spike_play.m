clc 
clear all
close all


% Load spike data
load spike_data_new_new.mat

% Total number of neurons
num_neurons = length(neuron_spk_times);
fsz = 20;

%-------------------------------------------
%              VISUALIZE DATA
%-------------------------------------------
figure(1)
for j = 1:num_neurons

    % Find last time of spike for each neuron
    if (isempty(neuron_spk_times{j}) == 0)
        max_t(j) = max(neuron_spk_times{j});
    end
    
    % For plotting
    l = ones(1,length(neuron_spk_times{j}));
    plot(neuron_spk_times{j}, j*l, '.', 'MarkerSize', 12)
    hold on
end
xlabel('Time','FontSize',fsz)
ylabel('Neuron', 'FontSize', fsz)
xline(0, 'k')
set(gca, 'FontSize', 15)


%-------------------------------------------
%              STATE SPACE MODEL
%-------------------------------------------
dt = 0.05;
time = -0.2 : dt : max(max_t);
T = length(time);
num_groups = 3;
n=1;
i = 0.2;
step = i;
count = 1;

% For each time step
for t = 1:T

    % For each neuron, get number of spikes within time bin
    for j = 1 : num_neurons

        if (t == 1)
            y(j, t) = sum(neuron_spk_times{j} < time(t));
        else
            y(j, t) = sum( (neuron_spk_times{j} < time(t) ))  - sum(y(j,1:t-1));
        end
        y_tot(j,t) = sum(y(j,1:t));

    end

    % Find max rate
    if (time(t) > i)
        
        
        final_x = sum(y(:,n:t),2);
        [idx, clust] = group(num_groups, final_x, time(t));
        for k = 1:num_groups
            y_avg(k,n:t) = mean(y(clust{k},n:t),1);
        end
        clust_store{count} = clust;
        i = i + step;
        n = t;
        count = count+1;
    end
end

final_x = y(:,T);
[idx, clust] = group(num_groups, final_x, time(t));
for k = 1:num_groups
    y_avg(k,n:T) = mean(y(clust{k}, n:T),1);
end

%-------------------------------------------
%              NUMBER OF SPIKES
%-------------------------------------------

%-------------------------------------------
%               FIRING RATES
%-------------------------------------------

% Plot a few neurons
plot_neurons = 58;
lwd = 1.5;


figure(2)
for j = 1:plot_neurons
    %j = datasample(1:num_neurons, 1);
    plot(time, y(j,:), 'Linewidth', lwd)
    hold on
end
xline(0, 'k')
set(gca, 'FontSize', 15)
xlabel('Time','FontSize',fsz)
ylabel('Number of times the neuron spiked', 'FontSize', fsz)
title('dt = ',num2str(dt), 'FontSize', 20)




%-------------------------------------------
%               GROUPS
%-------------------------------------------



figure(3)
% Cumulative spikes
for j = 1:plot_neurons
    plot(time(1:t-1), y_tot(j,1:t-1), 'Linewidth', lwd-0.5)
    hold on
end
xline(0, 'k')
xlabel('Time','FontSize',fsz)
ylabel('Cumulative number of spikes', 'FontSize', fsz)
title('dt = ',num2str(dt), 'FontSize', 20)
set(gca, 'FontSize', 15)


% Average firing rate of groups
figure(4)
for k = 1:num_groups
    %kk = datasample(1:num_groups, 2);
    plot(time, y_avg(k,:), 'Linewidth', lwd)
    hold on
end
xline(0, 'k')
xlabel('Time','FontSize',fsz)
ylabel('Average firing rate of Groups', 'FontSize', fsz)
title('dt = ',num2str(dt), 'FontSize', 20)
set(gca, 'FontSize', 15)


%save('firing_rate_data.mat', 'y', 'y_avg', 'time', 'T', 'clust_store', 'dt', 'num_neurons', 'num_groups')


