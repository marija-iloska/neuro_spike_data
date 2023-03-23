%clear all
close all
clc

load spike_all.mat

num_neurons = length(spikes_all_trials);

%-------------------------------------------
%              STATE SPACE MODEL
%-------------------------------------------

for j = 1:num_neurons
    max_temp(j) = max(spikes_all_trials{j});
end
max_t = max(max_temp);

dt = 0.1;
time = -0.1 : dt : max(max_t);
T = length(time);

% For each time step
for t = 1:T

    % For each neuron, get number of spikes within time bin
    for j = 1 : num_neurons

        if (t == 1)
            x(j, t) = sum(spikes_all_trials{j} < time(t));
        else
            x(j, t) = sum( (spikes_all_trials{j} < time(t) ))  - sum(x(j,1:t-1));
        end
        x_tot(j,t) = sum(x(j,1:t));
    end
end


% Plot a few neurons
plot_neurons = 58;
lwd = 1.5;
fsz = 20;
idx = 1:150;

figure(2)
for j = 1:plot_neurons
    %j = datasample(1:num_neurons, 1);
    plot(time(idx), x_tot(j,idx), 'Linewidth', lwd)
    hold on
end
xline(0, 'k')
set(gca, 'FontSize', 15)
xlabel('Time','FontSize',fsz)
ylabel('Number of times the neuron spiked', 'FontSize', fsz)
title('dt = ',num2str(dt), 'FontSize', 20)


