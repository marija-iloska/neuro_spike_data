clear all
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

dt = 0.05;
time = -0.1 : dt : max(max_t);
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

        % Get number of spikes of neuron per time step
        if (t == 1)
            y(j, t) = sum(spikes_all_trials{j} < time(t));
        else
            y(j, t) = sum( (spikes_all_trials{j} < time(t) ))  - sum(y(j,1:t-1));
        end

        % Get cumulative
        y_tot(j,t) = sum(y(j,1:t));

    end

    % Clustering
    if (time(t) > i)
        % Get neuron activity within last step interval        
        final_x = sum(y(:,n:t),2);

        % Get clusterings
        [idx, clust] = group(num_groups, final_x, time(t));

        % Take the average firing rates in each clusters
        for k = 1:num_groups       
            y_avg(k,n:t) = mean(y(clust{k},n:t),1);
        end

        % Store indices
        clust_store{count} = clust;

        % Hop a step
        i = i + step;

        % Save step
        n = t;

        % Count how many steps
        count = count+1;
    end
end

% Repeat for final group
final_x = y(:,T);
[idx, clust] = group(num_groups, final_x, time(t));
for k = 1:num_groups
    y_avg(k,n:T) = mean(y(clust{k}, n:T),1);
end


% Plot a few neurons
plot_neurons = 58;
lwd = 1.5;
fsz = 20;
idx = 1:40;

figure(2)
for k = 1:num_groups
    %j = datasample(1:num_neurons, 1);
    plot(time(idx), y_avg(k,idx), 'Linewidth', lwd)
    hold on
end
xline(0, 'k')
set(gca, 'FontSize', 15)
xlabel('Time','FontSize',fsz)
ylabel('Number of times the neuron spiked', 'FontSize', fsz)
title('dt = ',num2str(dt), 'FontSize', 20)



%save('all_dt005.mat')
