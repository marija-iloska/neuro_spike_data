clear all
close all
clc

%load firing_rate_data.mat
load all_dt005.mat

% Number of particles
M = 200;
beta = 0.2;

% PF estimation until t0
t0 = T;

% Prediction till t0 + tau
tau = 20;

[x, w, x_particles] = tpf_neurons(y_avg, num_groups, t0, M, beta, dt);

% Start predictions
for t = t0+1 : t0+tau

    for k = 1:num_groups

        % Propose particles
        x_particles(k, :) = exprnd( x_particles(k,:));

        % Get predictive distributions
        x_pred(k,t-t0) = squeeze(w(k,t0, :))'* x_particles(k, :)';
    end

end

k = 3;
for t = 1:t0
    y_test(t) = mean(poissrnd(x(k,t)*dt, 1,1));
end

for t = t0+1:t0+tau
    y_pred(t-t0) = mean(poissrnd(x_pred(k,t-t0)*dt, 1,M));
end


idx = 1:t0;

idxp = t0+1:t0+tau;

% lwd = 1.5;
% fsz = 20;
% figure(1)
% plot(time(idx), y_avg(k,idx), 'k', 'linewidth', lwd)
% hold on
% plot(time(idx), y_test(idx), 'b', 'linewidth', lwd-0.5)
% ylabel('Number of Spikes in dt', 'FontSize', fsz)
% xlabel('Time', 'FontSize', fsz)
% set(gca, 'FontSize', fsz)
% legend('Truth', 'Estimated', 'FontSize', fsz)
% ylim([0, max(max(y_avg))])

lwd = 1.5;
fsz = 20;
figure(2)
plot(time(idxp), y_avg(k,idxp), 'k', 'linewidth', lwd)
hold on
plot(time(idxp), y_pred(idxp-t0), 'b', 'linewidth', lwd-0.5)
ylabel('Number of Spikes in dt', 'FontSize', fsz)
xlabel('Time', 'FontSize', fsz)
set(gca, 'FontSize', fsz)
legend('Data', 'Predicted', 'FontSize', fsz)
ylim([0, max(max(y_avg))])


