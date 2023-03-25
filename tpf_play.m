clear all
close all
clc

%load firing_rate_data.mat
load all_dt005.mat

% Number of particles
M = 50;


% Initialize particles
for k = 1:num_groups
    x0(k,1:M) = exprnd(1, 1,M);
end

x_particles = x0;
beta = 0.2;

for t = 2:T

    % Propose particles
    for k = 1:num_groups
        mu_k = x_particles(k,:);
        x(k,:) = exprnd(mu_k);

        % Find best particles    
        log_probs = y_avg(k,t)*log( x(k,:)*dt ) - x(k,:)*dt; % - log(factorial(round(y_avg(k,t))));
        m_star = find(log_probs == max(log_probs));

        if (length(m_star > 1))
            m_star =datasample(1:length(m_star), 1);
        end

        % Modify proposal
        new_mu_k = beta*x(m_star) + (1-beta)*mu_k;

        % Propose new particles
        x_particles(k,:) = exprnd(new_mu_k);

        % Resample
        log_probs = y_avg(k,t)*log( x_particles(k,:)*dt ) - x_particles(k,:)*dt; % - log(factorial(round(y_avg(k,t))));     
        w  = exp(log_probs - max(log_probs));
        w = w./sum(w);
        idx = datasample(1:M, M, 'Weights',w);

        % Est
        x_particles(k, :) = x_particles(k,idx);
        x_est(k,t) = mean(x_particles(k,:),2);
    end

    

end



k = 2;
for t = 1:T
    y_test(t) = mean(poissrnd(x_est(k,t)*dt, 1,1));
end

idx = 5000:5100;

lwd = 1.5;
fsz = 20;
plot(time(idx), y_avg(k,idx), 'k', 'linewidth', lwd)
hold on
plot(time(idx), y_test(idx), 'b', 'linewidth', lwd-0.5)
ylabel('Number of Spikes in dt', 'FontSize', fsz)
xlabel('Time', 'FontSize', fsz)
set(gca, 'FontSize', fsz)
legend('Truth', 'Estimated', 'FontSize', fsz)
ylim([0, max(max(y_avg))])



figure(2)
for k = 1:num_groups
    plot(time(idx), x_est(k,idx), 'linewidth', lwd)
    hold on
end
ylabel('Inferred firing rates', 'FontSize', fsz)
xlabel('Time', 'FontSize', fsz)
set(gca, 'FontSize', fsz)
legend('Group 1', 'Group 2', 'Group 3', 'FontSize', fsz)



