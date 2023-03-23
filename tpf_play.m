clear all
close all
clc

load firing_rate_data.mat


% Number of particles
M = 20;


% Initialize particles
a = 10; b = 10; var = 0.1; var_k = 1;
for k = 1:num_groups
    x0(k,1:M) = gamrnd(a,b, 1,M);
end

x_particles = x0;
beta = 0.2;

for t = 2:T

    % Propose particles
    for k = 1:num_groups
        mu_k = x_particles(k,:);
        a = mu_k.^2./var_k;
        b = var_k./mu_k;
        %x(k,:) = gamrnd(a, b) + eps;
        x(k,:) = exprnd(mu_k);
        %x(k,:) = normrnd(mu_k, var_k);

        % Find best particles    
        log_probs = y_avg(k,t)*log( x(k,:)*dt ) - x(k,:)*dt; % - log(factorial(round(y_avg(k,t))));
        %probs = ((x(k,:)*dt).^y_avg(k,t)).*exp(-(x(k,:)*dt));
        m_star = find(log_probs == max(log_probs));

        if (length(m_star > 1))
            m_star =datasample(1:length(m_star), 1);
        end

        % Modify proposal
        new_mu_k = beta*x(m_star) + (1-beta)*mu_k;
        new_var_k = beta*var_k + (1-beta)*var;

        % Propose new particles
        a = new_mu_k.^2./new_var_k;
        b = new_var_k./new_mu_k;
        %x_particles(k,:) = gamrnd(a,b);
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

k = 1;
for t = 1:T
    y_test(t) = mean(poissrnd(x_est(k,t)*dt, 1,1));
end

lwd = 1.5;
fsz = 20;
plot(time, y_avg(k,:), 'k', 'linewidth', lwd)
hold on
plot(time, y_test, 'g', 'linewidth', lwd-0.5)
ylabel('Number of Spikes in dt', 'FontSize', fsz)
xlabel('Time', 'FontSize', fsz)
set(gca, 'FontSize', fsz)
legend('Truth', 'Estimated', 'FontSize', fsz)

