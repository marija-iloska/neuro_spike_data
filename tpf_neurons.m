function [x_est, w_store, x_particles] = tpf_neurons(y, num_groups, T, M, beta, dt)

% Initialize particles
for k = 1:num_groups
    x0(k,1:M) = exprnd(1, 1,M);
end

w_store = zeros(num_groups, T, M);
x_particles = x0;

for t = 2:T

    % Propose particles
    for k = 1:num_groups
        mu_k = x_particles(k,:);
        x(k,:) = exprnd(mu_k);

        % Find best particles    
        log_probs = y(k,t)*log( x(k,:)*dt ) - x(k,:)*dt; % - log(factorial(round(y_avg(k,t))));
        m_star = find(log_probs == max(log_probs));

        if (length(m_star > 1))
            m_star =datasample(1:length(m_star), 1);
        end

        % Modify proposal
        new_mu_k = beta*x(m_star) + (1-beta)*mu_k;

        % Propose new particles
        x_particles(k,:) = exprnd(new_mu_k);

        % Resample
        log_probs = y(k,t)*log( x_particles(k,:)*dt ) - x_particles(k,:)*dt;   
        w  = exp(log_probs - max(log_probs));
        w = w./sum(w);
        w_store(k, t, :) = w;
        idx = datasample(1:M, M, 'Weights',w);

        % Est
        x_particles(k, :) = x_particles(k,idx);
        x_est(k,t) = mean(x_particles(k,:),2);
    end

    

end




end


