function [idx, clust] = group(num_groups, x, time)

% Prep data for clustering
time = time*ones(1, length(x));
X = [x, time'];

% K means clustering
% The indices of the cluster that that point belongs to
idx = kmeans(X,num_groups);


% If I want to get cluster indices
for i = 1:num_groups
    % How many are in each group
    groupk(i) = sum(find(idx == i));
end

[~, idx_g] = sort(groupk);

for i = 1:num_groups
    clust{i} = find(idx == idx_g(i));
end


% Manual clustering
% div = rate_max./num_groups;
% div = div*[1:num_groups];
% 
% y = x;
% 
% for i = 1:length(div)-1
%     % Find top k elements
%     idx{i} = find(x == maxk(y, 3));
%     y = setdiff(y, x(idx{i}));
% 
% %     a = find( x < div(i+1));
% %     b = find( x > div(i));
% %     idx{i} = intersect(a,b);
% end
% %idx{i+1} = find(x > div(i+1));
% idx{i+1} = find(x == y);

end