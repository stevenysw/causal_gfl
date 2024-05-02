clear;
rng('default')
N = 500;
n = 10000;
K = 5;
d = 5;

sigma = 0.5;

e = repmat(0.5, [1, n]);
%lambda = logspace(-2,2,50);
lambda = 0.5;
mse = zeros(1,N);

for k = 1:N
    X = rand(d,n);
    
    m = sin(pi*X(1,:).*X(2,:)) + 2*(X(3,:)-0.5).^2 + X(4,:) + 0.5 * X(5,:);
    
    tau0 = cellfun(@d_tau, num2cell(X, 1));
    
    z = binornd(1, e);
    y = m + (z - 0.5) .* tau0 + sigma * normrnd(0, 1, [1, n]);
    
    Id1 = nearestneighbour(X, 'num', K)';
    Id2 = repmat(1:n, 1, K);
    edge = [Id1(:) Id2(:)];
    data = repmat([0 0], [length(edge),1]) - edge;
    dist = sqrt(data(:,1).^2 + data(:,2).^2);
    [c,ia,ib] = unique(dist);
    edgeunique = edge(ia,:);
    edgeunique = sort(edgeunique, 2);
    edge1 = edgeunique(:,1);
    edge2 = edgeunique(:,2);

    tau_est = causal_gfl(y, z, edge1, edge2, lambda, e);
    %[beta_est, tau_est] = causal_gfl_cov(y, X, z, edge1, edge2, lambda);
    mse(k) = sum((tau0 - tau_est).^2) / n;
end

mean(mse)