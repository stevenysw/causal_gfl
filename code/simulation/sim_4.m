clear;
rng('default')
N = 500;
n = 10000;
K = 5;
d = 10;

sigma = 0.5;

lambda = 0.1;
mse = zeros(1,N);

for k = 1:N
    X = rand(d,n);
    beta = rand(d, 1);
    
    m = 3 * double((X' * beta)' > (round(d/2) - 2));
    m(m == 0) = -3;
    
    e = double((X.^2' * beta)' > (floor(sqrt(d)) - 1));
    e(e == 0) = 0.25;
    e(e == 1) = 0.75;
    
    tau0 = 1;
    
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

    tau_est = causal_gfl(y, z, edge1, edge2, lambda);
    %[beta_est, tau_est] = causal_gfl_cov(y, X, z, edge1, edge2, lambda);
    mse(k) = sum((tau0 - tau_est).^2) / n;
end

mean(mse)
