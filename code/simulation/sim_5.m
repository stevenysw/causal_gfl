clear;
rng('default')
N = 500;
n = 1000;
K = 5;
d = 10;

sigma = 0.5;

%lambda = logspace(-2,2,50);
lambda = 0.1;
mse = zeros(1,N);
beta_iter = zeros(d,N);

for k = 1:N
    X = rand(d,n);
    beta = [0.1, 0.7, 0.2, 0.4, 0.6, -0.3, -0.5, -0.9, 0.8, -1];
    
    m = beta * X;
    
    e = cellfun(@d05, num2cell(X, 1));
    
    tau0 = cellfun(@d10, num2cell(X, 1));
    
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
    %[beta_est, tau_est] = causal_gfl_cov(y, X, z, edge1, edge2, lambda);;
    mse(k) = sum((tau0 - tau_est).^2) / n;
    beta_iter(:,k) = beta_est;
end

mean(mse)

%for i = 1:10
%    mean(beta_iter(i,:))
%end