clear;
rng('default')
N = 500;
k = 100;
n = k*k;

x1 = repelem((1:k)/k, k);
x2 = repmat((1:k)/k,[1,k]);
X = [x1; x2];

sigma = 1;

edge1 = [];
edge2 = [];

for i = 1:k-1
    for j = 1:k
        edge1(end + 1) = (j - 1) * k + i;
        edge2(end + 1) = (j - 1) * k + i + 1;
        edge1(end + 1) = (i - 1) * k + j;
        edge2(end + 1) = (i - 1) * k + j + k;
    end
end

edge1 = edge1';
edge2 = edge2';

%lambda = logspace(-2,2,50);
lambda = 1;
mse = zeros(1,N);

e = double((X(1,:) - 0.25).^2 + (X(2,:) - 0.25).^2 < (X(1,:) - 0.75).^2 + (X(2,:) - 0.75).^2);
e(e == 0) = 0.1;
e(e == 1) = 0.9;
m = double(0.25 * X(1,:) + 0.75 * X(2,:) < 0.5);
tau0 = X(1,:) < 0.6 & X(2,:) < 0.6;

for k = 1:N
    z = binornd(1, e);

    y = (m + (z - 0.5) .* tau0 + sigma * normrnd(0, 1, [1, n]));

    tau_est = causal_gfl(y, z, edge1, edge2, lambda);
    %[beta_est, tau_est] = causal_gfl_cov(y, X, z, edge1, edge2, lambda);
    mse(k) = sum((tau0 - tau_est).^2) / n;
end

mean(mse)

subplot(1,2,1)
scatter(X(1,:), X(2,:), 5, tau0)
xlabel('X1') 
ylabel('X2') 

subplot(1,2,2)
scatter(X(1,:), X(2,:), 5, tau_est)
xlabel('X1') 
ylabel('X2') 
colorbar()
caxis([0 1])
