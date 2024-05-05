function [tau_opt, beta_opt, lambda_opt] = causal_gfl_covariate(y, X, z, edge1, edge2, lambda, varargin)
    n = length(y);
