function [mu_x_ret,mu_y_ret] = optMu(p_x_H, p_y_H, ...
    lambda_1_arg, lambda_2_arg)
%OPTMU ¡¾Optimize mu¡¿
%   Given x,y-direction gradients of H_star,
%   return mu_x and mu_y.
%   ==========================================
%   Mathematic model is:
%       Phi(t) = { -k|t|,        |t|<l_t
%                { -(ax^2+b),    |t|>=l_t
%   Optimization goal is:
%       lambda_1 * [|Phi(mu_x)|+|Phi(mu_y)|] +
%       lambda_2 * [(mu_x - p_x_H)^2+(mu_y - p_y_H)^2]

%% Initialize parameters and variables
% parameters:
l_t = 20;
k = 0.16624;
a = 3.2e-5;  % 3.2 * 10^(-5)
g_max = 300;
% variables:
lambda_1 = lambda_1_arg;
lambda_2 = lambda_2_arg;

%% Normalize p_x_H and p_y_H
x_max = max(p_x_H(:));
x_ratio = g_max / x_max;
p_x_H = p_x_H .* x_ratio;
y_max = max(p_y_H(:));
y_ratio = g_max / y_max;
p_y_H = p_y_H .* y_ratio;

%% Calculate mu_x and mu_y
mu_x_ret = calcMu(p_x_H, lambda_1, lambda_2, ...
    l_t, k, a);
mu_y_ret = calcMu(p_y_H, lambda_1, lambda_2, ...
    l_t, k, a);
mu_x_ret = mu_x_ret ./ x_ratio;  % denormalize
mu_y_ret = mu_y_ret ./ y_ratio;  % denormalize

end


%% Function on call
function [mu_ret] = calcMu(p_H, lam_1, lam_2, ...
    l_t_arg, k_arg, a_arg)
%CALCMU  ¡¾Calculate mu matrix¡¿
%   Given gradient matrix, return mu matrix.
%   There to 3 cases to discuss.

g_size = size(p_H);
mu_ret = p_H;
mask = false(g_size);
% 1) 
mask(:,:) = p_H(:,:) >= (lam_1*k_arg)/(2*lam_2) & ...
            p_H(:,:) < l_t_arg+(lam_1*k_arg)/(2*lam_2);
mu_ret(mask) = p_H(mask)-(lam_1*k_arg)/(2*lam_2);
% 2)
mask(:,:) = p_H(:,:) > -l_t_arg-(lam_1*k_arg)/(2*lam_2) & ...
            p_H(:,:) < -(lam_1*k_arg)/(2*lam_2);
mu_ret(mask) = p_H(mask)+(lam_1*k_arg)/(2*lam_2);
% 3)
mask(:,:) = p_H(:,:) >= (lam_1*a_arg/lam_2+1)*l_t_arg | ...
            p_H(:,:) <= -(lam_1*a_arg/lam_2+1)*l_t_arg;
mu_ret(mask) = p_H(mask)*lam_2/(lam_1*a_arg+lam_2);
end
