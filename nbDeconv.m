function [H_star_ret] = nbDeconv(H_tilde_arg, ...
    deconv, gau)
% NBDECONV ¡¾Execute non-blind deconvolution.¡¿
%   Accept:
%       H_tilde_arg:    degraded image
%       lam_1:          lambda_1
%       lam_2:          lambda_2
%       iter_max_arg;   iteration times of
%                       minimization loop
%   Variables list:
%       g_size:         the size of image passed in,
%                       also be the size everywhere.
%       mu_x, mu_y:     variables used in minimization
%       f_PSF:          convolution kernel
%       partial_x, partial_y: gradient matrix
%   Return:
%       H_star_ret
%       

%% Initialize parameters
gau_size = gau.size;
gau_var = gau.var;
cvg_crt = 1e-3;  % converge criterion
iter_max = deconv.iter_max;
lam_1 = deconv.lambda_1;
lam_2 = deconv.lambda_2;

%% Initialize variables
g_size = size(H_tilde_arg);
mu_x = zeros(g_size);
mu_y = zeros(g_size);
% [mu_x, mu_y] = imgradientxy(H_tilde_arg);
f_PSF = fspecial('gaussian', gau_size, gau_var);
partial_x = [1 -1];
partial_y = [1; -1];
H_tilde_arg = edgetaper(H_tilde_arg, f_PSF);  % process border
F_H_star = 0;  % Initialize F_H_star

%% Convert variables to freq-domain-vars
% First move PSF center to origin and pad
shift_PSF = zeros(size(f_PSF));
shift_PSF(floor(size(f_PSF,1)/2)+1, floor(size(f_PSF,2)/2)+1) = 1;
F_f_PSF = conj(fft2(shift_PSF, g_size(1), g_size(2))) .* ...
    fft2(f_PSF, g_size(1), g_size(2));
% Then transform partial operator
F_partial_x = fft2(partial_x, g_size(1), g_size(2));
F_partial_x_conj = conj(F_partial_x);
F_partial_y = fft2(partial_y, g_size(1), g_size(2));
F_partial_y_conj = conj(F_partial_y);
% Some fourier variables used in H_star estimation
F_f_PSF_conj = conj(F_f_PSF);
F_H_tilde = fft2(H_tilde_arg);
tmp_down = F_f_PSF_conj.*F_f_PSF + lam_2*(...
    F_partial_x_conj.*F_partial_x + ...
    F_partial_y_conj.*F_partial_y);

%% Iteratively minimize the energy
for i = 1:iter_max
    F_H_star_last = F_H_star;
    mu_x_last = mu_x;
    %% Estimate H_star
    % Do the multiplication in freq domain, note that
    % mu_x and mu_y are updated in iterations.
    tmp_up = F_f_PSF_conj.*F_H_tilde + lam_2*(...
        F_partial_x_conj.*fft2(mu_x) + ...
        F_partial_y_conj.*fft2(mu_y));
    F_H_star = tmp_up./tmp_down;
    
    %% Estimate mu
    p_x_H_star = real(ifft2(F_partial_x.*F_H_star));
    p_y_H_star = real(ifft2(F_partial_y.*F_H_star));
    [mu_x, mu_y] = optMu(p_x_H_star, p_y_H_star, ...
        lam_1, lam_2);
%     FIG_HANDLE_p = figure('Name', 'mu_x');
%     movegui(FIG_HANDLE_p, 'center');
%     imshow(mu_x);
    
    %% Judge if converge
    rel_diff = ...
        sqrt(sum(abs(F_H_star_last(:)-F_H_star(:)).^2))/ ...
        sqrt(sum(abs(F_H_star(:)).^2));
    disp([newline, ' F_H_star_diff: ', num2str(rel_diff)]);
    rel_diff_mu_x = ...
        sqrt(sum(abs(mu_x_last(:)-mu_x(:)).^2))/ ...
        sqrt(sum(abs(mu_x(:)).^2));
    disp([' mu_x_diff: ', num2str(rel_diff_mu_x)]);
    disp([' ------------>¡¾Minimization ', ...
        num2str(i), ' done.¡¿']);
    if rel_diff < cvg_crt
        disp(' ------------>¡¾Converged!¡¿');
        break;
    end
    
    %% Increase lambda_2
    lam_2 = 3 * lam_2;
end

%% Return H_star
H_star_ret = real(ifft2(F_H_star));
end

