% ¡¾Reproduce "Fast Upsampling"¡¿
% This is the entrance of the whole project.
% An implementation of "Fast Image Upsampling".

%% Initialize parameters
clear; close all; clc;
FACTOR = 3;  % magnification factor
ITERATION = 5;  % times of iteration
DEC.lambda_1 = 0.01;  % parameters
DEC.lambda_2 = 50;
DEC.iter_max =9;  % max iteration times for deconvolution

%% Load original image and split channels
% Stuff about image L (input low resolution image).
L_RGB = im2double(imread('images/ceremony_small.jpg'));
L_YCbCr = rgb2ycbcr(L_RGB);
L_Y = L_YCbCr(:,:,1);
L_Cb = L_YCbCr(:,:,2);
L_Cr = L_YCbCr(:,:,3);
clear L_RGB L_YCbCr;

%% Show image L
% Image L duplicatedly resized to show.
L_show_Y = pixeldup(L_Y, FACTOR);
L_show_Cb = pixeldup(L_Cb, FACTOR);
L_show_Cr = pixeldup(L_Cr, FACTOR);
L_show = cat(3, L_show_Y, L_show_Cb, L_show_Cr);
L_show = ycbcr2rgb(L_show);
clear L_show_Y L_show_Cb L_show_Cr;
close all;
FIG_HANDLE_L = figure('Name', 'IMG - Low Resolution');
movegui(FIG_HANDLE_L, 'west');
imshow(L_show);
title('Original image has been resized to display.');

%% Upsample initially and normalize
% Y channel is upsampled using feed-back control, while
% Cr and Cb channel simply bicubicly interpolation.
H_tilde_Y = im2double(imresize(L_Y, FACTOR, 'bicubic'));
H_Cb = im2double(imresize(L_Cb, FACTOR, 'bicubic'));
H_Cr = im2double(imresize(L_Cr, FACTOR, 'bicubic'));
[~, low, gap] = simpnormimg(H_tilde_Y);
clear L_Cb L_Cr;

%% Feed-back control loop (**essential step**)
FIG_HANDLE_H = figure('Name', 'IMG - High Resolution');
movegui(FIG_HANDLE_H, 'east');
for i = 1:ITERATION
%     [H_tilde_Y, ~, ~] = simpnormimg(H_tilde_Y);
    %¡¾Non-blind deconvolution¡¿
    H_star_Y = nbDeconv(H_tilde_Y, DEC.lambda_1, ...
        DEC.lambda_2, DEC.iter_max);
    disp([newline,' ============>¡¾Iteration ', num2str(i), ' done!¡¿',newline]);
    [H_star_Y, ~, ~] = simpnormimg(H_star_Y);
    H_temp = showH(FIG_HANDLE_H, i, H_star_Y, H_Cb, H_Cr);
    if(i == ITERATION),break;end
    
    %¡¾Reconvolution¡¿
    PSF = fspecial('gaussian', 11, 1.05);
    H_s_Y = imfilter(H_star_Y, PSF, 'same', 'conv');
    
    %¡¾Pixel substitution¡¿
    H_tilde_Y = H_s_Y;
    H_tilde_Y(1:FACTOR:end,1:FACTOR:end) = L_Y(:,:);
end
clear i L_Y H_Cb H_Cr;


%% Function on call
function [I low gap] = simpnormimg(G)

lb = min(G(:));
ub = max(G(:));

gap = ub-lb;
low = lb;
I = (G - low) ./ gap;
end
   