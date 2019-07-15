% ¡¾Reproduce "Fast Upsampling"¡¿
% This is the entrance of the whole project.
% An implementation of "Fast Image Upsampling".

%% Initialize parameters
UPSAMP.factor = 4;  % magnification factor
UPSAMP.iter = 4;  % upsampling feed-back iterations
DECONV.lambda_1 = 0.01;  % lambda_1
DECONV.lambda_2 = 5;  % lambda_2
DECONV.iter_max = 10;  % max deconv iterations
GAU.size = 11;  % size of gaussian kernel
GAU.var = 1.5;  % variance of gaussian kernel
IMG.name = 'street';
IMG.read = ['images/', IMG.name, '_small.jpg'];

%% Load original image and split channels
% Stuff about image L (input low resolution image).
L_RGB = im2double(imread(IMG.read));
L_YCbCr = rgb2ycbcr(L_RGB);

%% Show image L
% Image L duplicatedly resized to show.
L_show_Y = pixeldup(L_YCbCr(:,:,1), UPSAMP.factor);
L_show_Cb = pixeldup(L_YCbCr(:,:,2), UPSAMP.factor);
L_show_Cr = pixeldup(L_YCbCr(:,:,3), UPSAMP.factor);
L_show = cat(3, L_show_Y, L_show_Cb, L_show_Cr);
L_show = ycbcr2rgb(L_show);
clear L_show_Y L_show_Cb L_show_Cr;
close all;
FIG_HANDLE.L = figure('Name', 'IMG - Low Resolution');
movegui(FIG_HANDLE.L, 'west');
imshow(L_show);
title('Original image has been resized to display.');

%% Upsample initially
H_tilde_Y = im2double(imresize(L_YCbCr(:,:,1), UPSAMP.factor, 'bicubic'));
H_Cb = im2double(imresize(L_YCbCr(:,:,2), UPSAMP.factor, 'bicubic'));
H_Cr = im2double(imresize(L_YCbCr(:,:,3), UPSAMP.factor, 'bicubic'));
H_bicubic = cat(3, H_tilde_Y, H_Cb, H_Cr);
H_bicubic = ycbcr2rgb(H_bicubic);

%% Feed-back control loop (**essential step**)
for i = 1:UPSAMP.iter
    %¡¾Normalize¡¿
    [H_tilde_Y_norm, low, gap] = simpnormimg(H_tilde_Y);
    
    %¡¾Non-blind deconvolution¡¿
    H_star_Y = nbDeconv(H_tilde_Y_norm, DECONV, GAU);
    IMG.appro = 'mine';
    IMG.other = [num2str(DECONV.lambda_2), '_lam2',...
                 '_', num2str(DECONV.iter_max), '_inner',...
                 '_', 'color'];

    %¡¾Print and denormalize¡¿
    disp([newline,' ============>¡¾Iteration ',...
    num2str(i), ' done!¡¿',newline]);
    H_star_Y = H_star_Y*gap + low;
    
    if(i == UPSAMP.iter),break;end
    
    %¡¾Reconvolution¡¿
    PSF = fspecial('gaussian', 3, 1);
    H_s_Y = imfilter(H_star_Y, PSF, 'same', 'conv');
    
    %¡¾Pixel substitution¡¿
    H_tilde_Y = H_s_Y;
    H_tilde_Y(1:UPSAMP.factor:end,...
              1:UPSAMP.factor:end) = L_YCbCr(:,:,1);
end

%% Show and save result
H_star = cat(3, H_star_Y, H_Cb, H_Cr);
H_star = ycbcr2rgb(H_star);
FIG_HANDLE.H = figure('Name', 'IMG - High Resolution');
movegui(FIG_HANDLE.H, 'east');
figure(FIG_HANDLE.H);
imshow(H_star);
title(['After ', num2str(i), ' iterations']);

% Write in the disk
IMG.folder = ['results/', IMG.name, '/'];
IMG.src = [IMG.folder, IMG.name, '_lr.jpg'];
IMG.bicubic = [IMG.folder, IMG.name, '_bicubic',...
    ['_', num2str(UPSAMP.factor), 'x'], '.jpg'];
IMG.write = [IMG.folder, IMG.name,...
    ['_', IMG.appro,],...
    ['_', num2str(UPSAMP.factor), 'x'],...
    ['_', num2str(UPSAMP.iter), '_outer'],...
    ['_', num2str(GAU.size),'-', num2str(GAU.var), '_gau'],...
    ['_', IMG.other],...
    '.jpg'];
% Create folder
if ~isfolder(IMG.folder)
    mkdir(IMG.folder)
end
% Save origin image
if ~isfile(IMG.src)
    imwrite(L_show, IMG.src);
    disp([newline,' ============>¡¾Image "',...
        IMG.src, '" saved!¡¿',newline]);
end
% Save bicubic image
if ~isfile(IMG.bicubic)
    imwrite(H_bicubic, IMG.bicubic);
    disp([newline,' ============>¡¾Image "',...
        IMG.bicubic, '" saved!¡¿',newline]);
end
% Save result image
if ~isfile(IMG.write)
    imwrite(H_star, IMG.write);
    disp([newline,' ============>¡¾Image "',...
        IMG.write, '" saved!¡¿',newline]);
end

%% Function on call
function [I, low, gap] = simpnormimg(G)

lb = min(G(:));
ub = max(G(:));

gap = ub-lb;
low = lb;
I = (G - low) ./ gap;
end

   