function DEMO_gray_func(factor_arg, iter_outer,...
    lam_2_arg, iter_inner, gau_var)

%% Initialize parameters
UPSAMP.factor = factor_arg;  % magnification factor
UPSAMP.iter = iter_outer;  % upsampling feed-back iterations
DECONV.lambda_1 = 0.01;  % lambda_1
DECONV.lambda_2 = lam_2_arg;  % lambda_2
DECONV.iter_max = iter_inner;  % max deconv iterations
GAU.size = 11;  % size of gaussian kernel
GAU.var = gau_var;  % variance of gaussian kernel
IMG.name = 'ceremony';
IMG.read = ['images/', IMG.name, '_small.jpg'];

%% Load original image and binarize
% Stuff about image L (input low resolution image).
L = im2double(imread(IMG.read));
L = rgb2gray(L);
L = imbinarize(L);

%% Show image L
% Image L duplicatedly resized to show.
L_show = pixeldup(L, UPSAMP.factor);
close all;
FIG_HANDLE.L = figure('Name', 'IMG - Low Resolution');
movegui(FIG_HANDLE.L, 'west');
imshow(L_show);
title('Original image has been resized to display.');

%% Upsample initially
H_tilde = im2double(imresize(L, UPSAMP.factor, 'bicubic'));
H_bicubic = H_tilde;

%% Feed-back control loop (**essential step**)
for i = 1:UPSAMP.iter
    %��Normalize��
    [H_tilde_norm, low, gap] = simpnormimg(H_tilde);
    
    %��Non-blind deconvolution��
    H_star = nbDeconv(H_tilde_norm, DECONV, GAU);
    IMG.appro = 'mine';
    IMG.other = [num2str(DECONV.lambda_2), '_lam2',...
                 '_', num2str(DECONV.iter_max), '_inner'];

    %��Print and denormalize��
    disp([newline,' ============>��Iteration ',...
    num2str(i), ' done!��',newline]);
    H_star = H_star*gap + low;
    
    if(i == UPSAMP.iter),break;end
    
    %��Reconvolution��
    PSF = fspecial('gaussian', 3, 1);
    H_s = imfilter(H_star, PSF, 'same', 'conv');
    
    %��Pixel substitution��
    H_tilde = H_s;
    H_tilde(1:UPSAMP.factor:end,...
            1:UPSAMP.factor:end) = L(:,:);
end

%% Show and save result
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
    disp([newline,' ============>��Image "',...
        IMG.src, '" saved!��',newline]);
end
% Save bicubic image
if ~isfile(IMG.bicubic)
    imwrite(H_bicubic, IMG.bicubic);
    disp([newline,' ============>��Image "',...
        IMG.bicubic, '" saved!��',newline]);
end
% Save result image
if ~isfile(IMG.write)
    imwrite(H_star, IMG.write);
    disp([newline,' ============>��Image "',...
        IMG.write, '" saved!��',newline]);
end

%% Function on call
function [I, low, gap] = simpnormimg(G)

lb = min(G(:));
ub = max(G(:));

gap = ub-lb;
low = lb;
I = (G - low) ./ gap;
end

end
   