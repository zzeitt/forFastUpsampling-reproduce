function [H_Y] = showH(fig_handle_arg, iter_num, H_Y, H_Cb, H_Cr)
%SHOWDSTIMAGE ¡¾Show image H on specific figure.¡¿
%   Pass in handel of figure and the
%   number of iteration, function will
%   concantenate and display image H.

%% Constraint intensity of H_Y
mapminmax(H_Y(:),-1, 1);

%% Concantenate channels
H_YCrCb = cat(3, H_Y, H_Cb, H_Cr);
H_RGB = ycbcr2rgb(H_YCrCb);

%% Display image
figure(fig_handle_arg);
imshow(H_RGB);
title(['Iteration ', num2str(iter_num)]);
end

