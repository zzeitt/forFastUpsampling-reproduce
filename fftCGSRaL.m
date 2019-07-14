function U = fftCGSRaL(G, H)
% fast nonblind deconv using augmented Lagrangian approach (split-Bregman iterations)
%
% Problem formulation: g = H*u + n, 
% where g is a vector of blurred image
% H is a blur matrix, u is the reconstructed image, and n is noise.
%
% Solving:  min_u  { gamma/2*|| g - H*u ||^2 + alpha*||grad(u)||_p^p },
% where 0<=p<=1.
% 
% G ... input LR images in a cell array
% h ... input blurs in a cell array, empty or PSF size 
% gamma ... weight of the fidelity term


maxiter = 10;  % number of iterations
ccreltol = 1e-3;  % convergence criterion
gamma = 1e2;
alpha = 1e-2*gamma;
beta = 1e0*gamma;
Lp = 0.3;  % which Lp norm to use
usize = size(G);  % size of image U
usize(3) = size(G,3);

% vrange ... range of intensity values in each color channel
vrange = zeros(usize(3),2);
for c=1:usize(3)
    vrange(c,:) = [min(reshape(G(:,:,c),[],1)), max(reshape(G(:,:,c),[],1))];
end

% If we work with FFT, we have to move H center into the origin
hshift = zeros(size(H));
hshift(floor(size(H,1)/2)+1, floor(size(H,2)/2)+1) = 1;

% FU ... FFT of u
FU = 0;

% FDx, FDx ... FFT of x and y derivative operators
FDx = repmat(fft2([1 -1],usize(1),usize(2)),[1 1 usize(3)]);
FDy = repmat(fft2([1; -1],usize(1),usize(2)),[1 1 usize(3)]);

% FH ... FFT of PSF
FH = repmat(conj(fft2(hshift,usize(1),usize(2))) .* ...
            fft2(H,usize(1),usize(2)), [1 1 usize(3)]); % FT of H (RGB)
FHTH = conj(FH).*FH; % FFT of (H^T)H (RGB)

% FGs ... FFT of H^T*g
% Note that we use edgetaper to reduce border effect
eG = edgetaper(G,H);
FGu = fft2(eG);
FGs = conj(FH).*FGu;

DTD = conj(FDx).*FDx + conj(FDy).*FDy;

% extra variables for Bregman iterations
Bx = zeros(usize);
By = zeros(usize);
Vx = zeros(usize);
Vy = zeros(usize);

% main iteration loop, do everything in the FT domain
for i = 1:maxiter
    
    disp(['nonblind deconv step ',num2str(i)]);
    
    FUp = FU;
    b = FGs + beta/gamma*(conj(FDx).*fft2(Vx+Bx) + conj(FDy).*fft2(Vy+By));
	FU = b./(FHTH + beta/gamma*DTD);
    
    % Prepare my Lp prior
    Pr = asetupLnormPrior(Lp,alpha,beta);
    % get a new estimation of the auxiliary variable v
    % see eq. 2) in help above 
    xD = real(ifft2(FDx.*FU));
    yD = real(ifft2(FDy.*FU));
    xDm = xD - Bx;
    yDm = yD - By;
    nDm = repmat(sqrt(sum(xDm.^2,3) + sum(yDm.^2,3)),[1 1 usize(3)]);
    Vy = Pr.fh(yDm,nDm);
    Vx = Pr.fh(xDm,nDm);
    
    % update Bregman variables
    Bx = Bx + Vx - xD;
    By = By + Vy - yD;
    
    % we can increase beta after every iteration
    % it should help converegence but probably not necessary
    %beta = 2*beta;
    % Calculate relative convergence criterion
    relcon = sqrt(sum(abs(FUp(:)-FU(:)).^2))/sqrt(sum(abs(FU(:)).^2));
    
    disp(['relcon:',num2str([relcon])]);
    if relcon < ccreltol
        break;
    end
end

U = real(ifft2(FU));
% impose constraints on U 
U = uConstr(U,vrange);

end

function newU = uConstr(U,vrange)
%
% newU = uConstr(U,vrange)
%
% impose constraints on the image U
%
% The only constraint currently implemented is the intensity value range.
% We want to have our estimated HR image inside the itensity range of 
% input images.

newU = U;
for c = 1:size(U,3)
	m = false(size(U));
	m(:,:,c) = U(:,:,c)<vrange(c,1);
	newU(m) = vrange(c,1);
	m(:,:,c) = U(:,:,c)>vrange(c,2);
	newU(m) = vrange(c,2);
end
end