function P = asetupLnormPrior(q,alpha,beta)
switch q
    case 1
        v_star = 0;
        u_star = alpha/beta; 
    case 0
        v_star = sqrt(2*alpha/beta);
        u_star = sqrt(2*alpha/beta);
    otherwise % for 0<q<1
        leftmarker = fzero( @(v) -v+alpha/beta*v^(q-1)*(1-q)*q, [eps 10]);
        v_star = fzero( @(v) -0.5*v^2+alpha/beta*v^q*(1-q), [leftmarker 10]);
        u_star = v_star + alpha/beta*q*v_star^(q-1);
end
P.fh = @(x,nx) aLn(x,nx,u_star,u_star-v_star);
end

function V = aLn(DU,normDU,u_star,k)
%
% V = aLn(DU,normDU,u_star,k)
%
% half-quadratic algorithm 
% additive version
%
% For prior function phi(s) = alpha*|s-t|^q, |s|>u_star 
% phi(s) = beta/2*s^2, |s|<=u_star 
%
% Shrinkage formula:
% min_v { (u-v)^2 + lambda*|v| } 
% v_min = sign(u)*max(|u|-lambda,0);
%
% This is generalized for problems with |v|^q, 0<=q<=1

V = zeros(size(DU));

DUp = DU(normDU>u_star);
normDUp = normDU(normDU>u_star);
V(normDU>u_star) = DUp.*(normDUp-k)./normDUp; 
end