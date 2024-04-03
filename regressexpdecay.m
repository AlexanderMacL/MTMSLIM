% function to perform an explicit regression to
% y = exp(px)(a sin(wx) + b cos(wx)) + c
% coeffs = [a,b,c,p,w];

function [coeffs] = regressexpdecay(x,y,shortway,expexp)
    if nargin<4
        expexp = false; % fits to y = exp(exp(px)(a sin(wx) + b cos(wx)) + c) instead
        % to control exactly how peaky this is you would instead have to fit to y = d + exp(exp(px)(a sin(wx) + b cos(wx)) + c)
        % which is not linearisable so that's sad
        if nargin<3
            shortway = false;
        end
    end
    y = y(~isnan(x));
    x = x(~isnan(x));
    n = length(x);
    
    if (expexp)
        Y = y;
        y = log(y);
    end
    
    S(1) = 0;
    SS(1) = 0;
    for k=2:n
        S(k) = S(k-1) + 0.5*(y(k)+y(k-1))*(x(k)-x(k-1));
        SS(k) = SS(k-1) + 0.5*(S(k)+S(k-1))*(x(k)-x(k-1));
    end
    AA = [sum(SS.^2), sum(SS.*S), sum(SS.*(x.^2)), sum(SS.*x), sum(SS); ...
        sum(SS.*S), sum(S.^2), sum(S.*(x.^2)), sum(S.*x), sum(S); ...
        sum(SS.*(x.^2)), sum(S.*(x.^2)), sum(x.^4), sum(x.^3), sum(x.^2); ...
        sum(SS.*x), sum(S.*x), sum(x.^3), sum(x.^2), sum(x); ...
        sum(SS), sum(S), sum(x.^2), sum(x), n];
    BB = [sum(SS.*y), sum(S.*y), sum((x.^2).*y), sum(x.*y), sum(y)]';
    X = linsolve(AA,BB);
    A = X(1);
    B = X(2);
    p = B/2;
    w = sqrt(-1*(A+p^2));
    alpha = exp(p.*x).*sin(w.*x);
    beta = exp(p.*x).*cos(w.*x);
    AAA = [sum(alpha.^2), sum(alpha.*beta), sum(alpha); ...
        sum(alpha.*beta), sum(beta.^2), sum(beta); ...
        sum(alpha), sum(beta), n];
    BBB = [sum(y.*alpha), sum(y.*beta), sum(y)]';
    XX = linsolve(AAA,BBB);
    a = XX(1);
    b = XX(2);
    c = XX(3);
    if (shortway)
        coeffs = [a,b,c,p,w];
        return;
    end
    rho = sqrt(a^2+b^2);
    phi = atan(b/a);
    if (a<0); phi = phi + pi; end
    K = round((w.*x+phi)./pi);
    r = (y-c).*exp(-p.*x);
    theta = zeros(1,n);
    for k = 1:n
        if (rho^2>r(k)^2)
            theta(k) = (-1)^K(k)*atan(r(k)/sqrt(rho^2-r(k)^2))+pi*K(k);
        elseif (r(k)>0)
            theta(k) = pi/2*(-1)^K(k)+pi*K(k);
        else
            theta(k) = -pi/2*(-1)^K(k)+pi*K(k);
        end
    end
    T = linsolve([sum(x.^2),sum(x);sum(x),n],[sum(theta.*x),sum(theta)]')';
    w3 = T(1);
    phi2 = T(2);
    alpha = exp(p.*x).*sin(w3.*x);
    beta = exp(p.*x).*cos(w3.*x);
    AAA = [sum(alpha.^2), sum(alpha.*beta), sum(alpha); ...
        sum(alpha.*beta), sum(beta.^2), sum(beta); ...
        sum(alpha), sum(beta), n];
    BBB = [sum(y.*alpha), sum(y.*beta), sum(y)]';
    XX = linsolve(AAA,BBB);
    a3 = XX(1);
    b3 = XX(2);
    c3 = XX(3);
    coeffs = [a3,b3,c3,p,w3];
end