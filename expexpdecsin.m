function [Y] = expexpdecsin(beta,X)
    A = beta(1);
    B = beta(2);
    C = beta(3);
    P = beta(4);
    W = beta(5);
    Y = exp(exp(P.*X).*(A.*sin(W.*X)+B.*cos(W.*X))+C);
end