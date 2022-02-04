function    S = sigmoidin(x,Rz,B)
S     = (1./(1 + exp(-Rz(:).*x(:) + B))) - 1/(1 + exp(B));

if any(S>0.45)
    a = 1;
end