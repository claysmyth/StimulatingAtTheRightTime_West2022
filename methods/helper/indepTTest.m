function p = indepTTest(Xbar1,Xbar2,S1,S2,N1,N2)
S1 = S1.^2; S2 = S2.^2; % convert SD to Var
% Compute T-statistic
    numer = Xbar1-Xbar2;
    sd = sqrt((S1./N1)+(S2./N2));
    tstat = numer./sd;
% Compute degrees of freedom
    df_num = ((S1./N1)+(S2./N2)).^2;
    df_denom = ((S1./N1).^2)./(N1-1) +...
        ((S2./N2).^2)./(N2-1);
    df = df_num./df_denom;
    
    p =  arrayfun(@(x,y) 2*tcdf(-abs(x),y),tstat,df);