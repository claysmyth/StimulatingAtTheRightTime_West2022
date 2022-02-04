function z = logisticfit(x,y,stpnt)
if nargin<3
    stpnt = [1 1 1];
end
% options = fitoptions;
options.MaxIter = 600
z = fit(x,y,'a./(1+exp(b+c*x)) + d','MaxIter',1e8,'StartPoint',stpnt);
% logtype =fittype('(a./(1+exp(b+c*x)))','dependent', {'y'}, 'independent',{'x'},...
%                  'coefficients',{'a','b','c'});
% z=fit(x,y,logtype,'StartPoint',stpnt);