function z = logfit(x,y,stpnt)
if nargin<3
    stpnt = [1 1];
end

logtype =fittype('a +b*log(x)','dependent', {'y'}, 'independent',{'x'},...
                 'coefficients',{'a','b'});
z=fit(x,y,logtype,'StartPoint',stpnt);