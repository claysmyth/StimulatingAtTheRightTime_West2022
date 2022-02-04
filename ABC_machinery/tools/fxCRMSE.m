function CRMSE = fxCRMSE(y,yhat)
E = angle(exp(1i*y)./exp(1i*yhat));   % Errors
SE = E.^2; % Squared Error
CRMSE = sqrt(sum(SE)./numel(yhat)); %root mean squared error

