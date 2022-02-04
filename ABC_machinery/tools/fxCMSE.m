function CMSE = fxCMSE(y,yhat)
E = angle(exp(1i*y)./exp(1i*yhat));   % Errors
SE = E.^2; % Squared Error
CMSE = sum(SE)./numel(yhat); %root mean squared error

