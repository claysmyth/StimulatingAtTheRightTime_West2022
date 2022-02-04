function RMSE = fxRMSE(y,yhat)
E = (y - yhat);    % Errors
SE = E.^2; % Squared Error
RMSE = sqrt(sum(SE)./numel(yhat)); %root mean squared error

