function MSE = fxMSE(y,yhat)
E = (y - yhat);    % Errors
SE = E.^2; % Squared Error
MSE = sum(SE)./numel(yhat); %mean squared error

