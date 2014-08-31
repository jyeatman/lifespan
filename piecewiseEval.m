function yhat = piecewiseEval(params,X)
% Evaluate a piecewise linear model
%
% yhat = piecewiseEval(params,X)
%
% Returns the prediction of the piecewise linear model. See piecewiseFit.m
%
% Inputs:
% params(1) = Slope 1
% params(2) = Intercept 1
% params(3) = Slope 2
% params(4) = Cutpoint 1
% params(5) = Slope 3
% params(6) = Cutpoint 2
% X         = Vector of points to evaluate model
%
% Outputs:
% yhat      = Model prediction at each X
%
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.

%% If there are 2 parameters then it is a single line
if length(params) == 2
    % If there are only 2 parameters then it is just a line and can be
    % evaluated with polyval. param(1) is slope param(2) is intercept
    yhat = polyval(params,X);
    
    %% if there ar 4 params it is 2 lines joined at a cutpoint
elseif length(params) == 4
   
    % The polynomial form of the line below the cutpoint
    p1 = params(1:2);
    % The cutpoint
    c = params(4);
    % The intercept of the second line is calulated by evaluating the first
    % line at the cutpoint and subtracting the slope of line 2* the
    % cutpoint
    i2 = polyval(p1,c) - params(3)*c;
    p2 = [params(3) i2];
    % Evaluate the values of each polynomial
    yhat = X; % allocate yhat to be the same size as X
    yhat(X<=c) = polyval(p1,X(X<=c));
    yhat(X>c) = polyval(p2,X(X>c));
    
    %% if there are 6 params then it is 3 lines joined at 2 cutpoints
elseif length(params) == 6
    
    % The polynomial form of the line below the cutpoint
    p1 = params(1:2);
    % The cutpoint
    c = params(4);
    % The intercept of the second line is calulated by evaluating the first
    % line at the cutpoint and subtracting the slope of line 2* the
    % cutpoint
    i2 = polyval(p1,c) - params(3)*c;
    p2 = [params(3) i2];
    
    % The intercept of the third line is calulated by evaluating the second
    % line at the cutpoint and subtracting the slope of line 3* the
    % cutpoint
    c2 = params(6);
    i3 = polyval(p2,c2) - params(5)*c2;
    p3 = [params(5) i3];
    
    % Evaluate the values of each polynomial
    yhat = X; % allocate yhat to be the same size as X
    yhat(X<=c) = polyval(p1,X(X<=c));
    yhat(X>c & X<=c2) = polyval(p2,X(X>c & X<=c2));
    yhat(X>c2) = polyval(p3,X(X>c2));

end

