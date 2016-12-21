function [ fitresult ] = fitFaraic( toBefitted, realT )

ft = fittype( '(a / (x+t)^0.5) + b/c * exp( -(x+t)/d )', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower =      [-10000 -100 10     1   0];
opts.StartPoint = [ 10     30  10000  20  0];
opts.Upper =      [ 10000  100 100000 100 2];

    % Fit model to data.
    imax = (length(toBefitted)/realT);
    for ( i= 1:1:imax)
        disp( [ num2str(i) ' / ' num2str(imax)]);
        Y = toBefitted( ((i*realT)-realT+1) : (i*realT));
        [ fitresult{i} ] = fit( [1:realT]', Y, ft, opts );
    end
    
end

%FITFARAIC Summary of this function goes here
%   Detailed explanation goes here