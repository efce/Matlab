function [cf_, gd_, Sx0]=fitPoly1( x, y )
    x=x(:);
    x_1 = x;
    y=y(:);
    [ cf_, gd_ ] = fit(x_1,y,'poly1');
    tSyx=0;
    for i=1:1:numel(y)
        tSyx=tSyx+( y(i) - (cf_.p1*x(i)+cf_.p2) )^2;
    end
    if ( numel(y) > 2 )
        Syx=sqrt(tSyx/(numel(y)-2));
    else
        Sx0=NaN;
    end
    
    sumX=0;
    for i=1:1:numel(y)
        sumX = sumX + ( x(i) - mean(x) )^2;
    end
    Sx0 = ( Syx/cf_.p1 ) * sqrt ( 1/numel(y) + 1/numel(y) + ( ( y(1)-mean(y) )^2 / ( cf_.p1^2 * sumX ) ) );
end

