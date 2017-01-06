% src = prepareStructFromRawData( PbASV_od_0ppb_do_60ppb(:,4:end), [ 0 10 20 30 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
% range = [ 20:165];
% peak = [ 70:120 ];
% abc_degree =3;
% abc_looptimes = 20;

src = prepareStructFromRawData( Tl_120s_vavg_(:,4:end), [ 0 .10 .20 .30 .40 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
range = [ 50:150];
peak = [ 80:120 ];
abc_degree =2;
abc_looptimes = 20;

% src = prepareStructFromRawData( CdASV_od_0ppb_do_60ppb(:,4:end), [ 0 10 20 30 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
% range = [ 20:90];
% peak = [ 30:50 ];
% abc_degree =1;
% abc_looptimes = 20;


nrofsens= length(unique(src.SENS));
clear s;
for i=1:nrofsens;
    s(:,i) = (src.SENS == i);
end

clear res_linfit;
clear res_finval
clear res_s_r; 
clear res_s_x0;
clear res_bezbkg; 
clear res_dokal;
clear bkg;

[ bkg, res_bezbkg ] = bkgautomatic([ [1:length(range)]' src.Y(range,:)*-1], abc_degree,abc_looptimes);
res_bezbkg = res_bezbkg * -1;
res_bezbkg = res_bezbkg(:,2:end);
plot(res_bezbkg)
res_dokal = max(res_bezbkg(peak-range(1),:));


for i=1:nrofsens
    n=sum(s(:,i));
    res_dokals = res_dokal(s(:,i));
    concs = src.CONC(s(:,i));
    %res_linfit(:,i) = polyfit(src.CONC(s(:,i)), res_dokal(s(:,i)), 1);
    [res_cf_{i}, res_gd_{i}, res_Sx0(i)]=fitPoly1( concs, res_dokals );
    res_finval(i)=res_cf_{i}.p2/res_cf_{i}.p1;
%     slope = res_linfit(1,i);
%     intercept = res_linfit(2,i);
%     res_finval(i) = intercept / slope;
%     res_s_r(i) = sqrt( (1/(n-2) ) * sum( (res_dokals - (slope*concs + intercept) ).^2 ) );
%     res_s_x0(i) = (res_s_r(i)/intercept) * sqrt( 1 + 1/n + (res_dokals(1)-mean(res_dokals))^2 / (intercept^2*sum( (concs-mean(concs)).^2 )));
end
disp([ 'Wynik: ' num2str(res_finval(1))]);
disp(['Przeidzal ufnosci: ' num2str(res_Sx0(1)*tinv(1-0.05/2,n-1)/sqrt(n)) ])
disp(['R^2: ' num2str(res_gd_{1}.rsquare) ])

