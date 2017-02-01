% src = prepareStructFromRawData( PbASV_od_0ppb_do_60ppb(:,4:end), [ 0 10 20 30 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
% bkgstruct = prepareStructFromRawData( PbASV_od_0ppb_do_60ppb(:,1), [ 0 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
% for i=1:4:16
%    src.Y(:,i) =  src.Y(:,i) - bkgstruct.Y(:,1);
%    src.Y(:,i+1) =  src.Y(:,i+1) - bkgstruct.Y(:,2);
%    src.Y(:,i+2) =  src.Y(:,i+2) - bkgstruct.Y(:,3);
%    src.Y(:,i+3) =  src.Y(:,i+3) - bkgstruct.Y(:,4);
% end
%== res=48.8283; ci=2.1873; r2=0.99251
% range = [ 20:165];
% peak = [ 70:120 ];
% abc_degree =3;
% abc_looptimes = 20;
%== res=44.5845; ci=1.457; r2=0.9966
% range = [ 20:165];
% peak = [ 70:120 ];
% abc_degree =1;
% abc_looptimes = 20;
%== res=35.8325; ci=0.9434; r2=0.99856
% range = [ 30:170];
% peak = [ 70:120 ];
% abc_degree =1;
% abc_looptimes = 100;

% 
src = prepareStructFromRawData( Tl_120s_vavg_(:,4:8), [ 0 .10 .20 .30 .40 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
% bkgstruct = prepareStructFromRawData( Tl_120s_vavg_(:,1), [ 0 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
% for i=1:4:16
%    src.Y(:,i) =  src.Y(:,i) - bkgstruct.Y(:,1);
%    src.Y(:,i+1) =  src.Y(:,i+1) - bkgstruct.Y(:,2);
%    src.Y(:,i+2) =  src.Y(:,i+2) - bkgstruct.Y(:,3);
%    src.Y(:,i+3) =  src.Y(:,i+3) - bkgstruct.Y(:,4);
% end
%== res=0.32828; ci=0.019698; r2=0.99111
% range = [ 50:150]; 
% peak = [ 80:120 ];
% abc_degree =2;
% abc_looptimes = 20;
%== res=0.34926; ci=0.019769; r2=0.99099
range = [ 35:160];
peak = [ 80:120 ];
abc_degree =3;
abc_looptimes = 20;
%== res= 0.28192; ci=0.013132; r2=0.99596
% range = [ 35:155];
% peak = [ 80:120 ];
% abc_degree =2;
% abc_looptimes = 100;


% src = prepareStructFromRawData( CdASV_od_0ppb_do_60ppb(:,4:end), [ 0 10 20 30 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
% bkgstruct = prepareStructFromRawData( CdASV_od_0ppb_do_60ppb(:,1), [ 0 ], 40, 3, [ 20 23 26 29 ], 'dpasv');
% for i=1:4:16
%    src.Y(:,i) =  src.Y(:,i) - bkgstruct.Y(:,1);
%    src.Y(:,i+1) =  src.Y(:,i+1) - bkgstruct.Y(:,2);
%    src.Y(:,i+2) =  src.Y(:,i+2) - bkgstruct.Y(:,3);
%    src.Y(:,i+3) =  src.Y(:,i+3) - bkgstruct.Y(:,4);
% end
%== res=31.9736; ci=0.53123; r2=0.99953
% range = [ 10:90];
% peak = [ 30:55 ];
% abc_degree =4;
% abc_looptimes = 20;
%== res=31.0435; ci=1.3388; r2=0.99691
% range = [ 10:70];
% peak = [ 30:55 ];
% abc_degree = 2;
% abc_looptimes = 100;


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
for i=1:4
    disp([ 'Wynik #' num2str(i) ': ' num2str(res_finval(i))]);
    disp(['Przeidzal ufnosci #' num2str(i) ': ' num2str(res_Sx0(i)*tinv(1-0.05/2,n-1)/sqrt(n)) ])
    disp(['R^2 #' num2str(i) ': ' num2str(res_gd_{i}.rsquare) ])
end
subplot(221); plot(res_bezbkg(:,1:4:end));
subplot(222); plot(res_bezbkg(:,2:4:end));
subplot(223); plot(res_bezbkg(:,3:4:end));
subplot(224); plot(res_bezbkg(:,4:4:end));