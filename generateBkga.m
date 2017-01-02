src = prepareStructFromRawData( CdASV_od_0ppb_do_60ppb(:,4:end), [ 0 10 20 30 ], 40, 3, [ 20 23 26 29 ], 'dpasv');

range = [ 1:130];
peak = [ 21:71 ];

nrofsens= length(unique(src.SENS));
clear s;
for i=1:nrofsens;
    s(:,i) = (src.SENS == i);
end
[ bkg, res_bezbkg ] = bkgautomatic([ [1:length(range)]' src.Y(range,:)*-1], 3,20);
res_bezbkg = res_bezbkg * -1;
res_bezbkg = res_bezbkg(:,2:end);
plot(res_bezbkg)
res_dokal = max(res_bezbkg(peak,:));
clear linfit;
clear finval;

for i=1:nrofsens
    n=sum(s(:,i));
    res_dokals = res_dokal(s(:,i));
    concs = src.CONC(s(:,i));
    res_linfit(:,i) = polyfit(src.CONC(s(:,i)), res_dokal(s(:,i)), 1);
    res_finval(i) = res_linfit(2,i) / res_linfit(1,i);
    res_s_r(i) = sqrt( (1/(n-2)) * sum( (res_dokals - (res_linfit(1,i)*concs+res_linfit(2,i)) ).^2 ) );
    res_s_x0(i) = (res_s_r(i)/res_linfit(2,i)) * sqrt( 1 + 1/n + (res_dokals(1)-mean(res_dokals))^2 / (res_linfit(2,i)^2*sum( (concs-mean(concs)).^2 )));
end