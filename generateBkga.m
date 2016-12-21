src = dataC.Y;
range = [ 100:400];
peak = [ 175:275 ];
s(:,1) = [ 1 2 3 4 5 ];
s(:,2) = [ 6 7 8 9 10 ];
s(:,3) = [ 11 12 13 14 15 ];
s(:,4) = [ 16 17 18 19 20 ];
s(:,5) = [ 21 22 23 24 25 ];

[ bkg, bezbkg ] = bkgautomatic([ [1:length(range)]' src(range,:)*-1], 5, 30);
bezbkg = bezbkg * -1;
dokal = max(bezbkg(peak-range(1),:));
for i=1:5
    linfit(:,i) = polyfit([ 0: 4 ], dokal(s(:,i)+1), 1);
    finval(i) = linfit(2,i) / linfit(1,i);
end
