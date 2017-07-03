function [bkg,bufo]=bkgautomatic(buf, deg, nn)
%
% buf - w kolumnach potencja�, pr�d
% deg - stopie� wielomianu
% nn - kroki
%
[m,n]=size(buf);

x=buf(:,1);
bufo=x;
bkg=x;

for (k=2:n)
y=buf(:,k);
y=smooth(y,11,'sgolay',3);

y1=y;
for (i=1:nn)
    a=polyfit(x,y1,deg);
    y2=polyval(a,x);
    for (j=1:m)
        if (y1(j) < y2(j)) y1(j)=y2(j); end;
    end;
 end;
 bufo(:,k)=y-y1;
 bkg(:,k)=y1;
end;

% figure;
% subplot(3,1,1); plot(buf(:,2:n),'color', [0.5 0.5 0.5]); % xlim([-1248 -700 ]); set(gca,'XDir','reverse');
% %ylabel('Pr�d / \muA');
% %ylabel('Pr�d / j.a.');
% subplot(3,1,2); plot(buf(:,2:n),'color', [0.5 0.5 0.5]); % xlim([-1248 -700 ]); set(gca,'XDir','reverse');
% %ylabel('Pr�d / \muA');
% %ylabel('Pr�d / j.a.');
% hold on; plot(x,bkg(:,2:n),'k', 'LineWidth', 2); hold off;
% subplot(3,1,3); plot(bufo(:,2:n),'color', [0.5 0.5 0.5]); % xlim([-1248 -700 ]); set(gca,'XDir','reverse');
% ylabel('Pr�d / \muA');
% xlabel('Potencja� / mV vs. Ag/AgCl');
%ylabel('Pr�d / j.a.');
%xlabel('Numer punktu');


% figure;
% plot(bufo(:,2:n),'r');


% bkg=y1;
% bufo=y-y1;