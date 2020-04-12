g=input('Enter Dir : ');
c=g; c=[c '\'];
g=['dir /b ' g '\*.mat'];
[s w1]=system(g);
fnd=find(double(w1)==10);
start_index=1;

for ii=1:length(fnd)
    flname=w1(start_index:fnd(ii)-1);
    flname1=w1(start_index:fnd(ii)-5)
    start_index=fnd(ii)+1;
  figure(ii);    
    load([c flname]);
   d2= [flname1,'-ms'];
  
h=figure;
x=imagesc(dwt,dwf*1000,log(stim_stat_zero));axis square;
colormap(jet);
caxis([-8 11.5])
axis xy;
hold on;

colorbar;
set(findobj('Tag','Colorbar'),'FontSize',14,'FontName','Times New Roman');
hold on;

power_tol=0.10;
ranktol=power_tol*max(max(abs(stim_stat_zero)));
contour(dwt,dwf*1000,log(abs(stim_stat_zero)),log([ranktol,-ranktol]),'k-');
hold on;

power_tol=0.010;
ranktol=power_tol*max(max(abs(stim_stat_zero)));
contour(dwt,dwf*1000,log(abs(stim_stat_zero)),log([ranktol,-ranktol]),'k-');
hold on;

power_tol=0.0010;
ranktol=power_tol*max(max(abs(stim_stat_zero)));
contour(dwt,dwf*1000,log(abs(stim_stat_zero)),log([ranktol,-ranktol]),'k-');
axis([-50 50 -1.3  1.3]);
title(flname1,'FontSize',16);
xlabel('Temporal Modulations (Hz)','FontName','Times New Roman','FontSize',16);
ylabel('Spectral Modulations (Cycles/kHz)','FontName','Times New Roman','FontSize',16);
set(gca,'FontSize',16);

saveas(h,flname1,'fig');
saveas(h,flname1,'jpg');
close all;
end