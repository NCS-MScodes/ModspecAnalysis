clear
g=input('Enter Dir : ');
c=g; c=[c '\'];
g=['dir /b ' g '\*.mat'];
[s w1]=system(g);
fnd=find(double(w1)==10);
start_index=1;


fid=fopen('Quantifiers01.xls','a');
fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n','Name', 'Separability', 'Shammasep', 'Asymmetry', 'Moddepth', 'area(0.1)', 'area(0.01)', 'area(0.001)', 'shapesep');

for ii=1:length(fnd)
    flname=w1(start_index:fnd(ii)-1);
    flname1=w1(start_index:fnd(ii)-5)
    d=flname1
    start_index=fnd(ii)+1;
    
load([c flname]);


% load matrix and mparameters required
z=(stim_stat_zero);

% z(fp2:nwf,:) = 0;
power_tol=[0.1 0.01 0.001];
%tolerance
% tolerance=0.001;
    
% --------------calculate separability and shamma separability---------------
u=svd(z);
Separability=u(1)/sum(u);
Shammasep=1-((u(1)^2)/sum(u.^2));

%----------------- Asymmetry--------------------------------
p1_bs=0;
p2_bs=0;
% power in (-+) quadrant
    for j=1:nb-1-1
        for i=ntt+1+1:nt
            p1_bs=p1_bs+z(j,i);
        end
    end   
% power in (++) quadrant
    for j=1:nb-1-1
        for i=1:ntt-1-1
            p2_bs=p2_bs+z(j,i);
        end
    end
    
Asymmetry=(p2_bs-p1_bs)/(p2_bs+p1_bs)
 
%--------------- Modulation Depth
totalpower=sum(sum(abs(stim_stat_zero)));
dcpower=stim_stat_zero(nb,ntt);
Moddepth=realsqrt((totalpower-dcpower)/dcpower) ;

%---------------------contour area-----------------------

for jj=1:length(power_tol)
  dstart=find(dwf==0);
  dend=length(dwf);
  stim_selected=stim_stat_zero(dstart:dend,:);
  powerlimit=power_tol(jj)*max(max(abs(stim_stat_zero)))
  area(jj)=length(find(stim_selected>=powerlimit));
  
end
%shape separability
maxstim=max(max(abs(stim_stat_zero)));
stim_stat_zero(find(stim_stat_zero<power_tol(2)*maxstim))=0;
stim_stat_zero(find(stim_stat_zero>=power_tol(2)*maxstim))=1;
p=(stim_stat_zero);
u1=svd(p);
shapesep=u1(1)/sum(u1);


 fprintf(fid,'%s\t %4.4f\t %4.4f\t %4.4f\t %4.4f\t %4.4f\t %4.4f\t %4.4f\t %4.4f\n',flname1, Separability, Shammasep, Asymmetry, Moddepth, area(1), area(2), area(3), shapesep);

end 
fclose(fid);