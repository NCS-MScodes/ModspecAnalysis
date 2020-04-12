clear
g=input('Enter Dir : ');
c=g; c=[c '\'];
g=['dir /b ' g '\*.mat'];
[s w1]=system(g);
fnd=find(double(w1)==10);
start_index=1;
fid=fopen('SpecTemporalLimit_01.xls','a');
fprintf(fid,'%s\t %s\t %s\t %s\n','FileName', 'Power', 'SpecLimit', 'TempLimit');


for ii=1:length(fnd)
     flname=w1(start_index:fnd(ii)-1);
    flname1=w1(start_index:fnd(ii)-5)
    d=flname1
    start_index=fnd(ii)+1;
  power_tol=0.01;  
load([c flname]);

      ranktol=power_tol*max(max(abs(stim_stat_zero)));
     contour(dwt,dwf*1000,log(abs(stim_stat_zero)),log([ranktol,-ranktol]),'k-');
        a = contour(dwt(301:end),dwf(nb:end)*1000,log(abs(stim_stat_zero([nb:end],[301:end]))),log([ranktol]),'k');
        for i = 1:length(a)
        if a(1,i) ~=  log(ranktol)
        a1(i) = a(1,i); 
        a2(i)=a(2,i);
        end
        end

        for i=1:length(a1)
         if a2(i)==0
           x(i)=a1(i);
         end
        end
        [x2]=find(x~=0);
        tempLimit=min(x(x2));

       [ro col]=find(stim_stat_zero>=ranktol);
       specLimit=dwf(max(ro))*1000
 fprintf(fid,'%s\t %4.4f\t %4.4f\t %4.4f\n',flname1,power_tol,specLimit, tempLimit);

end 

fclose(fid);