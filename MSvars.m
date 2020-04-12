g=input('Enter Dir : ');
c=g; c=[c '\'];
g=['dir /b ' g '\*.wav'];
[s w1]=system(g);
fnd=find(double(w1)==10);
start_index=1;

for i=1:length(fnd)
    flname=w1(start_index:fnd(i)-1);
    flname1=w1(start_index:(fnd(i)-5));
    d=flname1;
   
    start_index=fnd(i)+1;
    % read the sound file
    [x fs nbits]=wavread([c flname]);
    %reading the wav file for modspec
input=x;
  
  % Set the parameters
%subplot(1,2,2);
fwidthHz=45; % in Hz (ask this parameter)
fstep=fwidthHz;
fwidthms=1000/(2*pi*fwidthHz);
ampsamprate=1000;  % or commonly used value 10*fwidthms
samprate=fs;
increment=floor(fs/ampsamprate);
winLength=floor(6*fs/(2*pi*fwidthHz));
doFFTshift=1;
if size(input,1)>1
    input=input';
end
inputLength=length(input);
if inputLength < winLength
    input(winLength)=0;
    inputLength=winLength;
end
frameCount=floor((inputLength-winLength)/increment)+1;
fftLen=winLength;
wx2=((1:winLength)-winLength/2).^2;
wvar=(winLength/6)^2;
ws=exp(-0.5*(wx2./wvar)); % Gaussian Window
if rem(fftLen,2)
    s=zeros((fftLen+1)/2+1,frameCount); % winLength is odd
else
    s=zeros(fftLen/2+1,frameCount); % winLength is even
end
pg=zeros(1,frameCount);
for i=1:frameCount
    start=(i-1)*increment+1 ;
    last=start+winLength-1 ;
    
    f=zeros(fftLen,1);
    f(1:winLength)=ws.*input(start:last);
    pg(i)=std(f(1:winLength));
    
    if doFFTshift
        f=[f(round(winLength/2)+1:winLength); zeros(fftLen-winLength,1); f(1:round(winLength/2))];
    end

    specslice=fft(f);
    if rem(fftLen,2)
        s(:,i)=specslice(1:((fftLen+1)/2+1));
    else
        s(:,i)=specslice(1:(fftLen/2+1));
    end
end
tmp = log(abs(s) +1); % chk if v need to take log10
if rem(size(tmp,1), 2)    
   cutoff = (size(tmp,1) +1)/2;
else
   cutoff = size(tmp,1)/2+1;
end
outSpectrum = 10*tmp(2:cutoff,:);
%imagesc(20*log10(abs(outSpectrum)+eps)); axis xy; colormap(jet);
s=outSpectrum;
%save s s;
%clear ;
%load s.mat ;
t=s; % where y is the spectrogram
clear s;
nband=size(t,1);   % number of frequency bands
nlen=size(t,2);	      %  number of time points or windows
ntrails=size(t,1)     %  number of frequency bands
stim_env=t ; 	      % specgram
TimeLag=300            % check the value of timelag
twindow=[-TimeLag TimeLag];

stim_avg=0;
count_avg=0;
%Step 7.2 	take log of spec
%stim_env=log(stim_env+1.0); % spec is already a log chk this one
%Step 7.3	total energies in each band
stim_avg=sum(stim_env*ntrails,2);
%Step 7.4     count
count_avg=nlen*ntrails;
%Step 7.5	average
stim_avg=stim_avg/count_avg;

tot_corr=diff(twindow)+1;
spa_corr=(nband*(nband-1))/2+nband;
CS=zeros(spa_corr,tot_corr);
%stim_env=t;   % specgram
xb=1;
stimval=zeros(nband,nlen);

for tgood = 1:nlen
      stimval(:, tgood) = stim_env(1:nband, tgood) - stim_avg(1:nband);
end

for ib1=1:nband
	     for ib2=ib1:nband
			CS(xb,:)=CS(xb,:)+xcorr(stimval(ib2,:),stimval(ib1,:),twindow(2))*ntrails;
			       xb=xb+1;
	     end
end

lengthVec=ones(1,nlen);
CS_ns=zeros(1,tot_corr);
CS_ns(1,:)=CS_ns(1,:)+ntrails*xcorr(lengthVec,lengthVec,twindow(2));
nonzero_ns=isinf(1./CS_ns)+CS_ns;
% normalize CS matrix 
for i=1:spa_corr
    CS(i, :) = CS(i, :) ./ nonzero_ns;
end

% Get dimension size of modulation spectrum from CS
ncorr = size(CS, 1);
nt = size(CS, 2);
nb = (-1 + sqrt(1+8*ncorr))/2;

stim_modspec =zeros(2*nb-1,nt);
% Take average of stim along space
stimdiag = zeros(nb, nt);

for xboff=1:nb
    xb=xboff;
    for i=1:nb-xboff+1     
	    %stim_mean(xboff,i) = sqrt(mean(stim(xb,:)));
        stimdiag(xboff,:) = stimdiag(xboff,:)+CS(xb,:);
	    xb = xb + (nb-i+1);
    end
    stimdiag(xboff, :) = stimdiag(xboff, :)/(nb-xboff+1);
end

% Stuff a symmetric correlation by repeating the frequency.
for xb=1:nb
    stim_modspec(nb-1+xb,:) = stimdiag(xb,:); 
    
    if ( xb ~= 1)
        stim_modspec(nb+1-xb,:)=fliplr(stimdiag(xb,:));
        
    end   
end

% Smoothing purpose, define hanning window for windowing
% Get a hanning window for windowing in time and frequency;
wht = hanning(nt);
whf = hanning(2*nb-1);

% Make a 2-d window
whtf = zeros(size(stim_modspec));
for it=1:nt
    whtf(:,it)=whf;
end
for ib=1:2*nb-1
    whtf(ib,:) = whtf(ib,:).*wht';
end

% Inserted on Nov 11, 2005
dt=0:nt-1; % timescale in msec
df=0:2*nb-2; % frequency bands
% Nov 11 insertion over

% Subtract DC value and window
stim_modspec0 = stim_modspec-mean(mean(stim_modspec));
stim_modspec0w = stim_modspec0.*whtf;
stim_modspecw = stim_modspec.*whtf;

% Next take 2d fft of stim_modspec
stim_modspec_f   =fft2(stim_modspec);
stim_modspec0w_f =fft2(stim_modspec0w);
stim_modspec0_f  =fft2(stim_modspec0);
stim_modspecw_f  =fft2(stim_modspecw);

% Prepare for plotting
stim_modspec_f = fftshift(abs(stim_modspec_f));
stim_modspec0w_f = fftshift(abs(stim_modspec0w_f));
stim_modspec0_f = fftshift(abs(stim_modspec0_f));
stim_modspecw_f = fftshift(abs(stim_modspecw_f));

% Find labels for x and y axis
% fstep is the separation between frequency bands
for i=1:2*nb-1
   dwf(i)= (i-nb)*(1/(2*fstep*(nb-1)));
end
% 1 ms (1000 on the numerator) is the sampling rate
ntt=(nt+1)/2;
for i=1:2*ntt-1
   dwt(i) = (i-ntt)*(ampsamprate/(2*(ntt-1)));
end

% Zero outs the modulation spectrum outside the rectangle.
% k is the constant that relates f_width to the upper frequency of the
% sampled rectangle
k = 2.57;
wt_max = k*fwidthHz;

% tp1 and tp2 are the index points for dwt boundaries
nwt = length(dwt);
tp1 = 1;
tp2 = nwt;
for i=1:nwt-1
    if (dwt(i) <= -wt_max & dwt(i+1) > -wt_max )
        tp1 = i;
    end
    if (dwt(i) < wt_max & dwt(i+1)>= wt_max )
        tp2 = i+1;
    end
end

% fp1 and fp2 are the index points for the dwf boundaries
wf_max = 1./wt_max;
nwf = length(dwf);
fp1 = 1;
fp2 = nwf;
for i=1:nwf-1
    if (dwf(i) <= -wf_max & dwf(i+1) > -wf_max )
        fp1 = i;
    end
    if (dwf(i) < wf_max & dwf(i+1)>= wf_max )
        fp2 = i+1;
    end
end

stim_stat_zero = stim_modspec0w_f;
stim_stat_zero(:,1:tp1) = 0;
stim_stat_zero(:,tp2:nwt) = 0;
stim_stat_zero(1:fp1,:) = 0;
stim_stat_zero(fp2:nwf,:) = 0;

save (d,'dwt','dwf','stim_stat_zero','nb','nt','ntt');


end   
    