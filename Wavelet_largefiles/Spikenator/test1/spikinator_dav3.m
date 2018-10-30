function [y,spk] = spikinator_dav3(x,windw,thr,amp,smth,cleanthr,frame,fade)
%
%
%Remove large amplitude spikes from a time series.
%
%[y,spk] = spikinator(x,windw,thr,amp,smth,cleanthr,frame);
%
%x: spikey signal
%windw: spike window size in samples (... not time units); set equal to 1
%       to plot detection threshold before running program in standard 
%       mode        
%thresh: threshold for spike detection (greater or equal to 1); default=1.4
%amp: spike removal filter amplitude; default=2.5e-4
%    ****If amp = 0, THIS SETS PROGRAM TO EXCISION OF SPIKES***
%smth: magnitude of extra smoothing of spike stumps (0 to 1); default=0
%cleanthr: degree of spike cleaning (between 0 and 1); default=0.5
%          higher = greater cleaning, but may incur signal distortion
%frame: window frame margin (fraction of window size); default=0.2
%fade: when in excision mode (amp=0) specifies the region on either side of the
%           excision in which to apply weighted average fading
%           (units are in samples); default=windw/4
%y: cleaned signal
%spk: time series containing removed spikes (or spike markers in excision
%           mode)
%
%
% Osbert Zalay, Oct 2007
%       Fade & excision mode spikes added by David Stanley 2009
%
%*** For academic/research use. Please acknowledge author. ***

% hold off;
%load spikinatorJPG;
%subplot(1,2,1),image(spikinatorJPG)
%axis off

if nargin < 2
    disp('Need to input x and windw to work...exiting');
    y=[]; spk=[];
    return
end
if nargin < 3
    thr=1.4;
    amp=0.00025;
    smth=0;
    cleanthr=0.5; %0.2;
    frame=0.2;
    fade = round(windw/4);
end
if nargin < 4
    amp=0.00025;
    smth=0;
    cleanthr=0.5; %0.2;
    frame=0.2;
    fade = round(windw/4);    
end
if nargin < 5
    smth=0;
    cleanthr=0.5; %0.2;
    frame=0.2;
    fade = round(windw/4);
end
if nargin < 6
    cleanthr=0.5; %0.2;
    frame=0.2;
    fade = round(windw/4);
end
if nargin < 7
    frame=0.2;
    fade = round(windw/4);
end
if nargin < 8
    fade = round(windw/4);
end


testing_mode = 1;       % Set=1 to, when in excision mode, leave gaps where spikes originally
                        % were located, so that excised trace will be
                        % aligned with original.

lenx=length(x);
hwindw=round(windw/2);
cleanthr=2*(1-cleanthr);
spk=x;
x=x.';
xm=x-min(x);
[mv,indv]=lmaxv(xm,1);
%mvbar=mean(mv); %global threshold (mean of peak-to-peak values)
mvbar=rmsv(xm); %global threshold (RMS)
mvcut=thr*mvbar; mvc1=mvcut;
indmv=find(mv > mvcut);
indspk=indv(indmv);
if windw == 1
%     plot([1:lenx],x,[1:lenx],mvcut+min(x),'r',indspk,x(indspk),'*');
    figure;plot([1:lenx],x);
    hold on; plot([1:lenx], (mvcut+min(x))*ones(lenx,1),'r')
    hold on; plot(indspk,x(indspk),'*');
    y=[]; spk=[];
    return;
end


%return
lenspk=length(indspk);
hwdwsize=hwindw;
count=0;
for i=1:lenspk
    indi=indspk(i);
    while (hwdwsize >= indi) | (hwdwsize >= (lenx-indi))
        hwdwsize=round(hwdwsize/2);
    end
    lb=indi-hwdwsize;
    ub=indi+hwdwsize;
    xs=x(lb:ub);
    minxs=min(xs);
    xsm=xs-minxs;
    mvbar=rmsv(xsm(find(~isnan(xsm)))); %local threshold (RMS) - disregarding other nearby spikes (NaNs)
    mvcut=cleanthr*mvbar;
    mvx=lmaxv(xsm,1);
    indmvx=find(mvx > mvcut);
    %subplot(1,2,2),plot(xs);
%     plot(xs);
    pause(0.0001);
    
    
    if ~isempty(indmvx)
        if amp > 0
            lenxs=length(xs);
            mxsm=max(xsm);
            [env,lf,uf]=envelope(lenxs,1,frame);
            xflt=1-env.*amp.*xsm./mxsm;
            while (mxsm > mvcut) & (count < 10000)
            %while (mxsm > (cleanthr*mvbar+a0)) & (count < 10000)
                count=count+1;
                xsm=xflt.*xsm;
                %cm=CofM(xsm,lenxs);
                %indi=indi+(cm-round(lenxs/2));
                if smth~=0
                    lenflt=max([3 ceil(smth*(uf-lf))]);
                    xsm(lf:uf) = smooth(xsm(lf:uf),lenflt);
                end
                xsm=xsm+minxs;
                x(lb:ub)=xsm;
                hwdwsize=hwindw;
                while (hwdwsize >= indi) | (hwdwsize >= (lenx-indi))
                    hwdwsize=round(hwdwsize/2);
                end
                lb=indi-hwdwsize;
                ub=indi+hwdwsize;
                xs=x(lb:ub);
                minxs=min(xs);
                xsm=xs-minxs;
                %a0=mean(xsm);
                mxsm=max(xsm);
                lenxs=length(xs);
                [env,lf,uf]=envelope(lenxs,1,frame);
                xflt=1-env.*amp.*xsm./mxsm;
            end
        end
    end
    if amp <= 0
%                   If in excision mode, we should remove the spike
%                   regardless of whether or not the local thresholding
%                   came back positive, since this could be distorted by
%                   previous spike removal.
        %subplot(1,2,2),plot(xs);
%             plot(xs);
        pause(0.0001);
        x(lb:ub)=NaN;
    end
end



if amp>0                        % daveadded
    y=x(find(~isnan(x))).';
    spk=spk-y;    
elseif (amp <= 0)               % Apply fade
%     Fade operates by overlapping two regions of size fade/2 taken from either
%     side of the excision. Then, it performs a weighted average
%     of these two signals. The weights for the signal taken from the region to the left of
%     excision start at 100% and decrease linearly to zero (when traveling
%     from left to right across the overlapping region). Vivce versa is true
%     for the weights of the right region's signal.

    
    xnan = isnan(x);
    nanstarts = find ( (xnan(1:(end-1))==0 & xnan(2:end)==1) );
    nanstops = find ( (xnan(1:(end-1))==1 & xnan(2:end)==0) ) + 1;
    if ( ~isempty(nanstarts) && ~isempty(nanstops))
        tocheck = 1;
        while (tocheck)     % Incase the spike resides on the edge of the dataset, disregard it
            if nanstarts(1)<(fade)
                x = x(nanstops(1):end); xnan = xnan(nanstops(1):end);
                nanstarts = find ( (xnan(1:(end-1))==0 & xnan(2:end)==1) );
                nanstops = find ( (xnan(1:(end-1))==1 & xnan(2:end)==0) ) + 1;
                tocheck = 1;
                fprintf ('Removing 1 spike from the left\n');
            else
                tocheck = 0;
            end    
            if (nanstops(end) + (fade)) > length(x)
                x = x(1:nanstarts(end)); xnan = xnan(1:nanstarts(end));
                nanstarts = find ( (xnan(1:(end-1))==0 & xnan(2:end)==1) );
                nanstops = find ( (xnan(1:(end-1))==1 & xnan(2:end)==0) ) + 1;
                tocheck = 1;
                fprintf ('Removing 1 spike from the right\n');
            else
                tocheck = 0;
            end
        end
        spk = zeros(1,length(x)); spk(find(isnan(x)))=NaN; spk(nanstarts-fade-1)=1;    % Generate array to mark spike locations
        numgaps = length(nanstarts);
        nanstops = repmat (nanstops',1,fade+1);
        nanstarts = repmat (nanstarts',1,fade+1);
        fademat_ind = repmat(0:fade, numgaps, 1);
        nanstops = nanstops+fademat_ind;
        nanstarts = nanstarts+fliplr(-fademat_ind);
        fleft = x(nanstarts);
        fright = x(nanstops);
        faded = fleft .* repmat(fliplr(0:(1/fade):1), numgaps, 1) + fright .* repmat (0:(1/fade):1, numgaps, 1);

        x(nanstops) = faded;
        x(nanstarts) = NaN;
        spk(nanstarts) = NaN;
    else
        spk = 0;
    end

    xcond=x(find(~isnan(x))).';
    spk = spk*std(xcond)*thr + mean(xcond);
    
    y = x;
    if ~testing_mode
        y = xcond;
        spk = spk(find(~isnan(spk))).';
    end

else
    fprintf ('Debugging time. (You should not see this!)');
end                                     %end_daveadded
    


function [y] = Linflt(x,SR,wn,filt_type,order)
[mx nx]=size(x);
if nx>1 
    x=x.';
    mx=nx;
end
normCutoff=2*wn/SR;
len=mx;
h=fir1(order,normCutoff,filt_type);
H=fft(h,len);
X=fft(x,len).';
Y=H.*X;
y=ifft(Y);
y=shft(SR,order,y);


function [uShift] = shft(SR,order,u)
lenu=length(u);
uShift=zeros(1,lenu);
shiftindex=order/2;
tdelay=(shiftindex)/SR;
k=0;
for i=shiftindex:lenu
    k=i-shiftindex+1;
    uShift(k)=u(i);
end
for i=1:shiftindex-1
    uShift(k+i)=u(i);   
end


function [lmval,indd]=lmaxv(xx,filtr)
x=xx;
len_x = length(x);
	fltr=[1 1 1]/3;
  if nargin <2, filtr=0; 
	else
x1=x(1); x2=x(len_x); 
	for jj=1:filtr,
	c=conv(fltr,x);
	x=c(2:len_x+1);
	x(1)=x1;  
        x(len_x)=x2; 
	end
  end
lmval=[]; indd=[];
i=2;		
    while i < len_x-1,
	if x(i) > x(i-1)
	   if x(i) > x(i+1)	
lmval =[lmval x(i)];
indd = [ indd i];
	   elseif x(i)==x(i+1)&x(i)==x(i+2)	
i = i + 2;  		
	   elseif x(i)==x(i+1)
i = i + 1;		
	   end
	end
	i = i + 1;
    end
if filtr>0 & ~isempty(indd),
	if (indd(1)<= 3)|(indd(length(indd))+2>length(xx)), 
	   rng=1;	
	else rng=2;
	end
	  for ii=1:length(indd), 
	    [val(ii) iind(ii)] = max(xx(indd(ii) -rng:indd(ii) +rng));
	    iind(ii)=indd(ii) + iind(ii)  -rng-1;
	  end
  indd=iind; lmval=val;
else
end


function [env,lf,uf] = envelope(len,divar,frame)
k=100./len;
ivar=[divar:divar:len];
lf=round(frame*len);
uf=round((1-frame)*len);
env=(1+exp(-k.*(ivar-lf))).^(-1)-(1+exp(-k.*(ivar-uf))).^(-1);


function [yrms]=rmsv(y);
len=length(y);
yrms=0;
for i=1:len
    yrms = yrms+y(i)^2;
end
yrms=(yrms/len)^0.5;


function cm = CofM(xx,len)
r=[1:len];
cm=round(sum(xx.*r)/sum(xx));

