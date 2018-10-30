function [Y,tau] = embed(y,d,typ,plt,prnt)
%
%Embed time series in phase space using delay reconstruction method.
%
%[yr,tau] = embed(y,emb,delay,plt);
%
%emb: embedding dimension
%delay: enter a numeric value, or one of the following estimation methods:
%    {'m'  for first minimum of auto mutual information;
%     'm3' for 1/3 of first minimum of auto mutual information (default)
%     'a'  (1-1/e) of autocorrelation (similar to 'm3')
%     'z'  for first zero of autocorrelation}
%plt: plot flag = {0 (default) for no plot of delay estimation function;
%                  1 for a plot} 
%
%yr: time-delay phase space reconstruction of y 
%    (i-th time series = i-th row of yr)
%tau: time delay used
%
%
%Osbert Zalay, March 2008  
%
%*** For academic/research use. Please acknowledge author. ***


if nargin < 3
    typ = 'm3';
    plt=0;
    prnt=0;
end
if nargin < 4
    plt=0;
    prnt=0;
end
if nargin < 5
    prnt=0;
end

[len,m]=size(y);
if m > len
    y=y.';
    len=m;
end

if prnt
    fprintf('\n Embedding...');
end

%estimate time delay, then create time-delay signals and embed
Y=zeros(d,len);
Y(1,:)=y.';
if isnumeric(typ)
    tau=typ;
else
    if (strcmpi(typ,'m') | strcmpi(typ,'m3'))
        tau=optdelayMI(y,typ,plt,prnt);
    else
        tau=optdelayXCorr(y,typ,plt);
    end
end
for i=2:d
    Y(i,:)=circshift(y,-(i-1)*tau).';
end
endind=len-(d-1)*tau;
if endind > 100
    Y=Y(1:d,1:(len-(d-1)*tau));
end
if prnt
    fprintf('Done. \n');
end


function [tau] = optdelayMI(x,typ,plt,prnt)
%estimate the optimal time delay based on first minimum of auto mutual info
tauD=3;
len=length(x);
endpt=round(0.2*len);
hendpt=round(endpt/2);
mi=zeros(hendpt,1);
x=x(1:endpt);
x=x./max(abs(x));
hx = shannonEnt(x,[-1 1],20); %estimate the shannon entropy
cnt=1;
for i=1:hendpt
    mi(i) = hx - conditionalInfo(x,circshift(x,i)',[-1 1],20); %auto mutual info (AMI)
    if (prnt)
        if (cnt == 50)
            cnt=1;
            fprintf('.');
        end
        cnt=cnt+1;
    end
end
mi=mvavg(mi,3,10);
i=1;
while (i<(endpt-1)) & (mi(i+1)<mi(i)) & (i < hendpt-1)
    i=i+1; %find first minimum of AMI
end
if strcmpi(typ,'m')
    tau=i;
else
    tau=round(i/3);
end
if tau < 1
    tau=1;
end
if plt
    ii=[1:length(mi)];
    plot(ii,mi,'b',tau,mi(tau),'r*'), title('Auto Mutual Info');
end


function [tau] = optdelayXCorr(x,typ,plt)
%estimate the optimal time delay based on 'decorrelation time'
tauD=3;
if nargin < 2
    plt = 0;
end
x=x-mean(x);
lx=length(x);
xx=xcorr(x,'coeff');
lxx=length(xx);
xxh=xx(lx:lxx);
lxxh=length(xxh);
xxh=mvavg(xxh,3,10);
i=1;
if strcmpi(typ,'z')
    while (xxh(i) > 0) & (i < lxxh)
        i=i+1; %find first zero-crossing of auto-correlation
    end
    if (i >= lxxh) %if no zero-crossing, use (1-1/e) of autocorrelation
        val=(1-exp(-1))*xxh(1); 
        i=1;
        while (xxh(i) > val) & (i < lxxh)
            i=i+1;
        end
    end
else
    val=(1-exp(-1))*xxh(1); %(1-1/e) of initial value
    i=1;
    while (xxh(i) > val) & (i < lxxh)
        i=i+1;
    end
end
tau=i;
if tau < 1
    tau=1;
end
if plt
    ii=[1:length(xxh)];
    plot(ii,xxh,'b',tau,xxh(tau),'r*'), title('Autocorrelation');
end


function [s,p]=shannonEnt(x,bounds,bincount)
%Estimate the Shannon entropy (information) of a signal.
db=(bounds(2)-bounds(1))/bincount;
edges=[bounds(1):db:bounds(2)];
n=histc(x,edges);
p=n./length(x);
p=p(find(p)); %get non-zero values
s=-sum(p.*log2(p));


function h = conditionalInfo(x,y,bounds,bincount)
%conditional information between two signals 
db=(bounds(2)-bounds(1))/bincount;
edges=[bounds(1):db:bounds(2)];
numedges=length(edges);
numVal = length(x);   
h = 0; 
for i = 2:numedges %for each different value found in x 
    iOcur = find(edges(i-1) <= x < edges(i));   
    numOcur = numel(iOcur); 
    px = numOcur/numVal;
    for j = 2:numedges   %for each different value found in y
        jOcur = numel(find(edges(i-1) <= y(iOcur) <edges(i)));
        if jOcur > 0
            py_x = jOcur/numOcur;  %conditional entropy for specific values of y and given x
            h = h - px*py_x*log2(py_x);   %evaluate and sum conditional entropies
        end
    end
end


function [x]=mvavg(x,taps,flt)
lenx=length(x);
fltr=ones(1,taps)/taps; 
x1=x(1); x2=x(lenx); 
for j=1:flt %smooth out quantization errors
	c=conv(fltr,x);
	x=c(2:lenx+1);
	x(1)=x1;  
    x(lenx)=x2; 
end