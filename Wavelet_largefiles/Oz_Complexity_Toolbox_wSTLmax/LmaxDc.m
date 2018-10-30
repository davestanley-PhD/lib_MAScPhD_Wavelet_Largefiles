function [Lmax,Dc,lnDivergence,CorrSum,yL,yD,dg] = LmaxDc(x,m,tau,SR,dim,divergefrc,W,flt,prnt)
%
%Short time maximum lyapunov exponent estimator (Lmax) for a time series. 
%Optional computation of correlation dimension (Dc). Derived from algorithm 
%developed by Michael T. Rosenstein et al. (Physica D, 65: 117-134, 1993)
%
%[Lmax,Dc,lnD,CI,yL,yD] = LmaxDc(x,emb,tau,SR,dim,diverge,W,flt);
%
%x: time series (minimum length ~500; recommended length: 1000-10000)
%emb: embedding dimension
%tau: time delay for phase space reconstruction   
%SR: sampling rate (optional); Lmax ~ SR*slope<ln(divergence)>; default: 1 
%dim: calculate correlation dimension (optional); default: 0 (no); 1 (yes)
%diverge: fraction of time series length to compute <ln(divergence)> 
%         (optional);  range: >0 to 1; default: set by internal algorithm
%W: Theiler window (optional); default: zero of autocorrelation
%flt: 3-point moving average filter for <ln(divergence)> (optional); 
%     default=1; use higher integer value for noisier datasets
%
%Lmax: maximum lyapunov exponent estimate (base: natural logarithm); if SR
%      not specified, Lmax is dimensionless and expressed per unit sample
%Dc:   Correlation dimension estimate (Grassberger-Procaccia method)
%lnD:  <ln(divergence)> (column 2) vs. divergence time (column 1)
%CI:   correlation integral (column 2) vs. ln(r) (column 1)
%yL:   line of best fit (least squares optimization) to lnD
%yD:   line of best fit to CI
%
%
%Osbert Zalay, March 2008  
%
%*** For academic/research use. Please acknowledge author. ***


%input argument checking
if nargin<4
    SR=1;
    dim=0;
    divergefrc=0;
    W=0;
    flt=1;
    prnt=1;
end
if nargin<5
    dim=0;
    divergefrc=0;
    W=0;
    flt=1;
    prnt=1;
end
if nargin<6
    divergefrc=0;
    W=0;
    flt=1;
    prnt=1;    
end
if nargin<7
    W=0;
    flt=1;
    prnt=1;
end
if nargin<8
    flt=1;
    prnt=1;
end
if nargin<9
    prnt=1;
end
if SR==0
    SR=1;
end

%initialize and allocate memory
femto=1e-15;
maxLnR=12;
minLnR=-12;
nLnR=600;
neighborIndex=1; 
maxIndex=0;
i=0; 
j=0;
k=0; 
T=0; 
nPts=0;
r=0;
CSumIndex=0;
numPairsDone=0;
nVectors=0;
distnce=0;
d=0;
ds=0; 
kNorm=0;
expn1=exp(-1);

k1 = nLnR/(maxLnR-minLnR);
k1 = round(k1*0.5);
k2 = nLnR/2;
	
if ~W
    %set Theiler window
    W = TheilerTime(x);
end

nPts = length(x);

%set divergence time
if ~divergefrc
    if nPts <= 500
        divergefrc=1;
    else
        divergefrc=max(500/nPts,0.1);
    end
end
divergeT = ceil(divergefrc*nPts);

nVectors = nPts-tau*(m-1);	
isNeighbor = zeros(nVectors, 1);
numDivergence = zeros(divergeT,1); 
lnDivergence=numDivergence;
CorrSum=zeros(nLnR,1);

if (dim && prnt==1)
    fprintf('\n Processing Lmax and Dc...');
elseif prnt==1
    fprintf('\n Processing Lmax...');
else
end

if (nVectors < 10 | (W >= nVectors))
    if prnt<2
        fprintf('\n\n *** Warning: switching to sub-optimal delay and Theiler window. \n'); 
        fprintf('     Longer time series needed for optimal results. \n\n');
    end
    while nVectors < 10
        tau=tau-1;
        nVectors = nPts-tau*(m-1);
    end
    while (W > 0.9*nVectors && W > tau)
        W=W/2;
    end
    W=round(W);
    isNeighbor = zeros(nVectors, 1);
end

% iterate through phase space vectors
i = 1;
cnt=1;

if dim
    
while(i<=nVectors)
   if(~isNeighbor(i))
      distnce = 1e10;  
      % find the nearest neighbor for the vector at i,
      % and ignore neighbors on same orbit
      if(i>W)
          for j=1:(i-W)
          % calculate distance squared
              d=0;
              for k=1:m
                  T = (k-1)*tau;
                  ds = x(i+T)-x(j+T);
                  ds = ds^2;
                  d = d + ds;
              end
              d = d + femto;
              % map squared distance to series index
              CSumIndex = round(k1*log(d)+k2);
              if(CSumIndex<=0)
                CSumIndex = 1;
              end
              if(CSumIndex>=nLnR)
                CSumIndex = nLnR-1;
              end
              CorrSum(CSumIndex)=CorrSum(CSumIndex)+1;
              if (d<distnce)
                  distnce=d;
                  neighborIndex=j;
              end
          end
      end
      if (i<(nVectors-W))
          for j=(i+W):nVectors
              d=0;
              for k=1:m
                  T = (k-1)*tau;
                  ds = x(i+T)-x(j+T);
                  ds = ds^2;
                  d = d + ds;
              end
              d = d + femto;
              CSumIndex = round(k1*log(d)+k2);
              if(CSumIndex<=0)
                CSumIndex = 1;
              end
              if(CSumIndex>=nLnR)
                CSumIndex = nLnR-1;
              end
              CorrSum(CSumIndex)=CorrSum(CSumIndex)+1;
              if(d<distnce)
                  distnce=d;
                  neighborIndex=j;
              end
          end
      end
      isNeighbor(neighborIndex) = 1;
      
      % record divergence
      for j=1:(divergeT)
         maxIndex = nPts-m*tau-j-2;
         if ((i < maxIndex) && (neighborIndex < maxIndex))
             d=0;
             for k=1:m
                 T = (k-1)*tau+j;
                 ds = x(i+T)-x(neighborIndex+T);
                 ds = ds^2;
                 d = d + ds;
             end
             d = d + femto;
             numDivergence(j)=numDivergence(j)+1;
             ds = 0.5*log(d);
             lnDivergence(j) =  lnDivergence(j)+ds;
         end
      end
      numPairsDone=numPairsDone+1;
   end
   i=i+1;
   cnt=cnt+1;
   if (cnt == prnt*200)
       cnt=1;
       fprintf('.');
   end
end

else
    
while(i<=nVectors)  
   if(~isNeighbor(i))
      distnce = 1e10;  
      % find the nearest neighbor for the vector at i,
      % and ignore neighbors on same orbit
      if(i>W)
          for j=1:(i-W)
          % calculate distance squared
              d=0;
              for k=1:m
                  T = (k-1)*tau;
                  ds = x(i+T)-x(j+T);
                  ds = ds^2;
                  d = d + ds;
              end
              d = d + femto; 
              if (d<distnce)
                  distnce=d;
                  neighborIndex=j;
              end
          end
      end
      if (i<(nVectors-W))
          for j=(i+W):nVectors
              d=0;
              for k=1:m
                  T = (k-1)*tau;
                  ds = x(i+T)-x(j+T);
                  ds = ds^2;
                  d = d + ds;
              end
              d = d + femto;
              if(d<distnce)
                  distnce=d;
                  neighborIndex=j;
              end
          end
      end
      isNeighbor(neighborIndex) = 1;
      
      % record divergence
      for j=1:(divergeT)
         maxIndex = nPts-m*tau-j-2;
         if ((i < maxIndex) && (neighborIndex < maxIndex))
             d=0;
             for k=1:m
                 T = (k-1)*tau+j;
                 ds = x(i+T)-x(neighborIndex+T);
                 ds = ds^2;
                 d = d + ds;
             end
             d = d + femto;
             numDivergence(j)=numDivergence(j)+1;
             ds = 0.5*log(d);
             lnDivergence(j) =  lnDivergence(j)+ds;
         end
      end
      numPairsDone=numPairsDone+1;
   end
   i=i+1;
   cnt=cnt+1;
   if (cnt == prnt*200)
       cnt=1;
       fprintf('.');
   end
end

end
    
% Lyapunov sum average
for i=1:divergeT
  if(numDivergence(i)>0)
    lnDivergence(i) = lnDivergence(i)/numDivergence(i);
  end
end

if dim
    % calculate correlation integral
    for i=2:(nLnR-1)
        CorrSum(i) = CorrSum(i)+CorrSum(i-1);
    end
    kNorm = 1/CorrSum(nLnR-1);
    for i=1:(nLnR-1)
        CorrSum(i) = CorrSum(i)*kNorm;
    end
    for i=1:(nLnR-1)
        temp = CorrSum(i);
        if (temp<0.000045) || (temp>0.990050) 
            CorrSum(i) = 0;
        else
            CorrSum(i) = log(temp);
        end
    end
  
    % estimate correlation dimenstion
    lnR=find(CorrSum);
    CorrSum=CorrSum(lnR);
    lnR=k1.*log(lnR)+k2;
    dg=diff(CorrSum);
    dg=find(dg);
    lnR=lnR(dg);
    CorrSum=CorrSum(dg);
    c=polyfit(lnR,CorrSum,1);
    Dc=c(1);
    yD=polyval(c,lnR);
end

% estimate the maximum lyapunov exponent from scaling of log(divergence)
errthr=0.05;
k=1.1;
x=[1:divergeT]'/SR;
lnDivergence=mvavg(lnDivergence,3,flt*100); %smooth if necessary
dg=diff(lnDivergence); %determine slope information
dg1=dg(1);
ldg=numel(dg);
mindg=min(dg);
dg=abs(dg);
maxdg=max(dg);
meddg=median(dg);
imin=find(dg==maxdg); 
c=polyfit(x,lnDivergence,1);
y=polyval(c,x);
ind=divergeT;
se0=(sum((y-lnDivergence).^2)/ind).^0.5; %determine initial fit error
se=se0;
se_old=1e10;
i=1;
stp=0;
if (maxdg < 3*meddg) | (mindg == -maxdg) | (numel(imin) > 1) | dg1 < 0
    while (ind > 0) && (i > 0) && ~stp
        if  (se < se_old) && (se > errthr*se0)
            se_old=se;
            ind=divergeT-i;
            c=polyfit(x(1:ind),lnDivergence(1:ind),1);
            y=polyval(c,x(1:ind));
            se=(sum((y-lnDivergence(1:ind)).^2)/ind).^0.5;
            stepsz=ceil(k*i);
            while (stepsz > divergeT)
                k=0.9*k;
                stepsz=ceil(k*i);
            end
            i=stepsz;
        elseif (se > se_old) && (se > errthr*se0)
            se_old=se;
            i=i-1;
            ind=divergeT-i;
            c=polyfit(x(1:ind),lnDivergence(1:ind),1);
            y=polyval(c,x(1:ind));
            se=(sum((y-lnDivergence(1:ind)).^2)/ind).^0.5;
        else
            stp=1;
        end
    end
else
    i=divergeT;
end

%if there exists saturation region, fine tune estimate
%bypass initial transient, and find where lnDivergence levels off
errthr=0.1;
i=divergeT-i;
dg=dg/maxdg;
i1=0; i2=i; ibench=0;
limit=round(0.5*divergeT);
if (i < limit)
    if isempty(imin)
        imin=1;
    end
    %istop=find(dg<=(1-0.5*expn1^2)*dg(imin));
    istop=find(dg(imin:end)<=(1-0.5*expn1^2)*dg(imin))+imin-1;
    itrans=find(dg(imin:end)<0.1)+imin-1;
    if isempty(itrans)
        itrans=W;
    end
    if isempty(istop)
        istop=divergeT;
    end
    istop=istop(1);
    if imin > 1 && imin < ldg
        i1=imin;
        ibench=imin;
        if ibench < limit
            while (dg(ibench+1) < dg(ibench) && ibench < limit)
                ibench=ibench+1;
            end
        end
        if istop < ibench
            i1=istop;
            if dg(ibench) < 0.15
                if ibench < 1.2*W
                    istop=divergeT-1;
                    val=1e10;
                else
                    val=(1-dg(ibench))*(1-expn1)+dg(ibench);  
                end
            else
                istop=ibench;
                val=(1-expn1^2)*dg(istop);
            end
            while (dg(istop) > val && istop < limit)
                istop=istop+1;
            end
        end
        i2=istop;
    else
        ibench=imin;
        if ibench < limit
            while (dg(ibench+1) < dg(ibench) && ibench < limit)
                ibench=ibench+1;
            end
        end
        if (ibench < istop | ibench >= limit | istop > 0.2*divergeT)
            i1=1;
            istop=find(dg<=(1-expn1));
            istop=istop(1);
        elseif (istop < ibench && dg(ibench) < 0.05 && ibench < 1.2*W)
            i1=ibench;
            istop=divergeT;
        else
            i1=istop;
            istop=ibench;
            val=expn1*dg(istop);
            while (dg(istop) > val && istop < limit)
                istop=istop+1;
            end
        end
        i2=istop;
    end
    if i1==i2 
        i2=i1+W;
        if i2 > divergeT
            i2=divergeT;
        end
    end
    c=polyfit(x(i1:i2),lnDivergence(i1:i2),1);
elseif (i >= limit) && c(1) > 0
    if imin==1 & numel(imin)==1
         ibench=find(dg<=(1-0.5*expn1^2)*dg(imin));
         if isempty(ibench)
             ibench=divergeT;
         end
         ibench=ibench(1);
         while (dg(imin+1) < dg(imin) && imin < limit)
            imin=imin+1;
         end
         if ibench < imin
             i1=ibench;
             i2=ind;
             c=polyfit(x(i1:i2),lnDivergence(i1:i2),1);
         end      
    end
end            
Lmax=c(1);
yL=polyval(c,x);

%prepare outputs
lnDivergence=[x,lnDivergence];
if dim
    CorrSum=[lnR,CorrSum];
else
    clear CorrSum Dc yD;
    CorrSum=NaN;
    Dc=NaN;
    yD=NaN;
end

if prnt==1
    fprintf('Done. \n');
end


function [tau] = TheilerTime(x)
%estimate of window based on 'decorrelation time'
x=x-mean(x);
lx=length(x);
xx=xcorr(x,'coeff');
lxx=length(xx);
xxh=xx(lx:lxx);
lxxh=length(xxh);
i=1;
while (xxh(i) > 0) & (i < lxxh)
    i=i+1; %find first zero-crossing of auto-correlation
end
if (i >= lxxh) %if no zero-crossing, use (1/e) of autocorrelation
    %val=0.1*xxh(1);
    val=exp(-1)*xxh(1); 
    i=1;
    while (xxh(i) > val) & (i < lxxh)
        i=i+1;
    end
end
tau=i;


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
