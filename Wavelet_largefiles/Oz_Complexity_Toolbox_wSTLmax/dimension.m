function [Dc,DGP,DT,indGP,indT,tau] = dimension(x,delay)
%
%Fast correlation dimension calculation for a time series. Adapted from
%algorithm developed by M. Hein and J-Y. Audibert (Proceedings of the 22nd 
%ICML, pp. 289-296, Eds. L. de Raedt and S. Wrobel, 2005).
%
%[Dc,DG,DT] = dimension(x,delay);
%
%x:     time series
%delay: specify value, or method of estimating embedding time delay. 
%       *See embed() function. Do not enter anything for default method*
%
%Dc:
%   -the first column gives the values of correlation dimension as
%    estimated from the Grassberger-Procaccia (GP) algorithm versus
%    embedding dimension
%   -the second column gives Takens' maximum likelihood estimate of
%    the correlation dimension versus embedding dimension
%   (the 'x-axis' is the embedding dimension, which corresponds
%    directly to the row index)
%
%DG: Estimated correlation dimension according to GP
%DT: Estimated correlation dimension according to Takens' method
%    *DG/DT returns NaN (not a number) if no suitable estimate is found*
%
%
%Osbert Zalay, March 2008  
%
%*** For academic/research use. Please acknowledge author. ***

if nargin < 2
    delay='m3';
end

[len,m]=size(x);
if m > len
    x=x.';
    len=m;
end

%Eckmann-Ruelle limits:

%maximum correlation dimension for a time series of length 'len'
rho=0.1;
Dc_max=2*log(len)/log(1/rho);

%minimum time series size needed to estimate Lmax: N = exp(Dc*log(1/rho))

MaxEmbDim=floor(2*Dc_max+1); %Maximum embedding dimension based on Dc_max

%Calculate correlation dimension for different embedding dimensions
continu=1;
D=zeros(MaxEmbDim,3);
[X,tau]=embed(x,1,delay,0,1);
D(1,:)=corrdim(X);
fprintf('\n Processing dimension...');
if MaxEmbDim > 1
    for i=2:MaxEmbDim
        fprintf('.');
        X=embed(x,i,tau);
        D(i,:)=corrdim(X);
        if (D(i,2)==0 | D(i,2) > 10*Dc_max)
            D(i,2)=D(i,1);
        end
        if (D(i,3)==0 | D(i,3) > 10*Dc_max)
            D(i,3)=D(i,1);
        end
        %if exceeds Dc_max, terminate
        if (D(i,2)>Dc_max | D(i,3)>Dc_max)
            if i>1
                D=D(1:(i-1),:);
            end
            continu=0;
            break
        end
    end
end
Dc=D(:,2:3);
if numel(Dc)<=2
    Dc=NaN(2,2);
end

%Estimate correlation dimension from saturation region of Dc
%provided termination conditions are satisfied
slopethr0=0.2;
slopethr=slopethr0;
minSaturate=ceil(0.1*MaxEmbDim); %define minimum plateau size

%find the plateau
dDc=abs(diff(Dc));
indGP=find(dDc(:,1) < slopethr); 
indT=find(dDc(:,2) < slopethr);
if length(indGP) > 1
    lenGP=indGP(end)-indGP(1)+1;
else
    lenGP=0;
end
if length(indT) > 1
    lenT=indT(end)-indT(1)+1;
else
    lenT=0;
end

%keep going if no plateau is found, and within Eckmann-Ruelle limits
while ((lenGP < minSaturate || lenT < minSaturate)...
        && (MaxEmbDim <= 50) && continu)
    D=repmat(D,2,1);
    newMaxEmbDim=MaxEmbDim*2;
    for i=(MaxEmbDim+1):newMaxEmbDim
        fprintf('.');
        X=embed(x,i,tau);
        D(i,:)=corrdim(X);
        if (D(i,2)==0 | D(i,2) > 10*Dc_max)
            D(i,2)=D(i,1);
        end
        if (D(i,3)==0 | D(i,3) > 10*Dc_max)
            D(i,3)=D(i,1);
        end
        if (D(i,2)>Dc_max | D(i,3)>Dc_max)
            if i>1
                D=D(1:(i-1),:);
            end
            continu=0;
            break
        end
    end
    Dc=D(:,2:3);
    dDc=abs(diff(Dc));
    indGP=find(dDc(:,1) < slopethr);
    indT=find(dDc(:,2) < slopethr);
    if length(indGP) > 1
        lenGP=indGP(end)-indGP(1)+1;
    else
        lenGP=0;
    end
    if length(indT) > 1
        lenT=indT(end)-indT(1)+1;
    else
        lenT=0;
    end
    MaxEmbDim=newMaxEmbDim;
end

if ~lenGP | ~lenT
    dDc=abs(diff(Dc));
    slopethr=slopethr0+0.1;
    while (~lenGP & slopethr <= 0.5)
        indGP=find(dDc(:,1) < slopethr);
        slopethr=slopethr+0.1;
        if length(indGP) > 1
            lenGP=indGP(end)-indGP(1)+1;
        else
            lenGP=0;
        end
    end
    slopethr=slopethr0+0.1;
    while (~lenT & slopethr <= 0.5)
        indT=find(dDc(:,2) < slopethr);
        slopethr=slopethr+0.1;
        if length(indT) > 1
            lenT=indT(end)-indT(1)+1;
        else
            lenT=0;
        end
    end
end

%Compute plateau estimates
if ~isempty(indGP)
    DGP=D(indGP(1):(indGP(end)+1),2);
    DGP=DGP(find(~isnan(DGP)));
    DGP=DGP(find(DGP > 1e-10));
    DGP=mean(DGP);
    indGP=[indGP(1) (indGP(end)+1)];
else
    DGP=NaN;
end
if ~isempty(indT) 
    DT=D(indT(1):(indT(end)+1),3);
    DT=DT(find(~isnan(DT)));
    DT=DT(find(DT > 1e-10));
    DT=mean(DT);
    indT=[indT(1) (indT(end)+1)];
else
    DT=NaN;
end

fprintf('Done. \n\n');