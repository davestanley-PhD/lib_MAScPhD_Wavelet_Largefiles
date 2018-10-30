function [d,l,ind] = wcomplexity(x,typ,cap,SR,diverge,delay,emb,LsampleSz,flt)
%
%Sliding window complexity analysis. Type 'help complexity' for more 
%information.
%
%[Dc,Lmax,ind]=wcomplexity(x,typ,cap,SR,diverge,delay,emb,numL,flt);
%
%cap: window size and sliding increment
%     format: cap = [increment size, window size]
%     e.g. cap = [50 1000] is increment 50 samples per slide, 1000 samples
%                for the window size (i.e. overlapping windows)
%          cap = 2000 is increment 2000 samples with 2000 sample window
%                (i.e. no overlap)
%          *** Make sure long enough segments are used to satisfy
%              Eckmann-Ruelle limits, otherwise estimates may be inaccurate
%
%Dc: correlation dimension estimate: returns Not a Number (NaN) if no 
%   saturation of embedding dimension vs. fractal dimension
%Lmax: maximum lyapunov estimate; returns NaN if not calculated 
%      (as specified by operating mode 'typ'), or due to error condition
%ind: time (sample) indices for denoting each window analyzed
%
%
%Osbert Zalay, April 2008
%
%*** For academic/research use. Please acknowledge author. ***

if nargin < 2
    typ='lgp';
    cap=[1000 1000];
    SR=1;
    diverge=0;
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 3
    cap=[1000 1000];
    SR=1;
    diverge=0;
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 4
    SR=1;
    diverge=0;
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 5
    diverge=0;
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 6
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 7
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 8
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 9
    flt=1;
end
if SR==0
    SR=1;
end
if length(cap)==1
    cap=[cap cap];
end

[len,m]=size(x);
if m > len
    x=x.';
    len=m;
end

dt=1/SR;

%set divergence time
if ~diverge
    if len <= 500
        diverge=1;
    else
        diverge=max(500/len,0.2);
    end
end

len=length(x);
numSegments=floor((len-cap(2))/cap(1));

l=zeros(numSegments,2);
d=zeros(numSegments,2);
ind=zeros(numSegments,1);
i1=0; i2=0; ii=0;
mdpt=round(cap(2)/2);

for j=0:numSegments
    i1=j*cap(1)+1;
    i2=i1+cap(2)-1;
    ii=j+1;
    clc;
    fprintf('\n\n *** Analyzing points %d to %d of %d (%d/%d): \n',i1,i2,len,ii,numSegments+1);
    [d(ii,:),l(ii,:)]= complexity(x(i1:i2),typ,SR,diverge,delay,emb,LsampleSz,flt);
    ind(ii)=dt*(mdpt+i1);
end
