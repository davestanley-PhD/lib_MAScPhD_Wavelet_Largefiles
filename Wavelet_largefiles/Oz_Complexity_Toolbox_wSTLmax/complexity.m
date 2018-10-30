function [Dc,Lmax,D,L,X]=complexity(x,typ,SR,diverge,delay,emb,LsampleSz,flt)
%
%Automated complexity analysis for a time series. Returns the maximum
%lyapunov exponent and correlation dimension estimates.
%
%[Dc,Lmax,D,L,X]=complexity(x,typ,SR,diverge,delay,emb,numL,flt)
%
%x: time series
%typ: operation mode; type any of the following ('...'):
%        gp = Grassberger Procaccia estimate of correlation dimension
%        takens = Takens' estimate of correlation dimension
%        l = maximum lyapunov estimate (Rosenstein Method)
%        stl = maximum lyapunov estimate (STLmax, Iasemidas et al.)
%        ! = fast mode (Rosenstein Method only; accuracy not guaranteed)
%        * = error checking on dimension estimate (uses slow GP algorithm
%            in case of poor convergence)
%        Modes can be run independently or in combination; To combine 
%        modes, just type string in series (default='lgp*'); user cannot 
%        combine both dimension measures or both lyapunov measures in one 
%        run, as only one version will be used to produce the estimate.         
%SR: sampling rate (optional)
%diverge: (optional for Rosenstein) type 'help LmaxDc' for information
%delay: (optional) embedding time delay; type 'help embed' for information
%emb: embedding dimension (optional)
%numL: specify embedding dimension range for calculating Lmax (optional); 
%      format: [numL(1) numL(2)]  OR  numL;  numL = -1 for full range;
%              with dimension calculation: 2*Dc+numL(1) to numL(2);
%              without:                    numL(1) to numL(2)  OR
%                                          2 to numL
%flt: smoothing filter for <ln(divergence)>; aids with noisy signals in
%     obtaining convergence for estimate of Lmax (only implemented for
%     Rosenstein algorithm) Range: positive integers (default=1)        
%Dc: correlation dimension estimate [mean +- std]
%Lmax: maximum lyapunov exponent estimate [mean +- std]
%D: individual Dc values computed over range of embedding dimensions
%L: individual values of Lmax over embedding range specified
%X: time-delay phase space reconstruction
%
%
%Osbert Zalay, April 2008
%
%*** For academic/research use. Please acknowledge author. ***


%input argument checking
if nargin < 2
    typ='lgp*';
    SR=1;
    diverge=0;
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 3
    SR=1;
    diverge=0;
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 4
    diverge=0;
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 5
    delay='m3';
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 6
    emb=0;
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 7
    LsampleSz=[0 0];
    flt=1;
end
if nargin < 8
    flt=1;
end
if SR==0
    SR=1;
end

[len,m]=size(x);
if m > len
    x=x.';
    len=m;
end

%set divergence time
if ~diverge
    if len <= 500
        diverge=1;
    else
        diverge=max(500/len,0.2);
    end
end

%input checking for mode of operation
typ=lower(typ);
lflg=0; dim=0; err=0;
s1=~isempty(findstr(typ,'l'));
s2=~isempty(findstr(typ,'gp'));
s3=~isempty(findstr(typ,'takens'));
s4=~isempty(findstr(typ,'*'));
s5=~isempty(findstr(typ,'!'));
s6=~isempty(findstr(typ,'st'));
if s1 && ~s5
    lflg=1;
end
if s1 && s5
    lflg=2;
end
if s6
    lflg=1.5;
end
if s2
    dim=1;
end
if s3
    dim=2;
end
if s4
    err=1;
end

%check range on embedding dimension for Lmax calculation
if length(LsampleSz)==1
    LsampleSz=[2 LsampleSz];
end

%initialize some variables
cap=1000; %size of segment for fast algorithm, and dimension recalculation 
Dc=NaN;
embflg=0;
lenL=0;
mdL=0;
Lm=0;
aL=0;
if lflg==1.5
    x=x.';
    p=STLmax_getConfig();
    p.dt=(1/SR);
end

%Eckmann-Ruelle limits:

%maximum correlation dimension for a time series of length 'len'
rho=0.1;
Dc_max=2*log(len)/log(1/rho);
limDc=floor(Dc_max-1);
X=[]; L=[];

%minimum time series size needed to estimate Lmax: N = exp(Dc*log(1/rho))

MaxEmbDim=ceil(2*Dc_max+1); %Maximum embedding dimension based on Dc_max
minSaturate=ceil(0.1*MaxEmbDim); %define minimum plateau size
slopethr=0.2; %threshold for maximum slope magnitude of plateau

if dim
if ~emb
    %Start with fast estimate correlation dimension
    [D,DG,DT,ig,it,tau]=dimension(x,delay);
    if dim==2
        D=D(:,2);
        ind=it;
        if ~isempty(ind)
            Dc=[DT std(D(ind(1):ind(2)))];
        else
            Dc=NaN;
        end
    elseif dim==1
        D=D(:,1);
        ind=ig;
        if ~isempty(ind)
            Dc=[DG std(D(ind(1):ind(2)))];
        else
            Dc=NaN;
        end
    end
    %If fast estimate fails, resort to alternative method
    if (Dc < 0.5 | (isnan(Dc) & D(end) < limDc) |...
            ~isempty(find(isnan(D)))) & err                
        fprintf(' *** Recalculating using slow GP algorithm...');
        if len >= cap
            numSegments=floor(len/cap);
            %Limit Dc_max (cap sample segments) for speed
            Dc_max=2*log(cap)/log(1/rho); 
            MaxEmbDim=ceil(2*Dc_max+1);
            xs=zeros(cap,numSegments);
            taus=zeros(numSegments,1);
            L=0;
            Ds=zeros(numSegments,1);
            Ls=Ds;
            D=zeros(MaxEmbDim,1);
            for j=0:(numSegments-1)
                xs(:,(j+1))=x((j*cap+1):(j+1)*cap);
                [Xs,taus(j+1)]=embed(xs(:,(j+1)),1);
            end
            for j=1:MaxEmbDim
                for i=1:numSegments
                    %compute segment correlation dimension
                    [Ls(i),Ds(i)]= LmaxDc(xs(:,i),j,taus(i),SR,1,0.5,0,flt,4);
                end
                %get the average over all segments for given embedding
                D(j)=mean(Ds);
                L(j)=mean(Ls);
                %fprintf('.');
            end
        else
            L=zeros(MaxEmbDim,1);
            D=zeros(MaxEmbDim,1);
            for i=1:MaxEmbDim
                %compute Lmax and correlation dimension
                [L(i),D(i)]= LmaxDc(x,i,tau,SR,1,diverge,0,flt,4);%
            end
        end
            
        %Obey Ruelle-Eckmann limit
        D=D(find(D <= Dc_max));
            
        %Estimate correlation dimension from saturation region of Dc
        %provided termination conditions are satisfied
        minSaturate=ceil(0.1*MaxEmbDim); %define minimum plateau size
            
        %find the plateau
        lenD=0;
        dD=abs(diff(D));
        ind=find(dD(:,1) < slopethr); 
        if ~isempty(ind)
            lenD=ind(end)-ind(1)+1;
        end
        if (lenD > minSaturate)
            %compute the correlation dimension estimate
            Dc=[mean(D(ind(1):ind(end))) std(D(ind(1):ind(end)))];
        else
            Dc=NaN;
        end
        fprintf('Done. \n\n');
    end
else
    %For user-specified embedding, find correlation dimension
    [X,tau]=embed(x,emb,delay,0,1); embflg=1;
    fprintf('\n Processing dimension...');
    D=corrdim(X);
    if dim==2
        Dc=D(3);
    elseif dim==1
        Dc=D(2);
    end
    fprintf('Done. \n');
    if (isnan(D) & err)                
        fprintf(' *** Recalculating using slow GP algorithm...');
        %compute Lmax and correlation dimension
        [L,Dc]= LmaxDc(x,emb,tau,SR,1,diverge,0,flt,4);%
        fprintf('Done. \n');
    end
    D=Dc;
end
end

if err
    %reset MaxEmbDim
    Dc_max=2*log(len)/log(1/rho);
    MaxEmbDim=ceil(2*Dc_max+1);
end

if ~LsampleSz
    %Set automatic sample size if not user-specified
    if dim
        if (lflg < 2 | len > 2000)
            %LsampleSz = [2 (ceil(log(len)/log(1/rho))+1)];
            LsampleSz = [1 (ceil(log(len)/log(1/rho))+1)];
        elseif (lflg == 2 && len <= 2000)
            %LsampleSz = [2 MaxEmbDim];
            LsampleSz = [1 MaxEmbDim];
        end
    else
        LsampleSz = [2 MaxEmbDim];
    end
end

%Set minimum embedding dimension
if dim
    MinEmbDim=ceil(2*Dc(1)+LsampleSz(1));
    if isnan(MinEmbDim) 
        MinEmbDim = 7;
    end
    if (MinEmbDim < 1)
        MinEmbDim = 1;
    end
else
    MinEmbDim=LsampleSz(1);
end

%Use EndDim = MaxEmbDim if user specifies this
if LsampleSz(2)==-1
    if dim
        LsampleSz(2)=MaxEmbDim-MinEmbDim+1;
    else
        LsampleSz(2)=MaxEmbDim;
    end
end

%Process the maximum lyapunov exponent
if (lflg==1 | lflg==1.5) %regular mode
    if ~emb
        if dim      
            %if Dc is known
            fprintf(' Processing Lmax...');
            EndDim=MinEmbDim+LsampleSz(2)-1;
            if (MinEmbDim >= EndDim)
                %check dimension precedence
                MinEmbDim=max(1,round(EndDim/2));
            end        
            L=NaN(EndDim,1);
            for i=MinEmbDim:EndDim
                %compute Lmax over selected range of embedding dimensions
                if lflg~=1.5
                    L(i)= LmaxDc(x,i,tau,SR,0,diverge,0,flt,4);%
                else
                    p.tau=tau; p.DIMM=i;
                    L(i)=STLmax(x,p);
                end
            end
            L=L(find(~isnan(L)));
            if length(L)>1
                %Lm=min(abs(L(2:end)));
                %L=L(find(L > -100*Lm & L < 100*Lm));
                Lmax=[mean(L) std(L)];
            else
                Lmax=[L 0];
            end
            fprintf('Done. \n');
        else
            %if Dc not known, estimate Lmax by convergence plateau
            if isnumeric(delay)
                tau=delay;
            else
                [X,tau]=embed(x,3,delay,0,1);
            end
            fprintf('\n Processing Lmax...');
            EndDim=LsampleSz(2);
            if (MinEmbDim >= EndDim)
                %check dimension precedence
                MinEmbDim=max(1,round(EndDim/2));
            end        
            L=NaN(EndDim,1);
            for i=MinEmbDim:EndDim
                if lflg~=1.5 
                    L(i)= LmaxDc(x,i,tau,SR,0,diverge,0,flt,4);%
                else
                    p.tau=tau; p.DIMM=i;
                    L(i)=STLmax(x,p);
                end
                %fprintf('.');
            end
            L=L(find(~isnan(L)));
            if length(L)>1
                aL=abs(L(2:end));
                Lm=min(aL);
                mdL=mean(aL);
                %L=L(find(L > -100*Lm & L < 100*Lm));
                lenL=length(L);
            else
                mdL=abs(L);
                lenL=1;
            end
            ind=[];
            if lenL > 1
                dL=abs(diff(L));
                ind=find(dL < slopethr*mdL);
                if isempty(ind)
                    dind=0;
                else
                    dind=ind(end)-ind(1)+1;
                end
                while (dind < minSaturate) && (dind < lenL)...
                        && (slopethr <= 1)
                    clear ind; ind=[];
                    slopethr=slopethr+0.1;
                    ind=find(dL < slopethr*mdL);
                    if isempty(ind)
                        dind=0;
                    else
                        dind=ind(end)-ind(1)+1;
                    end
                end
            else
                Lmax=[L 0];
            end
            if ~isempty(ind) && lenL
                Lmax=[mean(L(ind)) std(L(ind))];
                %Lmax=[mean(L(ind(1):ind(end))) std(L(ind(1):ind(end)))];
            end
            fprintf('Done. \n');
        end
    else
        %single case Lmax for user-specified embedding
        if lflg==1.5
            fprintf('\n Processing Lmax...');
        end
        if isnumeric(delay)
            tau=delay;
        elseif ~embflg
            [X,tau]=embed(x,3,delay,0,1);
        end
        if lflg~=1.5
            L=LmaxDc(x,emb,tau,SR,0,diverge,0,flt);
        else
            p.tau=tau; p.DIMM=emb;
            L=STLmax(x,p);
        end        
        Lmax=L;
        if lflg==1.5
            fprintf('Done. \n');
        end
    end
elseif lflg==2 %fast mode
    if ~emb
        Dc_max=2*log(cap)/log(1/rho); 
        MaxEmbDim=ceil(2*Dc_max+1);
        if dim
            fprintf(' Processing Lmax...');
            EndDim=MinEmbDim+LsampleSz(2)-1;
        else         
            if isnumeric(delay)
                tau=delay;
            elseif ~embflg
                [X,tau]=embed(x,3,delay,0,1);
            end
            fprintf('\n Processing Lmax...');
            EndDim=LsampleSz(2);
        end
        if (MinEmbDim >= EndDim)
            %check dimension precedence
            MinEmbDim=max(1,round(EndDim/2));
        end        
        if length(L)<2
            L=NaN(EndDim,1);
            if len >= cap
                numSegments=floor(len/cap);
                xs=zeros(cap,numSegments);
                taus=zeros(numSegments,1);
                Ls=zeros(numSegments,1);
                for j=0:(numSegments-1)
                    xs(:,(j+1))=x((j*cap+1):(j+1)*cap);
                    [Xs,taus(j+1)]=embed(xs(:,(j+1)),1);
                end
                for j=MinEmbDim:EndDim
                    for i=1:numSegments
                        %compute segment Lmax
                        Ls(i)= LmaxDc(xs(:,i),j,taus(i),SR,0,0.5,0,flt,4);
                    end
                    %get the average over all segments for given embedding
                L(j)=mean(Ls);
                %fprintf('.');
                end
            else
                for i=MinEmbDim:EndDim
                    %compute Lmax
                    L(i)= LmaxDc(x,i,tau,SR,0,diverge,0,flt,4);%
                end
            end
            L=L(find(~isnan(L)));
        else
            L=L(MinEmbDim:EndDim);
        end
        if length(L)>1
            aL=abs(L(2:end));
            Lm=min(aL);
            mdL=mean(aL);
            %L=L(find(L > -100*Lm & L < 100*Lm));
            lenL=length(L);
        else
            mdL=abs(L);
            lenL=1;
        end
        if dim 
            if lenL>1
                Lmax=[mean(L) std(L)];
            else
                Lmax=[L 0];
            end
        else
            ind=[];
            if lenL > 1
                dL=abs(diff(L));
                ind=find(dL < slopethr*mdL);
                if isempty(ind)
                    dind=0;
                else
                    dind=ind(end)-ind(1)+1;
                end
                while (dind < minSaturate) && (dind < lenL)...
                        && (slopethr <= 1)
                    clear ind; ind=[];
                    slopethr=slopethr+0.1;
                    ind=find(dL < slopethr*mdL);
                    if isempty(ind)
                        dind=0;
                    else
                        dind=ind(end)-ind(1)+1;
                    end
                end
            else
                Lmax=[L 0];
            end
            if ~isempty(ind) && lenL
                Lmax=[mean(L(ind)) std(L(ind))];
                %Lmax=[mean(L(ind(1):ind(end))) std(L(ind(1):ind(end)))];
            end
        end
        fprintf('Done. \n');
    else
        %single case Lmax for user-specified embedding
        fprintf('\n Processing Lmax...');
        if isnumeric(delay)
            tau=delay;
        elseif ~embflg
            [X,tau]=embed(x,3,delay,0,1);
        end
        L=LmaxDc(x,emb,tau,SR,0,diverge,0,flt);
        Lmax=L;
        fprintf('Done. \n');
    end
end
    
%output checking
if ~lflg && dim
    Lmax=NaN;
    L=NaN;
elseif ~dim && lflg
    Dc=NaN;
    D=NaN;
end
if isempty(Lmax)
    Lmax=NaN;
end
if isempty(L)
    L=NaN;
end
if isempty(Dc)
    Dc=NaN;
end
if isempty(D)
    D=NaN;
end

if isempty(X)
    X=embed(x,3,tau,0,0);
end

fprintf('\n');