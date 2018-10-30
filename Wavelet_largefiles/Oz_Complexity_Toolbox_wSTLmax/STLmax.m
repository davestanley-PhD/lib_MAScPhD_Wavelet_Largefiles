function L = STLmax(x,params)
%STLmax algorithm by tne Iasemidis Group; see L.D., Sackellares, J.C., 
%Zaveri, %H.P., Williams, W.J., Brain Topogr. 2, 187–201, 1990.
%
%L = STLmax(x,params);
%
%L = short-time estimate of maximum lyapunov exponent
%
%x = time series row vector
%params = struct of parameters for algorithm; 
%   Usage: type 'params = STLmax_getConfig()'; modify parameters as needed
%   and call function as displayed above.
%
%Note: function returns Not-A-Number (NaN) if time series segment (or 
%      window size if using in the context of wcomplexity.m) is not 
%      sufficiently long for the embedding and time-delay chosen.

lenx=numel(x); rho=0.1; mxd=2*log(lenx)/log(1/rho); 
if nargin == 1
    params = STLmax_getConfig();
elseif nargin == 2
    if ~isstruct(params)
       L=NaN;
       fprintf('\n Parameters must be given as a struct. Exiting... \n\n');
       return;
    end 
end
p=params;
%error-checking: if delay-reconstruction exceeds time-series length, exit
%if (p.DIMM-1)*p.tau > lenx
if (p.tau*(3*p.DIMM-1) >= lenx) | (p.DIMM > ceil(2*mxd+2))
    L=NaN;
    return;
end
L = lmaxw(x,p.DIMM,p.tau,p.dt,p.bb,p.cc,p.mult,p.angmx1,p.fiduc1);


