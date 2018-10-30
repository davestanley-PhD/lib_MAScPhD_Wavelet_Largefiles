% Modification of parameters
% 1) Get the default parameter struct using p = lmax_getConfig()
% 2) Modify parameters in struct
% 3) Run STLmax with modified parameter struct as second parameter

function [L] = lmax(x,params)
if nargin == 1
    p = lmax_getConfig();
    L = lmaxw(x,p.DIMM,p.tau,p.dt,p.bb,p.cc,p.mult,p.angmx1,p.fiduc1);
else
    if nargin == 2
        p = params;
        L = lmaxw(x,p.DIMM,p.tau,p.dt,p.bb,p.cc,p.mult,p.angmx1,p.fiduc1);
    else
        fprintf('Usage: L = lmax(x,params)\n');
        fprintf('x - is row vector of input data\n');
        fprintf('params - a struct of input parameters that can\n');
        fprintf('be accessed using params = lmax_getConfig()\n');
    end
end

