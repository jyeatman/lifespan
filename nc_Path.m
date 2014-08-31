function path = nc_Path
% Returns path to "lifespan" code directory
% 
% path = nc_Path
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications

path = fileparts(which('nc_Path'));
