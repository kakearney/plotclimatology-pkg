function plotclimatology(varargin)
%PLOTCLIMATOLOGY Plot climatology with bounds
%
% plotclimatology(t1, x1, t2, x2, ...);
% plotclimatology(..., dt);
% plotclimatology(..., funhandles)
% plotclimatology(..., cmap)
%
% Input variables:
%
%   t#:         time values
%
%   x#:         x values
%
%   dt:         time interval tolerance used in consolidation
%
%   funhandles: 1 x 3 cell array of functions handles for center line,
%               lower limit, and upper limit (default is mean, 25th
%               percentile, and 75th percentile)
%
%   cmap:       n x 3 array, colormap to use for line colors
%
% Copyright 2010 Kelly Kearney

% Parse input

isfun = cellfun(@(x) iscell(x) && length(x)==3, varargin);
                 
if ~any(isfun)
    fun = {@mean, ...
           @(x) prctile(x,25), ...
           @(x) prctile(x,75)};
else
    fun = varargin{isfun};
    isfh = cellfun(@(x) isa(x, 'function_handle'), fun);
    if ~all(isfh)
        error('Unrecognized input');
    end
end

isdt = cellfun(@(x) isnumeric(x) && isscalar(x), varargin);
if any(isdt)
    dt = varargin{isdt};
else
    dt = [];
end

iscmap = cellfun(@(x) isnumeric(x) && size(x,2)==3 && all(x(:)>=0 & x(:)<=1), varargin);
if any(iscmap)
    cmap = varargin{iscmap};
    modifycol = true;
else
    modifycol = false;
end
    
tx = varargin(~isfun & ~isdt & ~iscmap);

t = tx(1:2:end);
x = tx(2:2:end);

nx = length(x);


if ~all(cellfun(@isvector, x))
    error('all x must be vectors');
end

% Calculate climatologies

for ix = 1:nx
    [tclim{ix}, xclim{ix}, xlo, xhi] = climatology(t{ix}, x{ix}, fun{:}, dt);
    rng{ix} = [xclim{ix}-xlo xhi-xclim{ix}];
end

dv = cellfun(@datevec, tclim, 'uni', 0);
yr = cellfun(@(x) unique(x(:,1)), dv);
yr = min(yr);
for ix = 1:nx
    dv{ix}(:,1) = yr;
    tclim{ix} = datenum(dv{ix});
end

% Check for NaNs; these will mess up the plotting

for ii = 1:nx
    isbad = isnan(xclim{ii}) | any(isnan(rng{ii}),2);
    tclim{ii} = tclim{ii}(~isbad);
    xclim{ii} = xclim{ii}(~isbad);
    rng{ii} = rng{ii}(~isbad);
end


% Plot

% col = {'b', 'g', 'r', 'c', 'm', 'y', 'k'};
% col = repmat(col, ceil(nx/length(col)), 1);
% col = col(1:length(tclim));

data = [tclim; xclim; rng];

if ~modifycol
    cmap = get(gca, 'colororder');
end

[hl, hp] = boundedline(data{:}, 'cmap', cmap);

    
