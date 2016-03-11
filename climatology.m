function varargout = climatology(time, data, varargin)
%CLIMATOLOGY Creates climatology from timeseries
%
% [timec, datac] = climatology(time, data)
% [timec, datac] = climatology(time, data, tbin)
% [timec, datac] = climatology(time, data, 'noleap')
% [timec, datac] = climatology(time, data, 'expand')
% [timec, datac, day] = climatology(...)
%
% This function creates a climatology by calculating the average yearly
% cycle of data. 
%
% Input variables:
%
%   time:       time data, in any Matlab date format (datenumbers, date
%               vectors, or date strings)
%
%   data:       nt x m array of data timeseries values.
%
%   tbin:       bin edges (similar to histc) used to group time data.
%               Recommended if data is at all irregular (i.e. data not at
%               exact same month/day/time each year).
%                   
%   'noleap':   if included, indicates that time data is in days since a
%               pivot year, with no leap year (i.e. all 365-day years).
%
%   'expand':   if included, the climatology is expanded to all input time
%               points, rather than just including the one-year cycle 
%
%   fun:        Function to apply to consolidated data.  Can be either a
%               single function handle or a cell array of function handles
%               [@nanmean]
%
% Output variables:
%
%   timec:      climatology time values, in either date numbers or days
%               since pivot year, no leap, depending on input format. 
%
%   datac:      climatology data, one per function applied to data
%
%   day:        day of year corresponding to each value of input time

% Copyright 2010-2013 Kelly Kearney

% Determine input and output options

noleap = any(strcmp('noleap', varargin));
expand = any(strcmp('expand', varargin));

isbin = cellfun(@isnumeric, varargin);
if any(isbin)
    tbin = varargin{isbin};
    binflag = true;
else
    binflag = false;
end

isfun = cellfun(@(x) isa(x, 'function_handle'), varargin);

if isvector(time)
    time = time(:);
end
if isvector(data)
    data = data(:);
end

if any(isfun)
    fun = varargin(isfun);
    nfun = length(fun);
else
    fun = {@nanmean};
    nfun = 1;
end

    
% Calculate time aggregator corresponding to day of year for each data
% point

if noleap % time should already be in days since some reference
    day = mod(time, 365);
else
    if ~isnumeric(time) || size(time,2) ~= 6
        time = datevec(time);
    end 
    yr1 = min(time(:,1));
    day = time;
    day(:,1) = yr1;
    day = datenum(day) - datenum(yr1,1,1);
end


% % Climatology will try to match time sampling of input data.
% 
% if ~noleap && (ischar(time) || ~isvector(time))
%     time = datenum(time);
%     dt = mean(diff(time));
% end

% Consolidate data, based on either unique day-fractions or by bin

if binflag
    [blah, bin] = histc(day, tbin);
    for ifun = 1:nfun
        [binc, tmp] = consolidator(bin, data, fun{ifun});
        datac{ifun} = nan(length(tbin), size(data,2));
        datac{ifun}(binc,:) = tmp;
    end
    timec = tbin;
else
    for ifun = 1:nfun
        [timec, datac{ifun}] = consolidator(day, data, fun{ifun});
    end
end

% Concatenate climatological data

if isvector(data)
    datac = cat(2, datac{:});
else
    ndim = ndims(data);
    datac = cat(ndim+1, datac{:});
end

% Expand if necessary

if expand

    error('Haven''t updated this option yet');
    timec = time;
    
    if nfun > 1
        for id = 1:length(datac)
            datac{id} = datac{id}(ind,:);
        end
    else
        datac = datac(ind,:);
    end
end
    
% Assign output

out = {timec datac day};
varargout = out(1:nargout);

% if nargout == 2
%     varargout{1} = timec;
%     varargout{2} = datac;
% elseif nargout == (1 + nfun)
%     varargout = [timec datac];
% else
%     error('Wrong number of output arguments');
% end

% if nargout == 3
%     varargout{1} = lower;
% elseif nargout == 4
%     varargout{1} = lower;
%     varargout{2} = upper;
% end


    
    