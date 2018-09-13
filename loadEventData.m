function [ eventData, fitData, traceMap ] = loadEventData( fileT, fileIndices, colIndices )
%LOADEVENTDATA loads event data from state files
%
% Input parameters:
%	fileT			file information table as returned by `getFiles`
%	fileIndices		vector of indices of files in `fileT` to be loaded
%	colIndices		vector of indices of columns of files to be loaded
%					(leave `colIndices` empty to load all columns)
%
% Output parameters:
%	eventData		n x m array of data
%						n: number of traces from all files
%						m: number of columns as indicated by `colIndices`
%	fitData			n x p array of data
%						n: number of traces from all files
%						p: number of fit parameters
%	traceMap		n x 3 array of indices
%						n: number of files given by `fileIndices`
%						1st column: index of file in `fileT`
%						2nd column: index of first trace of file in `eventData`
%						3rd column: index of last trace of file in `eventData`
%
% Copyright © 2018 Daniel Woschée <daniel.woschee@physik.lmu.de>
% Faculty of Physics / Ludwig-Maximilians-Universität München
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

%% Test input values
if nargin < 3 || ~isnumeric(colIndices)
	colIndices = [];
else
	if ~isvector(colIndices)
		colIndices = colIndices(:);
	end
	if any(colIndices < 1) || any(mod(colIndices, 1) ~= 0)
		error('`colIndices` must be positive integers.')
	end
end

%% Allocate arrays
eventData = [];
fitData = [];
traceMap = uint32.empty(0,3);

%% Populate event data array
% Get files
idx_paths = ismember(fileT.index_f, fileIndices);
state_paths = fileT.path_state(idx_paths);
params_paths = fileT.path_params(idx_paths);
idx_paths = fileT.index_f(idx_paths);

for iPath = 1:numel(idx_paths)
	% Test path
	if ismissing(state_paths(iPath))
		warning('Ignore cell with missing STATE file at index %d', ...
			idx_paths(iPath));
		continue
	end
	if ismissing(params_paths(iPath))
		warning('Ignore cell with missing PARAMS file at index %d', ...
			idx_paths(iPath));
		continue
	end

	% Read file content
	state_tab = csvread(state_paths{iPath});
	params_tab = csvread(params_paths{iPath});

	% Ensure data consistency
	if size(state_tab, 1) ~= size(params_tab, 1)
		warning('Inconsistent sizes for STATE and PARAMS in file %d', ...
			idx_paths(iPath));
		continue
	end

	% Get columns to be copied
	if isempty(state_tab) || isempty(params_tab)
		continue
	elseif isempty(colIndices)
		idcs = 1:size(state_tab, 2);
	elseif size(state_tab, 2) < max(colIndices)
		error('The file %s has less columns than necessary.', state_paths{iPath})
	else
		idcs = colIndices;
	end

	% Get first and last index of new traces
	idx_start = size(eventData, 1) + 1;
	idx_end = idx_start - 1 + size(state_tab, 1);

	% Append file content to data arrays
	eventData = [ eventData; state_tab(:,idcs) ];
	fitData = [ fitData; params_tab ];
	traceMap(size(traceMap,1)+1,:) = [ idx_paths(iPath) idx_start idx_end ];
end

end
