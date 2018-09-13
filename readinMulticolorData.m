function [ tracesArray, metaStruct, dataFiles ] = readinMulticolorData( fileStruct )
%readinMulticolorData reads given CSV files into memory.
%
% Parameters:
% ===========
%	fileStruct		struct array containing the fields:
%						path_state	path to file
%						marker		marker token (optional)
%						condition	condition token (optional)
%
% Returns:
% ========
%	tracesArray		array of traces (rows) and values (columns)
%	metaStruct		struct array containing the fields:
%						index_a		index of trace in tracesArray
%						data		numerical array
%						marker		marker token (may be empty)
%						condition	condition token (may be empty)
%						measurement	index of measurement (may be emtpty)
%						index_f		index of file containing trace (optional)
%						index_t		index of trace in file (optional)
%	dataFiles		cell array of file paths
%
% If three output parameters are required, dataFiles will contain the
% paths of the files from which the data in tracesArray and metaStruct come.
% metaStruct will then have the additional fields “index_f”, which
% indicates the index of the corresponding file path in dataFiles, and
% “index_t”, which indicates the index of the trace in the file.
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

%% Initialize output variables
tracesArray = [];
metaStruct = struct( ...
		'index_a', {}, ...
		'marker', {}, ...
		'condition', {}, ...
		'measurement', {}, ...
		'index_f', {}, ...
		'index_t', {} ...
	);
dataFiles = cell(0);

%% Populate output variables
for i_in = 1:length(fileStruct)

	% Get index and path of file
	index_f = length(dataFiles) + 1;
	path_f = fileStruct(i_in).path_state;
	measurement_token = getMeasurementToken(path_f);

	% Read in datafile
	data = csvread(path_f);
	n0_traces = uint32(size(tracesArray, 1));
	tracesArray = [ tracesArray; data ];

	% Write data to output struct
	for index_t = 1:size(data, 1)
		index_struct = length(metaStruct) + 1;
		metaStruct(index_struct).index_a = n0_traces + index_t;
		metaStruct(index_struct).marker = fileStruct(i_in).marker;
		metaStruct(index_struct).condition = fileStruct(i_in).condition;
		metaStruct(index_struct).measurement = measurement_token;
		metaStruct(index_struct).index_f = uint16(index_f);
		metaStruct(index_struct).index_t = uint32(index_t);
	end

	% Write path to path cell, if required
	if nargout > 1
		dataFiles{index_f} = path_f;
	end
end

%% Cleanup
if nargout == 2
	metaStruct = rmfield(metaStruct, {'index_f', 'index_t'});
end

end

function measurement = getMeasurementToken(filename)
%% This function extracts a measurement token out of the file path
%
% Parameter:
%	filename		path to the file
%
% Returns:
%	measurement		measurement token (part of filename)
%%

% Define default position of measurement token
token_position = 3;

% Split file name
[~, filename, ~] = fileparts(filename);
fn_split = textscan(filename, '%s', 'Delimiter','_');
fn_split = fn_split{1};

% Find measurement token
if length(fn_split) >= 3
	measurement = fn_split{token_position};
elseif length(fn_split) == 1
	measurement = fn_split{1};
else
	measurement = [];
end

end