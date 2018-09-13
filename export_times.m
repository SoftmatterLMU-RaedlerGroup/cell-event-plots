% Exports the event times to a matfile to be read with other software.
%
% Instructions:
% * Call `loadCorrelationData` and load the conditions and combinations
%   from which you want to export the event times
% * Call this script to create a matfile in the "out" directory
% * Call the python script "convert_matfile" in the "py" directory,
%   giving the path to the created matfile as a string. This will generate
%   in the "out" directory a python-readable pickled `pandas.DataFrame`
%   with a similar structure as `exportT` in this script.
%
% Note: the last step requires a python3.6 installation with the modules
% "scipy", "numpy" and "pandas" available. It is up to you how to ensure
% that. The recommended way is to create a virtual python environment in
% the "py" directory that contains these modules.
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

% Load metadata, if not done yet
if ~exist('fileT', 'var')
	loadCorrelationData;
end

% Apply time correction
% (delay between beginning of measurement and drug addition)
% Prevent error upon multiple applications
lockString = 'Delay applied';
measRow = find(ismember(corrT.Properties.VariableNames, 't_event'));
if numel(corrT.Properties.VariableDescriptions) >= measRow && ...
		strcmpi(corrT.Properties.VariableDescriptions{measRow}, lockString)
	warning('Delay times have already been applied.')
else
	applyTimeDelays;
end
clear lockString measRow

% Set up table
exportT = table_template(0);

% Populate exportT
nCmb = height(corrT);
for iCmb = 1:nCmb
	% Get combination data
	markers = corrT.markers(iCmb,:);
	t_event = corrT.t_event{iCmb};
	nTraces = size(t_event, 1);

	% Prepare table for this combination
	ctab = table_template(nTraces);
	ctab.t_event = t_event;
	ctab.markers(:,:) = repmat(markers, nTraces, 1);

	% Get file map and check for consistence
	fm = corrT.file_map(iCmb,:);
	fm_check = find(any(~all(fm{1}(:,2:3) == fm{2}(:,2:3), 2)));
	fm = fm{1};
	if numel(fm_check) > 0 || fm(end) ~= nTraces
		error('Cannot proceed. Corrupt corrT.');
	end
	clear fm_check

	% Insert file information to table
	for iF = 1:size(fm, 1)
		file_idx = fm(iF, 1);
		tr_idx_start = fm(iF, 2);
		tr_idx_end = fm(iF, 3);

		file_tab = fileT(fileT.index_f == file_idx,:);
		ctab.condition(tr_idx_start:tr_idx_end) = file_tab.condition;
		ctab.measurement(tr_idx_start:tr_idx_end) = file_tab.measurement;
		ctab.position(tr_idx_start:tr_idx_end) = file_tab.position;

		ctab.trace(tr_idx_start:tr_idx_end) = 1:(tr_idx_end - tr_idx_start + 1);
	end

	% Append this combination to global table
	exportT = [exportT; ctab];
end
clear ctab file_idx file_tab fm iCmb markers nTraces
clear t_event tr_idx_start tr_idx_end

% Insert global trace id
exportT.id(:) = 1:height(exportT);

% Calculate delay of onset times
exportT.delay = exportT.t_event(:,2) - exportT.t_event(:,1);
exportT.delay(~isfinite(exportT.delay)) = NaN;

% Save table columns separately for compatibility with python
id = exportT.id;
condition = cellstr(exportT.condition);
measurement = cellstr(exportT.measurement);
position = cellstr(exportT.position);
trace = exportT.trace;
markers = cellstr(exportT.markers);
t_event = exportT.t_event;
delay = exportT.delay;

% Export table to matfile version 7
% (because newest version 7.3 cannot be read with scipy)
matfilename = fullfile(out_dir, [getTime, '_event_times.mat']);
save(matfilename, 'exportT', 'id', 'condition', 'measurement', ...
	'position', 'trace', 'markers', 't_event', 'delay', '-v7');
fprintf('Saved event times to: %s\n', matfilename)

% Table template function
function tab = table_template(nRows)
col_u32 = zeros(nRows, 1, 'uint32');
col_f64 = zeros(nRows, 1);
col_str = strings(nRows, 1);

tab = table( ...
	col_u32, ... cell id
	col_str, ... condition
	col_str, ... measurement
	col_str, ... position
	col_u32, ... trace index in files
	[col_str col_str], ... markers
	[col_f64 col_f64], ... event times
	col_f64, ... delay
	'VariableNames', ...
	{'id', 'condition', 'measurement', 'position', 'trace', 'markers', 't_event', 'delay'} ...
	);
end