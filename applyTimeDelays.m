%applyTimeDelays adds delay times to the event times in `corrT`
%
% The delay times are read by calling `timeDelays`.
% The description of `corrT.t_event` is set to "Delay applied" in order to
% prevent unintentional re-application of the delay times.
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


% Prevent multiple applications of delay table
lockString = 'Delay applied';
measRow = find(ismember(corrT.Properties.VariableNames, 't_event'));
if numel(corrT.Properties.VariableDescriptions) >= measRow && ...
		strcmpi(corrT.Properties.VariableDescriptions{measRow}, lockString)
	clear lockString measRow
	error('Delay table has been applied already.');
end

% Import delays
delayT = timeDelays;

% Loop over `corrT`
for iCmb = 1:height(corrT)
	idxFiles = corrT.file_map{iCmb,1}(:,1);
	idxTraces = corrT.file_map{iCmb,1}(:,2:3);
	nFiles = length(idxFiles);
	nTraces = max(idxTraces(:));
	measurements = repmat("", nFiles, 1);
	delays = zeros(nTraces, 1);

	% Get measurements for all file indices
	for iFile = 1:nFiles
		measurements(iFile) = fileT.measurement( ...
			ismember(fileT.index_f, idxFiles(iFile)) );
	end

	% Prepare delay vector
	[measurements, ~, idxMeas] = unique(measurements);
	for iMeas = 1:numel(measurements)

		% Find row in delay table corresponding to measurement
		rowDelay = find(ismember(delayT.Date, measurements{iMeas}));

		if isempty(rowDelay)
			warning('No delay time found for measurement: %s', ...
				measurements{iMeas});
			continue

		elseif numel(rowDelay) > 1
			warning(['Multiple occurrences of measurement "%s" found ' ...
				'in measurement table; using first'], measurements{iMeas});
			rowDelay = rowDelay(1);
		end

		% Write delays to delay vector
		rowMap = find(idxMeas == iMeas)';
		for iRow = rowMap
			delays(idxTraces(iRow,1):idxTraces(iRow,2)) = delayT.Delay(rowDelay);
		end
	end

	% Apply delay time vector to event times
	if length(delays) ~= size(corrT.t_event{iCmb}, 1)
		error(['Length of delay vector (%d) does not match number of ' ...
			'traces (%d) in `corrT.t_event{%d}`.'], length(delays), ...
			size(corrT.t_event{iCmb}, 1), iCmb);
	end
	corrT.t_event{iCmb} = corrT.t_event{iCmb} + delays;
end

% Prevent multiple applications of delay table
corrT.Properties.VariableDescriptions{measRow} = lockString;

% Clean up
clear delayT delays iCmb idxFiles idxMeas idxTraces iFile iMeas iRow
clear lockString measRow measurements nFiles nTraces rowDelay rowMap