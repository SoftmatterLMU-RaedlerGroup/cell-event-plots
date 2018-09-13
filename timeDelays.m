function varargout = timeDelays(cmd)
%timeDelays imports and exports the time delay tables for the measurements.
%
% Input options:
%	cmd			command to be executed [optional]
%					Possible commands are:
%					'read'		read the delay table from the default file
%								and return it (no change to base workspace)
%								This is the default.
%
%					'import'	read the delay table from the default file
%								and write it into the base workspace
%
%					'export'	read the delay table from base workspace
%								and write it to the default file
%
% Return value:
%	varargout{1}	delay table; empty array in case of error
%						varargout{1} is only set if `nargout` >= 1
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

%% Define constants
% Define filename of delay table file
fn = 'delays.json';

% Define commands
cmd_read = 'read';
cmd_import = 'import';
cmd_export = 'export';

% Define delay table name in base workspace
name_tab = 'delayT';

% Define default command
if nargin == 0
	cmd = cmd_read;
end

% Define default output
switch nargout
	case 0
		varargout = {};
	case 1
		varargout = {[]};
end

% Construct full path of delay table file
path = mfilename('fullpath');
[path,~,~] = fileparts(path);
fn = fullfile(path, fn);

%% Execute depending on command
if strcmpi(cmd, cmd_import) || strcmpi(cmd, cmd_read)
	% Import delay table from file

	%% Read JSON from file
	[fh,err] = fopen(fn, 'r');
	if fh == -1
		warning('Cannot open file %s: %s', fn, err);
		return
	end
	delayJ = fscanf(fh, '%s');
	fclose(fh);

	%% Convert delay values from JSON to table
	delayJ = jsondecode(delayJ);

	delayTab = table(string.empty(0,1), [], 'VariableNames', {'Date', 'Delay'});

	for iRow = 1:numel(delayJ)
		delayTab(iRow,:) = {delayJ(iRow).Date, delayJ(iRow).Delay};
	end

	%% Export table to base workspace
	if strcmpi(cmd, cmd_import)
		isTab = evalin('base', sprintf('exist(''%s'', ''var'')', name_tab));

		if isTab
			warning('Overwrite existing `delayT` table');
		end

		assignin('base', name_tab, delayTab);
	end

elseif strcmpi(cmd, cmd_export)
	% Export delay table to file

	%% Import table from base workspace
	isTab = evalin('base', sprintf('exist(''%s'', ''var'')', name_tab));

	if ~isTab
		warning('No delay table found; do nothing.');
	end

	delayTab = evalin('base', name_tab);

	if ~istable(delayTab)
		error('Delay table is not a table.');
	end

	%% Convert table to JSON format
% 	delayJson = jsonencode(delayTab);
	delayJson = toJSON(delayTab);

	%% Write converted table to file
	[fh,err] = fopen(fn, 'w');

	if fh == -1
		warning('Cannot open file %s: %s', fn, err);
		return
	end

	fprintf(fh, '%s', delayJson);
	fclose(fh);

else
	% Undefined command; throw error
	error('Unknown option: %s', cmd);
end

if nargout > 0
	varargout = {delayTab};
end
end

%% Convert to formatted JSON
function json = toJSON(tab)
json = sprintf('[');
comma = '';

for iRow = 1:height(tab)

	json = [json sprintf('%s\n\t{\n\t\t"%s": "%s",\n\t\t"%s": %f\n\t}', ...
		comma, 'Date', tab{iRow,'Date'}, 'Delay', tab{iRow,'Delay'})];

	if iRow == 1
		comma = ',';
	end

end

json = sprintf('%s\n]', json);
end