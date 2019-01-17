function [corrT, combT, fileT, out_dir, restrict_to_conditions] = ...
	loadCorrelationData_batch(inS)
% loadCorrelationData_batch loads correlation data in an automated fashion
%
% Input arguments:
%	inS					optional structure with optional fields:
%		.fileT				file table; loaded if missing
%		.out_dir			target directory; use default if missing
%		.combT				combination table; loaded if missing
%		.comb_sel			vector of indices of selected combinations;
%								prompted if missing
%		.restrict_to_conditions		Condition filter; prompted if missing
%
% Return values:
%	corrT			correlation data table
%	combT			combination table
%	fileT			file table
%	out_dir			target directory
%	restrict_to_conditions		condition filter
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

%% Parse input structure
hasConditionRestrictions = false;
hasFileInfo = false;
hasCombinationInfo = false;
hasCombinationSelection = false;

if nargin >= 1 && isstruct(inS)
	if isfield(inS, 'restrict_to_conditions')
		hasConditionRestrictions = true;
	end
	if isfield(inS, 'fileT') && isfield(inS, 'out_dir')
		hasFileInfo = true;
	end
	if isfield(inS, 'raw_dir') && isfield(inS, 'state_dir')
		raw_dir = inS.raw_dir;
		state_dir = inS.state_dir;
	else
		raw_dir = '';
		state_dir = '';
	end
	if isfield(inS, 'combT')
		hasCombinationInfo = true;
	end
	if isfield(inS, 'comb_sel')
		hasCombinationSelection = true;
	end
end

%% Import data
% Get file information
if hasFileInfo
	fileT = inS.fileT;
	out_dir = inS.out_dir;

else
	if ~ischar(raw_dir) || isempty(raw_dir) || ~ischar(state_dir) || isempty(state_dir)
		useDefault = questdlg(['Use the default input directories ' ...
			'(datadir_raw and datadir_state) or define custom input directories?'], ...
			'Define input directories', ...
			'Default', 'Custom', 'Default');

		if strcmpi(useDefault, 'Custom')
			raw_dir = uigetdir('', 'Raw measurement files');
			state_dir = uigetdir('', 'State files');
		end
	end

	if ischar(raw_dir) && ~isempty(raw_dir) && ischar(state_dir) && ~isempty(state_dir)
		[fileT, out_dir] = getFiles(raw_dir, state_dir);
	else
		[fileT, out_dir] = getFiles;
	end
end

% Get condition restrictions
if hasConditionRestrictions
	restrict_to_conditions = inS.restrict_to_conditions;
else
	all_conds = unique(fileT.condition);
	[cond_sel, ok] = listdlg('ListString', all_conds, ...
		'Name', 'Condition selection', ...
		'PromptString', 'Select the conditions you want to analyze:');

	if ok && ~isempty(cond_sel)
		restrict_to_conditions = all_conds(cond_sel);
	else
		restrict_to_conditions = {};
	end
end

% Get combination info
findArgs = {fileT, restrict_to_conditions};
if hasCombinationInfo
	findArgs{end+1} = inS.combT;
else
	findArgs{end+1} = initCombinations;
end

combT = findIdenticalCells(findArgs{:});

% Get combinations to select
if hasCombinationSelection
	comb_sel = inS.comb_sel;
else
	% Prompt user for combinations to analyze
	text_arrow = [ ' ' char(8594) ' ' ];
	[comb_sel, ~] = listdlg('ListString', join(combT.markers, text_arrow), ...
		'Name', 'Combination selection', ...
		'PromptString', 'Select the combinations you want to load:');
end

%% Load data for selected combinations into `corrT`
% Initialize table with correlation data
corrT = table(string.empty(0,2), cell(0,1), cell(0,1), cell(0,2), cell(0,2), cell(0,2), double.empty(0,3), ...
	'VariableNames', {'markers', 't_event', 'found', 'state', 'params', 'file_map', 'color'});

for i_cmb = 1:length(comb_sel)
	i_sel = comb_sel(i_cmb);	% indices of selected files

	% Load STATE and PARAMS files of selected files
	[state1, params1, file_map1] = loadEventData(fileT, combT.files{i_sel}(:,1));
	[state2, params2, file_map2] = loadEventData(fileT, combT.files{i_sel}(:,2));

	% Extract event times from state
	if ~isempty(state1) && ~isempty(state2)
		t_ev1 = state1(:,2);
		t_ev2 = state2(:,2);
		t_ev = [t_ev1, t_ev2];
	else
		% No traces found for this combination
		warning('No traces found for combination: %s-%s', combT.markers(i_cmb,:));
		continue
	end

% 	% Delete unneeded state information (trace index and event time)
% 	state1(:,1:2) = [];
% 	state2(:,1:2) = [];

	% Get entries with wholly and partially identified events
	% found(:,1)	number of cells with events found for both markers
	% found(:,2)	number of cells with event found only for first marker
	% found(:,3)	number of cells with event found only for second marker
	found = false(size(t_ev, 1), 3);
	found(:,1) = isfinite(t_ev1) & isfinite(t_ev2);
	found(:,2) = isfinite(t_ev1) & ~isfinite(t_ev2);
	found(:,3) = ~isfinite(t_ev1) & isfinite(t_ev2);

	% Write selected content of selected files into `corrT`
	i_tab = height(corrT) + 1;	% index of new entry in `corrT`
	corrT(i_tab,:) = table( combT.markers(i_sel,:), ...	marker names
		{t_ev}, ...										event times
		{found}, ...									found indices
		{state1, state2}, ...							trace state
		{params1, params2}, ...							fit parameters
		{file_map1, file_map2}, ...						file map
		combT{i_sel, 'color'} ...						color info
		);
end

end