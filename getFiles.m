function [ fileT, out_dir ] = getFiles( raw_dir, state_dir, out_dir )
%GETFILES Get a list of files to read from
%
% Parameters:
% ===========
%	in_dir		Directory to read from (optional)
%	out_dir		Directory to write to (optional)
%
%
% Returns:
% ========
%	fileT			Table of file names
%	out_dir			target directory for results
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

%% Get constants

% Get markers to exclude
blacklist_markers = "hoechst33342";

% Default input directories
default_rawdir = 'datadir_raw';
default_statedir = 'datadir_state';
% default_paramdir = default_statedir;

% Default output directory
default_outdir = 'out';

% Skript directory
mfilepath = fileparts(mfilename('fullpath'));

%% Evaluate parameters
% Ensure that raw_dir is a non_empty cell array of valid directory names
if nargin >= 1
	if ~isstring(raw_dir)
		raw_dir = string(raw_dir);
	end
else
	raw_dir = string(fullfile(mfilepath, default_rawdir));
	if ~exist(raw_dir{1}, 'dir')
		error('No input directory for raw files found.')
	end
end

% Ensure that state_dir is a non-empty cell array of valid directory names
if nargin >= 2 && ~isempty(state_dir)
	if ~iscell(state_dir)
		state_dir = { state_dir };
	end

	for i_dir = numel(state_dir):-1:1
		if ~exist(state_dir{i_dir}, 'dir')
			state_dir(i_dir) = [];
		end
	end
end

if nargin < 2 || isempty(state_dir)
	state_dir = { fullfile(mfilepath, default_statedir) };
	if ~exist(state_dir{1}, 'dir')
		error('No input directory for state files found.')
	end
end

% Set `param_dir` to `state_dir`
param_dir = state_dir;

% Ensure that out_dir is a valid directory name
if nargin < 3
	out_dir = fullfile(mfilepath, default_outdir);
end

if ~exist(out_dir, 'dir')
	[s,~,~] = mkdir(out_dir);
	if ~s
		error('Could not create output directory.')
	end
end

%% Analyze state directory
typeTokens = ["STATE", "PARAMS"];

% Initialize output struct
dirS = struct('condition',{}, 'marker',{}, 'path',{});
fileT_state = table(uint32.empty(0,1), string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
	string.empty(0,1), string.empty(0,1), string.empty(0,1), 'VariableNames', {'index_f', 'condition', ...
	'marker', 'position', 'measurement', 'path_state', 'path_params'});

% Populate output struct with input directories
for i_stdir = numel(state_dir):-1:1
	dirS(i_stdir).path = state_dir{i_stdir};
end

% Read in first directory depth level
for i_ds = length(dirS):-1:1
	temp_path = dirS(i_ds).path;
	[fileTab, dirList] = analyzeDir(temp_path, typeTokens);
	dirS(i_ds) = [];

	if ~isempty(fileTab)
		idxStart = height(fileT_state) + 1;
		idxEnd = height(fileT_state) + height(fileTab);
		fileT_state(idxStart:idxEnd, ...
			{'position', 'measurement', 'path_state', 'path_params'}) = ...
			fileTab(:, {'position', 'measurement', typeTokens{:}});

	elseif ~isempty(dirList)
		fulldirList = fullfile(temp_path, dirList);
		for i_dl = 1:numel(dirList)
			index = length(dirS) + 1;
			dirS(index).path = fulldirList{i_dl};
			dirS(index).marker = dirList{i_dl};
		end
	end
end

% Read in second directory depth level
if isempty(fileT_state)
	for i_ds = length(dirS):-1:1
		temp_path = dirS(i_ds).path;
		marker = dirS(i_ds).marker;
		[fileTab, dirList] = analyzeDir(temp_path, typeTokens);
		dirS(i_ds) = [];

		if ~isempty(fileTab)
			idxStart = height(fileT_state) + 1;
			idxEnd = height(fileT_state) + height(fileTab);
			fileT_state(idxStart:idxEnd, ...
				{'position', 'measurement', 'path_state', 'path_params'}) = ...
				fileTab(:, {'position', 'measurement', typeTokens{:}});
			fileT_state.marker(idxStart:idxEnd) = marker;

		elseif ~isempty(dirList)
			fulldirList = fullfile(temp_path, dirList);
			for i_dl = 1:numel(dirList)
				index = length(dirS) + 1;
				dirS(index).path = fulldirList{i_dl};
				dirS(index).marker = dirList{i_dl};
				dirS(index).condition = marker;
			end
		end
	end
end

% Read in third directory depth level
if isempty(fileT_state)
	for i_ds = length(dirS):-1:1
		temp_path = dirS(i_ds).path;
		marker = dirS(i_ds).marker;
		condition = dirS(i_ds).condition;
		[fileTab, dirList] = analyzeDir(temp_path, typeTokens);
		dirS(i_ds) = [];

		if ~isempty(fileTab)
			idxStart = height(fileT_state) + 1;
			idxEnd = height(fileT_state) + height(fileTab);
			fileT_state(idxStart:idxEnd, ...
				{'position', 'measurement', 'path_state', 'path_params'}) = ...
				fileTab(:, {'position', 'measurement', typeTokens{:}});
			fileT_state.marker(idxStart:idxEnd) = marker;
			fileT_state.condition(idxStart:idxEnd) = condition;

		elseif ~isempty(dirList)
			fulldirList = fullfile(temp_path, dirList);
			for i_dl = 1:numel(dirList)
				index = length(dirS) + 1;
				dirS(index).path = fulldirList{i_dl};
				dirS(index).marker = dirList{i_dl};
				dirS(index).condition = marker;
			end
		end
	end
end

% Read in fourth directory depth level
if isempty(fileT_state)
	for i_ds = length(dirS):-1:1
		temp_path = dirS(i_ds).path;
		marker = dirS(i_ds).marker;
		condition = dirS(i_ds).condition;
		[fileTab, ~] = analyzeDir(temp_path, typeTokens);
		dirS(i_ds) = [];

		if ~isempty(fileTab)
			nFiles = height(fileTab);
% 			fileT_state(height(fileT_state)+(1:nFiles), ...
% 				{'position', 'measurement', 'path_state', 'path_params'}) = ...
% 				fileTab(:, {'position', 'measurement', typeTokens{:}});
% 			fileT_state.marker(idxStart:idxEnd) = marker;
% 			fileT_state.condition(idxStart:idxEnd) = condition;
			
			fileTab.Properties.VariableNames{'STATE'} = 'path_state';
			fileTab.Properties.VariableNames{'PARAMS'} = 'path_params';
			fileT_state = [ ...
				fileT_state; [ table(zeros(nFiles, 1, 'uint32'), ...
				repmat(string(condition), nFiles, 1), ...
				repmat(string(marker), nFiles, 1), ...
				'VariableNames', {'index_f', 'condition', 'marker'}), ...
				fileTab ] ];
		end
	end
end

% % If there are no files in the directory, throw an error
% if isempty(fileT_state)
% 	error('There were no files found.')
% end

%% Find raw traces
fileT_raw = get_raw_files(raw_dir);

%% Merge tables of raw and state data
fileT = outerjoin(fileT_state, fileT_raw, 'Type','full', 'MergeKeys',true);
if isempty(fileT)
	error('There were no files found.');
end
index_f = 1:height(fileT);
fileT.index_f = index_f';

% End of function

%% Auxiliary function
function fileT_raw = get_raw_files(raw_dir)

%% Get content of raw directories (-> measurement directories)
tab_mes = table(string({}), string({}), 'VariableNames', {'measurement', 'path'});
for i_raw = 1:numel(raw_dir)
	cont_raw = dir(raw_dir{i_raw});
	
	for i_cont = 1:length(cont_raw)
		% Filter out entries that are no directories
		if ~cont_raw(i_cont).isdir
			continue
		end

		% Filter out entries beginning with “.”
		if any(strfind(cont_raw(i_cont).name, '.') == 1)
			continue
		end

		% Filter out entries with an underscore
		if ~isempty(strfind(cont_raw(i_cont).name, '_'))
			continue
		end

		% Add valid entry to table `tab_mes`
		tab_mes(end+1,:) = { cont_raw(i_cont).name, ...
			fullfile(cont_raw(i_cont).folder, cont_raw(i_cont).name) };
	end
end

%% Get content of measurement directories (-> condition directories)
tab_cond = table(string({}), string({}), string({}), ...
	'VariableNames', {'condition', 'measurement', 'path'});
for i_mes = 1:height(tab_mes)
	cont_mes = dir(tab_mes.path{i_mes});

	for i_cont = 1:length(cont_mes)
		% Filter out entries that are no directories
		if ~cont_mes(i_cont).isdir
			continue
		end

		% Filter out entries beginning with “.”
		if any(strfind(cont_mes(i_cont).name, '.') == 1)
			continue
		end

		% Filter out entries with an underscore
		if ~isempty(strfind(cont_mes(i_cont).name, '_'))
			continue
		end

		% Add valid entry to table `tab_cond`
		tab_cond(end+1,:) = { cont_mes(i_cont).name, ...
			tab_mes.measurement{i_mes}, ...
			fullfile(cont_mes(i_cont).folder, cont_mes(i_cont).name) };
	end
end

%% Get content of condition directories (-> markers)
tab_mark = table(string({}), string({}), string({}), string({}), ...
	'VariableNames', {'marker', 'condition', 'measurement', 'path'});
for i_cond = 1:height(tab_cond)
	cont_cond = dir(tab_cond.path{i_cond});

	for i_cont = 1:length(cont_cond)
		% Filter out entries that are no directories
		if ~cont_cond(i_cont).isdir
			continue
		end

		% Filter out entries beginning with “.”
		if any(strfind(cont_cond(i_cont).name, '.') == 1)
			continue
		end

		% Filter out entries with an underscore
		if ~isempty(strfind(cont_cond(i_cont).name, '_'))
			continue
		end

		% Filter out blacklisted markers
		if ismember(cont_cond(i_cont).name, blacklist_markers)
			continue
		end

		% Add valid entry to table `tab_mark`
		tab_mark(end+1,:) = { cont_cond(i_cont).name, ...
			tab_cond.condition{i_cond}, ...
			tab_cond.measurement{i_cond}, ...
			fullfile(cont_cond(i_cont).folder, cont_cond(i_cont).name) };
	end
end

%% Get content of marker directories (-> raw files)
fileT_raw = table(string([]), string([]), string([]), ...
	string([]), string([]), 'VariableNames', {'condition', ...
	'marker', 'position', 'measurement', 'path_raw'});
for i_mark = 1:height(tab_mark)
	cont_mark = dir(tab_mark.path{i_mark});

	for i_cont = 1:length(cont_mark)
		% Filter out entries that are no directories
		if cont_mark(i_cont).isdir
			continue
		end

		% Filter out entries beginning with “.”
		if any(strfind(cont_mark(i_cont).name, '.') == 1)
			continue
		end

		% Find measurement position
		pos = regexpi(cont_mark(i_cont).name, [ '^(.+)_' ...
			regexptranslate('escape', tab_mark.measurement{i_mark}) ], 'tokens');

		% Add valid entry to table `tab_mark`
		fileT_raw(end+1,:) = { tab_mark.condition(i_mark), ...
			tab_mark.marker(i_mark), ...
			pos, ...
			tab_mark.measurement(i_mark), ...
			fullfile(cont_mark(i_cont).folder, cont_mark(i_cont).name) };
	end
end

end %end of function `fileT_raw`


end % end of function `getFiles`

function [pathTab, dirList] = analyzeDir(query_dir, typeTokens)
%% analyzeDir extracts a list of valid files and of directories in a
% directory.
%
% Parameter:
% ==========
%	query_dir		the directory to be tested
%	typeTokens		[optional] token required in file names
%
% Returns:
% ========
%	fileList		a cell array of valid files
%	dirList			a cell array of directories

% Validate input
if ~exist(query_dir, 'dir')
	return
end

if nargin < 2
	typeTokens = ["STATE", "PARAMS"];
end

% Prepare output data
pathTab = table(string.empty, string.empty, ...
	'VariableNames', {'position', 'measurement'});
for iTok = 1:numel(typeTokens)
	pathTab = [pathTab, table(string.empty, 'VariableNames', {typeTokens{iTok}})];
end
fileList = string.empty;
dirList = cell(0);

% Process directory content
dirCont = dir(query_dir);

for i_cont = 1:length(dirCont)
	if strncmp(dirCont(i_cont).name, '.', 1)
		continue
	elseif dirCont(i_cont).isdir
		dirList = [ dirList ; dirCont(i_cont).name ];
	else
		fileList = [ fileList ; dirCont(i_cont).name ];
	end
end

% Perform matching
keywords = strjoin(typeTokens, '|');

regex = join(['^(?<position>.*)_(?<measurement>[^_]*)_ALL_(?<type>' ...
	keywords ')(?:_[^_]*)?\.(?:csv|txt)$'], '');
regexRes = regexpi(fileList, regex, 'once', 'names');

% Enter match results into table
for iRes = 1:numel(regexRes)
	res = regexRes{iRes};
	if isempty(res)
		continue
	end

	pos = res.position;
	meas = res.measurement;
	type = res.type{1};

	idx = find(ismember(pathTab.position, pos) & ismember(pathTab.measurement, meas));
	if isempty(idx)
		idx = height(pathTab) + 1;
% 		pathTab(idx, {'position', 'measurement'}) = {pos, meas};
		pathTab = [pathTab; pathTabTemplate(pos, meas)];
	end
	pathTab{idx, type} = string(fullfile(query_dir, fileList{iRes}));
end

	function pathTabRow = pathTabTemplate(pos, meas)
		%pathTabTemplate creates a template row for `pathTab`
		pathTabRow = table(pos, meas, ...
			'VariableNames', {'position', 'measurement'});
		for jTok = 1:numel(typeTokens)
			pathTabRow = [pathTabRow, table(missing, 'VariableNames', {typeTokens{jTok}})];
		end
	end

end % end of function `analyzeDir`

