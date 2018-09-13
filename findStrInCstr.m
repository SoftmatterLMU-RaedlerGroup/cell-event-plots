function [idx] = findStrInCstr(cstr, pattern, options)
% This function finds the string pattern in the cellstring cstr
% and returns the index idx of its first occurrence.
% If pattern is not found in cstr, [] is returned.
% The argument options is not required. If options is a char array and
% contains an 'i', a case insensitive search will be performed. If options
% is a char array and contains an 'm', idx will be a vector containing the
% indices of all occurrences of pattern in cstr. If options is a char array
% and contains a '0', idx will be 0 if no occurrences were found.
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

%% Initialize variables
idx = [];
ignoreCase = false;
multipleOccurrences = false;
zeroNegative = false;

if nargin > 2 && ischar(options)

	% Assess whether to perform a case-insensitive search
	if ~isempty(strfind(options, 'i'))
		ignoreCase = true;
	end

	% Assess whether to search for multiple occurrences
	if ~isempty(strfind(options, 'm'))
		multipleOccurrences = true;
	end

	% Assess whether to indicate a negative result by zero
	if ~isempty(strfind(options, '0'))
		zeroNegative = true;
	end
end

%% Perform search
for i_str = 1:numel(cstr)
	
	if ignoreCase
		occurrenceFound = strcmpi(cstr{i_str}, pattern);
	else
		occurrenceFound = strcmp(cstr{i_str}, pattern);
	end

	if occurrenceFound
		idx = [ idx i_str ];

		if ~multipleOccurrences
			break
		end
	end
end

%% Postprocess
if zeroNegative && isempty(idx)
	idx = 0;
end