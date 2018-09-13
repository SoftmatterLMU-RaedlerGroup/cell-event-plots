function fullline = labelcat(segments)
%labelcat concatenates segments of labels
%
% Input argument:
%	segments	vector of labels grouped as contiguous segments,
%				the segments separated by NaN
%
% Return value:
%	fullline	vector of label sequences that are combined as long as
%				possible; colliding segments are separated by NaN.
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


% Build contiguous line
nl = {};
sgm = [];
p_start = [];
p_end = [];

for iP = length(segments):-1:1

	% Read next element
	p = segments(iP);
	if ~isfinite(p)

		if ~isempty(p_start) && ~isempty(p_end)
			sgm = segments(p_start:p_end);
			if ~iscolumn(sgm)
				sgm = sgm';
			end
			segments(p_start:p_end) = [];
		end

		p_start = [];
		p_end = [];

	elseif isempty(p_end)
		p_end = iP;

	else
		p_start = iP;
	end

	if isempty(sgm)
		continue
	end

	% Segment detected; search for matching point
	m1 = searchMatchPoint(sgm(1), nl);
	m2 = searchMatchPoint(sgm(end), nl);

	if isempty(m1) && isempty(m2)
		nl{numel(nl)+1} = sgm;

	elseif isempty(m1)
		if m2(2) == 1
			nl{m2(1)} = [ sgm(1:end-1); nl{m2(1)} ];
		else
			nl{m2(1)} = [ nl{m2(1)}; flip(sgm(1:end-1)) ];
		end

	elseif isempty(m2)
		if m1(2) == 1
			nl{m1(1)} = [ flip(sgm(2:end)); nl{m1(1)} ];
		else
			nl{m1(1)} = [ nl{m1(1)}; sgm(2:end) ];
		end

	elseif ~isempty(m1) && ~isempty(m2)
		if length(sgm) > 2
			sgm = sgm(2:end-1);
		else
			sgm = [];
		end

		if m1(2) == 1
			nl{m1(1)} = flip(nl{m1(1)});
		end
		if m2(2) ~= 1
			nl{m2(1)} = flip(nl{m2(1)});
		end
		nl{m1(1)} = [ nl{m1(1)}; sgm; nl{m2(1)} ];
		nl(m2(1)) = [];
	end
	sgm = [];
end

if isempty(nl)
	fullline = [];
	return
end

% Concatenate segments to one vector
fullline = nl{1};
for i = 2:numel(nl)
	fullline = [ fullline; NaN; nl{i} ];
end

end % end of `labelcat`

function idx = searchMatchPoint(p, Pts)
%searchMatchPoint returns the first free occurrence of `p` in `Pts`.
% Free occurrences are occurrences at the first or last index of a vector.
% `p` is a scalar point label
% `Pts` is a cell array of nx1 vectors of n point labels
% `idx` is a mx2 matrix of m matches, where the first colum gives the
% element of `Pts` and the second column the index of `Pts{idx(1)}`
	idx = zeros(0, 2);

	for iC = 1:numel(Pts)
		isThere = find(ismember(Pts{iC}, p), 1);
		if ~isempty(isThere) && ismember(isThere, [1 length(Pts{iC})])
			idx = [iC isThere];
			break
		end
	end
end
