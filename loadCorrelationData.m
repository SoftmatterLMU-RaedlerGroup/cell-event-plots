% This script calls `loadCorrelationData_batch` and saves its return values
% in the current workspace.
% If variables to be used by `loadCorrelationData_batch` are found in the
% workspace, they are given to it.
%
% This makes the script convenient for manual use of
% `loadCorrelationData_batch`.
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
callS = struct;

if exist('raw_dir', 'var')
	callS.raw_dir = raw_dir;
end
if exist('state_dir', 'var')
	callS.state_dir, state_dir;
end
if exist('fileT', 'var')
	callS.fileT = fileT;
end
if exist('out_dir', 'var')
	callS.out_dir = out_dir;
end
if exist('combT', 'var')
	callS.combT = combT;
end
if exist('restrict_to_conditions', 'var')
	callS.restrict_to_conditions = restrict_to_conditions;
end

[corrT, combT, fileT, out_dir, restrict_to_conditions] = ...
	loadCorrelationData_batch(callS);

clear callS