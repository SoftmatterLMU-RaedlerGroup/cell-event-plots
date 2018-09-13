function timeNow = getTime
%GETTIME builds a time string for use in time-dependent file names
%
% Returns:
%	timeNow		character vector of current time
%				format: YYYY-MM-DD–hhmmss
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

text_dash = char(8211);
timeNow = clock;
timeNow = sprintf('%4d-%02d-%02d%s%02d%02d%02.0f', ...
	timeNow(1), timeNow(2), timeNow(3), ...
	text_dash, timeNow(4), timeNow(5), timeNow(6));

end

