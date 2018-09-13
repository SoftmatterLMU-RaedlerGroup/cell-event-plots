function s = center_errors(data)
%center_errors calculates the errors of the ellipse coordinates
% The errors are calculated as standard deviation of the mean
% according to Demtröder Bd. 1, Kap. 1.8.3, Gl. 1.14
%
% `data` is a nxm matrix of n datasets in m dimensions
% `s` is a 1xm vector of the errors for each dimension
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
	n = size(data, 1);
	if n < 2
		s = NaN(1, size(data, 2));
		return
	end
	s = sqrt(sum((mean(data) - data).^2) / n / (n - 1));
end
