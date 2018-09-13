%% Tests clustering algorithms
%  Requires data d
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

nClusters = 3;

d_clean = d(all(isfinite(d), 2), :);

idx = kmeans(d_clean, nClusters, 'Distance', 'cosine');

fh = figure;
clr = lines(nClusters);

for iD = 1:length(idx)
	if isfinite(idx(iD))
		plot(d_clean(iD,1), d_clean(iD,2), '.', 'Color', clr(idx(iD),:));
		hold on
	end
end
