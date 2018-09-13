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

% Time scale (convert to hours)
time_scale = 1 ;%/ 6;

%% Define clustering parameters
% %% For combination tmrm-ros
% eps = 2.7;
% minPts = 50;

% %% For combination lyso-tmrm
% % Large clusters
% eps = 3.9;
% minPts = 15;
% % Small clusters
% eps = 2.5;
% minPts = 25;

% %% For combination lyso-ros
% eps = 2.5;
% minPts = 15;

% %% For combination casp-toto
eps = 2.5;
minPts = 30;

%% Logical part
d_scaled = d * time_scale;

idx = dbscan(d_scaled, eps, minPts, 'squaredeuclidean');
nClusters = max(idx);
fprintf('Found %d clusters.\n', nClusters)

% Make noise (idx == 0) black
clusterSizes = zeros(nClusters, 1);
idx_unq = unique(idx);
idx_unq = idx_unq(isfinite(idx_unq));
for i_idx = 1:length(idx_unq)
	clusterSizes(i_idx) = length( idx( idx == idx_unq(i_idx) ) );
end
if idx_unq ~= 0
	clusterSizes = [0;clusterSizes];
end
idx = idx + 1;
clr = log_color(clusterSizes);

% Prepare plot
fh = figure;
for iD = 1:length(idx)
	if isfinite(idx(iD))
		plot(d_scaled(iD,1), d_scaled(iD,2), '.', 'Color', clr(idx(iD),:));
		hold on
	end
end

% Format plot
ax = gca;
ax.YDir = 'reverse';

l = line(xlim, ylim, 'Color', 'k');
uistack(l, 'bottom')

xlabel('$t_\mathrm{event}^{(n)}$ [h]', 'interpreter','latex')
ylabel('$t_\mathrm{event}^{(n+1)}$ [h]', 'interpreter','latex')

%% Make a colormap for clusters weighted by size
function cm = log_color(cs)
	% c0: default colors: c0(1,:) noise, c0(2:ncc+1,:) largest ncc clusters
	c0 = [ ...
			0 0 0 ... % noise: black
			; 1 0 0 ...
			; 0 .45 .75 ...
% 			; 0 1 0 ...
% 			; 0 1 1 ...
		];
	ncc = size(c0,1) - 1;

	% nc: number of colors/clusters (incl. noise)
	nc = length(cs);

	% np: number of points
	np = sum(cs);


	% cm: final colormap
	cm = zeros(nc, 3);

	% Assign to each cluster a color
	for ic = 2:length(cs)

		% Populate final colormap
% 		cm(ic, 1:3) = cs(ic) / np;	% size-dependent grayscane
		cm(ic, 1:3) = .6;			% constant value
	end

	% i_lrg: indices of three largest clusters
	[~, csi] = sort(cs(2:end), 'descend');
	i_lrg = [];
	for i_csi = 1:length(csi)
		i_lrg(i_csi) = csi(i_csi) + 1;
		if i_csi >= ncc
			break
		end
	end

	% Make three largest clusters have default colors
	for il = 1:length(i_lrg)
		cm(i_lrg(il),:) = c0(il+1,:);
	end

end