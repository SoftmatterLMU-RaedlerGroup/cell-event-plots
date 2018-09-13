function [ counts, centers, edges ] = hist3count( data_x, data_y, nbins_x, nbins_y, bounds_x, bounds_y )
%HIST3COUNT counts data points in two-dimensional bins
%
% Inupt values:
%	data_x			vector of the data in x-direction
%	data_y			vector of the data in y-direction
%	n_bins_x		number of bins in x-direction
%	n_bins_y		number of bins in y-direction
%	bounds_x		bounds of x-axis
%	bounds_y		bounds of y-axis
%
% `n_bins_x` and `n_bins_y` are optional. If omitted, they will be assigned
% a default value of 100.
% `bounds_x` and `bounds_y` are optional. If omitted, they will be set to
% the minimum and maximum, respectively, of the data vectors. If a minimum
% or maximum is not finite, it will be replaced by 0 or 1, respectively.
%
% Output values:
%	counts				structure with bin count data
%		.matrix				`n_bins_y`-by-`n_bins_x` matrix with bin counts
%		.linearized			`n_bins_y`*`n_bins_x` vector with bin counts
%		.normfac			normalization factor
%	centers				structure with bin center position data
%		.x					`n_bins_x` vector of bin center x-positions
%		.y					`n_bins_y` vector of bin center y-positions
%		.linearized			matrix for indexing `counts.linearized`
%	edges				structure of bin edge position data
%		.x					`n_bins_x`+1 vector of edges x-positions
%		.y					`n_bins_y`+1 vector of edges y-positions
%
% Divide `counts.matrix` or `counts.linearized` by `counts.normfac` to
% normalize them.
% `centers.linearized` has as many rows as `counts.linearized` and two
% columns. The first column are the x-positions and the second column the
% y-positions of the corresponding elements in `counts.linearized`.
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
nbins_default = 100;

%% Validate input data
data_x = data_x(:);
data_y = data_y(:);

% nbins_x
if nargin < 3 || ~isscalar(nbins_x) || ~isnumeric(nbins_x)
	nbins_x = nbins_default;
end

% nbins_y
if nargin < 4 || ~isscalar(nbins_y) || ~isnumeric(nbins_y)
	nbins_y = nbins_default;
end

% bounds_x
if nargin < 5
	bounds_x = [];
end
[min_x, max_x] = iGetBounds(bounds_x, data_x);

% bounds_y
if nargin < 6
	bounds_y = [];
end
[min_y, max_y] = iGetBounds(bounds_y, data_y);

%% Create histogram information
edges.x = linspace(min_x, max_x, nbins_x + 1);
edges.y = linspace(min_y, max_y, nbins_y + 1);

width_x = (max_x - min_x) / nbins_x;
width_y = (max_y - min_y) / nbins_y;

centers.x = edges.x(1:end-1) + width_x/2;
centers.y = edges.y(1:end-1) + width_y/2;

% Create and populate array of bin values
counts.matrix = zeros(nbins_y, nbins_x);

for iy = 1:nbins_y
	idx_y = data_y >= edges.y(iy) & data_y < edges.y(iy + 1);
	for ix = 1:nbins_x
		idx_x = data_x >= edges.x(ix) & data_x < edges.x(ix + 1);
		counts.matrix(iy, ix) = counts.matrix(iy, ix) + sum(idx_y & idx_x);
	end
end

% Create linearized representations for fitting (linearized input data)
counts.linearized = zeros(nbins_y * nbins_x, 1);
centers.linearized = zeros(nbins_y * nbins_x, 2);

il = 0;
for iy = 1:nbins_y
	for ix = 1:nbins_x
		il = il + 1;
		counts.linearized(il) = counts.matrix(iy, ix);
		centers.linearized(il, 1:2) = [centers.x(ix) centers.y(iy)];
	end
end

% Calculate normalization factor `counts.normfac` such that
% sum(counts.matrix(:)) / counts.normfac == 1
counts.normfac = sum(counts.matrix(:)) * width_x * width_y;

end

function [min_val, max_val] = iGetBounds(bounds, data)
% Retrieve bounds for histogram from input parameters or from data values
min_val = NaN;
max_val = NaN;

if isnumeric(bounds) && length(bounds) >= 1
	min_val = bounds(1);
	
	if length(bounds) >= 2
		max_val = bounds(2);
	end
end

if ~isfinite(min_val)
	min_val = min(data(isfinite(data)));
	if isempty(min_val)
		min_val = 0;
	end
end

if ~isfinite(max_val)
	max_val = max(data(isfinite(data)));
	if isempty(max_val)
		max_val = 1;
	end
end

end
