function [ allIDs ] = dbscan( setPts, eps, minPts, dist_fn )
%DBSCAN Performs a DBSCAN clustering
%
% According to:
% M. Ester, H.-P. Kriegel et al.: A density-based algorithm for discovering
% clusters in large spatial databases with noise. KDD-96, 1996.
%
% Input:
% ======
%	setPts			n-to-m array of n datapoints that are m-dimensional
%	eps				maximum neighborhood distance; scalar or vector of m
%	minPts			minimum number of points in a cluster
%	dist_fn			distribution function; function handle or string that
%					can be used to indicate distance measure in `pdist`
%
% Returns:
% ========
%	allIDs			vector of length n; indicates cluster labels of points,
%					where 0 indicates noise (no cluster) and NaN indicates
%					nonfinite input values
%
% eps, minPts, and dist_fn can be left out or empty. In this case, they
% will be set to default values.
% The cluster labels are positive integers. Noise (point without cluster)
% is indicated by the cluster label 0. Any rows in setPts that contain NaN
% or Inf values are assigned the cluster label NaN.
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

%% Filter out non-finite values
fnt_idx = all(isfinite(setPts), 2);
setPts = setPts(fnt_idx,:);

%% Initialize values
nPts = size(setPts, 1);
clusterID = NaN(nPts, 1);
currID = 1;

if nargin < 2 || isempty(eps) || ~isfinite(eps)
	% Mark eps for further processing
	eps = NaN;
end

if nargin < 3 || isempty(minPts) || minPts <= 0
	minPts = floor(log10(nPts));
end

if nargin >= 4 && isa(dist_fn, 'function_handle')
elseif nargin >= 4 && (ischar(dist_fn) || isstring(dist_fn))
		dist_fn = @(X) pdist(X, char(dist_fn));
else
		dist_fn = @(X) pdist(X);
end

%% Get matrix of distances
distances = zeros(nPts);
dist_tri = dist_fn(setPts);
i_dist = 1;

for i_1 = 1 : nPts-1
	for i_2 = i_1+1 : nPts
		distances(i_2, i_1) = dist_tri(i_dist);
		distances(i_1, i_2) = dist_tri(i_dist);
		i_dist = i_dist + 1;
	end
end

% Get distance, if not given by user
if isnan(eps)
	q = quantile(dist_tri, [.25 .75]);
	eps = q(2) - q(1);
	eps = 0.5 * eps;
end

%% Main loop
% Loop through all non-visited points
for i_pt = 1:length(clusterID)

	if isnan(clusterID(i_pt)) && expandCluster(i_pt)
		incrementID();
	end
end

%% Write final clusterID vector
% Original points with nonfinite values are assigned NaN
allIDs = nan(length(fnt_idx), 1);
allIDs(fnt_idx) = clusterID;

%% Auxiliary functions

	function isInCluster = expandCluster(i_pt)
		%% Expands the cluster

		% Get neighborhood of current point
		ngbrhood = getNeighbors(i_pt);

		if length(ngbrhood) < minPts
			% Neighborhood too small; mark point as noise
			clusterID(i_pt) = 0;
			isInCluster = false;

		else
			% Neighborhood is large enough; mark as cluster
			clusterID(ngbrhood) = currID;

			% Delete current point from neighborhood list
			ngbrhood(ngbrhood == i_pt) = [];

			% Iterate through neighborhood 
			while ~isempty(ngbrhood)

				% Get indirect neighborhood
				result = getNeighbors(ngbrhood(1));

				% Enough indirect neighbors:
				% Add them to neighborhood, if they still aren’t neighbors
				if length(result) >= minPts
					for i_res = 1:length(result)
						res = result(i_res);

						if isnan(clusterID(res)) || clusterID(res) == 0
							% Add indirectly connected point to neighborhood
							if isnan(clusterID(res))
								ngbrhood(end+1) = res;
							end

							% Mark as cluster
							clusterID(res) = currID;
						end
					end
				end
				ngbrhood(1) = [];

			end
			isInCluster = true;

		end
	end

	function neighbors = getNeighbors(i_pt)
		%% Returns a vector of indices of neighbors in clusterID

		% Find indices with distance < eps
		neighbors = find(distances(i_pt,:) < eps);
	end

	function id = incrementID
		%% Returns the next available cluster ID
		currID = currID + 1;
		id = currID;
	end
end
