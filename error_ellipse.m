function [ellipse_rotated, x0, y0, eigvec_max, eigval_max, ...
	eigvec_min, eigval_min, isInEllipse] = error_ellipse(data)
%error_ellipse calculates an error ellipse
%
% Input:
%	data	n-by-2 matrix of data points who are to be indicated by the
%			ellipse. `n` is the number of points. The first column holds
%			the x-values, the second column holds the y-values. `data` may
%			contain NaN and +/-Inf, but must contain at least 2 points with
%			finite coordinates.
%
% Returns:
%	ellipse_rotated		m-by-2 matrix of coordinates of the error ellipse,
%						with `m` points and the x- (y-) values in the first
%						(second) column. `m` is the default number of
%						values returned by `linspace`.
%						If only one output value is requested, the
%						coordinates are shifted to the data location. Else,
%						the error ellipse is centered at the origin and can
%						be shifted to the data location by adding
%						`[x0,y0]`.
%
%	x0					The x-coordinate of the data location (where you
%						might want the ellipse to be centered).
%
%	y0					The y-coordinate of the data location (where you
%						might want the ellipse to be centered).
%
%	eigvec_max			The eigenvector of the covariance matrix of the
%						data corresponding to the larger eigenvalue.
%
%	eigval_max			The larger eigenvalue of the covariance matrix of
%						the data.
%
%	eigvec_min			The eigenvector of the covariance matrix of the
%						data corresponding to the smaller eigenvalue.
%
%	eigval_min			The smaller eigenvalue of the covariance matrix of
%						the data.
%
%	isInEllipse			logical vector of length `m`, whose i-th component
%						indicates whether the point of the i-th line in
%						`data` is inside the error ellipse (true) or not
%						(false). Points on the ellipse are counted as
%						inside.
%
% The algorithm is based on this Github gist by Piyush3dB:
% https://gist.github.com/Piyush3dB/bf2c83a8eb7344798644
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

% Filter out non-finite values (inf, NaN)
idx_finite = all(isfinite(data), 2);
data = data(idx_finite, :);

% Allocate inside-ellipse-index
isInEllipse = false(length(idx_finite), 1);

% Catch error: data empty
if isempty(data) || size(data, 1) < 2 || size(data, 2) < 2
	ellipse_rotated = NaN(1, 2);
	x0 = NaN;
	y0 = NaN;
	eigvec_max = NaN;
	eigval_max = NaN;
	eigvec_min = NaN;
	eigval_min = NaN;
	return
end

% Calculate the eigenvectors and eigenvalues
covar_data = cov(data);
[eigvecs, eigvals] = eig(covar_data);

% Check which is the smaller/larger eigenvalue
[eigvals_sorted, eigvals_sorted_idx] = sort(max(eigvals));
eigval_min = eigvals_sorted(1);
eigval_max = eigvals_sorted(2);
eigvec_min = eigvecs(:,eigvals_sorted_idx(1));
eigvec_max = eigvecs(:,eigvals_sorted_idx(2));

% Calculate larger and smaller semi-axis
interval = 1; % 1 for no scaling, 2.3 for 1sigma, 2.4477 for 2sigma
sx_max = interval * sqrt(eigval_max);
sx_min = interval * sqrt(eigval_min);

% Construct x- and y-coordinates of ellipse
thetas = linspace(0, 2 * pi);
ellipse_x = sx_max * cos(thetas);
ellipse_y = sx_min * sin(thetas);

% Rotate ellipse to right direction if tilted
angle = atan2(eigvec_max(2), eigvec_max(1));
cos_ang = cos(angle);
sin_ang = sin(angle);
R = [cos_ang, sin_ang; -sin_ang, cos_ang];
ellipse_rotated = [ellipse_x; ellipse_y]' * R;

% Calculate center of mass of data
avg = median(data);
x0 = avg(1);
y0 = avg(2);

% If only one output argument is used, move ellipse to center of mass
% If eight output arguments are used, check which points are in ellipse
if nargout == 1
	ellipse_rotated(:,1) = ellipse_rotated(:,1) + x0;
	ellipse_rotated(:,2) = ellipse_rotated(:,2) + y0;
elseif nargout >= 8
	isInEllipse(idx_finite) = ...
	( ((data(:,1)-x0)*cos_ang + (data(:,2)-y0)*sin_ang) / sx_max ).^2 + ...
	( ((data(:,2)-y0)*cos_ang - (data(:,1)-x0)*sin_ang) / sx_min ).^2 <= 1;
end