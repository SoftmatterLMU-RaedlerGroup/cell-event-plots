function [ellipse_rotated, x0, y0, eigvec_max, eigval_max, ...
	eigvec_min, eigval_min, isInEllipse, info] = error_potato(data, direction)
%error_potato calculates an asymmetric (“potato-shaped”) error ellipse
%
% Input:
%	data	n-by-2 matrix of data points who are to be indicated by the
%			ellipse. `n` is the number of points. The first column holds
%			the x-values, the second column holds the y-values. `data` may
%			contain NaN and +/-Inf, but must contain at least 2 points with
%			finite coordinates.
%
%	direction	one of 'both', 'large', 'small'. Determines in direction of
%				which semi-axis the ellipse shall be asymmetric.
%				'both': in direction of both semi-axes
%				'large': only in direction of larger semi-axis
%				'small': only in direction of smaller semi-axis
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
% The algorithm for the symmetric ellipse is based on this Github
% gist by Piyush3dB:
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

% Get directions of asymmetry
dir_asym_large = true;
dir_asym_small = true;
if exist('direction', 'var')
	if ischar(direction)
		switch lower(direction)
			case 'both'
			case 'large'
				dir_asym_small = false;
			case 'small'
				dir_asym_large = false;
			otherwise
				error('Unknown direction: %s', direction)
		end
	else
		error('`direction` must be a char vector.')
	end
end

% Filter out non-finite values (inf, NaN)
idx_finite = all(isfinite(data), 2);
data = data(idx_finite, :);

% Catch error: data empty
if isempty(data) || size(data, 1) < 2 || size(data, 2) < 2
	ellipse_rotated = NaN(1, 2);
	x0 = NaN;
	y0 = NaN;
	eigvec_max = NaN;
	eigval_max = NaN;
	eigvec_min = NaN;
	eigval_min = NaN;
	isInEllipse = [];
	info = struct;
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

% Rotate ellipse to right direction if tilted
angle = atan2(eigvec_max(2), eigvec_max(1));
cos_ang = cos(angle);
sin_ang = sin(angle);
R = [cos_ang, sin_ang; -sin_ang, cos_ang];

% Rotate data for better analysis
data_rot = data * R';

% Get the coordinates of the data mean
avg = median(data);
x0 = avg(1);
y0 = avg(2);
avg_rot = avg * R';

% Get the standard deviations in different directions
if dir_asym_large
	stdXpos = custom_std(avg_rot(1), data_rot(data_rot(:,1) >= avg_rot(1), 1));
	stdXneg = custom_std(avg_rot(1), data_rot(data_rot(:,1) <= avg_rot(1), 1));
else
	stdXpos = custom_std(avg_rot(1), data_rot(:,1));
	stdXneg = stdXpos;
end
if dir_asym_small
	stdYpos = custom_std(avg_rot(2), data_rot(data_rot(:,2) >= avg_rot(2), 2));
	stdYneg = custom_std(avg_rot(2), data_rot(data_rot(:,2) <= avg_rot(2), 2));
else
	stdYpos = custom_std(avg_rot(2), data_rot(:,2));
	stdYneg = stdYpos;
end

% Build the ellipse
interval = 1; % 1 for no scaling, 2.3 for 1sigma, 2.4477 for 2sigma
thetas = linspace(0, 2 * pi);
ellipse_rotated = NaN(length(thetas), 2);
for iT = 1:length(thetas)
	theta = thetas(iT);
	r = semiaxes(theta);
	ellipse_rotated(iT,1) = r(1) * cos(theta);
	ellipse_rotated(iT,2) = r(2) * sin(theta);
end
ellipse_rotated = ellipse_rotated * R;

% Prepare output values
if nargout == 1
	% If only the ellipse is selected, move it to right position
	ellipse_rotated(:,1) = ellipse_rotated(:,1) + x0;
	ellipse_rotated(:,2) = ellipse_rotated(:,2) + y0;

elseif nargout >= 8
	% Identify points inside of ellipse

	% Initialize index vector
	isInEllipse = false(length(idx_finite), 1);

	% Get vectors of points with origin as center and their norms
	orig_vec = data_rot(idx_finite,:) - avg_rot;
	orig_len = sqrt(sum(orig_vec.^2, 2));

	% Get angles of origin vectors in [0;2*pi]
	orig_ang = atan2(orig_vec(:,2) ./ orig_len, orig_vec(:,1) ./ orig_len);
	idx_negAng = orig_ang < 0;
	orig_ang(idx_negAng) = 2 * pi + orig_ang(idx_negAng);

	% Test for each point if it fulfills the ellipse equation
	for iP = find(idx_finite)'
		r = semiaxes(orig_ang(iP));
		isInEllipse(iP) = ...
			(orig_vec(iP,1) / r(1)).^2 + (orig_vec(iP,2) / r(2)).^2 <= 1;
	end
end
if nargout >= 9
	info.angle = angle;
	info.stdXpos = stdXpos;
	info.stdXneg = stdXneg;
	info.stdYpos = stdYpos;
	info.stdYneg = stdYneg;
end

	function r = semiaxes(theta)
	%semiaxes returns the semi-axes of the error_potato at a given angle
	%	r = [ <major semi-axis>, <minor semi-axis> ]

		% Initialize radius
		r = NaN(1,2);

		% Confine `theta` to [0; 2*pi]
		pi2 = 2 * pi;
		while theta < 0
			theta = theta + pi2;
		end
		while theta > pi2
			theta = theta - pi2;
		end

		% Select major (`r(1)`) and minor (`r(2)`) semi-axes
		if theta <= pi/2 || theta >= 3/2 * pi
			r(1) = stdXpos;
		else
			r(1) = stdXneg;
		end
		if theta <= pi
			r(2) = stdYpos;
		else
			r(2) = stdYneg;
		end

		% Scale semi-axes to chisquare value
		r = r * interval;
	end % end of `semiaxes`

end % end of `error_potato`


function vrnc = custom_std(ctr, pts)
%custom_variance calculates a standard deviation based on a given center
vrnc = sqrt( sum((pts - ctr).^2) / size(pts, 1) );
end % end of `custom_var`
