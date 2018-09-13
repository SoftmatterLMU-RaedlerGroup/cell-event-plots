function pdf = logNormBiN(N, coord, params)
%logNormBiN computes a multimodal bivariate log-normal distribution
%
%	The bivariate log-normal density function is taken from:
%	Yue, Sheng: “The bivariate lognormal distribution for describing 
%	joint statistical properties of a multivariate storm event”,
%	2002, Environmetrics 13:811–819
%	http://onlinelibrary.wiley.com/doi/10.1002/env.483/abstract
%
% Input arguments:
%	N			number of distributions (multi-modes)
%
%	coord		coordinates at which to calculate the distribution; can be
%				given as numerical array or as cell array
%				* Numerical: give X-values in `coord(:,:,1)` and Y-values
%				in `coord(:,:,2)`
%				* Cell: give X-values as vector `coord{1}` and Y-values as
%				vector `coord{2}`
%
%	params		numerical vector of parameters to be used
%				Meaning of components: rho, muX, muY, sigmaX, sigmaY
%				1. a; scaling factor of distribution
%				2. rho; correlation coefficient with -1 <= rho <= 1
%				3. muX; logarithm of mean X value
%				4. muY; logarithm of mean Y value
%				5. sigmaX; logarithm of X standard deviation
%				6. sigmaY; logarithm of Y standard deviation
%				Repeat this scheme for multiple modes.
%
% Return value:
%	pdf			computed distribution function; numerical array of size
%				`size(length(Y), length(X))`
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
%%


%% Define constants
nParams = 6;

%% Validate input
if ~isscalar(N) || ~isnumeric(N) || mod(N,1) ~= 0 || N <= 0
	error('N must be a positive integer scalar.')
end
if ~isnumeric(params) || ~isvector(params) || any(~isfinite(params(:)))
	error('params must be a finite numeric vector.')
end
if isnumeric(coord)
	if any(~isfinite(coord(:)))
		error('No argument must be NaN or Infinite.')
	elseif size(coord, 3) ~= 2
		error('When giving coordinates as numerical array, X and Y must be given in third dimension.')
	end
	X = coord(:,:,1);
	Y = coord(:,:,2);

elseif iscell(coord)
	if numel(coord) ~= 2
		error('When giving coordinates as cell array, it must have two elements.')
	end

	X = coord{1};
	Y = coord{2};

	if ~isnumeric(X) || ~isnumeric(Y)
		error('Elements of coordinate cell array must be numeric.')
	elseif ~all(isfinite(X(:))) || ~all(isfinite(Y(:)))
		error('Only finite values are allowed for coordinates.')
	end
else
	error('Coordinates must be given either as numerical or as cell array.')
end
if length(params) ~= N * nParams
	error('For every partial distribution, %d parameters are needed.', nParams)
end

%% Prepare coordinates
X = unique(X);
Y = unique(Y);

if ~isrow(X)
	X = X.';
end
if ~iscolumn(Y)
	Y = Y.';
end

validX = X > 0;
validY = Y > 0;

%% Calculate bivariate lognormal distribution
pdfL = zeros(sum(validY), sum(validX));

for iN = 1:N
	tp = params( (iN - 1) * nParams + (1:nParams) );
	pdfL = pdfL + logNormBiv(X(validX), Y(validY), tp(1), tp(2), tp(3), tp(4), tp(5), tp(6));
end

pdf = zeros(length(Y), length(X));
pdf(validY, validX) = pdfL;


end % end of `logNormBiN`


function pdf = logNormBiv(x, y, a, rho, muX, muY, sigmaX, sigmaY)

	%% Test input for validity
	if rho < -1 || rho > 1
		error('`rho` must be between 1- and 1.')
	end

	%% Compute subterms
	% Compute factors for exponents
	expfacX = (reallog(x) - muX) / sigmaX;
	expfacY = (reallog(y) - muY) / sigmaY;

	% Compute exponent
	q = (expfacX.^2 + expfacY.^2 - 2 * rho .* expfacX .* expfacY) / (1-rho^2) / 2;

	%% Assemble probability density function
	pdf = exp(-q) ./ (2 * pi * sigmaX * sigmaY * sqrt(1-rho^2) * (x .* y)) * a;

end
