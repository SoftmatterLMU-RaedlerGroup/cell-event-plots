function pdf = logNormBi( rho, muXlog, sigmaXlog, muYlog, sigmaYlog, X, Y )
%logNormBi computes the bivariate lognormal probability distribution
%
%	For more information, see:
%	Yue, Sheng: “The bivariate lognormal distribution for describing 
%	joint statistical properties of a multivariate storm event”,
%	2002, Environmetrics 13:811–819
%	http://onlinelibrary.wiley.com/doi/10.1002/env.483/abstract
%
%
% Input values:
%	rho				correlation coefficient between x- and y-distribution
%	muXlog			mean x-value of distribution
%	sigmaXlog		standard deviation of distribution in x-direction
%	muYlog			mean y-value of distribution
%	sigmaYlog		standard deviation of distribution in y-direction
%	X				x-values for which to calculate distribution
%	Y				y-values for which to calculate distribution
%
% Output values:
%	pdf				log-normal probability density function
%
% `pdf` has the same size as `X .* Y`.
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

%% Validate input
if ~isscalar(rho) || ~isnumeric(rho)
	error('`rho` must be a scalar numeric.')
end

if ~isscalar(muXlog) || ~isnumeric(muXlog)
	error('`muX` must be a scalar numeric.')
end

if ~isscalar(sigmaXlog) || ~isnumeric(sigmaXlog)
	error('`sigmaX` must be a scalar numeric.')
end

if ~isscalar(muYlog) || ~isnumeric(muYlog)
	error('`muY` must be a scalar numeric.')
end

if ~isscalar(sigmaYlog) || ~isnumeric(sigmaYlog)
	error('`sigmaY` must be a scalar numeric.')
end

if rho < -1 || rho > 1
	error('`rho` must be between 1- and 1.')
end

if any(X(:) < 0) || any(Y(:) < 0)
	error('`X` and `Y` must be larger than 0.')
end

% Eliminate zeros from input vectors because of logarithm
idx = X == 0;
if any(idx)
	X(idx) = eps;
end

idx = Y == 0;
if any(idx)
	Y(idx) = eps;
end

%% Compute subterms
% Compute factors for exponents
expfacX = (reallog(X) - muXlog) / sigmaXlog;
expfacY = (reallog(Y) - muYlog) / sigmaYlog;

% Compute exponent
q = (expfacX.^2 + expfacY.^2 - 2 * rho .* expfacX .* expfacY) / (1-rho^2) / 2;

%% Assemble probability density function
pdf = exp(-q) ./ (2 * pi * sigmaXlog * sigmaYlog * sqrt(1-rho^2) * (X .* Y));

end
