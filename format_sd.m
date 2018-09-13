function s = format_sd(val, err, width)
%FORMAT_SD formats a value with error with significant digits
%
% The returned string has the format:
%    12.34±0.05
%
% The error is rounded up to the first significant digit. If the first
% significant digit is 1 or 2, it is rounded up to the second significant
% digit. The value is rounded to the same number of decimal digits as the
% error.
%
% Arguments:
%    val     The value (left part)
%    err     The error (right part)
%    width   (optional) The minimum width of the returned string.
%            If width is a scalar, the returned string is padded with
%            spaces if it would be shorter than width.
%            If width is a vector of two elements, the first element is
%            the minimum width of the value, and the second element is
%            the minimum width of the error. Note that the returned string
%            has at least length `width(1)+width(2)+1` due to the ± sign.
%
% Returns:
%    s       Formatted string
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
pm_sign = char(177);

% Get optional widths
w_val = 0;
w_err = 0;
w_all = 0;
if nargin < 3
elseif length(width) == 1
	w_all = width;
elseif length(width) == 2
	w_val = width(1);
	w_err = width(2);
end

% Get order of magnitude of error
magn_err = floor(log10(err));
if err < 3 * 10^magn_err
	magn_err = magn_err - 1;
end

% Get numbers of decimal digits
n_dec = max(0, -magn_err);

% Round value and error (error is always rounded up)
val = round(val, n_dec);
err_round = round(err, n_dec);
if err_round < err
	err = err_round + 10^magn_err;
else
	err = err_round;
end

% Build formatted string
s = sprintf('%*.*f%s%*.*f', w_val, n_dec, val, pm_sign, w_err, n_dec, err);
if w_all > length(s)
	s = sprintf('%*s', w_all, s);
end
