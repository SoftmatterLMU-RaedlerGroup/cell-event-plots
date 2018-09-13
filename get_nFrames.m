function [ nFrames ] = get_nFrames( token, scale, specfile )
%GET_NFRAMES returns the number of frames in a measurement.
%
% Parameters:
%	token		the measurement token as string
%	scale		factor for time unit conversion
%	specfile	the file to be used for lookup (optional)
%
% Returns:
%	nFrames		the number of frames in the measurement
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
default_specfile = 'measurements_Tmax.csv';

if isnumeric(token)
	token = num2str(token);
elseif ~ischar(token)
	error('Measurement token must be numeric!')
end

if nargin < 2 || ~isnumeric(scale)
	scale = 1;
end

if nargin < 3 || ~ischar(specfile)
	specfile = default_specfile;
end

fid = fopen(specfile);
specs = textscan(fid, '%s %d %d', 'Delimiter',' \t\b', 'CommentStyle','%');
fclose(fid);

index = findStrInCstr(specs{1}, token, 'i');

if isempty(index)
	warning(['No frame number found for measurement token: ' token])
	nFrames = [];
else
	nFrames = scale * (specs{3}(index) - specs{2}(index) + 1);
end

end

