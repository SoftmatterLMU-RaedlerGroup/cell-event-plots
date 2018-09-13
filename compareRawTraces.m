% This script expects:
%	mkr_dir1
%	mkr_dir2
%	file_list
%	mkr_name1
%	mkr_name2
%	out_dir
%
% To convert the resulting multi-page PS-files to PDF, use:
% for i in *.ps; do ps2pdf -dAutoRotatePages=/None $i; done
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
time_scale = 1/6;
color1 = 'blue';
color2 = 'red';

if ~exist('out_dir', 'var')
	out_dir = 'out';
end
if ~exist(out_dir, 'dir')
	[~,~,~] = mkdir(out_dir);
end

fh = figure('Visible', 'off');

for i_file = 1:length(file_list)
	path1 = fullfile(mkr_dir1, file_list{i_file});
	path2 = fullfile(mkr_dir2, file_list{i_file});
	[~, filename, ~] = fileparts(file_list{i_file});

	if ~exist(path1, 'file') || ~exist(path2, 'file')
		warning('File not found: %s', filename);
		continue
	end

	data1 = csvread(path1);
	data2 = csvread(path2);

	if size(data1, 1) ~= size(data2, 1)
		warning('Data lengths differ in %s.', file_list{i_file});
		continue
	end
	if size(data1, 2) ~= size(data2, 2)
		warning('Trace numbers differ in %s.', file_list{i_file});
		continue
	end

	t_vec = data1(:,1) * time_scale;

	for i_trace = 2:size(data1, 2)
		idx_trace = i_trace - 1;

		clf(fh)
		ax = axes(fh);

		yyaxis(ax, 'left');
		plot(ax, t_vec, data1(:,i_trace), '-', 'Color', color1)
		ax.YColor = color1;
		ylabel(ax, [ mkr_name1 ' Fluorescence Intensity [a.u.]' ]);

		yyaxis(ax, 'right');
		plot(ax, t_vec, data2(:,i_trace), '-', 'Color', color2)
		ax.YColor = color2;
		ylabel(ax, [ mkr_name2 ' Fluorescence Intensity [a.u.]' ]);

		xlabel(ax, 'Time [h]')
		title(ax, sprintf('File: %s, Cell: %03d', filename, idx_trace), 'Interpreter','none')
% 		legend(ax, {mkr_name1, mkr_name2}, 'Location', 'best');

		fh.PaperPositionMode = 'manual';
		fh.PaperUnits = 'centimeters';
		fh.PaperPosition = [0 0 9 9];
		fh.PaperSize = [9 9];
		
		print(fh, fullfile(out_dir, sprintf('%s_compare_%s-%s.ps', ...
			filename, mkr_name1, mkr_name2)), '-dpsc', '-append')
	end
end

close(fh)
