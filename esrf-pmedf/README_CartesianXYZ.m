% Test of drawing a CARTESIAN MATRIX by various output utilities.
%
% Cartesian matrix z means a matrix
%	z = rows(z) x cols(z) = length(x) x length(y)
% and thus
%	z(i,j) = f(x(i),y(j))
%
% Then, in the usual output of Octave/Matlab, typing:
%        z        % position your monitor on its left-hand edge
%        rot90(z) % view z normally as you would expect from a cartesian system
%
% Result of this script: It presents the way how to work with cartesian
% systems.
%
% Author: Petr Mikulik
%
% History:
%	10. 8. 2010: Minor update.
%	 1. 6. 2009: Updated for Octave 3
%	13. 6. 2004: Original version for Octave 2

% Create a cartesian matrix.
x = linspace(0,100,101);
y = linspace(-20,70,51);
[yy, xx] = meshgrid(y,x); % NOTE THE REVERSED ORDER!
z = xx + yy/2;
% => now z is cartesian, try to have a look to elements z(ix,iy)

if exist('pm3d')
    colormap(pm3d(128));
end

% Display the cartesian matrix by "imagesc" command.
imagesc(x,y,z'); axis xy; colorbar
'Press ENTER to continue...'
pause;

% Note save it to gnuplot binary file, and let gnuplot draw it:
if exist('savegpbin')
    savegpbin(x, y, z, 'a.gpbin');
    if 0
	% gset 'xlabel "x_axis"'
	% gset 'ylabel "y_axis"'
	% gset autoscale fix
	% gset pm3d map
	% graw('splot "a.gpbin" binary\n');
    end
end

% Now with output via edf files:
h = pmedf_emptyHeader;
h = pmedf_putInHeader(h, 'DataType', 'FloatValue');
h = pmedf_putInHeader(h, 'PSize_1', x(2)-x(1));
h = pmedf_putInHeader(h, 'PSize_2', y(2)-y(1));
h = pmedf_write('a.edf', h, fliplr(z));
h = pmedf_writeWithRasterAxes('b.edf', h, fliplr(z), 'XrightYdown');
h = pmedf_writeWithRasterAxes('c.edf', h, z, 'XrightYup');

% eof
