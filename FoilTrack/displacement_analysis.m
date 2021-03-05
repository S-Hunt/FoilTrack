%% displacement_analysis
%   Calcualates the displacements for the ImageAnalysis suit of scripts by
%   finding the minimum in the Sum Squared Differences of the selected
%   areas in the X-radiographs. 
%
%   Syntax:
%       [out, SSD_values, all_y] = displacement_analysis(type, Iref, Icomp,...
%                                          boxX, boxY, options...)
%       where:  'type'   - 'precrocess' or ??
%               'Iref'   - reference image 
%               'Icomp'  - comparison image 
%               'boxX'   - X position of comparison boxes
%               'boxY'   - Y position of comparison boxes
%       options:
%               'k', 'number'  - distance to look for comparison
%               'spot_removal' - remoce bright spots? on if present off if not)
%               'NaN_bg'       - make bg NaN? (on if present off if not)
%               'plot'         - display comparison figure (on if present off if not)
%               'y offset'     - move the search box by y (in the search direction).
%               'x offset'     - move the search box by x (perpendicual to search direction).
%               'warning'      - turn the warnings on. They are set to be off by default
%
%   See also:  ImageAnalysis, SelectBoxes

%   $ version: 4.1 $ 10th April 2019 $
%
%    intially developed by Li Li, 2002
%    further developed by Simon Hunt 2006 - 2008
%    made generic by Andrew Walker 1/4/2009
%    further additions Simon Hunt 2009 - 2019
%       - for details of changes between versions see end of file.

function [out, SSD_values, all_y] = displacement_analysis_dev(type, I1, varargin)% I2, boxX, boxY, k, spot_removal, NaN_bg, varargin)

%% initial setup

%check there are enough input arguments.
if strcmpi(type, 'preprocess') && nargin < 6
    error('Not enough input arguments')
elseif nargin < 6
    error('Not enough input arguments')
end

%extract image 2 from varargin
% the position of image2 is fixed in the syntax.
if nargin >= 3
    I2 = varargin{1};
    boxX = varargin{2};
    boxY = varargin{3};
end

% Set defaults for options.
k = 5; %search distance
spot_removal = 0;
NaN_bg = 0;
scale = 0;
plots_on = 0; %plotting on or off.
Y_offset = 0;
X_offset = 0;

%turn off warnings so that mypolyfit does not does not keep spitting them out.
ws = warning('off','all');

% Process the optional arguments overiding defaults
iarg = 4 ;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case {'k', 'search'}
            k = varargin{iarg+1};
            iarg = iarg + 2;
        case 'scale'
            scale = 1;
            iarg = iarg + 1;
        case 'spot removal'
            spot_removal = 1;
            iarg = iarg + 1;
        case 'nan bg'
            NaN_bg = 1;
            iarg = iarg + 1;
        case 'y offset'
             Y_offset = varargin{iarg+1};
            iarg = iarg + 2;
        case 'x offset'
            X_offset = varargin{iarg+1};
            iarg = iarg + 2;
        case 'plot'
            plots_on = 1;%varargin{iarg+1};
            iarg = iarg + 1;
        case 'warning'
            ws = warning(varargin{iarg+1},'all');
            iarg = iarg + 1;
        otherwise
            error(['Unknown option: ' varargin{iarg}]) ;
    end
end

%set variables that are used here.
[h1 ~] = size(I1);
number_boxes = size(boxX,1);

%%set variables that used to be needed and might be needed again...
%Y_offset = 0;
%X_offset = 0;

%% changes string variable into numbers which speeds up program
%if strcmpi(spot_removal, 'yes') == 1
%    spot_removal = 1;
%elseif strcmpi(spot_removal, 'no') == 1
%    spot_removal = 0;
%end

%array setup
SSD_values = NaN(number_boxes,2*k+1);
minima = zeros(1,number_boxes);
all_y = zeros(number_boxes,2*k+1);


%% find global minimum (find offsets)
if strcmpi(type,'preprocess') ~= 1 && (strcmpi(Y_offset, 'find') == 1 || strcmpi(X_offset, 'find') == 1)
    
    %define area which to match between the images
    Xmatch = [min(boxX(:)), max(boxX(:))];
    Ymatch = [min(boxY(:)), max(boxY(:))];
    
    %calculate the offsets
    [X_offset , Y_offset , ~ , ~] = displacement_analysis_PosCorr(I1 , I2 , Xmatch , Ymatch);
    
    disp( [num2str(X_offset), ' ', num2str(Y_offset)] )
end



%% misfit minimum
%finds the minimum misfit for each of the foils in the y direction using a
%sum of squares difference (ssd) or cross correlation method.

% position change - box by box.
for n = 1 : number_boxes

    shift_y = Y_offset;
    shift_yold = .1;
    shift_yold2 = -.1;

    % get all data might ever need for area1.
    k_down = k;
    k_up = k;
    if boxY(n,1) - k <= 0
        k_down = boxY(n,1) - 1;
    elseif boxY(n,2) + k > h1
        k_up = h1 - boxY(n,2) - 1;
    end
    area1extended = (I1(boxY(n,1)-k_down:boxY(n,2)+k_up, boxX(n,1):boxX(n,2)));

    if strcmpi(type,'preprocess') == 1; %preprocesses the selected data from image1.
        if spot_removal == 1
            area1extended = spot_remove(area1extended);
        end
        if NaN_bg == 1
            area1extended = bg_makeNaN(area1extended); %original function
        end
        I1(boxY(n,1)-k_down:boxY(n,2)+k_up, boxX(n,1):boxX(n,2)) = area1extended;

        out = I1;
        
        %end of preprocessing
        
    else %calculates the displacement
        area1 = area1extended(k_down+1:end-k_up, :, ones(1,k_up+k_down+1)); %cuts area1 to size and replicates area1 in third dimension

        %calcuates the offsets
        area2 = zeros(size(area1,1), size(area1,2), k_up+k_down+1);
        vals = zeros(size(area1,1), size(area1,2), k_up+k_down+1);
        
        while shift_y ~= shift_yold;
            if boxY(n,1)-k_down+shift_y <= 0
                shift_y = k_down - boxY(n,1) + 1;
            elseif boxY(n,2)+k_up+shift_y > h1
                shift_y = h1 - (boxY(n,2) + k_up);
            end

            area2all = I2(round(boxY(n,1)-k_down+shift_y):round(boxY(n,2)+k_up+shift_y), boxX(n,1)+X_offset:boxX(n,2)+X_offset);
            if spot_removal == 1
                area2all = spot_remove(area2all);
                % NaN background creation is not needed for area2 because of NaN mask to area1
            end
            if NaN_bg == 1
                % area2all = stripe_remove(area2all); %this is part of the bg_removal and is needed here for crapy data.
                
                area1extended = bg_makeNaN(area1extended); %original function
            end
            if scale == 1
                area2all = scale_image(area1extended, area2all);
            end

            %populate image 2 array for comparison.
            extent = size(area1);
            for i = 1 : k_up+k_down+1
                area2(:,:,i) = area2all(i:i+extent(1)-1, :);
                % imagesc(area2(:,:,i)), colorbar, pause
 
                %These lines were an attempt to use the Pearson correlation coefficient instead of the SSD to match the foils. 
                %It did not work because the variation in the Pearson value is very small compared to that of the SSD. 
                % area1linear = area1extended(i:i+extent(1)-1, :);
                % area1linear = area1linear(:);
                % area2linear = area2all(i:i+extent(1)-1, :);
                % area2linear = area2linear(:);
                % ij = find(area1linear == NaN)
                % area1linear(ij) = NaN;
                % area2linear(ij) = NaN;
                % ij = find(area2linear == NaN)
                % area1linear(ij) = NaN;
                % area2linear(ij) = NaN;
                % pearson(n,i) = corr(area1linear,area2linear,'rows','pairwise', 'type','Kendall');

            end

            % these lines are a vectorised version of those above - but are significantly slower than the looped version.
            % There were sets of zero arrays for each of the new elements which made this a bit faster but it did not help. 
            % box_w = boxX(n,2) - boxX(n,1) + 1;
            % box_h = boxY(n,2) - boxY(n,1) + 1;
            % [xm, ym, zm] = meshgrid(1 : box_h+2*k : box_w*(box_h+2*k), 0:box_h-1 , 0:2*k);
            % area2 = area2all(xm + ym + zm);
            
            vals = (area1-area2).^2;
            if spot_removal == 1 | strcmpi(NaN_bg,'yes') == 1  %checks criteria to see if NaN values are likely to be in the matrices.
                nan_totals = sum(sum(~isnan(vals),1),2);
                vals(isnan(vals)) = 0;
                SSD_values(n,k-k_down+1: k+k_up+1) = sum(sum(vals,1),2) ./ nan_totals; %divide by number of entries in each column (normalised for NaN removal
            else
                SSD_values(n,k-k_down+1: k+k_up+1) = sum(sum(vals,1),2); %this is faster if no NaN values are present
            end
            %         figure(2), plot(values(n,:)), pause

            shift_yold2 = shift_yold;
            shift_yold = shift_y;

            miny = find(SSD_values(n,:) == min(min(SSD_values(n,:))));
            if miny == 1 | miny == 2 | miny == 3
                shift_y = shift_y - k;
            elseif miny == 2*k+1 | miny == 2*k | miny == 2*k-1
                shift_y = shift_y + k;
                %         else
                %             shift_yold = shift_y;
            end
            %shift_y
            %checks to keep offsets inside the image
            if boxY(n,1) + shift_y - k <= 0
                break
            elseif boxY(n,2) + shift_y +k >= h1
                break
            elseif shift_yold2 == shift_y
                break
            end

        end
        ij = -k + shift_y : k + shift_y;

        %keep copy of y values for plotting and position correlations.
        all_y(n,:) = (boxY(n,1) + boxY(n,2))/2 + ij;

        switch type

            case 'poly' %polynomial fit to the missfits

                miss_poly6 = mypolyfit(ij,SSD_values(n,:),6);

                %minimum in polynomial by differentiation.
                poly_deriv = mypolyder(miss_poly6);
                no_slope = roots(poly_deriv);
                roots_of_interest = [];
                j=1;
                for i = 1 : length(no_slope)
                    if isreal(no_slope(i)) == 1 && no_slope(i) <= k+shift_y && no_slope(i) >= -k+shift_y
                        roots_of_interest(j) = no_slope(i);
                        j=j+1;
                    end
                end
                if length(roots_of_interest) > 1
                    depth = polyval(miss_poly6,roots_of_interest);
                    actual = find(min(depth)==depth,1,'last');
                    minima(n) = roots_of_interest(actual);
                elseif isempty(roots_of_interest) == 1
                    minima(n) = 10;
                else
                    minima(n) = roots_of_interest;
                end

                %calcualte profile for plotting
                if plots_on == 1
                    fit_plot(1,:) = ij;%SSD_values(n,:);
                    fit_plot(2,:) = polyval(miss_poly6, ij);
                end
                
            case 'spline'
                resolution = 0.01;

                yi = zeros(1,2*k/resolution+1);
                
                if exist('fitting','var') == 0
                    fitting = zeros(1, length(yi));
                
                    fitvalues = zeros(number_boxes,length(yi));
                    allxi = zeros(number_boxes,2*k/resolution+1);
                end
                yi = -k + shift_y : resolution : k + shift_y;

                allxi(n,:) = yi;
                %spline minimum to the misfits.

                %             x = zeros(2*k+1);
                cols = length(yi);

                % 1 dimsonsional data interploation to find the minimum misfit. Fits a cubic spline to the points and then
                % calculates the minimum misfit of the 'foils' to a maximum reolution as
                % defined by the variable 'resolution', which is the point spacing of the
                % points caclulated for the cubic fit.
                fitting = spline(ij,SSD_values(n,:),yi);
                fitvalues(n,1:cols) = fitting(1:cols);

                % minimum of the fitted polynomial function displacement of the line and its displacement from the original
                min_displacement = find(fitvalues(n,:) == min(fitvalues(n,:)));

                minima(n) = yi(min_displacement(1));

                %calcualte profile for plotting
                if plots_on == 1
                    fit_plot(1,:) = yi;%SSD_values(n,:);
                    fit_plot(2,:) = fitting;
                end
                
                
                
            case 'gauss'
                error('This option needs fixing')
                %gaussian munimum to the misfits

                init = [1,1,1,1,1]; %this is needed because of a flaw in the fitgaussian program which was not initialising properly without it.
                SSD_values = double(SSD_values); %this is also needed for the same reason as above
                [f, eX, err, it] = fitgaussian(ij, SSD_values(n,:));
                all_gauss(n,:) = f;

                %eX(3) is divided by 2 first because the Gaussian fitting routine
                %returns 2*sigma as the third value.
                % displacement_error(n) = eX(3)/2/sqrt(2*k+1);
                
                
            case 'minimum'
                
                [~, b] = min(SSD_values(n,:));
                minima(n) = ij(b);

                %profile for plotting
                if plots_on == 1
                    fit_plot(1,:) = ij;
                    fit_plot(2,:) = SSD_values;
                end

        end

        if plots_on == 1 %plot the output.... if required.
                title_str = ['Box ',num2str(n),'; ', type, ' minimum = ', num2str(minima(n))];
                if size(area1extended,2) ~= 1
                    subplot(1,3,1); imagesc(area1extended);
                    subplot(1,3,2); imagesc(area2all); %set(gca, 'clim', [0 256])
                    subplot(1,3,3); plot(ij,SSD_values(n,:),'.',fit_plot(1,:), fit_plot(2,:),'r:');
                else
                    subplot(2,2,1); plot(area1extended);xlim([0,length(area1extended)]);
                    rectangle('Position',[k_down,min(area1extended),boxY(2)-boxY(1),max(area1extended)-min(area1extended)],'EdgeColor','m')
                    subplot(2,2,3); plot(area2all); hold on, xlim([0,length(area1extended)]);%set(gca, 'clim', [0 256])
                    rectangle('Position',[minima(n)-Y_offset+k_down,min(area2(:)),boxY(2)-boxY(1),max(area2(:))-min(area2(:))],'EdgeColor','m')
                    plot((1:length(area1extended))+minima(n)-Y_offset, area1extended,'r')
                    
                    subplot(2,2,[2,4]); plot(ij,SSD_values(n,:),'.',fit_plot(1,:), fit_plot(2,:),'r:');
                    xlim([min(ij),max(ij)]);
                    
                end
                suptitle(title_str);
%                 keyboard; 
                pause;% clf
        end
        
        out = minima;
    end

end %loop over number of boxes


warning(ws);


% %plots the normalised least squares between the data sets. and the fitted
% %polynomial.
% % if fig_offsets == 'on'
%     fig_col = 10;
%     fig_row = length(boxX) / fig_col;
% %    figure(fig); fig = fig + 1; clf;
%     for n = 1 : number_boxes
%         title_str = strcat('box ',num2str(n));
%
%         subplot(fig_row, fig_col, n);
%         plot(all_x(n,:),values(n,:),'.',allxi(n,:),poly_calc6(n,:),'r-', [minima(n), minima(n)], [1,0],'b'); axis([all_x(n,1) all_x(n,2*k+1) 0 1]); title(title_str, 'FontWeight', 'bold'); hold on;
%         %subplot(fig_row, fig_col, n);
%         %plot(allxi(n,:),poly_calc6(n,:),'r-',allxi(n,:),poly_calc4(n,:),'r:'); %plots the gaussian fit to the data.
%         %pause; clf
%     end
%     pause; clf
% % end


%
%         subplot(1,2,1); imagesc(I1); colormap(gray)
%         set(gca, 'clim', [0 256])
%         for n = 1 : number_boxes
%             x = boxX(n,1);
%             y = boxY(n,1);
%             h = boxX(n,2)-boxX(n,1);
%             w = boxY(n,2)-boxY(n,1);
%             rectangle('position',[x,y,h,w],'EdgeColor','r');
%         end
%
%         subplot(1,2,2); imagesc(I2); set(gca, 'clim', [0 256])
%         for n = 1 : number_boxes
%             x = boxX(n,1);
%             y = boxY(n,1)+minima(n);
%             h = boxX(n,2)-boxX(n,1);
%             w = boxY(n,2)-boxY(n,1);
%             rectangle('position',[x,y,h,w],'EdgeColor','r');
%         end
%         pause (.2);


end %displacement analysis.


%% SUBFUNCTIONS

function p = mypolyfit(x,y,n) %[p,S,mu] = mypolyfit(x,y,n)
%POLYFIT Fit polynomial to data.
% This the an edited version of the Matlab version of polyfit: Revision: 5.17.4.10

% Any warnings that the polynomial is poor has been removed --
% because it slows the function down significantly.
% The warnings are now at begining and end of displacement_analysis file
% Extra code has been removed.

x = x(:);
y = y(:);

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1);%,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
% ws = warning('off','all');
p = R\(Q'*y);    % Same as p = V\y;
% warning(ws);

% r = y - V*p;
p = p.';          % Polynomial coefficients are row vectors by convention.

end


function out = mypolyder(in)
powers = length(in)-1:-1:1;
out = in(1:end-1).*powers;
% out = out(1:end-1);
end



function Im = spot_remove(Im)

Im2 = sort(Im(:),1,'descend');
diffs = Im2(1:10) - Im2(2:11);
large = find(diffs > 3,10,'last');
if isempty(large) ~= 1;
    remove = find(Im >= min(Im2(large)));
    Im(remove) = NaN;
end

end %spot_remove

% function I_clean = bg_makeNaN(Im)
% 
% alpha = .5;
% [I_bar I_std] = nanstats(Im,1);
% bg_coeffs = splinefit(double(1:length(I_bar)), double(I_bar-alpha*I_std),7,3);
% bg_spline = ppval(bg_coeffs,1:length(I_bar));
% Im_bg = repmat(round(bg_spline),size(Im,1),1);
% Cut = find(Im > Im_bg);
% I_clean = Im;
% I_clean(Cut) = NaN;%Im_bg(Cut);
% 
% subplot(1,2,1),imagesc(Im),
% subplot(1,2,2), imagesc(I_clean)
% 
% keyboard
% end %bg_makeNaN


function I_clean = bg_makeNaN(Im)

Im_new1 = stripe_remove(Im);
Im_new = smoothn(Im_new1);

alpha = .5;
[I_bar I_std] = nanstats(Im_new,1);
bg_coeffs = splinefit(double(1:length(I_bar)), double(I_bar-alpha*I_std),7,3);
bg_spline = ppval(bg_coeffs,1:length(I_bar));
Im_bg = repmat(round(bg_spline),size(Im_new,1),1);
Cut = find(Im_new > Im_bg);
I_clean = Im_new1;
I_clean(Cut) = NaN;%Im_bg(Cut);

if (0)
    % range = [min(Im_new1(:)), max(Im_new1(:))] ;
    subplot(1,3,1),surf(Im_new1), %clim(range)
    subplot(1,3,2),surf(Im_new),% clim(range)
    subplot(1,3,3), surf(I_clean),  %clim(range)
    pause
end
end %bg_makeNaN

function part2scaled = scale_image(part1, part2)

% hist(part1), pause
% hist(part2), pause

%scale using means and standard deviations
std1 = std(part1(:));
mean1 = mean(part1(:));
std2 = std(part2(:));
mean2 = mean(part2(:));

part2scaled = (part2-mean2)./std2.*std1 + mean1;


end %scale image

function [x_bar x_std] = nanstats(x,dim)
%
% Combined copy of nanmean and nanstd after the scripts of:
%    Jan Gl?scher, Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
% Optimsed for use here by Simon Hunt

%NANMEAN
%numerator
nans = isnan(x);

% denominator
count = size(x,dim) - sum(nans,dim);

xb = x;
xb(isnan(x)) = 0;

% denominator
% count = size(x,dim) - sum(nans,dim);
% Protect against all NaNs in one dimension
i = find(count==0);
count(i) = ones(size(i));

x_bar = sum(xb,dim)./count;
x_bar(i) = i + NaN;

% FORMAT: Y = NANSTD(X,DIM,FLAG)

% Find NaNs in x and nanmean(x)
% nans = isnan(x);
% avg = x_bar; %nanmean(x,dim);

% create array indicating number of element
% of x in dimension DIM (needed for subtraction of mean)
tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

% remove mean
xs = x - repmat(x_bar,tile);

% count = size(xs,dim) - sum(nans,dim);

% Replace NaNs with zeros.
xs(isnan(xs)) = 0;


% Protect against a  all NaNs in one dimension
i = find(count==0);

% if flag == 0
	x_std = sqrt(sum(xs.*xs,dim)./max(count-1,1));
% else
% 	x_std = sqrt(sum(xs.*xs,dim)./max(count,1));
% end
x_std(i) = i + NaN;


end %nanstats



%% versions
%   - 4.1   -- 10th April 2019
%       - Added code to revive offset in image position
%   - 4.0   -- 10th April 2019
%       - Revised the syntax for using the script.
%       - Fixed a serious bug that changes the box of interest if it is too near the ends of the data. Fixed it. 
%   - 3.9.2 -- 17 May 2016
%       - Fixed bug regarding calculating number of boxes.
%   - 3.9.1 -- 9 May 2016
%       - Updated help documentation at the top of the file.
%       - Added the void lines for pearson correlation coefficient (for record keeping).
%       - Added some void lines for vectorising the code which are slower than those used (for record keeping).
%       - Added a line making vals as a set of zeros before it is used which sped the code up quite abit.
%   - 3.9 -- 2 Sept 2014
%       - Added output for SSD and the position at which the SSD was calculated. The script does this 
%       by default and decides what to save in the ImageAnalysis script.
%   - 3.8 -- 25 Jan 2012
%       - threshold the images as a way of reducing the SSD and error in displacement - area outside
%       foil becomes NaN.
%       - moved the image1 spot removal and NaN background into a separate part of the loop so that it is only done once per file.
%   - v 3.7.4 -- 25 Jan 2012
%       - removed more lines from mypolyfit to make it run faster.
%   - v 3.7.3 -- 9 Jan 2012
%       - added lines to make sure area1extended is always within the image.
%   - v 3.7.2 -- 13 Oct 2011
%       - changed polyder to mypolyder so runs faster.
%       - vectorised more of the displacement calculation
%       - changed spot removal so that it is done for both images once per box -- adds consistency to analysis
%       - changed loop finding minima to make sure stops -- added shift_yold2
%   - v3.7 -- 13 Feb 2011
%		- add ability to remove bright spots from individual sub areas -- replace with NaN and don't count them
%       - removed some of the redundant code
%   - v3.6a -- 29/Jan/2011
%       - reorder and optimisation of SSD routine to reduce loops and overhead time.
%   - v3.6 -- 27/Jan/2011
%       - added ability to reduce seleced box to reduce noise (hopefully) -rows, 2x2 and 3x3
%       - removed the header creation routine into ImageAnalysis (v1.4+)
%       - analysis variables are now fed to the code by ImageAnalysis (v1.4+)
%   - v 3.5
%       - 1st numbered version, 1 Jan 2011
%       - Added mypolyfit to speed up function. Mypolyfit differs from polyfit by removal of error
%   catching.
%   - v3
%       - Andrew Walker's generic version.
%   - v2
%       - Simon Hunt's developments (code in Hunt, thesis, 2009)
%   - v1
%       - Li's version of the code
