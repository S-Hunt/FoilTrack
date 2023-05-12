function [NewBoxesY, NewBoxesX] = SelectBoxes_position(I1, varargin)
% SelectBoxes_position. Position the boxes over the image by the selected routine. 
%   syntax:
%   [NewBoxesY, NewBoxesX] = SelectBoxes_position(Image,...
%           boxY, boxX, box_fit_type, search_distance, empty_boxes,
%           HeightBox);
%   ['options'] = SelectBoxes_position('posibilities')
%
%   See also: SelectBoxes, ImageAnalysis

% version 1.2   Simon Hunt 2017 - June 2018
% Used to be part of Select Boxes but was removed when SelectBoxes was
% rewritten (version 2).


%Return the list of options if called
if strcmpi(I1, 'posibilities')
    %Lists the fit options avaliable in this file.
    % If I1 is 'possibilities' then return the list.
    
	% call SelectBoxes_getCases on this file and list all options in the switch case.
    %caller = [mfilename, '.m'];
    fit_opts = getCases;%SelectBoxes_getCases(caller, 115);
            
    NewBoxesY = fit_opts;
    NewBoxesX = [];
    return
    
end
    
%process varargin
boxY = varargin{1};
boxX = varargin{2};
box_fit_type = varargin{3};
search = varargin{4};
empty_boxes = varargin{5};
HeightBox = varargin{6};

%FIX ME: replace this input list with a 'proper' parseable list which
%could include options for plotting the fitted boxes.

%FIX ME: Need to be able to define the number of halfwidths to use.


%make sure HeightBox has one element for each box
if length(HeightBox) ~= length(boxX)
   HeightBox = ones(length(boxX),1)*HeightBox;
end

%make sure box_fit_type has one element for each box
if length(box_fit_type) ~= length(boxX)
   box_fit_type = repmat({box_fit_type},length(boxX),1);
end

%gets sizes of the image
[h1 ~] = size(I1);

number_boxes = size(boxY,1);

sigmas = 2.5; %number of standard deviations to include (if required)


%% Position the box
% centres the chosen areas over the foil - in the manner descrided in the switch setting
% goes through one box at a time.

for n = 1 : number_boxes
%     if boxY(n,2) < boxY(n,1)
%         
%         [boxY(n,2), boxY(n,1)] = swap(boxY(n,2), boxY(n,1));
%         disp(["swap", boxY(n,1), boxY(n,2)])
%     end
    extra_space = search - round(HeightBox(n)/2);
    top_of_area = round(boxY(n,1)) - extra_space;
    bot_of_area = round(boxY(n,2)) + extra_space;
    
%     if bot_of_area < top_of_area
%         
%         [bot_of_area, top_of_area] = swap(bot_of_area, top_of_area);
%         disp(["swap", top_of_area, bot_of_area])
%     end

    if ismember(n, empty_boxes) == 1 %is the box empty?
        position_min(n) = mean(boxY(n,:));
        box_range(n) = (boxY(n,2)-boxY(n,1))/2;
        
    else       
        
        %keep search area with in the image. 
        if top_of_area <= 0
            top_of_area = 1;
        end
        if bot_of_area >= h1
            bot_of_area = h1;
        end

        area = I1(top_of_area : bot_of_area, boxX(n,1):boxX(n,2));

        %setup numbers used in switch cases
        profile = double(nanmean(area,2));
        X =  top_of_area:bot_of_area;
        posminimum = find(min(profile) == profile,1); %used in a number of the options so called here once. 

        profile_range  = max(profile) - min(profile);
        if profile_range == 0 %there is nothing to fit and the fitting routines will spit out an error or two.
            box_fit_type{n} = 'None';
        end
        
        %position box based on fitting type chosen. 
        switch box_fit_type{n}  %    fit_opts = {'Minimum intensity'; 'Maximum slope'; 'Fit polynomial'; 'Fit Gaussian'; 'None'};
            case 'Minimum intensity'
                position_min(n) = X(posminimum(1));
                box_range(n) = HeightBox(n)/2;
                box_width_type = 1;
                
            case 'Minimum intensity interpolated (spline)'
                splin = interpolant_spline(X, profile, 4);
                position_min(n) = spline_values(splin,'min');
                
                box_range(n) = HeightBox(n)/2;
                box_width_type = 1;
                               
            case 'Minimum intensity interpolated (poly)'
                range = 7;
                x_range = posminimum-range : posminimum+range;
                x_range(x_range<1) = [];
                x_range(x_range>length(X)) = [];
%                 x_vals = X(x_range)';
%                 y_vals = profile(x_range);
%                 x_vals = x_range';
                x_vals = x_range'-posminimum;
                y_vals = profile(x_range);
                p = polyfit(x_vals,y_vals,3);
                
                if 1 %plot the fitted polnomial
                    if ~exist('another')
                        another.fig = figure('Position', [ 884   362   560   420],'Name','Minima Locations');
                        xlabel('Position')
                        ylabel('Mean Intensity')
                        title(['Box ',num2str(n)])
                    else
%                         figure(another);
                        delete(another.plot);
                    end
                    pval = polyval(p,x_vals(1):.1:x_vals(end));
                    another.plot = plot(X,profile,'b.', X(x_vals(1)+posminimum):.1:X(x_vals(end)+posminimum),pval,'r-');
%                     pause(.5)
                    
                    if n == max(setdiff(1:number_boxes, empty_boxes))
                        close(findobj('type','figure','Name','Minima Locations'))
                    end
                end
                
                k = polyder(p); %differentiate polynomial
                r = roots(k);   %find the roots.
                r = r( r>x_vals(1) & r<x_vals(end)); %discard the roots outside of the search range.
                mins = polyval(p, r); %calculate y values for roots in range -- needed incase length(r) > 1
                position_min(n) = X(posminimum) + r(min(mins)==mins); 
%                 position_min(n) = r( r>x_vals(1) & r<x_vals(end) ); %using logical indexing which is faster then find 
%                 position_min(n) = X(posminimum) + r( r>x_vals(1) & r<x_vals(end) ); %using logical indexing which is faster then find 
                box_range(n) = HeightBox(n)/2;
                box_width_type = 1;
 
            case {'Maximum slope'}
                gradient = profile(1:end-1) - profile(2:end);
                if boxY(n,2) < h1/2
                    max_abs_grad = find(max(gradient) == gradient);
                else
                    max_abs_grad = find(min(gradient) == gradient);
                end
                position_min(n) = X(max_abs_grad(1));
                box_range(n) = HeightBox(n)/2;
                box_width_type = 1;
            
            case {'Maximum slope interpolated (spline)'}
                
                [splin] = interpolant_spline(X, profile, 10);  %was 4
                sp_grad = ppdiff(splin);
                
                if boxY(n,2) < h1/2 %+100 %FIXME This was here to made sure that Cu_02 processed well.
                    direction = 'min';
                else
                    direction = 'max';
                end
                position_min(n) = spline_values(sp_grad,direction);
                box_range(n) = HeightBox(n)/2;
                box_width_type = 1;
                
            case {'Minimum slope'}
                gradient = profile(1:end-1) - profile(2:end);
                if boxY(n,2) < h1/2
                    max_abs_grad = find(min(gradient) == gradient);
                else
                    max_abs_grad = find(max(gradient) == gradient);
                end
                position_min(n) = X(max_abs_grad(1));
                box_range(n) = HeightBox(n)/2;
                box_width_type = 1;
            
            
            case {'Minimum slope interpolated (spline)'}
                
                [splin] = interpolant_spline(X, profile, 6);  %was 4
                sp_grad = ppdiff(splin);
                
                if boxY(n,2) < h1/2
                    direction = 'max';
                else
                    direction = 'min';
                end
                position_min(n) = spline_values(sp_grad,direction);
                box_range(n) = HeightBox(n)/2;
                box_width_type = 1;
                
            case {'Fit Gaussian on gradient'}
                warning('The type of Gaussian needs to be specified; it is currently a modified gaussian peak')
                slope = (profile(end)-profile(1))/length(X);
                offset = profile(1) - X(1)*slope;
                coeffs = double([offset slope max(profile)-min(profile) X(find(profile == min(profile),1)) 10 1]);
                coeffs = fminsearch(@fit_gaussian_and_bg,coeffs,[], X, profile);

                position_min(n) = coeffs(4);
                if abs(coeffs(5)) >= 1.5*HeightBox(n)
                    box_range(n) = 1.5*HeightBox(n);
                else
                    box_range(n) = abs(coeffs(5));
                end
                coeffs_keep(n,:) = coeffs;
                box_width_type = 2;

                
            case {'Fit Gaussian over step'}
                warning('The type of Gaussian needs to be specified; it is currently a modified gaussian peak')
                bg1 = profile(1);
                bg2 = profile(end);
                coeffs = double([bg1 bg2 max(profile)-min(profile) X(find(profile == min(profile),1)) 10 1]);
                coeffs = fminsearch(@fit_gaussian_and_step,coeffs,[], X, profile);

                position_min(n) = coeffs(4);
                box_range(n) = abs(coeffs(5));
                coeffs_keep(n,:) = coeffs;
                box_width_type = 2;
                
            case {'Fit sigmoidal peak on step'}    
                coeffs = double([profile(1) profile(end) max(profile)-min(profile) X(find(profile == min(profile),1)) 10 10 10]);
                coeffs = fminsearch(@fit_sigmoids,coeffs,[], X, profile);

                position_min(n) = coeffs(4);
                box_range(n) = abs(coeffs(5));
%                 sigmas = 4; %number of standard deviations to include
                coeffs_keep(n,:) = coeffs;       
                box_width_type = 2;       
%      
%             case {'Fit spline on gradient'}
%                 coeffs = [X(find(profile == min(profile),1)) 15];
%                
%                 fitted_vals = fminsearch(@fit_spline_onbg, coeffs, [], X,profile);
%                 
%                 position_min(n) = fitted_vals(1);
%                 box_range = coeffs(2);
%                 box_width_type = 2;
%             
%             case {'Fit spline over step'}
%                 coeffs = [X(find(profile == min(profile),1)) 10 profile(1) profile(end)];
%                
%                 fitted_vals = fminsearch(@fit_spline_onstep, coeffs, [], X,profile);
%                 
%                 position_min(n) = fitted_vals(1);
%                 box_range = coeffs(2);
%                 box_width_type = 2;
%                 
            
%             case {'Fit polynomial'}
%                 warning('This is not fitting as chosen -- it is fitting gaussian rather than polynomially')
%                 break
%                 %FIX ME this is a gaussian not a polynomial
%                 coeffs = [profile(1) (profile(end)-profile(1))/length(profile) max(profile)-min(profile) find(profile == min(profile),1) 4];
%                 coeffs = fminsearch(@gauss_and_bg,coeffs,[]);
% 
%                 position_min(n) = round(coeffs(4)) + top_of_area-1;
%                 box_range(n) = coeffs(5);
%                 if width_foil(n) > 50;
%                     width_foil(n) = 50;
%                 end
%                 box_width_type = 2;
               
            case {'None'}
                position_min(n) = (bot_of_area + top_of_area)/2;
                box_range(n) = HeightBox(n)/2;
                box_width_type = 1;

            otherwise
                error('No known type of box positioning')
        end
    end
end

%add the width of the box onto the position to define the boxes
if box_width_type == 1;
        boxY(:,1) = position_min - round(box_range);
        boxY(:,2) = position_min + round(box_range);
elseif box_width_type == 2;
        boxY(:,1) = position_min - (sigmas/2.*box_range);
        boxY(:,2) = position_min + (sigmas/2.*box_range);
else
    error('Unknown option')
end

NewBoxesX = boxX;
NewBoxesY = boxY;
end

%% gaussian functions
function [SSD] = fit_gaussian_and_bg(coefs,X,Y)
    
    func = gauss_plus_bg(coefs, X);
    
    SSD = sum((Y - func').^2);

    %plot(X,func,'r',X,coefs(1)+coefs(2)*X,'g',X,Y,'b'), %pause(.05)
end

function func = gauss_plus_bg(coefs, loc)

    a = coefs(1); 
    b = coefs(2);
    c = coefs(3);
    d = coefs(4);
    e = coefs(5);
    f = abs(coefs(6));
    
    func = (a +b*loc) -  c*ngaussian(loc,d,e,f);%(c*exp(-(((loc'-d)./e).^2)));

end

function [SSD] = fit_gaussian_and_step(coefs,X,Y)
    
    func = gauss_plus_step(coefs, X);    
    SSD = sum((Y - func').^2);

    %figure
    %plot(X,func,'r',[X(1), X(end)],[func(1), func(end)],'g',X,Y,'b'), % pause(.05)
end

function func = gauss_plus_step(coefs, loc)

    a1 = coefs(1); 
    a2 = coefs(2);
    height = coefs(3);
    pos = coefs(4);
    hw = coefs(5);
    mod = abs(coefs(6));
    
    lx=length(loc);
    hx=find(loc<pos,1,'last');
    bg(1:lx) = a1;
    bg(hx+1:length(loc)) = a2;
    
    low = height*ngaussian(pos,pos,hw,mod);
    scale = low - (a1 - a2);
    
    peak = ngaussian(loc,pos,hw,mod);
    peak(1:hx) = peak(1:hx)*height;
    peak(hx+1:end) = peak(hx+1:end)*scale;
    
    func = bg - peak;

end

function g = gaussian(x,pos,wid)
g = exp(-((x-pos)./(wid)) .^2);
end

function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = flattened Gaussian centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
% Shape is Gaussian when n=1. Becomes more rectangular as n increases.
%  T. C. O'Haver, 1988, revised 2014
% modified S. Hunt 2016.
g = exp(-((abs((x-pos))./(wid)).^(2*n)));
end

function g = bingaussian(x,pos,wid,n,m)
% BiGaussian (different width on leading edge and trailing edge).
% pos=position; wid=width 
% m determines shape; symmetrical if m=50.
%  T. C. O'Haver, 2012
% modified S. Hunt 2016.
lx=length(x);
hx=find(x<=pos,1,'last');%val2ind(x,pos);
g(1:hx)=ngaussian(x(1:hx),pos,wid*(m/100),n);
g(hx+1:lx)=ngaussian(x(hx+1:lx),pos,wid*(1-m/100),n);
end


%% sigmoid functions
function [SSD] = fit_sigmoids(coefs,X,Y)
    
    func = sigmoidpeak_onstep(X,coefs);
    
    SSD = sum((Y - func').^2);

%     plot(X,func,'r',X,coefs(1)+coefs(2)*X,'g',X,Y,'b'), pause(.05)
    plot(X,func,'r',[X(1),X(end)],[coefs(1) coefs(2)],'g',X,Y,'b'), pause(.05)
end

function [sigmoidpeak] = sigmoidpeak_onbg(X, coefs)

    a = coefs(1); 
    b = coefs(2);
    c = coefs(3);  
    pos = coefs(4);
    t1 = coefs(5);
    t2a = coefs(6);
    t2b = coefs(7);
    
    lx=length(X);
    hx=find(X<pos,1,'last');%val2ind(x,pos);
    g(1:hx)=1/2 + 1/2* erf(real((X(1:hx)+t1-pos)/sqrt(2*t2a)));
    g1=.5-.5*erf(real((X(hx:lx)-t1-pos)/sqrt(2*t2b)));
    g(hx:lx) = g1/g1(1)*g(hx);
    
    sigmoidpeak = (a +b*X) - c* g;

end

function [sigmoidpeak] = sigmoidpeak_onstep(X, coefs)

    a = coefs(1); 
    b = coefs(2);
    c = coefs(3);  
    pos = coefs(4);
    t1 = coefs(5);
    t2a = coefs(6);
    t2b = coefs(7);
    

    lx=length(X);
    hx=find(X<pos,1,'last');%val2ind(x,pos);
    g(1:hx)=(1/2 + 1/2* erf(real((X(1:hx)+t1-pos)/sqrt(2*t2a)))) * c ;
    low = 1/2 + 1/2* erf(real((pos+t1-pos)/sqrt(2*t2a)));
    bg(1:hx) = a;
    
    g1=.5-.5*erf(real((X(hx+1:lx)-t1-pos)/sqrt(2*t2b)));
    high=.5-.5*erf(real((pos-t1-pos)/sqrt(2*t2b)));
    bg(hx+1:length(X)) = b;
    
    scale = low*c - (a - b);
    g(hx+1:lx) = g1.*(scale);
    
    sigmoidpeak = bg - g;

end


%% spline functions 
% function [SSD] = fit_spline_onbg(coefs,X,Y)
%     
%     func = spline_onbg(X, Y, coefs);
%     
%     SSD = sum((Y - func').^2);
% 
% %     plot(X,func,'r',X,coefs(1)+coefs(2)*X,'g',X,Y,'b'), pause(.05)
%     plot(X,func,'r',[X(1),X(end)],[func(1) func(end)],'g',X,Y,'b'), pause(.05)
% end

% function [SSD] = fit_spline_onstep(coefs,X,Y)
%     
%     func = spline_onstep(X, Y, coefs);
%     
%     SSD = sum((Y - func').^2);
% 
% %     plot(X,func,'r',X,coefs(1)+coefs(2)*X,'g',X,Y,'b'), pause(.05)
%     plot(X,func,'r',[X(1),X(end)],[func(1) func(end)],'g',X,Y,'b'), pause(.05)
% end

% function vals = spline_onbg(x,y, coefs)
% 
% pos = coefs(1);
% width = coefs(2);
%         
% % Breaks
% breaks = [x(1), pos+[-1.5*width:width:1.5*width], x(end)];
% 
% % Clamped endpoints, y' = 0 at pos
% xc = [x(1), x(end), pos];
% yc = [y(1) ,y(end), 0];
% cc = [1 1 0; 
%       0 0 1; 
%       0 0 0];
% con = struct('xc',xc,'yc',yc, 'cc',cc);
% pp1 = splinefit(x,y,breaks,con);
% 
% %calculate the line
% vals = ppval(pp1,x);
% 
% end

% function vals = spline_onstep(x,y, coefs)
% 
% pos = coefs(1);
% width = coefs(2);
% bg1 = coefs(3);
% bg2 = coefs(4);
%  
% %background
% lx=length(x);
% hx=find(x<pos,1,'last');
% bg(1:lx) = bg1;
% bg(hx+1:length(x)) = bg2;
%     
% % Breaks
% breaks = [x(1), pos+[-2*width:width:2*width], x(end)];
% breaks(breaks>x(end)) = [];
% % bg_breaks = ones(size(breaks)*bg1;
% % bg_breaks(breaks>pos) = bg2;
% 
% % Clamped endpoints, y' = 0 at pos
% xc = [(round(breaks(2))), (round(breaks(end-1))), pos];
% yc = [bg1, bg2, 0];
% cc = [1 1 0; 
%       0 0 1; 
%       0 0 0];
% con = struct('xc',xc,'yc',yc, 'cc',cc);
% pp1 = splinefit(x,y,breaks,con);
% 
% %calculate the line
% vals = ppval(pp1,x);
% 
% peaks_edges = round([breaks(2) breaks(end-1)]);
% edges = [find(peaks_edges(1)==x) find(peaks_edges(end)==x)] ;
% 
% vals(1:edges(1)) = bg(1:edges(1));
% vals(edges(2):end) = bg(edges(2):end);
% 
% 
% end

function [sp] = interpolant_spline(X,profile,varargin)
% fit a spline to the profile.

%set default section length for the spline.
if length(varargin) >=1
    length_section = varargin{1};
else
    length_section = 6;
end

%spline constraints
% xc = [X(1) X(end)];
% yc = [ 0 0];
% cc = [0 0; 1 1; 0 0];
% con = struct('xc',xc,'yc',yc, 'cc',cc);
%sp = splinefit(X, profile, round(length(X)/length_section), con);

%removed spline constrains because the fit is crap with them. Why they were
%even here I do not know.
sp = splinefit(X, profile, round(length(X)/length_section));

if (0)
    figure(5),
    plot(X,profile, 'k.', X, ppval(sp,X),'r-'),
    hold on
    plot(sp.breaks, ppval(sp,sp.breaks),'or')
    xlabel('Position')
    ylabel('profile')
    legend('Profile', 'spline interpolant', 'breaks')
    keyboard
end

end


function pos = spline_values(splin,maxmin)
%find the minimum or maximum value of the presented spline.

%swtich between max or min function
if strcmpi(maxmin, 'min') == 1
    point = @min;
elseif strcmpi(maxmin, 'max') == 1
    point = @max;    
else
    error('spline_tp:unknown type', 'The lookup type is unknown. the options are ''max'' or ''min''.')
end

X = splin.breaks(1):splin.breaks(end);
spline_profile = ppval(splin, X);
[~, min_loc] = point(spline_profile);

%find the segment of the spline which contains the minimum.
segment =  find(X(min_loc)>=splin.breaks,1,'last');
if isempty(segment) %error catching
    segment = 1;
end


% figure(6),
%  plot(X, spline_profile,'r-')

%find the location of the spline minimum.
combined_mins = [];
for x = 1:3
    seg_eval = segment-2+x;
    if seg_eval>=1 && seg_eval<=length(splin.breaks)-1
        
        segment_coefs =  splin.coefs(seg_eval,:);
        r = roots(polyder(segment_coefs));
        
        range = splin.breaks(seg_eval+1)-splin.breaks(seg_eval);
        r_licit = find(r>=0 & r<range & isreal(r));
        
        minima =  polyval(segment_coefs, r(r_licit));
        
        combined_mins = [combined_mins; minima, r(r_licit)+splin.breaks(seg_eval)];
    end
end

if ~isempty(combined_mins)
    [~,a] = point(combined_mins(:,1));
    pos = combined_mins(a,2);
elseif numel(segment) == 1
    pos = splin.breaks(segment);
else
    error('I should not have got to here. Something is very wrong.')
end





end

function allCases = getCases
%copied from https://stackoverflow.com/questions/17325614/programmatically-get-valid-switch-case-values
% on 4 june 2018. then edited.
st = dbstack('-completenames');
%line = st(1).line;
fLines = importdata(st(1).file, sprintf('\n'));
switchLine = find(~cellfun(@isempty, ...
    regexp(fLines(1:end), '^\s*switch\s', 'once')), 1, 'last');
otherwLine = find(~cellfun(@isempty, ...
    regexp(fLines(1:end), '^\s*otherwise\s*$', 'once')), 1, 'last');
caseLines = fLines(switchLine:otherwLine);
casesStr = regexprep(caseLines(~cellfun(@isempty, ...
    regexp(caseLines, '^\s*case\s', 'once'))), '^\s*case\s*', '');
casesCells = cell(size(casesStr));
for iCases = 1:numel(casesCells);
    casesCells{iCases} = evalin('caller', casesStr{iCases});
end
allCases = [casesCells{:}]';
end

function [b, a] = swap(a, b)
    % This function has no body!
end

%Versions
% v 1.2 4 June 2018
%   - Automatic listing of the posibilities included in the file because
%   SelectBoxes_getCases does not work on new mac/ in matlab 2017b.
%	- other bug fixes for random errors brought on by version change in matlab.
% v 1.1
%	- Automatic reading of the posibilities from SelectBoxes_??.m
%	- Edited so that each box can have its own positioning type.
% v 1
%	- as extracted from Select Boxes

