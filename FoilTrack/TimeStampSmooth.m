% TIME STAMP SMOOTH 
%
%   For smoothing the time stamps of tiff, bmp, etc images in which the 
%   time stamp precission is 1 second.
%
%   syntax:
%       new_time_stamps = TimeStampSmooth(time_stamps, options...)
%       new_time_stamps = TimeStampSmooth('*.lst file')
%       new_time_stamps = TimeStampSmooth('*.mat file')
%
%   options:
%       'check':
%       'cut off':
%       'discard':
%       'differences':
%       'plot':
%       'polynomial':
%       'replace':
%
%   $ version: 0.3.2 $ April 2019 $
%       - for details of changes between versions see end of file.
%
%   Simon Hunt 2015-2019
%
% FIX ME: the smoothing is set to work for time stamp differences <1. It
% will not smooth data with longer differences. This needs to be fixed so
% that interpolation can be done for the d-DAC experiments.

function [new_time_stamps, not_use] = TimeStampSmooth(in, varargin)

rotate = 0;

if ischar(in)
    [~, file_name, ext] = fileparts(in);
    if strcmpi(ext, '.txt') || strcmpi(ext, '.lst')
        file_id = in;
        data = importdata(file_id);
        %     keyboard
        time_stamps = data.data;
        time_stamps = time_stamps-time_stamps(1);
    elseif strcmpi(ext, '.mat')
        % Import the file
        newData1 = load('-mat', in);
        time_stamps = newData1.timestamps - newData1.timestamps(1);
%         keyboard
    else
        error('Unkown file type')
    end
    fig_title = file_name;
else
    time_stamps = in;
    fig_title = 'Time Stamps';
end
if size(time_stamps,1) ~= size(time_stamps(:),1)
    time_stamps = time_stamps(:);
    rotate = 1;
end

%set defaults 
diff_cut = 2; %difference in time steps to cut the fitting.
poly_order = 2; %order of the polynomial to fit the time stamps with.
cut_off = 0.8;  % maximum difference between input and output time stamps.
plot_on = 0;
check = 0;
discard = 0;
substitute = [];

%force to be ready for output
not_use = 0;

%parse varargin.
iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case 'plot'
            plot_on = 1;
            iarg = iarg + 1;
        case 'check'
            check = 1;
            iarg = iarg + 1;
        case 'discard'
            discard = 1;
            iarg = iarg + 1;
        case 'differences'
            diff_cut = varargin{iarg+1};
            iarg = iarg + 2;
        case 'polynomial'
            poly_order = varargin{iarg+1};
            iarg = iarg + 2;
        case 'replace' %this was inserted to attemp to fix the D-DAC expeirments. 
            % The option calculates the mean time step for the good data
            % and then uses this where the time stamps are really bad.
            substitute = varargin{iarg+1};
            iarg = iarg + 2;
        case 'cut off'
            cut_off = varargin{iarg+1};
            iarg = iarg + 2;
    end
end

%exit script if time stamps are already unique
if numel(unique(time_stamps)) == numel(time_stamps)
    new_time_stamps = time_stamps;
    return
end

%Error check to see if incoming data is processable.
% The script assumes that the time stamps always increase 
%script will error out if there is a negative time difference.
diff = time_stamps(2:end) - time_stamps(1:end-1);
if min(diff) < 0
    error('TimeStampSmooth:NotIncreasing','The time stamps do not always increase \n Check the input data and try again')
end



% keep the originals and cut the data we are fitting.
time_stamps_all = time_stamps;

% discard data with all the same time stamps at the end of the sequence.
% This is for the dDAC experiments at 16IDB that saved tiff radiographs.
if discard == 1
    allowed_diff = 0.01;
    sames = find(time_stamps(end)-time_stamps < allowed_diff);
    not_use = length(sames)-1;
else 
    not_use = 0;
end
time_stamps = time_stamps(1:end-not_use);


% if we are replaceing some of the time stamps 
ts_pos_all = 1:length(time_stamps-not_use);
ts_pos = ts_pos_all;
if ~isempty(substitute)
    
    time_stamp_key = logical(ones(size(time_stamps)));
    time_stamp_key(substitute) = 0;
    
    time_stamps = time_stamps(time_stamp_key);
    ts_pos = ts_pos_all(time_stamp_key);
    
    diff = time_stamps(2:end) - time_stamps(1:end-1);
end


%find any obvious gaps in the data (i.e. differences > 1)
loc = find(diff>diff_cut);

%fit polynomial to the time stamps and use to calcualte interpolated time
%stamps
% poly_order = polynomial_order;
% cut_off = 0.8;
loc = [0; loc(:); length(time_stamps)];

done = 0;
loc_old = [];
while done == 0
    
    for x = 1 : length(loc) - 1
        begin = loc(x)+1;
        stop = loc(x+1);
        X = begin:stop;
        Y = time_stamps(X);
        
        steps = round(length(X(:))/100*3);
        if steps == 0, steps = 1; end
        if length(X) < 10, order = 2; else order = poly_order; end
        pp1 = splinefit(X(:),Y(:),steps,order);
        new_time_stamps(X) = ppval(pp1,X);
        
    end
    
    loop = 0;
    while min(new_time_stamps(2:end) - new_time_stamps(1:end-1))<0 && loop < 5
        loop = loop+1
        
        grad = new_time_stamps(2:end) - new_time_stamps(1:end-1);
        locate = find(grad<0);
        
        for y = 1:length(locate)
            
            smooth_over = 3;
            J = locate(y) - smooth_over:locate(y) + smooth_over;
            J(J<1) = [];
            K = time_stamps(J);
        
            pp1 = splinefit(J(:),K(:),2,2);
            new_time_stamps(J) = ppval(pp1,J);
        end
    end
    
    residuals = time_stamps - new_time_stamps(:);
    
    %display time stamps and fitted time stamps. For debugging
    if plot_on
        figure(1)
        subplot(2,1,1)
        plot(1:length(time_stamps_all), time_stamps_all,'.-')
        hold on
        plot(ts_pos, new_time_stamps, 'r')
        hold off
        xlabel('Time stamp number in series')
        ylabel('Time (s)')
        legend({'Original time stamps','Smoothed time stamps'}, 'Location', 'NorthWest')
        title(fig_title, 'Interpreter', 'none')
        xlims = get(gca, 'Xlim');
        
        subplot(2,1,2)
        plot(ts_pos, time_stamps-new_time_stamps(:),'.-')
        xlabel('Time stamp number in series')
        ylabel('Residuals (s)')
        set(gca, 'Xlim', xlims)
        
%         pause(1)
%         return
    end
    
    if max(abs(residuals)) > cut_off %if there are large residuals
        
        big_resid = find(abs(residuals) > cut_off);
        unaccouted = big_resid(~ismember(big_resid,loc));
        vals = residuals(unaccouted);
        [~,biggest] = max(abs(vals));
        
        
        if ~isempty(biggest) && min(abs(loc - unaccouted(biggest))) > 5
            %is there a big residual which is not accounted for? and which is
            %not close to a previous break?
            
            % add residual location to list of breaks
            loc = [loc; unaccouted(biggest)];
            
        elseif isempty(unaccouted) && ~isempty(big_resid) %is there a big residual that is accounted for but still large
            
            if length(big_resid) ~= 1
                big_resid = big_resid(1);
            end
            
            loc(loc == big_resid) = [];
            add1 = min(big_resid+6, length(time_stamps));
            add2 = max(big_resid-6, 1);
            loc = [loc; add1; add2]

        elseif min(abs(loc - unaccouted(biggest))) < 5
            % is the residual close to a previous break
            
            if ismember(min(abs(loc - unaccouted(biggest))), (loc - unaccouted(biggest)))
            % is the residual below the break?
                loc = [loc; unaccouted(biggest)-11];
            else
                loc = [loc; unaccouted(biggest)+11];
            end
     
        end
                
    else
        done = 1;
    end
    
    if isequal(sort(unique(loc)),loc_old) && done ~= 1
        fprintf('TimeStampSmoothing: %s, The biggest residual is greater than %3.1f seconds.\n', fig_title, cut_off)
        disp('Release the process to keyboard control.')
        keyboard
    end
    
    
    loc = sort(unique(loc));
    loc_old = loc;
     
%     keyboard
end

% if we are replaceing some of the time stamps 
% keep the originals and cut the data we are fitting.
if ~isempty(substitute)
    
    step = (new_time_stamps(2:end) - new_time_stamps(1:end-1));
    
    time_stamps_step = ones([numel(time_stamps_all)-1,1]) * mean(step);
    
    time_stamp_key(end)=[];
    time_stamps_step(time_stamp_key==1) = step;
    time_stamps_step(~time_stamp_key) = mean(step);
    
    new_time_stamps = [0; cumsum(time_stamps_step)] + time_stamps_all(1);
end


%check to see if the output times are always increasing.
if min(new_time_stamps(2:end) - new_time_stamps(1:end-1)) < 0
    
    %the possible ways to deal with this are:
    % 1. Find the median time step between images and replace the non-unique time stamps with these gaps.
    % 2. REPLACE THE TIME STAMPS with those from XRD (for Fe anelasticity experiments)
    %       a. start at the begining and match the times (assume any missing images are at the end)
    %       b. start at the end and match the times (assume any missing images are at the start)
    %       c. stretch the times to match the range of the xrd data (assume the missing images are in the middle somewhere).
    %       
    warning on
    warning('TSS:OhNo', 'The smoothed time stamps are not continuously increasing \nRelease to keyboard')
    warning off
    keyboard
    
     
     
%     diffs = (new_time_stamps(2:end) - new_time_stamps(1:end-1));
%     diffs_non = diffs(diffs~=0);
%     
%     scal=1000;
%     step = unique(round(diffs_non*scal)/scal);
%     
%     if numel(step) == 1
%        new_time_stamps = (0:length(new_time_stamps)-1) * step;
%     else
%         step = mean(diffs_non);
%         diffs(diffs==0) = step;
%         new_time_stamps = [0; cumsum(diffs(:))];
% %         keyboard
%     end
  
end

new_time_stamps = new_time_stamps(:);

%set first time stamp to 0 if the input time stamps started at 0.
if time_stamps(1) == 0
    new_time_stamps = new_time_stamps - new_time_stamps(1);
end


%display time stamps and final fitted time stamps.
if plot_on || check
    %figure(2)
    subplot(2,1,1)
        plot(1:length(time_stamps_all), time_stamps_all,'.-')
        hold on
        plot(ts_pos_all, new_time_stamps, 'r')
        
%     plot(time_stamps,'.-')
%     hold on
%     plot(new_time_stamps, 'r')
    hold off
    xlabel('Time stamp number in series')
    ylabel('Time (s)')
    legend({'Original time stamps','Smoothed time stamps'}, 'Location', 'NorthWest')
    title(fig_title, 'Interpreter', 'none')
    
    subplot(2,1,2)
    plot(time_stamps_all(1:length(new_time_stamps))-new_time_stamps(:),'.-')
    xlabel('Time stamp number in series')
    ylabel('Residuals (s)')
    
    if ~check
        pause
    end
end


%% error checking of the output time stamps.
if numel(new_time_stamps) ~= numel(time_stamps_all) && not_use ~= 0
   disp('Some time stamps have been discarded because the time stamps at the end are identical')
end

if numel(new_time_stamps) ~= numel(time_stamps_all)-not_use
    error('The number of elements in the modified time stamps is not the same as in the original')
end

if min(new_time_stamps(2:end) - new_time_stamps(1:end-1))<0
%     new_time_stamps(2:end) - new_time_stamps(1:end-1)
   error('The new time stamps change the order of the images.') 
end

%if size(new_time_stamps) ~= size(in) & ~ischar(in)
if rotate == 1
    new_time_stamps = new_time_stamps';
end


%% VERSION info:
%  - v 0.4 -- April 2019
%       - added option to override discarding of frames
%  - v 0.3 -- Dec 2018
%       - added optin to discard some of the frames 
%		- added check options.
%  - v 0.2.1 -- May 2017
%       - Added option so the script can read *.lst files directly.
%  - v 0.2 -- Feb 2017
%       - Input changed so can feed varibles to the process.
%       - Added option to read time stamps from .mat files. 
%  - v 0.1 -- November 2015
%     - Tme stamp processing moved into a separate file