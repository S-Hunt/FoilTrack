% SelectBoxes
%   Selects the boxes to be used in the position displcaement analysis by
%   ImageAnalysis. 
%   The script is called by ImageAnalysis.m It should be fed the image that the boxes 
%   are required for and a of file names and it will lead one throught the creation 
%   of the boxes for image sequence.
%
%   syntax:
%     SelectBoxes(Image, name_outfile, root_name, ref_id, num_images)
%       Image        - image for boxes to be selectd from
%       name_outfile - name of the file to save the boxes in
%       root_name    - root name of the experiment
%       ref_id       - image number in the sequence
%       num_images   - number of images in the sequence
%          
% see also: ImageAnalysis

% $ version: 2.2.3 $ 5th April 2019 $
% Simon Hunt 2008 - 2019
%       - for details of changes between versions see end of file. 

function SelectBoxes(I1, name_outfile, root_name, ref_id, num_images, fast_boxes)

% persistent Im

%% define defaults and read previous boxes
%Defines the height of the boxes which are selected.
search = 50;    %looks for the criteria +/-'search' pixels from the middle of the box.
HeightBox = 20; %-defines the height of the box which is being looked for.

%Figure annotations.
% how many numbers to skip when adding them to the images.
text_skip = 5;

%defines where the output files go
file_dir = fileparts(name_outfile);
[~,title_name] = fileparts(name_outfile);

box_file_name = strrep(title_name,'position_change','boxes');

title_name = strrep(title_name,'_position_change_',' ');
title_name = strrep(title_name,'_position_change',' ');

%define how to plot images 
%options and 'log', 'lin'
Im_plot_type = 'lin';
%define where to start in the processing loop.
proc = 0;

%read in selected areas if they exist.
prev_boxes = exist([file_dir,filesep,'data_prev.mat'],'file');
if prev_boxes == 2
    previous_boxes = load('data_prev');
    number_boxes = length(previous_boxes.boxX); 
    if isfield(previous_boxes, 'search')
        search = previous_boxes.search;
    end
    if isfield(previous_boxes, 'height')
        HeightBox = previous_boxes.HeightBox;
    end
    proc = 1;
else
    proc = 0;
end

%reads sizes of the image
[h1 w1] = size(I1);

%% Image Figure

%set colour scale limits in image.
if 0%exist('prctile') == 2 %=2two is function in matlab search directories
    
    lis  =  PlotI(I1, Im_plot_type, 'GetLims', 5);
    limitbot = lis(1);
    limittop = lis(2);
    
else
    limitbot = min(min(I1));
    limittop = max(max(I1));
end

% variables for loop
if exist('previous_boxes', 'var') == 1
    boxXold = previous_boxes.boxX;
    boxYold = previous_boxes.boxY;
    boxX = previous_boxes.boxX;
    boxY = previous_boxes.boxY;
    box_fit_type = previous_boxes.box_fit_type;
    Image_rotation = previous_boxes.Image_rotation;
    try 
        empty_boxes = previous_boxes.empty_boxes;
    catch
        empty_boxes = [];
    end
    number_boxes = size(boxX,1);
    box_height = boxY(:,2) - boxY(:,1);
    processing = 'justposition'; %just position the boxes on the new foils
else
    empty_boxes = [];
    Image_rotation = 0; 
    processing = 'select'; %need to select the boxes
end

%create the image to be passsed around.
if Image_rotation ~= 0
    Irot = imrotate(I1,Image_rotation,'bilinear');
    [h1 w1] = size(Irot);  
else
    Irot = I1;
end

%figure set up
Tag_str = 'Select';
Tag_str2 = 'Previous';

%see if a previous figure exists. 
fig_list = findall(0,'type','figure');
tags = get(fig_list, 'Tag');
%             tagged = ~cellfun(@isempty,strfind(tags,Tag_str));
if ~isempty(tags)
    if ~iscell(tags)
        tags = {tags};
    end
    for x = 1:length(tags)    
        if strcmpi(tags{x}, Tag_str) == 1
            tag_select = x;
        elseif strcmpi(tags{x}, Tag_str2) == 1
            tag_prev = x;
        end
    end
end

if exist('tag_select','var') == 1 %if the Select figure exists. 
% if sum(tagged) == 1 %reuse the old figure
    Im.Fig = fig_list(tag_select);
    Im.Plot = findall(Im.Fig,'type','axes');
    Im.Title =  findall(Im.Plot,'Tag','title');%findobj(Im.Plot,'Tag','title');
    Im.Image = findall(Im.Plot,'type','image');
    Im.rect = findall(Im.Plot,'type','rectangle');
    Im.text = findall(Im.Plot,'type','text');
    Im.text(Im.text == Im.Title) = []; %remove title from text list.
    
    if fast_boxes == 0
        figure(Im.Fig) % -- skip if fast boxes is on so does not spend time making figure current.
    end
    
%     Im.Plot = imagesc(Irot);
    %set new image
    set(Im.Image, 'CData', PlotI(Irot, Im_plot_type))
    set(Im.Plot, 'Xlim', [0 size(Irot,2)])
    set(Im.Plot, 'Ylim', [0 size(Irot,1)])
    
    %remove old rectangles and text.
    delete(Im.rect)
    delete(Im.text)
    
    %set new title
    NewTitle = [title_name,' - reference image #', num2str(ref_id), '/', num2str(num_images)];
    set(Im.Title, 'String', NewTitle);
    %Im.Title = title(NewTitle, 'FontWeight','Bold', 'Interpreter', 'none', 'Tag', 'title'); %makes a new title here because delete(Im.text) removes all text from the figure.
    
else %make a new figure 
    
    %create the figure to disply the image
    Im.Fig = figure('Tag', Tag_str);
    set(Im.Fig,'Name',['SelectBoxes v', FileVersion([mfilename('fullpath'),'.m'])], 'NumberTitle', 'off')
    set(Im.Fig,'Units','normalized','Position',[.01 0.0467 0.48 0.85]);

    %plot the image
    Im.Plot = axes;
    Im.Image = imagesc(PlotI(Irot, Im_plot_type));
    colormap(gray);
    axis image;
    hold on
    
    %set title for image
    Im.Title = title([title_name,' - reference image #', num2str(ref_id), '/', num2str(num_images)],...
                        'FontWeight','Bold', 'Interpreter', 'none', 'Tag', 'title');
end
    
caxis([limitbot limittop]);
      
    
%% processing loop
% processing = NaN;
re_pos = 1;

while strcmpi(processing,'done') == 0
    
    processing = lower(processing);
    switch processing
        case 'select' %select new boxes
            
            %if boxes exist but one wants to change options pass the old
            %boxes along.
            if exist('boxX', 'var')
                additional = {boxX, boxY};
            else
                additional = {};
            end
            [boxXold, boxYold, box_fit_type, empty_boxes, Image_rotation, Irot, search] = SelectNewBoxes(Im, I1, Image_rotation, HeightBox, search, Im_plot_type, additional{:});
            
            boxX = boxXold;
            boxY = boxYold;
            previous_boxes.box_fit_type = box_fit_type;
            [h1 w1] = size(Irot);
            number_boxes = size(boxX,1);
            
            re_pos = 1;
            
            box_height = boxY(:,2) - boxY(:,1);
            
            if number_boxes <= 10;
                text_skip = 1;
            end
            
        case 'movexy' %Move all the boxes in X and Y
            [boxX, boxY] = SelectBoxes_MoveXY(Irot, boxXold, boxYold);
            
            %boxXold = boxX;
            %boxYold = boxY;
            re_pos = 1;          
            delete(Im.boxnew)
            if isfield(Im, 'text2')
                delete(Im.text2)
            end
            
        case 'moveindividual' %Move individual boxes in Y
            
            [boxX, boxY] = SelectBoxes_MoveSingleBoxes(Irot, boxX, boxY);
            re_pos = 0;            
            delete(Im.boxnew)
            if isfield(Im, 'text2')
                delete(Im.text2)
            end

        case 'defineempty' %Define the empty boxes
            
            prompt = {'Which Boxes are empty?'};
            dlg_title = 'Empty Box list';
            num_lines = 1;
            def = {empty_box_2_string(empty_boxes)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            empty_boxes = empty_box_2_num(answer{:});
            
        case 'justposition'
            re_pos = 1;
            
        otherwise
            warning('Unknown processing option')
    end
    
    %Reposition the boxes if required.
    if re_pos == 1
        boxY = SelectBoxes_position(Irot, boxY, boxX, box_fit_type, search, empty_boxes, box_height);
%         boxX = boxXold;
    end
    
    %plot boxes on image
    for n = 1 : number_boxes
        %plot boxes from previous data set.
        Im.boxold(n) = rectangle('position',[boxXold(n,1)-0.5,...
            boxYold(n,1)-0.5,...
            boxXold(n,2)-boxXold(n,1)+1,...
            boxYold(n,2)-boxYold(n,1)+1],'EdgeColor','g','LineStyle',':',...
            'Parent', Im.Plot);
        %plot new boxes
        Im.boxnew(n) = rectangle('position',[boxX(n,1)-0.5,...
            boxY(n,1)-0.5,...
            boxX(n,2)-boxX(n,1)+1,...
            boxY(n,2)-boxY(n,1)+1],'EdgeColor','r',...
            'Parent', Im.Plot);
        
        %puts indentifying number below each box
        if round((n)/text_skip) == ((n)/text_skip) || n == 1 || n == number_boxes
            textx = (boxX(n,2)-boxX(n,1)) /2 + boxX(n,1);
            texty = boxY(n,2);
            name = num2str(n);
            %writes the number of the box on the image.
            if number_boxes <= 100 || text_skip > 20
                Im.text(n) = text(textx,texty,name,'color','r','HorizontalAlignment','center','VerticalAlignment','top',...
                    'Parent', Im.Plot);
            else
                Im.text(n) = text(textx,texty+10,name,'color','r','HorizontalAlignment','left','VerticalAlignment','middle',...
                    'Parent', Im.Plot);
            end
        end
        %writes the distance moved by the box above it
        if mean(boxYold(n,:)) - mean(boxY(n,:)) ~= 0
            textx = (boxX(n,2)-boxX(n,1)) /2 + boxX(n,1);
            texty = boxY(n,1);
            name = sprintf('%3.1f',mean(boxYold(n,:)) - mean(boxY(n,:)));
            Im.text2(n) = text(textx,texty,name,'color','r','HorizontalAlignment','center','VerticalAlignment','bottom',...
                'Parent', Im.Plot);
        end
    end

    % Question if the new boxes are in the right place.
    if fast_boxes == 0 %skip this if everything is being done without question.
               
        options = {'Yes';'Move all X,Y';'Move Single foils';'Define Empty';'Select New Boxes';'Edit';'Quit'};
        choice = SelectBoxes_MultiButton('Correct Boxes',...
            'Are the boxes in the right place? If not how do you wish to adjust them?',...
            options);
        
    else
        choice = 'Yes';
    end
        
        
    % Handle response
    switch choice
        case 'Yes'
            processing = 'done';
            
        case 'Move all X,Y'
            processing = 'movexy';
            
        case 'Move Single foils'
            processing = 'moveindividual';
            
        case 'Define Empty'
            processing = 'defineempty';
            
        case 'Select New Boxes'
            processing = 'select';
            
        case 'Edit'
            keyboard;
            
        case {'No', 'Quit'}
            close all
            error('I''d rather be doing anything but this so ending program.')
                        
        otherwise
            error('The script doesn''t like where this is going. It is time for it to give up and have a nap.')
    end
end

%% saves the boxes and image of box position
if exist('Xfull') ~= 1
    Xfull = 0;
    Yfull = 0;
end

% box_file_name = strcat(root_name,'_boxes.mat');
var_Array = {'boxX', 'boxY', 'Xfull', 'Yfull', 'Image_rotation', 'box_fit_type', 'empty_boxes', 'HeightBox', 'search'};

if exist([box_file_name,'.mat'],'file') == 0 | fast_boxes == 1
    save(box_file_name, var_Array{:});
else
    % Gives dialogue box asking to overwrite the file...
    choice = questdlg('Box file already exists. Overwrite?', ...
        'Overwrite', ...
        'Yes','No','Yes');
    
    % Handle response
    switch choice
        case 'Yes'
            save(box_file_name, var_Array{:});
    end
end

%rewrite data_prev.mat
save('data_prev', var_Array{:});

%checks if an image of the boxes starting positions exists and if not generates one.
% image_boxes = exist(Im_starting_boxes,'file');
image_boxes = 0; % FIXME - we assume we won't clobber a file here!
if image_boxes == 0

    %if the old figure exists then use it
    if exist('tag_prev','var')==1 %if the Select figure exists. 
        Im2.Fig = fig_list(tag_prev);
        Im2.Plot = findall(Im2.Fig,'type','axes');
        Im2.Title = findobj(Im2.Plot,'Tag','title');
        Im2.Image = findall(Im2.Plot,'type','image');
        Im2.rect = findall(Im2.Plot,'type','rectangle');
        Im2.text = findall(Im2.Plot,'type','text');
        
        if fast_boxes ~= 1
            figure(Im2.Fig);
        end
        
        %set new image
        set(Im2.Image, 'CData', PlotI(Irot, Im_plot_type))
        
        %set new title
        ttl_str = [title_name,' - reference image #', num2str(ref_id), '/', num2str(num_images)];
        if strcmpi(Im_plot_type, 'log') == 1
            ttl_str = [ttl_str, ' -- colour scale = log(Intensity)']; 
        end
%         set(Im2.Title, 'String', ttl_str);
        
        %remove old rectangles and text.
        delete(Im2.rect)
        delete(Im2.text)
    else
        %create the figure to disply the image
        Im2.Fig = figure('Tag', Tag_str2);
        set(Im2.Fig,'Name',['Previous Boxes. SelectBoxes v', FileVersion([mfilename('fullpath'),'.m'])], 'NumberTitle', 'off')
        set(Im2.Fig,'Units','normalized','Position',[.51 0.0467 0.48 0.85]);
        
        %plot the image
        Im2.Image = imagesc(PlotI(Irot, Im_plot_type));
        Im2.Plot = findall(Im2.Fig,'type','axes');
        hold on;    
        colormap(gray);
        axis image;
        
        %set title for image 
        ttl_str = [title_name,' - reference image #', num2str(ref_id), '/', num2str(num_images)];
        if strcmpi(Im_plot_type, 'log') == 1
            ttl_str = [ttl_str, ' -- colour scale = log(Intensity)']; 
        end
%         Im2.Title = title(ttl_str,'FontWeight','Bold', 'Interpreter', 'none', 'Tag', 'title');
    end
    
    %plot the image
%     Im2.Plot = imagesc(Irot);
    %allows over writing of the image
    caxis([limitbot limittop]);
    title({ttl_str; 'Box Positions'},'FontWeight','Bold', 'Interpreter', 'none');

    for n = 1:number_boxes
        Im2.rect(n) = rectangle('position',[boxX(n,1)-0.5,boxY(n,1)-0.5,boxX(n,2)-boxX(n,1)+1,boxY(n,2)-boxY(n,1)+1],'EdgeColor','r',...
            'Parent', Im2.Plot);
        %puts indentifying number in the box
        if round((n-1)/text_skip) == ((n-1)/text_skip)
            textx = (boxX(n,2)-boxX(n,1)) /2 + boxX(n,1);
            texty = boxY(n,2);
            name = num2str(n);
            %writes the number of the box on the image.
            if number_boxes <= 100 || text_skip > 20
                Im2.text(n) = text(textx,texty,name,'color','r','HorizontalAlignment','center','VerticalAlignment','top',...
                    'Parent', Im2.Plot);
            else
                Im2.text(n) = text(textx,texty+10,name,'color','r','HorizontalAlignment','left','VerticalAlignment','middle',...
                    'Parent', Im2.Plot);
            end
        end
    end

    pause(0.01)
    
    if fast_boxes ~= 1 %
%     print('-dtiff', [box_file_name,'.tif']);
    saveas(Im2.Fig, [box_file_name,'.tif']);
    end
end


end %SelectBoxes




%% subfunctions
function [boxX, boxY, box_fit_type, empty_box, Image_Image_rotation, I1, search_new] = SelectNewBoxes(Im, I1, Image_Image_rotation, default_height, default_search, Im_plot_type, varargin)

    if ~isempty(varargin) 
       Xbox = varargin{1};
       Ybox = varargin{2};
    end

    % Check the Image_rotation of the image
    rot_final = -1;
    set(Im.Fig, 'Pointer', 'fullcrosshair');
    while rot_final ~= Image_Image_rotation;
        rot_final = Image_Image_rotation;
        
        %Image_rotation dialogue box
        prompt = {'Image Image_rotation'};
        dlg_title = 'Box setup';
        num_lines = 1;
        def = {num2str(Image_Image_rotation)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if  isempty(answer)==1
            answer = {0};
        end
        Image_Image_rotation = str2double(answer{1});
        
        %rotate the image
        Irot = imrotate(I1,Image_Image_rotation,'bilinear');
        
        %replot the image
%         Im.Plot = imagesc(Irot);
        set(Im.Image, 'CData', PlotI(Irot, Im_plot_type))
        set(Im.Plot, 'Xlim', [0 size(Irot,2)])
        set(Im.Plot, 'Ylim', [0 size(Irot,1)])
        
        
    end
    
    %if the image is rotated then redefine the image to work on and resize the limits.
    if Image_Image_rotation ~= 0
        I1 = Irot;
        clear Irot;
    end
    set(Im.Fig, 'Pointer', 'arrow');
    
    %get limits of I1
    [h1 w1] = size(I1);
    
    % select the boxes
    [Xpick, Ypick] = ginput;
    if ~isempty(Xpick) %if new boxes are picked with the pointer
        
        Xpick = round(Xpick);
        Ypick = round(Ypick);
        
        %FIX ME. Need to check the order of the inputs in case the data is
        %selected not top left and bottom right.
        
        %define the X and Y values of each box in terms of the edges instead of the corners
        for num = 1 : size(Xpick)
            if num/2 == round(num/2)
                col = 2;
                row = num/2;
            else
                col = 1;
                row = num/2 + 0.5;
            end
            boxXraw(row,col) = Xpick(num);
            boxY(row,col) = Ypick(num);
        end
        
        %     boxXraw = boxX;
        boxX = round(boxXraw);
        boxX(boxX < 1) = 1;
        boxX(boxX > w1) = w1;
        boxY(boxY < 1) = 1;
        boxY(boxY > h1) = h1;  
        
    else
        
        %copy old boxes into new boxes
        boxX = Xbox;
        boxY = Ybox;
        clear Xbox Ybox
    end
    
    %FOR THERMAL DIFFUSION
    %Split the selected boxes into the requested number of boxes.
    prompt = {'Align boxes vertically? (yes/no)'...
              'Height of the boxes (in pixels, either single number or one number per drawn box)'...
              'Search distance to find feature'...
              'Is this a thermal diffusion experiment? If so decide on the number of boxes per selected row'...
              'Or decide the width of each box (in pixels)'};
    dlg_title = 'Box setup';
    num_lines = 1;
    def = {'Yes', num2str(default_height), num2str(default_search), num2str(1), ''};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    %align the boxes vertically?
    if strcmpi(answer{1,1},'Yes') == 1 %align the boxes vertically
        boxX(:,1) = max(boxX(:,1));
        boxX(:,2) = min(boxX(:,2))-1;
    elseif strcmpi(answer{1,1},'No') == 1
        %there is nothing here deliberately 
        %i.e. do not align the boxes.
    else
        %if the answer is not yes or no!
        warning('SB:noalign','Answer is not recognised. The boxes have not been aligned')
    end
        
    %FIX ME. The cutting of the boxes by thermal_diff_boxes assumes that they are aligned. This
    %is not necessarily the case. 
    
    % set box heights. 
    heights = commasep_2_num(answer{2,1});
    if numel(heights) == 1
        heights = ones(size(boxY,1),1) * heights;
    end
    boxY = round([mean(boxY,2)-heights(:)/2, mean(boxY,2)+heights(:)/2 ])
    
    %set the search distance
    search_new = str2num(answer{3});
    
    %cut the boxes?
    if isempty(answer{5,1}) == 0 %the box width is defined.
        box_width = str2double(answer{5,1}); %box width
        total_width = boxX(1,2) -  boxX(1,1);
        cut_into = round(total_width/box_width); %number of boxes to cut into

        new_total_width = cut_into * box_width;
        adjust = (new_total_width - total_width) /2; %how far to move the array to align the boxes up
        boxX(:,1) = round(boxX(:,1) - adjust);
        boxX(:,2) = round(boxX(:,2) + adjust);
    else %number of boxes is defined.
        cut_into = str2double(answer{4,1});
    end

    empty_box = [];
    if cut_into ~= 1
        [boxXcut, boxYcut] = thermal_diff_boxes(boxX,boxY,cut_into);
        
        % defines which boxes do not have foil in them.
        % The boxes are labelled as empty if they are outside the selected areas.
        for x = 1 : size(boxX,1) %loop by row.
            row_empties = find(boxYcut(:,1) == boxY(x,1) & (boxXcut(:,2) <= boxXraw(x,1) | boxXcut(:,1) >= boxXraw(x,2) ));
            empty_box = [empty_box row_empties'];
        end
        if isempty(empty_box)
            clear empty_boxes;
        end
        boxX = boxXcut;
        boxY = boxYcut;
    end    
    
    % determine how to position the boxes in the images.
    fit_opts = SelectBoxes_position('posibilities'); %call possibilities from SelectBoxes_position
    
    selected = 0;
    while selected == 0
        [s,v] = listdlg('PromptString',{'Select method for positioning boxes:'},'SelectionMode','single','ListString',['INDIVIDUAL';fit_opts],'ListSize',[175 175],'Name','Box positioning');
        if v == 1
            if s == 1
                box_fit_type = SelectBoxes_tableUI(fit_opts, boxX);
                selected = 1;
            else
                box_fit_type = char(fit_opts(s-1));
                selected = 1;
            end
        else
            % Construct a questdlg
            choice = questdlg('You have not selected a Positioning type. Do you want to continue?', ...
                'No type selected', ...
                'Yes','No','No');
            % Handle response
            switch choice
                case 'Yes'
                    selected = 0;
                case 'No'
                    return %exit the sctipt.
            end
        end
    end %while selected == 0
    
end
    
%% box processing functions.

function [newboxesX, newboxesY] = thermal_diff_boxes(old_boxX,old_boxY,number)
% Cut the selected areas into the number of boxes desired 

%width and number of original selection
number_boxes = size(old_boxX);
total_width = old_boxX(1,2) - old_boxX(1,1);

%sets the size for each of the sub-boxes. The total box size might
%change slightly from original because of rounding.
each_width = round(total_width/number);

for x = 1 : number_boxes
    for num = 1 : number
        row = num+(x-1)*number;
        boxXnew(row,1) = old_boxX(1,1) + (num-1) * each_width;
        boxXnew(row,2) = old_boxX(1,1) + num * each_width - 1;
        boxYnew(row,1) = old_boxY(x,1);
        boxYnew(row,2) = old_boxY(x,2);
    end
end

newboxesX = boxXnew;
newboxesY = boxYnew;

end

function out = empty_box_2_string(box_list)
%turns list of boxes from numbers into string with commas and hyphens
if isempty(box_list) == 1
    empty = '';
elseif length(box_list) == 1
    empty = num2str(box_list);
else
    difs = [1 box_list(2:end)] - [0 box_list(1:end-1)];
    missing = find(difs~=1)-1;
    miss = [0 missing length(box_list)];
    empty = []; 
    for i = 1:length(miss)-1,
        if miss(i) ~= miss(i+1)-1
            empty = [empty, num2str(box_list(miss(i)+1)), '-', num2str(box_list(miss(i+1)))];
        else
            empty = [empty, num2str(box_list(miss(i)+1))];
        end
        if i ~= length(miss)-1
            empty = [empty, ', '];
        end
    end
end

out = empty;
end

function out = empty_box_2_num(box_string)

%turns list of comma separated numbers into sorted, unique list
out = commasep_2_num(box_string);
out = sort(unique(out));

end

function vals = commasep_2_num(in)
%turns string of input numbers into list of numbers
dividers = find(in == ',');

dividers = [1 dividers length(in)+1]; %adds leading zero and length to the variable so that the following for loop works.
vals = [];
for i = 1 : length(dividers)-1
    num = in(dividers(i):dividers(i+1)-1);
    hyphen_location = find(num == '-');
    if isempty(hyphen_location)
        vals = [vals str2num(num)];
    else
        start = str2num(num(1:hyphen_location-1));
        finish = str2num(num(hyphen_location+1:end));
        vals = [vals start:finish];
    end
end
end %commasep_2_num


function [version, content] = FileVersion(file)

%reads version of this file.
fid = fopen(file);
content = fscanf(fid,'%c');
dollars = find(content == '$');
version =  strtrim(content(dollars(1)+1:dollars(2)-1));
colon   = find(version == ':');
version =  strtrim(version(colon+1:end));

fclose(fid);

end %FileVersion


function choice = choosedialog(title, message, options)

%     options = {'Red';'Green';'Blue';'yellow'}

    box_position = [300 300 250 150];
    text_position = [20 80 210 40];
    button_positions = [89 20 70 25;
                        89 50 70 25;
                        89 80 70 25;
                        89 120 70 25];
                    
    
    d = dialog('Position',box_position,'Name',title);
    txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',text_position,...
           'String',message);
       
%     popup = uicontrol('Parent',d,...
%            'Style','popup',...
%            'Position',[75 70 100 25],...
%            'String',options,...  'String',{'Red';'Green';'Blue'},...
%            'Callback',@popup_callback);
    for x = 1:length(options)
        
        button(x) = uicontrol('Parent',d, 'Style', 'pushbutton',...
            'String', options(x), ...
            'Position', button_positions(x,:),...
            'Tag', options(x),...
            'Callback',button_callback);
        %[0.4 0.5 0.25 0.15],'Units', 'normalized') ;
%     btn.x = uicontrol('Parent',d,...
%            'Position',button_positions(x,:),...
%            'String',options(x),...
%            
    end
    
       
    % Wait for d to close before running to completion
    uiwait;
   
       function choice = button_callback(source,eventdata)
%           idx = popup.Value;
%           popup_items = popup.String;
          % This code uses dot notation to get properties.
          % Dot notation runs in R2014b and later.
          % For R2014a and earlier:
          get(source,'Value')
%           idx = get(popup,'Value');
          popup_items = get(source,'String');
          choice = char(popup_items(idx,:));
          delete(gcf);
       end
end

function Iout = PlotI(Iin, plot_type, varargin)

    if strcmpi(plot_type, 'Log') == 1
        Iout = log(Iin);
        
    elseif strcmpi(plot_type, 'Lin') == 1
        Iout = Iin;
        
    else
        error(' Unknown image ploting type')
    end

    if nargin > 2 && strcmpi(varargin{1}, 'GetLims') == 1
        
        if nargin > 3
            cut = varargin{2};
        else
            cut = 5;
        end
        
        limitbot = prctile(Iout(:), cut);
        limittop = prctile(Iout(:), 100-cut);
        
        Iout = [limitbot, limittop];
    end
end

%% Versions
% v2.2.3 - 5th April 2019
%   - added function to plot off figures as log or linear and with
%   percetile limits.
%   - Fixed move all x,y so that uses old boxes when doing the moving,
%   rather than the new (and misplaced) ones.
% v2.2.2 - 5th June 2018
%   - Fixed bug that meant the wrong type of fit was being selected after
%   the list dialoge.
% v2.2.1 - 4th June 2018
%   - Fixed bug that removed title from reused figure. The title is not
%   present.
% v2.2 - 30th June 2017
%   - Added search and BoxHeight to dataprev.mat so that these values are
%   carried through the box picking.
% v2.1 - 17th January 2017
%   - Bug fix. The script passed I1, rather than Irot to the sub functions
%   in the box positioning while loop. Now fixed and the loop passes Irot
%   around.
%   - Bug Fix. The X-Y repositing was not propagated properly through the
%   loop and got lost. Now fixed.
% v2.0 - June 2016
%   - rewritten so that it is possible to propagate the boxes for each data set
%   with only one click of the mouse. 
% v 1. The original version