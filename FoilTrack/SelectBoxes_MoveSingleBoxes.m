% SelectBoxes_MoveBoxes
%   Provides validation for the box position and allows individual boxes to be moved
%   around. The script is called by SelectBoxes.m. It should be fed the image that the boxes 
%   are required for and a of file names and it will lead one throught the creation 
%   of the boxes for image sequence.
%
%   syntax:
%     SelectBoxes_MoveBoxes(I1, boxX, boxY)
%          
% see also: SelectBoxes, ImageAnalysis

% Simon Hunt June 2016
% Used to be part of Select Boxes but was removed when SelectBoxes was
% rewritten (version 2).

function [boxX, boxY] = SelectBoxes_MoveBoxes(Im1, boxX, boxY)

%size of image.
[h1 ~] = size(Im1);

%presents options to choose wrong boxes and then correct them
prompt = {'Which boxes are in the wrong place?', 'Search for foil +/- ?? pixels'};
dlg_title = 'Right place?';
num_lines = 1;
def = {'box numbers separated by commas', num2str(20)};
answer = inputdlg(prompt,dlg_title,num_lines,def);

wrong = answer{1,1};
search = str2num(answer{2,1});  %#ok<ST2NM> %used to define search area later

%sorts the variable 'wrong' into a list of numbers.
wrong_boxes = empty_box_2_num(wrong);

%asks for input for correct places of boxes
extra_space = search;% - round(box_height/2);
if extra_space < 0, extra_space = 0; end

%plots the incorrect box and allows it to be fixed.
for m = 1 : length(wrong_boxes)
    done = 0;
    
    n = wrong_boxes(m);
    Y_height = boxY(n,2) - boxY(n,1);
    Y_pos = boxY(n,1) + Y_height/2;
    
    while done ~= 1
        
        top_of_area = round(Y_pos - Y_height/2 - extra_space);
        bot_of_area = round(Y_pos + Y_height/2 + extra_space);
        if top_of_area <= 0
            top_of_area = 1;
        end
        if bot_of_area > h1
            bot_of_area = h1;
        end
        
        area = Im1(top_of_area : bot_of_area, boxX(n,1):boxX(n,2));
        
        profile = mean(area,2);
%         profile_range  = max(profile) - min(profile);
        
        area_x = [max(profile) min(profile) min(profile) max(profile)];
        
        if Y_pos - Y_height/2 < top_of_area
            area_y = [1 1 Y_pos+Y_height/2 Y_pos+Y_height/2];
        elseif Y_pos + Y_height/2 > bot_of_area
            area_y = [Y_pos-Y_height/2 Y_pos-Y_height/2 bot_of_area bot_of_area];
        else
            area_y = [Y_pos-Y_height/2 Y_pos-Y_height/2 Y_pos+Y_height/2 Y_pos+Y_height/2];
        end
        
        Iplot = Im1(top_of_area : bot_of_area, boxX(n,1):boxX(n,2));
        
        %plot box and details.
        repos_figure_handle = findobj('type','figure','name','Reposition');
        if isempty(repos_figure_handle)
            repos_figure_handle = figure('Units','normalized','Position',[.01 0.0467 0.98 0.85], 'name','Reposition');
        end
        
        %plot entire image
        subplot(1,3,1);
        imagesc(Im1);
        if exist('prctile') == 2 %=2two is function in matlab search directories
            limitbot = prctile(Im1(:), 5);
            limittop = prctile(Im1(:), 95);
        else
            limitbot = min(min(Im1));
            limittop = max(max(Im1));
        end
        colormap('jet');
        rectangle('Position', [boxX(n,1) top_of_area boxX(n,2)-boxX(n,1) bot_of_area-top_of_area]);
        caxis([limitbot limittop]);
        title('Radiograph + position of box')
        
        %plot detail of box being investigated
        subplot(1,3,2);
        imagesc(Iplot);
%         colormap('default'); %hold on; %plot([min(profile),max(profile)],[position_min(n)+box_height/2,position_min(n)+box_height/2], 'k:', [min(profile),max(profile)],[position_min(n)-box_height/2,position_min(n)-box_height/2], 'k:');
        title('Detail of box and above/below it')
        
        %plot profile across box being investigated
        subplot(1,3,3); 
        fill(area_x, area_y, 'c'); %plot coloured area covering extent of box
        hold on;
        plot([min(profile),max(profile)],[Y_pos,Y_pos], 'b', 'LineWidth', 2); %plot position of center of box
        plot(profile, top_of_area:bot_of_area, 'k-o',  'MarkerFaceColor', 'w'); %plot intensity profile long box.
        set(gca,'YDir','reverse');
        axis tight;
        set(gca,'YTick',round(top_of_area:5:bot_of_area));
        title('Intensity profile through box and above/below it')
        hold off;
        
        try
            suptitle(['Box ',num2str(n)])%, '; range = ', num2str(profile_range)])
        end
            
        %Dialoge box to change position of Box.
        verbage = strcat('Y position of box ',num2str(n),': if no information input ''none''');
        prompt = {verbage, 'Box height'};
        dlg_title = 'Position of foil';
        num_lines = 1;
        def = {num2str(Y_pos), num2str(Y_height)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if isempty(answer) == 1
            % Construct a questdlg to catch the cancel option
            choice = questdlg('Did you mean to canel?', ...
                'Canelling?', ...
                'Yes','No','Quit box picking','Yes');
            % Handle response
            switch choice
                case 'Yes'
                    done = 1;
                    close(figure(3))
                case 'No'
                    done = 0;
                case 'Quit box picking'
                    close all
                    return
            end
        else
            position_out = str2num(answer{1,1});
            HeightBox = str2num(answer{2,1});
            
            pos_diff = position_out - Y_pos;
            height_diff = HeightBox - Y_height;
            display_rounding = 1e-3;
            
            if isempty(position_out) == 1
                done = 1;
                empty_boxes = [empty_boxes n];
                empty_boxes = sort(empty_boxes);
            elseif abs(pos_diff) < display_rounding && abs(height_diff) < display_rounding 
                % this condition is used to allow for the rounding to 1e-4 when the position is displayed by the script.
                % the rounding is such that it doesn't materially affect to box positions...
               
                done = 1;
                boxY(n,1) = position_out - HeightBox/2;
                boxY(n,2) = position_out + HeightBox/2;
            else
                Y_pos = position_out;
                Y_height = HeightBox;
            end
        end
    end
end
% end


close(repos_figure_handle)

end %of main function


%% box processing functions.

function out = empty_box_2_string(box_list)
%turns list of boxes from numbers into string with commas and hyphens
if length(box_list) == 1
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
%turns string of input numbers into list of numbers
dividers = find(box_string == ',');
if box_string(1) == 'n'  %no change to the option in the dialogue box returns 'none' and the string compare has to be the same length
    empty_boxes = 0;
    % elseif isempty(dividers) == 1
    %     empty_boxes = str2num(box_string);
else
    dividers = [1 dividers length(box_string)+1]; %adds leading zero and length to the variable so that the following for loop works.
    empty_boxes = [];
    for i = 1 : length(dividers)-1
        num = box_string(dividers(i):dividers(i+1)-1);
        hyphen_location = find(num == '-');
        if isempty(hyphen_location)
            empty_boxes = [empty_boxes str2num(num)];
        else
            start = str2num(num(1:hyphen_location-1));
            finish = str2num(num(hyphen_location+1:end));
            empty_boxes = [empty_boxes start:finish];
        end
    end
    empty_boxes = sort(unique(empty_boxes));
end
out = empty_boxes;
end

