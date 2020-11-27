% SelectBoxes_MoveXY
%   Provides validation for the box position and allows individual boxes to be moved
%   around. The script is called by SelectBoxes.m. It should be fed the image that the boxes 
%   are required for and a of file names and it will lead one throught the creation 
%   of the boxes for image sequence.
%
%   syntax:
%     SelectBoxes_MoveXY(I1, boxX, boxY)
%          
% see also: SelectBoxes, ImageAnalysis

% Simon Hunt June 2016
% Used to be part of Select Boxes but was removed when SelectBoxes was
% rewritten (version 2).

function [boxX, boxY] = SelectBoxes_MoveXY(Im1, boxX, boxY)

%size of image.
[h1 ~] = size(Im1);

number_boxes = numel(boxX)/2;

%plot box and details.
repos_figure_handle = figure('Units','normalized','Position',[.25 0.0467 0.5 0.85], 'name','Reposition');
        
imagesc(Im1);
colormap('gray')
   axis square
   
%plot boxes on image
for n = 1 : number_boxes
    rectangle('position',[boxX(n,1)-0.5,boxY(n,1)-0.5,boxX(n,2)-boxX(n,1)+1,boxY(n,2)-boxY(n,1)+1],'EdgeColor','r');
end



oldX = NaN;
oldY = NaN;
while or(~isequal(oldX, boxX) ,~isequal(oldY, boxY))
    
    oldX = boxX;
    oldY = boxY;
    
    %prompt for inputs to move all boxes
    prompt = {'Offset Y (up/down)?', 'Offset X (left/right)?'};
    dlg_title = 'Global positioning';
    num_lines = 1;
    def = {num2str(0), num2str(0)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        moveY = str2num(answer{1,1});
        moveX = str2num(answer{2,1});
    else
        moveY = 0;
        moveX = 0;
    end
    if moveY == 0 && moveX == 0 %if nothing has moved escape from this loop.
        oldX = boxX;
        oldY = boxY;
    else
        %reposition the boxes
        boxX = boxX + moveX;
        boxY = boxY + moveY;
        
        %plot the new boxes
        for n = 1 : number_boxes
            rectangle('position',[boxX(n,1)-0.5,boxY(n,1)-0.5,boxX(n,2)-boxX(n,1)+1,boxY(n,2)-boxY(n,1)+1],'EdgeColor','g');
        end
        
        % Gives dialogue box asking if the new boxes are close enough...
        choice = questdlg('Close enough??', ...
            'How near?', ...
            'New Boxes','Yes','No','Yes');
        % Handle response
        switch choice
            case 'Yes'
                oldX = boxX;
                oldY = boxY;
            case 'No'
                for n = 1 : number_boxes
                    rectangle('position',[boxX(n,1)-0.5,boxY(n,1)-0.5,boxX(n,2)-boxX(n,1)+1,boxY(n,2)-boxY(n,1)+1],'EdgeColor','y');
                end
        end
    end
end


close(repos_figure_handle)
