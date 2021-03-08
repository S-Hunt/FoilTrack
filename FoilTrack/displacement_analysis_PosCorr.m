
function [offset_x , offset_y , Xfull , Yfull] = displacement_analysis_PosCorr(image1 , image2 , Xfull , Yfull)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlates the relative position of the images. It does this by reducing1 imageanalysis.m
% the size of the image into k * k boxes to speed the calculation and then
% doing the same for the best part at full resolution
[h_in w_in ~] = size(image1);
shrink = 10;

width = floor(w_in/shrink);
height = floor(h_in/shrink);

Ishrunk1 = zeros(height ,width);
Ishrunk2 = zeros(height ,width);

%makes the images for the reduced resolution fit
for n = 1 : width
    for m = 1 : height
        Ishrunk1(m,n) = mean(mean(image1(shrink*m - shrink + 1 : shrink*m , shrink*n - shrink + 1 : shrink*n)));
        Ishrunk2(m,n) = mean(mean(image2(shrink*m - shrink + 1 : shrink*m , shrink*n - shrink + 1 : shrink*n)));
    end
end
if exist('Xfull','var') == 0
    figure (1); clf;
    subplot(1, 2, 1); imagesc(Ishrunk1); axis image; title('image1');
    subplot(1, 2, 2); imagesc(Ishrunk2); axis image; title('image2');
    if exist('Xfull','var') == 0;
        [Xshrunk , Yshrunk] = ginput(2);
        Xshrunk = round(Xshrunk);
        Yshrunk = round(Yshrunk);
        Xfull = Xshrunk * shrink;
        Yfull = Yshrunk * shrink;
        % save data Xfull Yfull
    end
else
    Xshrunk = round(Xfull / shrink);
    Yshrunk = round(Yfull / shrink);
end

%correlates the positions in the shrunken images
corr_width = Xshrunk(2) - Xshrunk(1);
corr_height = Yshrunk(2) - Yshrunk(1);
disp_y = height - corr_height;
disp_x = width - corr_width;

area1 = zeros(corr_height , corr_width);
area2 = zeros(corr_height , corr_width);
value1 = zeros(disp_y , disp_x);

area1 = Ishrunk1(Yshrunk(1):Yshrunk(2), Xshrunk(1):Xshrunk(2));

for i = 1 : disp_y
    for j = 1 : disp_x
        area2 = Ishrunk2(i : i + corr_height , j : j + corr_width);
        value1(i,j) = myxcorr(area1 , area2);
    end
end

[miny,minx] = find(value1 == min(min(value1)));
Xscale = zeros(width -1);
Yscale = zeros(height -1);
Xscale = -Xshrunk(1) + 1 : width - Xshrunk(2);
Yscale = -Yshrunk(1) + 1 : height - Yshrunk(2);
%correlates the positions in the full scale images using the same relative
%position as in the shrunk images
offset_x = Xscale(minx) * shrink;
offset_y = Yscale(miny) * shrink;

shift_y2 = offset_y;
shift_y2old = .1;
shift_x2 = offset_x;
shift_x2old = .1;

area1 = zeros(Yfull(2)-Yfull(1), Xfull(2)-Xfull(2));
area2 = zeros(Yfull(2)-Yfull(1), Xfull(2)-Xfull(2));

area1 = image1(Yfull(1):Yfull(2), Xfull(1):Xfull(2));

while or(shift_y2 ~= shift_y2old , shift_x2 ~=shift_x2old);
    for i = 1 : 2*shrink+1
        for j = 1 : 2*shrink+1
            area2 = image2((Yfull(1)+i-shrink -1+shift_y2):(Yfull(2)+i-shrink -1+shift_y2),...
                (Xfull(1)+j-shrink -1+shift_x2):(Xfull(2)+j-shrink -1+shift_x2));
            value2area(i,j) = myxcorr(area1 , area2);
        end
    end
    [miny2 ,minx2] = find(value2area == min(min(value2area)));
    if miny2 == 1 | miny2 == 2
        shift_y2 = shift_y2 - shrink;
    elseif miny2 == or(2*shrink+1, 2*shrink)
        shift_y2 = shift_y2 + shrink;
    elseif minx2 == 1 | minx2 == 2
        shift_x2 = shift_x2 - shrink;
    elseif minx2 == or(2*shrink+1, 2*shrink)
        shift_x2 = shift_x2 + shrink;
    else
        shift_y2old = shift_y2;
        shift_x2old = shift_x2;
    end
end
Xscale = -shrink+shift_x2 : shrink+shift_x2;
Yscale = -shrink+shift_y2 : shrink+shift_y2;

offset_x = Xscale(minx2);
offset_y = Yscale(miny2);

%generates a R2 plot of the misfit and shows where the minimum point is.
% figure(fig); fig = fig + 1; clf;
% contour(Xscale , Yscale , value1 , 15); hold on;
% plot(bestXfit , bestYfit , 'sk','MarkerEdgeColor ','k','MarkerFaceColor ','k','MarkerSize ',8);

%adds a rectangle to the plot of the images to show where the optimised fit
%is on the second of the two images
figure (1); %hold off;
subplot(1, 2, 1); imagesc(image1); rectangle('position',[Xfull(1), Yfull(1), Xfull(2)-Xfull(1), Yfull(2)-Yfull(1)],'EdgeColor','k'); axis image; title('image1');
subplot(1, 2, 2); imagesc(image2); rectangle('position',[Xfull(1)+offset_x , Yfull(1)+offset_y , Xfull(2)-Xfull(1), Yfull(2)-Yfull(1)],'EdgeColor','k'); axis image; title('image2');



% cross correlation (xcorr) subfunction
function value = myxcorr(area1 , area2)
value = sum(sum((area1 -area2).^2))/sum(sum(area1));