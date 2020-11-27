% Image scale code for use with 
% Thermal diffusivity analysis driver code for Image file (primarily NetCDF) data
%             Simon Hunt 2014
% 
% Takes intensity profiles from a suite of images in which sample is moved
% and calculates pixel scale for the experiment.
% 
%   Syntax:
%   called by ImageAnalysis('FILE_NAME','scale')
%
%   $ version: 0.1.DEV $ 13th Jan 2014 $
%       - for details of changes between versions see end of the file.


function ImageScale(I_profile)

% exclude regions which are static in the images
plot(I_profile,'.-')
[exclusion_x,~] = ginput;
exclusion_x = round(exclusion_x);

for x = 1: 2: length(exclusion_x)
    I_profile(exclusion_x(x):exclusion_x(x+1),:) = NaN;
end

%region(s) to match in profile 1
plot(I_profile(:,1),'.-')

[match_x,~] = ginput(2);
match_x = round(match_x);
ref = [];
for x = 1: 2: length(match_x)
    ref = [ref; match_x(x):match_x(x+1)];
end

%minimise SSD differences.
for x = 1 : size(I_profile,2)

     ref_profile = I_profile(ref);
     
     for y = 1 : size(I_profile,1)-length(ref_profile)
         
         SSD(y) = nansum((I_profile(y:y+length(ref_profile)-1 ,x) - ref_profile').^2);
         
     end
     
     [~, min_SSD(x)] = min(SSD);
     min_SSD(x) = ref(1) - min_SSD(x);
    
end
    
% get distance press moved between images. 
prompt = {'Distance press was moved between images (micrometers):'};
dlg_title = 'pixel scale factor';
num_lines = 1;
def = {'100'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
step_size = str2num(answer{1});
distances = (0:length(min_SSD)-1) * step_size;


fitobject = fit(distances',min_SSD','poly1','robust','bisquare')

plot(distances, min_SSD,'o'), hold on
plot(fitobject,'r'),
legend off
xlabel('Distance moved (microns)')
ylabel('Distance moved (pixels)')    , hold off
     


display(['Scale: ',num2str(fitobject.p1), ' pixels/micron']);

end