% Image_Manipulate_DESY_max200x
%   Performs the manipulating of the images (rotation inversion etc) on images from 
%   beamlines at the old DESY max200x beamline. 
%   ImageAnalysis script. 
% 
%   Syntax:
%     images_out = Image_Manipulate_DESY_max200x('type','images_in')
%       where:
%       - 'type':      image type (tiff, netcdf...) 
%       - 'images_in': images to be manipulated.
%
%       See also ImageAnalysis, ImageAnalysis_Functions_Tiff

%   Simon Hunt 2016
%   $ version: 0.1 $ 3rd June 2016 $
%     - This script was removed from ImageAnalysis version 1.9
%     - for changes see the end of the file. 

function Im_out = Image_Manipulate_DESY_max200x(im_type,Im_in)

switch im_type
    case {'tif', 'tiff'} %manipulate tiff images.
        
        % Matrix has to be rotated so axis parallel to pixels
        Im_out = imrotate(Im_in, -54.74, 'bilinear');
        
     %nocase for netcdf exists for this experiment location.
        
    otherwise
        error('DESY:notype', ['The image type ''', im_type, ''' is not recognised for DESY Max200x images.']) 

end

end