% Image_Manipulate_APS_13IDD
%   Performs the manipulating of the images (rotation inversion etc) on images from 
%   beamlines at the APS's 13-ID-D beamline. 
%   ImageAnalysis script. 
% 
%   Syntax:
%     images_out = Image_Manipulate_APS_13IDD('type','images_in')
%       where:
%       - 'type':      image type (tiff, netcdf...) 
%       - 'images_in': images to be manipulated.
%
%       See also ImageAnalysis, ImageAnalysis_Functions_Tiff

%   Simon Hunt 2016
%   $ version: 0.1 $ 3rd June 2016 $
%     - This script was removed from ImageAnalysis version 1.9
%     - for changes see the end of the file. 

function Im_out = Image_Manipulate_APS_16IDB(im_type,Im_in)

switch im_type
    case {'tif', 'tiff'} %manipulate tiff images.

        % Matrix has to be rotated
        for j = 1 :  size(Im_in,3)
            im_rot = Im_in(:,:,j);
            Im_out(:,:,j) = im_rot.';
        end
        
        
%         Im_out = Im_in';
        
        warning('The images could be upside down. The rotation is arbitary')

%     case {'nc' 'netcdf'} %manipulate netcdf images.
       
    otherwise
        error('16IDB:notype', ['The image type ''', im_type, ''' is not recognised for 16-ID-B images.']) 

end

end