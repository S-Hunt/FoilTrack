% Image_Manipulate_APS_16IDBxrd
%   Performs the manipulating of the images (rotation inversion etc) on the XRD patterns from 
%   the APS's 16-ID-D beamline. 
%   ImageAnalysis script. 
% 
%   Syntax:
%     images_out = Image_Manipulate_APS_13IDD('type','images_in')
%       where:
%       - 'type':      image type (tiff, netcdf...) 
%       - 'images_in': images to be manipulated.
%
%       See also ImageAnalysis, ImageAnalysis_Functions_Tiff

%   Simon Hunt 2017
%   $ version: 0.1 $ 7th June 2017 $
%     - This script was removed from ImageAnalysis version 1.9
%     - for changes see the end of the file. 

function Im_out = Image_Manipulate_APS_16IDBxrd(im_type,Im_in)

switch im_type
    case {'tif', 'tiff'} %manipulate tiff images.

        %Nothing to do here.
        % this is left as a comment
        
        Im_out = Im_in;

        % case {'nc' 'netcdf'} %manipulate netcdf images.
       
    otherwise
        error('16IDBxrd:notype', ['The image type ''', im_type, ''' is not recognised for 16-ID-B xrd images.']) 

end

end