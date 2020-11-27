% Image_Manipulate_ESRF_ID06LVP
%   Performs the manipulating of the images (rotation inversion etc) on images from 
%   beamlines at the old DESY max200x beamline. 
%   ImageAnalysis script. 
% 
%   Syntax:
%     images_out = Image_Manipulate_ESRF_ID06LVP('type','images_in')
%       where:
%       - 'type':      image type (tiff, netcdf...) 
%       - 'images_in': images to be manipulated.
%
%       See also ImageAnalysis, ImageAnalysis_Functions_Tiff

%   Simon Hunt 2016
%   $ version: 0.1 $ 3rd June 2016 $
%     - This script was removed from ImageAnalysis version 1.9
%     - for changes see the end of the file. 

function Im_out = Image_Manipulate_ESRF_ID06LVP(im_type,Im_in)

switch im_type
    case {'tif', 'tiff'} %manipulate tiff images.
        
        %no rotations or inversions are required for this location.         
        Im_out = Im_in;
        
    case {'edf'}
        
        Im_in = single(Im_in);
        
        % Matrix has to be transposed and fliped        
        for j = 1 :  size(Im_in,3)
            Im_out(:,:,j) = rot90(Im_in(:,:,j), -1 );
        end
        
    otherwise
        error('DESY:notype', ['The image type ''', im_type, ''' is not recognised for DESY Max200x images.']) 

end

end