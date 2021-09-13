% Image_Manipulate_generic
%   Performs the manipulating of the images on images from generic
%   beamlines. there is no rotation or flipping of the images.
% 
%   Syntax:
%     images_out = Image_Manipulate_generic('type','images_in')
%       where:
%       - 'type':      image type (tiff, netcdf...) 
%       - 'images_in': images to be manipulated.
%
%       See also ImageAnalysis, ImageAnalysis_Functions_Tiff

%   Simon Hunt 2021
%   $ version: 0.1 $ 13th Sept 2021 $
%     - This script was removed from ImageAnalysis version 1.9
%     - for changes see the end of the file. 

function Im_out = Image_Manipulate_generic(im_type,Im_in)

switch im_type
    case {'tif', 'tiff'} %manipulate tiff images.

        %there are apparently no rotations or inversions to do here. 
        Im_out = Im_in;

    case {'nc' 'netcdf'} %manipulate netcdf images.

        %get the size of the image array
        [m,n,o] = size(Im_in);

        % Data is unsigned int but NetCDF thinks it is signed -typecast changes type without
        % changing bit arrangement. This is much faster than the previous methods of doing it
        % 'manually' (v1.4a and before)

        % data is then converted to single precision number so can be used
        % in displacement analysis properly .. cannot have negative uint8 number.

        % single is significantly faster than double (around 40% in test I did) becuase the amount of
        % memory needed to store data is half the size smaller. When opening file for first time.

        if isa(Im_in,'int8') == 1
            Im_in = single(typecast(Im_in(:),'uint8'));
        elseif isa(Im_in,'int16') == 1
            Im_in = single(typecast(Im_in(:),'uint16'));
        else
            error('Unknown data type for images');
        end
            
        Im_in = reshape(Im_in,m,n,o);
        
        % Matrix has to be transposed
        if ndims(Im_in) == 2
            Im_in = Im_in';
        else
            for j = 1 :  size(Im_in,3)
                Im_rot = Im_in(:,:,j);
                Im_rot = Im_rot.';
                Im_out(:,:,j) = Im_rot;
            end
            Im_in = Im_out;
        end
        Im_out = Im_in;
       
    otherwise
        error('X17B2:notype', ['The image type ''', im_type, ''' is not recognised for generic-beamline images.']) 

end

end