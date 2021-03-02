% ImageAnalysis_Manipulate_X17B2
%   Performs the manipulating of the images (rotation inversion etc) on images from 
%   beamlines at the old X17B2 beamlines. 
%   ImageAnalysis script. 
% 
%   Syntax:
%     images_out = IA_Manipulate_X17B2('type','images_in')
%       where:
%       - 'type':      image type (tiff, netcdf...) 
%       - 'images_in': images to be manipulated.
%
%       See also ImageAnalysis, ImageAnalysis_Functions_Tiff

%   Simon Hunt 2016
%   $ version: 0.1 $ 6th May 2016 $
%     - This script was removed from ImageAnalysis version 1.9
%     - for changes see the end of the file. 

function Im_out = Image_Manipulate_NSLS_X17B2ss(im_type,Im_in)

switch im_type
    case {'tif', 'tiff'} %manipulate tiff images.
        
        %force the data format for the images.
        Im_in = single(Im_in);
        
        % truncate images to 8 bit. Mostly for NSLS side station camera.
        % FIX ME: This should be made into a separate option.
        %large = Im_in > 2^8;
        %Im_in(large) = 2^8;%2^8-1;

        % Image has to be transposed and fliped
        for j = 1 :  size(Im_in,3)
            im_rot = Im_in(:,:,j);
            [~,n] = size(im_rot);
            im_rot  = im_rot.';
            Im_out(:,:,j) = im_rot(n:-1:1,:);
        end

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

        % Matrix has to be transposed and fliped
        % Im_in = Im_in';
        % Im_in = fliplr(Im_in);
        % Im_in = flipud(Im_in);
        
        %this is a quicker version of the lines of code above
        if ndims(Im_in) == 2
            [m, n] = size(Im_in);
            Im_out = Im_in(m:-1:1,n:-1:1);
        else
            [m, n, ~] = size(Im_in);
            Im_out = Im_in(m:-1:1,n:-1:1,:);
        end
        
       
    otherwise
        error('X17B2:notype', ['The image type ''', im_type, ''' is not recognised for X17B2 images.']) 

end

end