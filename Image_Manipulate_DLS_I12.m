% Image_Manipulate_DLS_I12
%   Performs the manipulating of the images (rotation inversion etc) on images from 
%   beamline I12 at Diamond Light Source. 
%   ImageAnalysis script. 
% 
%   Syntax:
%     images_out = Image_Manipulate_DLS_I12('type','images_in')
%       where:
%       - 'type':      image type (tiff, netcdf...) 
%       - 'images_in': images to be manipulated.
%
%       See also ImageAnalysis, ImageAnalysis_Functions_Tiff

%   Simon Hunt 2018
%   $ version: 0.1 $ 15th October 2018 $
%     - This script was copied from Image_Manipulate_X17B2
%     - for changes see the end of the file. 

function Im_out = Image_Manipulate_DLS_I12(im_type,Im_in)

switch im_type
    case {'nxs'} %manipulate nxs images.
        
        
        cam = evalin('caller', 'cam');%Image_Functions_NXSlist('camera', data.textdata{1});
        
        if strcmpi(cam, 'camFloat1')
            % Matrix has to be transposed and fliped
            for j = 1 :  size(Im_in,3)
                im_rot = Im_in(:,:,j);
                [m,n] = size(im_rot);
                Im_out(:,:,j)  = im_rot.';
                %
                %Im_out(:,:,j) = im_rot(n:-1:1,1:1:m);
                %             im_rot  =
                %            Im_out(:,:,j) = im_rot.';
            end
        elseif strcmpi(cam, 'camFloat2')
            % Matrix has to be transposed and fliped
            for j = 1 :  size(Im_in,3)
                im_rot = Im_in(:,:,j);
                [m,n] = size(im_rot);
                im_rot  = im_rot.';
                Im_out(:,:,j) = im_rot(n:-1:1,1:1:m);
                %             im_rot  =
                %            Im_out(:,:,j) = im_rot.';
            end
            
        else
            
            error('The camera type is not recognised')
        end
        %Im_out = Im_in;
       
    otherwise
        error('I12:notype', ['The image type ''', im_type, ''' is not recognised for Diamond I12 images.']) 

end

end