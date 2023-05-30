% ImageAnalysis_Functions_Tiff
%   Performs the reading and organising of tiff images for the
%   ImageAnalysis script. 
% 
%   Syntax:
%     output = ImageAnalysis_Functions_Tiff('function','varargin')
%     The output of the function varies depending on the 'function'. The full set of options is a follows:
%       file_id     = ImageAnalysis_Functions_Tiff('open','FILE NAME')
%       timestamps  = ImageAnalysis_Functions_Tiff('times','file_id')
%       uids        = ImageAnalysis_Functions_Tiff('unique','file_id')
%       [max_frames image_size(1) image_size(2)]...
%                   = ImageAnalysis_Functions_Tiff('NumSize','file_id')
%       images      = ImageAnalysis_Functions_Tiff('getimage',file_id, expt_location, ref_id, num_images_wanted, image_size)
%       temperature = ImageAnalysis_Functions_Tiff('temperature','file_id')
%       status      = ImageAnalysis_Functions_Tiff('close','file_id')

%
%   For the single image / file formats the File Name opened has to be a
%   *.lst file made by MakeSingleLstFile
%       See also ImageAnalysis, MakeSingleLstFile, MakeManyTimes

%   Simon Hunt 2016
%   $ version: 0.1 $ 6th May 2016 $
%     - This script was removed from ImageAnalysis version 1.9.1
%     - for changes see the end of the file. 

function varargout = Image_Functions_Tiff(funct,varargin)

switch funct
    case 'open'
        % opens Tiff list file
        % this is not strictly needed but here for compatability with the other file
        % types
        
        if exist(varargin{1},'file') == 2
            varargout{1} = varargin{1};
        else
            error(['The file ', varargin{1}, ' does not exist.'])
        end
        
    case 'times'
        % pull timestamps from the list of tiff files
        file_id = varargin{1};

        data = importdata(file_id);
        varargout{1} = data.data(:,1);

    case 'unique'
        % pull timestamps and unique ID's from the list of tiff files
        file_id = varargin{1};

        data = importdata(file_id);
        varargout{1} = data.textdata;

    case 'NumSize'
        %get number of images and dimensions of the image array
        file_id = varargin{1};

        data = importdata(file_id);
        max_frames = length(data.data);
        
        if data.textdata{1}(end)==','
            data.textdata{1}=data.textdata{1}(1:end-1)
        end
        
        Im = imread(data.textdata{1});
        [image_size(1) image_size(2) ~] = size(Im);

        varargout(1:3) = {max_frames image_size(1) image_size(2)};

    case 'getimage'
        %get image from file
        file_id = varargin{1};
        manipulate_function = str2func(['Image_Manipulate_',varargin{2}]);
        start = varargin{3};
        num_frames = varargin{4};
        image_size = varargin{5};

        data = importdata(file_id);
       
        Im = zeros([image_size(1) image_size(2) num_frames]);
        %original code
%         for x = 1 : num_frames
%             if ~isempty(strfind(varargin{2},'ESRF'))
%                 Im_temp = importdata(data.textdata{start+x-1});
%                 if isstruct(Im_temp) %catch error incase image type is a bmp
%                     Im(:,:,x) = Im_temp.cdata;
%                 else
%                     Im(:,:,x) = Im_temp;
%                 end
%             else
%                 Im(:,:,x) = imread(data.textdata{start+x-1});
%             end
%         end

        %modified code for 16IDB images
        %modificed again so can read the dark field images.
        for x = 1 : num_frames
            if isnumeric(data)
                Im_temp = data;
            else
                if data.textdata{start+x-1}(end)==','
                    data.textdata{start+x-1}=data.textdata{start+x-1}(1:end-1);
                end
                try
                    Im_temp = imread(data.textdata{start+x-1});
                    % disp('imread')
                catch err_con
                    Im_temp = importdata(data.textdata{start+x-1});
                    % disp('importdata')
                end
            end
            if x==1 & size(Im_temp,1) == size(Im,2)
                clear Im
                Im = zeros([image_size(2) image_size(1) num_frames]);
            end
            if isstruct(Im_temp) %catch error incase image type is a bmp
                Im(:,:,x) = Im_temp.cdata;
            elseif ndims(Im_temp) == 3 %error catch for RGB images
                if isequal(Im_temp(:,:,1), Im_temp(:,:,2)) %error catch if the image in actually in colour
                    Im(:,:,x) = Im_temp(:,:,1);
                else
                    error('The tiff images appear to be in colour. This code cannot deal with that!')
                end
            else %'normal' 1 colour image
                Im(:,:,x) = Im_temp;
            end
        end

        varargout{1} = manipulate_function('tiff', Im);

%         if ~isempty(strfind(cd,'NSLS'))
%             varargout{1} = manupilate_image(Im);
%         elseif ~isempty(strfind(cd,'6BMB'))
%             varargout{1} = manupilate_image_6BMB(Im);
%         elseif ~isempty(strfind(cd,'DESY'))
%             varargout{1} = manupilate_image_DESY(Im);
%         elseif ~isempty(strfind(cd,'APS'))
%             varargout{1} = Im;
%         elseif ~isempty(strfind(cd,'ESRF'))
%             varargout{1} = Im;
% %            varargout{1} = manupilate_image_APS_PVT(Im);
%         else
%             error('Experiment location is not recognised')
%         end

    case 'temperature'
        % pull timestamps from the list of tiff files
        file_id = varargin{1};

        data = importdata(file_id);
        if size(data.data,2) >= 2
            varargout{1} = data.data(:,2);
        else
            varargout{1} = NaN(size(data.data,1),1);
        end
        
    case 'power'
        % pull timestamps from the list of tiff files
        file_id = varargin{1};

        data = importdata(file_id);
        if size(data.data,2) >= 3
            varargout{1} = data.data(:,3);
        else
            varargout{1} = NaN(size(data.data,1),1);
        end        
        
    case 'close'
        % close the image files.

        %not necessary but needs to be here for completeness

    otherwise
        error('This is not a recognised image property process')
end

%     function Im_out = manupilate_image(Im_in)
%         %manipulates the images -- written as a subfucntion so that the same thing
%         %has to be done to both images.
%         Im_in = single(Im_in);
%         
%         % truncate images to 8 bit. Mostly for NSLS side station camera.
%         % It does need to be made into a separate option.
%         large = Im_in > 2^8;
%         Im_in(large) = 2^8;%2^8-1;
% 
%         % Matrix has to be transposed and fliped
% 
%         for j = 1 :  size(Im_in,3)
%             im_rot = Im_in(:,:,j);
%             [~,n] = size(im_rot);
%             im_rot  = im_rot.';
%             Im_out(:,:,j) = im_rot(n:-1:1,:);
%         end
% 
%     end
% 
%     function Im_out = manupilate_image_DESY(Im_in)
%         %manipulates the images -- written as a subfucntion so that the same thing
%         %has to be done to both images.
% 
%         % Matrix has to be rotated so axis parallel to pixels
%         Im_out = imrotate(Im_in, -54.74, 'bilinear');
%     end
% 
%     function Im_out = manupilate_image_6BMB(Im_in)
% 
%         % Matrix has to be transposed and fliped
% 
%         for j = 1 :  size(Im_in,3)
%             im_rot = Im_in(:,:,j);
%             [m,n] = size(im_rot);
%             im_rot  = im_rot.';
%             Im_out(:,:,j) = im_rot(1:1:n,m:-1:1);
%         end
% 
%     end
% 
%     function Im_out = manupilate_image_APS_PVT(Im_in)
%         %manipulates the images -- written as a subfucntion so that the same thing
%         %has to be done to both images.
% 
%         % Matrix has to be rotated so axis parallel to pixels
%         for r = 1:size(Im_in,3);
%             Im_out(:,:,r) = flipud(rot90(Im_in(:,:,r)));
%             
%         end
%     end
end %tiff_func

%% versions
% v 1.2 - 2nd July 2018
%   - if no temperature in the list file then return NaN and added options to return power.
% v 1.1 - 16th June 2016
%   - Changed image loading code to deal with 16IDB images and to make it
%   non-location specific. i.e. it should just read the images.
%   The location specific functions are now in 'Image_Manipulate_...'
% v 1.0 - May 2016
%   - code removed from ImageAnalysis and image manipulation removed to
%   locaton specific files.