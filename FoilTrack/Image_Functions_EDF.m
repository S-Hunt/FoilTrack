% Image_Functions_EDF
%   Performs the reading and organising of EDF images from the ESRF for the
%   ImageAnalysis script. 
% 
%   Syntax:
%     output = Image_Functions_EDF('function','varargin')
%     The output of the function varies depending on the 'function'. The full set of options is a follows:
%       file_id    = Image_Functions_EDF('open','FILE NAME')
%       timestamps = Image_Functions_EDF('times','file_id')
%       uids       = Image_Functions_EDF('unique','file_id')
%       [max_frames image_size(1) image_size(2)]...
%                  = Image_Functions_EDF('NumSize','file_id')
%       images     = Image_Functions_EDF('getimage',file_id, expt_location, ref_id, num_images_wanted, image_size)
%       status     = Image_Functions_EDF('close','file_id')
%
%   For the single image / file formats the File Name opened has to be a
%   *.lst file made by MakeSingleLstFile
%       See also ImageAnalysis, MakeSingleLstFile, MakeManyTimes

%   Simon Hunt 2016
%   $ version: 0.1 $ 6th May 2016 $
%     - This script was removed from ImageAnalysis version 1.9.1
%     - for changes see the end of the file. 

function varargout = Image_Functions_EDF(funct,varargin)

switch funct
    case 'open'
        % opens Tiff list file
        varargout{1} = varargin{1};

    case 'times'
        % pull timestamps from the list of tiff files
        file_id = varargin{1};

        data = importdata(file_id);
        varargout{1} = data.data;

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

        
        [~, Im] = pmedf_read(data.textdata{1}); %imread(data.textdata{1});
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
        for x = 1 : num_frames
            [~, Im_temp] = pmedf_read(data.textdata{start+x-1});
            Im(:,:,x) = Im_temp;
        end

        varargout{1} = manipulate_function('edf', Im);
%         if ~isempty(strfind(cd,'ESRF'))
%             varargout{1} = manupilate_image(Im);
%         else
%             error('Experiment location is not recognised')
%         end

    case 'close'
        % close the image files.

        %not necessary but needs to be here for completeness

    otherwise
        error('This is not a recognised image property process')
end

%     function Im_out = manupilate_image(Im_in)
% 
%         %manipulates the images -- written as a subfucntion so that the same thing
%         %has to be done to both images.
%         Im_in = single(Im_in);
%         
%         % Matrix has to be transposed and fliped
% %         for j = 1:size(Im_in,3)
% %             Im_out = rot90(Im_in, -1 );
% %         end
% 
%         for j = 1 :  size(Im_in,3)
%             Im_out(:,:,j) = rot90(Im_in(:,:,j), -1 );
%         end
%         
%     end %manipulate image

end %edf_func
