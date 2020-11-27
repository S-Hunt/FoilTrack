% Image_Functions_NetCDF
%   Performs the reading and organising of NetCDF images for the
%   ImageAnalysis script. 
% 
%   Syntax:
%     output = Image_Functions_NetCDF('function','varargin')
%     The output of the function varies depending on the 'function'. The full set of options is a follows:
%       file_id    = Image_Functions_NetCDF('open','FILE NAME')
%       timestamps = Image_Functions_NetCDF('times','file_id')
%       uids       = Image_Functions_NetCDF('unique','file_id')
%       [max_frames image_size(1) image_size(2)]...
%                  = Image_Functions_NetCDF('NumSize','file_id')
%       images     = Image_Functions_NetCDF('getimage',file_id, expt_location, image_type, first_image, num_images_wanted, image_size)
%       status     = Image_Functions_NetCDF('close','file_id')
%
%   For the single image / file formats the File Name opened has to be a
%   *.lst file made by MakeSingleLstFile
%       See also ImageAnalysis, MakeSingleLstFile, MakeManyTimes

%   Simon Hunt 2016
%   $ version: 0.1 $ 6th May 2016 $
%     - This script was removed from ImageAnalysis version 1.9.1
%     - for changes see the end of the file. 

function varargout = Image_Functions_NetCDF(process,varargin)

switch process
    case 'open'
        %opens NetCDF file
        file_name = varargin{1};
        varargout{1} = netcdf.open(file_name,'NC_NOWRITE');

    case 'times'
        % pull timestamps from the NetCDF file
        file_id = varargin{1};
        varid = netcdf.inqVarID(file_id, 'timeStamp');
        varargout{1} = netcdf.getVar(file_id, varid);

    case 'unique'
        % pull unique ID's from the NetCDF file
        file_id = varargin{1};
        varid = netcdf.inqVarID(file_id,'uniqueId');
        varargout{1} = netcdf.getVar(file_id, varid);

    case 'NumSize'
        %get number of images and dimensions of the image array
        file_id = varargin{1};
        [~, max_frames] = netcdf.inqDim(file_id,0);
        [~, image_size(1)] = netcdf.inqDim(file_id,2);
        [~, image_size(2)] = netcdf.inqDim(file_id,1);

        varargout(1:3) = {max_frames image_size(1) image_size(2)};

    case 'getimage'
        %get image from file
        file_id = varargin{1};
        manipulate_function = str2func(['Image_Manipulate_',varargin{2}]);
        image_id = varargin{3};
        num_frames = varargin{4};
        image_size = varargin{5};

        imagesID = netcdf.inqVarID(file_id,'array_data');
        varargout{1} = netcdf.getVar(file_id,imagesID,[0,0,image_id-1],[image_size(1),image_size(2),num_frames]);

        
        varargout{1} = manipulate_function('netcdf', varargout{1});


    case 'close'
        % close the image files.
        file_id = varargin{1};
        netcdf.close(file_id);

    otherwise
        error('This is not a recognised image property process')
end

end %netcdf_func




%     function Im_in = manupilate_image(Im_in)
%         %manipulates the images -- written as a subfucntion so that the same thing
%         %has to be done to both images.
% 
%         [m,n,o] = size(Im_in);
% 
%         % Data is unsigned int but NetCDF thinks it is signed -typecast changes type without
%         % changing bit arrangement. This is much faster than the previous methods of doing it
%         % 'manually' (v1.4a and before)
% 
%         % data is then converted to single precision number so can be used
%         % in displacement analysis properly .. cannot have negative uint8 number.
% 
%         % single is significantly faster than double (around 40% in test I did) becuase the amount of
%         % memory needed to store data is half the size smaller. When
%         % opening file for first time.
% 
%         if isa(Im_in,'int8') == 1
%             Im_in = single(typecast(Im_in(:),'uint8'));
%         elseif isa(Im_in,'int16') == 1
%             Im_in = single(typecast(Im_in(:),'uint16'));
%         else
%             error('Unknown data type for images');
%         end
%             
%         Im_in = reshape(Im_in,m,n,o);
% 
%         if ~isempty(strfind(lower(cd),lower('NSLS')))
%             % Matrix has to be transposed and fliped
%             % Im_in = Im_in';
%             % Im_in = fliplr(Im_in);
%             % Im_in = flipud(Im_in);
%             
%             %this is a quicker version of the lines of code above
%             if ndims(Im_in) == 2
%                 [m, n] = size(Im_in);
%                 Im_in = Im_in(m:-1:1,n:-1:1);
%             else
%                 [m, n, ~] = size(Im_in);
%                 Im_in = Im_in(m:-1:1,n:-1:1,:);
%             end
%         elseif ~isempty(strfind(lower(cd),lower('6BMB')))
%             % Matrix has to be fliped. No transpose is required.
%             % Im_in = fliplr(Im_in);
%             
%             if ndims(Im_in) == 2
%                 Im_in = Im_in(:,end:-1:1);
%             else
%                 Im_in = Im_in(:,end:-1:1,:);
%             end
% 
%         elseif ~isempty(strfind(cd,'APS'))
%             % Matrix has to be transposed
%           
%             if ndims(Im_in) == 2
%                 Im_in = Im_in';
%             else
%                 for j = 1 :  size(Im_in,3)
%                     Im_rot = Im_in(:,:,j);
%                     Im_rot = Im_rot.';
%                     Im_out(:,:,j) = Im_rot;
%                 end
%                 Im_in = Im_out;
%             end
%             
%         else
%             error('Experiment location is not recognised')
%         end
%     end %manipulate_image