% Image_Functions_NXS
%   Performs the reading and organising of NSX images for the
%   ImageAnalysis script. 
% 
%   Syntax:
%     output = Image_Functions_NXS('function','varargin')
%     The output of the function varies depending on the 'function'. The full set of options is a follows:
%       file_id    = Image_Functions_NXS('open','FILE NAME')
%       timestamps = Image_Functions_NXS('times','file_id')
%       uids       = Image_Functions_NXS('unique','file_id')
%       [max_frames image_size(1) image_size(2)]...
%                  = Image_Functions_NXS('NumSize','file_id')
%       images     = Image_Functions_NXS('getimage',file_id, expt_location, image_type, first_image, num_images_wanted, image_size)
%       status     = Image_Functions_NXS('close','file_id')
%
%   The script assumes that the nxs file does not contain a time series of images.
%   If there is more than one image in the file they are added together (for intnsities)
%   and the time stamps are averaged.
%
%   The script currently only recognises CamFloat2.
%
%   To process *.nxs files as time series from a single file then call the
%   nsx file directly from ImageAnalysis. 
%
%   For the single image / file formats the File Name opened has to be a
%   *.nxslst file made by MakeSingleLstFile
%       See also ImageAnalysis, MakeSingleLstFile, MakeManyTimes

%   Simon Hunt 2018
%   $ version: 0.1 $ 30th September 2018 $
%     - This script was a copy of Image_Functions_NetCDF version 0.1
%     - for changes see the end of the file. 

function varargout = Image_Functions_NXSlist(process,varargin)

% 
% %check the camera exists if present as an option.
% if nargin > 2
%    if strcmpi(varargin{2}, 'camera') ~= 1
%        error('This is suposed to have a camera definition in here')
%    end
%    camera = varargin{3};
%    
%    try 
%        h5info(varargin{1}, ['/entry1/instrument/',camera]);
%    end
% end


switch process
    case 'open'
        % opens Tiff list file
        % this is not strictly needed but here for compatability with the other file
        % types
        
        if exist(varargin{1},'file') == 2
            varargout{1} = varargin{1};
        else
            error(['The file ', varargin{1}, ' does not exist.'])
        end
        
    case 'camera'
        
        %check a camera is defined in the call.
        %if strcmpi(varargin{2}, 'camera') ~= 1
        %    error('This is suposed to have a camera definition in here')
        %end
        %camera = varargin{3};
        
        
        % determine the camera type from the NXS file
        file_id = varargin{1};
        flds = h5info(file_id, '/entry1/instrument/');
        grps = flds.Groups;
        for x = 1:numel(grps) %list some of the bits of data in the file
            grp = grps(x).Name;
            l = strfind(grp, '/');
            item{x} = grp(l(end)+1:end);
        end
        
        %eliminate the bits that are not the camera.
        not_cam = {'actualTime';
            'earlyFramesIncluded';
            'ix';
            'ring';
            'source';
            't7_m3z'};
        ind = zeros(1,numel(item));
        for x = 1: numel(not_cam)
            index = cellfun(@(y) strcmp(y, not_cam{x}), item, 'UniformOutput', 1);
            ind = ind+index;
        end
        
        if sum(~ind) == 1
            varargout{1} = item{~ind};
        else
            error('Something is in the file that is not recognised')
        end
        
        
    case 'times'
        % pull timestamps from the NXS file
        file_id = varargin{1};
        
        [~,~,e] = fileparts(file_id);
        if strcmpi(e, '.nxslst')
            
            %if the list file is being read for times
            data = importdata(file_id);
            varargout{1} = data.data(:,1);

        else
            %if an individual frame be being read for times. 
           
            %need to get an average incase there is more than one frame in the 'image'
            varargout{1} = mean(h5read(file_id, '/entry1/instrument/actualTime/actualTime'));
        end
        
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

        cam = Image_Functions_NXSlist('camera', data.textdata{1});
        Im = h5read(data.textdata{1}, ['/entry1/instrument/',cam,'/data']);
        
        [image_size(1), image_size(2), ~] = size(Im);

        varargout(1:3) = {max_frames image_size(1) image_size(2)};


    case 'getimage'
        %get image from file
        file_id = varargin{1};
        manipulate_function = str2func(['Image_Manipulate_',varargin{2}]);
        image_id = varargin{3};
        num_frames = varargin{4};
        image_size = varargin{5};

        data = importdata(file_id);
        Im = zeros([image_size(1) image_size(2) num_frames]);
                
        
        cam = Image_Functions_NXSlist('camera', data.textdata{1});
        for x = 1 : num_frames
            %read images.
            Im(:,:,x) = sum(h5read(data.textdata{image_id+x-1}, ['/entry1/instrument/',cam,'/data']),3);
        end
        
        %imagesID = netcdf.inqVarID(file_id,'array_data');
        %varargout{1} = netcdf.getVar(file_id,imagesID,[0,0,image_id-1],[image_size(1),image_size(2),num_frames]);
        %varargout{1} = h5read(file_id, '/entry1/instrument/camFloat2/data');

        %manipulate the images.
        varargout{1} = manipulate_function('nxs', Im);


    case 'close'
        % close the image files.

        %not necessary but needs to be here for completeness

    otherwise
        error('This is not a recognised image property process')
end

end %Image_Functions_NXS

