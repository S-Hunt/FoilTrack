function [ out ] = ImageAnalysis_SSDio(proc, file_name, varargin)
%ImageAnalysis_SSDwrite. Reads/Writes SSD output files from ImageAnalysis.
%   Saves SSD values calculated by ImageAnalysis, along with the assiciated 
%   location and time values. 

switch proc
    case 'write'

        SSDvals = varargin{1};
        positionvals = varargin{2};
        timestamps = varargin{3};
        AnalysisInformation = varargin{4};
        
        %save the data in a mat file.
        save(file_name, 'SSDvals', 'positionvals', 'timestamps', 'AnalysisInformation')

        
    case 'read'
        %there are two formats for the SSD array files. The 'old' and 'new'
        %formats. 
        %The old format has the arrays: 'pos_array', 'SSD_array',
        %'image_time', 'boxX', 'boxY' and a full set of analysis options stored.
        %While new has arrays: 'SSDvals', 'positionvals', 'timestamps',
        %'Parameters' or 'headers'.
        
        data = load(file_name);        
        
        if numel(varargin) ~= 0
            warning('The listing of variables to read is not currently supported by ImageAnalysis_SSDio')
        end
        
        % FIX ME: The output arrays are static. It should be possible to
        % use varargin to force the output data to have the right fields.
        if isfield(data, 'SSDvals')
           
           out.SSD_array  = data.SSDvals;
           out.pos_array  = data.positionvals;
           out.image_time = data.timestamps;
           if isfield(data, 'Parameters')
               out.boxX       = data.Parameters{:}.boxes.boxX;
               out.boxY       = data.Parameters{:}.boxes.boxY;
           else
               out.boxX       = data.AnalysisInformation.boxes.boxX;
               out.boxY       = data.AnalysisInformation.boxes.boxY;
           end
           
        elseif isfield(data, 'SSD_array') %old data format. 
            
            out = data;
            
        else
            error('Ooops. This did not work');
        end

        
    case 'reformat'  
        
        
    otherwise
        error({['Unrecognsed process: ',proc,'.'];'Valid processes are ''read'' and ''write''.'});
end

