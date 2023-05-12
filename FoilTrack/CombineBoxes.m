% CombineBoxes
% Called by Phases or PhasesSSD to combine the boxes according to switches.
%
% syntax: [combined_data1, combined_data2, reference_length, horizontal_position_pixels, varargout]...
%                 = CombineBoxes(combine_type, data_type, data1, data2, data3, data4, options)
%           combine_type      -- 'anelastic', 'pair', 'single'
%           data_type         -- 'displacements', 'phases' or 'test'
%           data1             -- displacement array or phase array if present errors as second column
%           data2             -- time stamp array or amplitude array if present errors as second column
%           data3             -- box X positions 
%           data4             -- box Y positions 
%  options: 'getrid', value   -- remove unwanted row of data if the numner of rows in uneven.
%           'sym_type', value -- 'symmetric' or 'antisymmetric' -- how to combine the rows
%           
%   $ version: 0.2 $ 19th October 2018 $

% settings:
%   - combine type: anelastic, pair (symmetric, antisymmetric, remove_rows), single

%from phases:
%   - displacements, time stamps
%   - position of boxes

% from PhasesSSD
%   - phase, amplitude,
%   - position of boxes


function [combined_data1, combined_data2, reference_length, horizontal_position_pixels, varargout] = CombineBoxes(combine_type, data_type, data1, data2, Xbox, Ybox, varargin)



% parse the inputs.
% possible: get_rid symmetric, AmplitudesError, PhasesError.

%process varargin and set defaults if there is no option.
get_rid = 'none';
RefDisps = [];
symm_type = 'symmetric';
iarg = 1 ;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case 'get_rid'
            get_rid = varargin{iarg+1};
            iarg = iarg + 2;
        case 'symm_type'
            symm_type = varargin{iarg+1};
            iarg = iarg + 2;
        case 'amplitudeerror'
            AmplitudesError = varargin{iarg+1};
            iarg = iarg + 2;
        case 'phaseerror'
            PhasesError = varargin{iarg+1};
            iarg = iarg + 2;
        case 'refoffsets'
            RefDisps = varargin{iarg+1};
            iarg = iarg + 2;
        otherwise
            error(['Unknown option: ' varargin{iarg}]) ;
    end
end

% validate the inputs.

%check box array are the same size.
    if size(Xbox) ~= size(Ybox)
        error('The box arrays are not the same size')
    end
%check the box arrays are n by 2 in size
    if size(Ybox,2) ~= 2
        error('The box arrays are not correct')
    end
%check the combination type is licit
    if strcmpi(combine_type, 'phases')==1
        if size(Ybox,1) ~= size(data1,1)
            error('Input data is the wrong size')
        end
    elseif strcmpi(combine_type, 'displacements')==1
        if size(Ybox,1) ~= size(data1,2)
            error('Input data is the wrong size')
        end
    end

number_boxes_in = size(Ybox,1);

%check the number of rows in the data.
if number_boxes_in == 1
    combine_type = 'single';
    boxes_in_row = number_boxes_in;
    num_rows = 1;
else
    gaps = Xbox(2:end,1) - Xbox(1:end-1,2);
    row_ends = find(gaps<0);
    
    num_rows = length(row_ends) + 1;
    boxes_in_row =  size(Xbox,1)/num_rows;

    %error catching code
    if num_rows == 1 & strcmpi(combine_type, 'single')==0
        error('If there is only one row of data the combine type has to be ''single''.');
    end    
end

% make sure RefOffsets is the smae size as the data set.
if numel(RefDisps) ~= number_boxes_in
   RefDisps = zeros(number_boxes_in,1);
end
if size(RefDisps, 1) == 1
    RefDisps = RefDisps';
end

%% boxes to match
%make array of numbers for boxes corresponding to the numbering position in
%the images.
box_array = 1:number_boxes_in;
box_array = reshape(box_array,boxes_in_row, num_rows)';


switch combine_type
    %this switch makes an array of 2 by n which lists in each column the boxes to be combined.
    % if the array has one row then no combinations are necessary.
    
    
    % Phases -- values that need to be organised for each case are: phases, phases_err, amps, amp_err and SSD and reference length
    % displcements -- values to be organised are position, SSD? and refernce length.

    case {'anelastic', 'deformation'}
        % array to be made needs to be 
        % combine = [ 1 2 3 4 5 ... n-1;
        %             2 3 4 5 6 ... n  ]
        
        %assumes that there is only one column of boxes -- which is not
        %necessarily true
        if size(box_array,2) ~= 1
            error('The anelastic combination will not work for this data set')
        end
        
        combination_array = [1:number_boxes_in-1; 2:number_boxes_in]; %FIX ME: this code assumes that there is only one column of boxes to combine.
    
    case 'single'
        %array to be made is a single row of 1:length of data array.
        combination_array = 1:number_boxes_in;
        
    case 'pair'
        %array to be made depends on the symmetric/assymetric switch and the rows to be removed (if any)
        
        % checks that the number of rows is even (1 if even otherwise 0)
        % if the number of rows is not even remove a row of the data.
        even_rows = isequal(num_rows/2,round(num_rows/2));      
        if even_rows == 0 && ~strcmpi(get_rid, 'none')
           switch get_rid
               case 'end'
                   togo = size(box_array,1);
                   
               case 'middle'
                   togo = (size(box_array,1)+1)/2;
                   
               case {'1' '2' '3' '4' '5' '6' '7'}
                   togo = str2num(get_rid);

               otherwise
                   error('The uneven row removal mechanism is not defined')
           end %get_rid 
           box_array(togo,:) = [];
           
        end %even_rows == 0
           
        
        %section of array to reverse.
        rev = size(box_array,1)/2 +1;
        
        %if the data is to be combined asymmetrically revese the order of the bottom half of the array. 
        switch symm_type
            case 'symmetric'
                % this is left blank delibarately. nothing needs to be done to the array
                
            case 'asymmetric'
                box_array(rev:end,:) = fliplr(box_array(rev:end,:));
                
            otherwise
                   error('The symmetry to be used in the combination is not defined.')
        end %symm_type
        
        %make into a 2 x n array
        arr_top = box_array(1:rev-1,:)';
        arr_bot = fliplr(box_array(rev:end,:)');
        arr_top = arr_top(:);
        arr_bot = arr_bot(:);
        combination_array = [arr_top'; arr_bot'];
        
    otherwise
        error('The pairing type is not recognised')
        
end


%% combine the boxes.

number_pairs = size(combination_array,2);
midpointsY = mean(Ybox, 2);
midpointsX = mean(Xbox, 2);

if strcmpi(data_type, 'test') == 0 && strcmpi(combine_type, 'single') == 1
    combined_data1 = data1;
    combined_data2 = data2;
    reference_length = midpointsY;
    horizontal_position_pixels = midpointsX;
    if length(varargin) >= 2
        varargout = {PhasesError};
    end
    if length(varargin) >= 4
        varargout = [varargout, {AmplitudesError}];
    end
    return
end

switch data_type
    case 'test' %used for testing the combination array.
        combined_data1 = combination_array;
        combined_data2 = [];
        reference_length = [];
        horizontal_position_pixels = [];
        varargout = {[], []};
        return
        
    case 'displacements' %combine the displacement array. Called by Phases.m
        
        foil_change = data1;
        
        combined_data1 = foil_change(:,combination_array(2,:)) - foil_change(:,combination_array(1,:));
        combined_data2 = data2;
     
        
        
    case 'phases' %combine the phase array. Called by PhasesSSD.m
        
        if min(size(data1)) ~= 1 %if the array is not a column or row vector
            error('The data is not in the correct format. It is required to be a single row or column array')
        else %if the data array is a column or a row. -- it is as it should be.
            Amplitudes = data2;
            Phases = data1;
        end
        
        for x = 1: number_pairs
            
            %get the two data locations to be combined
            point1 = combination_array(1,x);
            point2 = combination_array(2,x);
            
            %combine the data
            if exist('AmplitudesError','var') == 1
                [combined_data2(x), combined_data1(x), combined_amplitude_error(x), combined_phase_error(x)] =...
                    Sinusoid_Difference(Amplitudes(point1), Phases(point1), Amplitudes(point2), Phases(point2), AmplitudesError(point1), PhasesError(point1), AmplitudesError(point2), PhasesError(point2));
            else
                [combined_data2(x), combined_data1(x)] =...
                    Sinusoid_Difference(Amplitudes(point1), Phases(point1), Amplitudes(point2), Phases(point2));
            end
            
        end
        
        if exist('combined_amplitude_error', 'var')
            varargout = {combined_phase_error, combined_amplitude_error};
        else 
            varargout = {[],[]};
        end
        
    otherwise
        error('Data combination type is unknown')
        
end

%% box position arrays
%reference_length = midpointsY(combination_array(2,:)) - midpointsY(combination_array(1,:));

%correction for phase angle of reference.
reference_length = (midpointsY(combination_array(2,:)) - RefDisps(combination_array(2,:)) ) - (midpointsY(combination_array(1,:)) - RefDisps(combination_array(1,:)) );

horizontal_position_pixels = mean([midpointsX(combination_array(1,:)), midpointsX(combination_array(2,:))],2);

reference_length = reference_length';
horizontal_position_pixels = horizontal_position_pixels';
 

%% version history.
%   - v0.2 - 19th October 2018
%       - Added correction to the lengths for the displacemnt at the
%       reference time of the driving wave.
%   - v0.1.1 - 3rd July 2017
%       - Bug fix so that single boxes (ie no cominbations) will run.
%   - v0.1 - 13th February 2016
%       - written function

