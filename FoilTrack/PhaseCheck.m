% PhaseCheck
% Checks that the periods in the sin fit files are similar at each nominal
% period and if not recalcualtes the erronious ones them to make the
% consistent.
%
% Syntax: 
%      PhaseCheck[Both]('option', 'setting'...)
%   possible options:
%       'type'        - 'disp' or 'ssd'. sets the default options for the
%       functions and files.
%       'file_str'   - override the default sine fit file ending
%       'to_process' - override the default data file ending
%       'funct'      - override the default Processing function (Phases).
%       'data dir'   - override the default data file location of the current directory.
%
% Simon Hunt 2010 - 2018.
%  Version 2.3.1, June 2018.

function PhaseCheck(varargin)

%process varargin and set defaults if there is no option.
type = 'unknown';
plots_on = 1;
directory = '.';
iarg = 1 ;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case 'type'
            type = varargin{iarg+1};
            iarg = iarg + 2;
        case 'to_process'
            to_process = varargin{iarg+1};
            iarg = iarg + 2;
        case 'file_str'
            file_str = varargin{iarg+1};
            iarg = iarg + 2;
        case 'funct'
            funct = varargin{iarg+1};
            iarg = iarg + 2;
        case 'data dir'
            directory = varargin{iarg+1};
            iarg = iarg + 2;
        case 'plotoff'
            plots_on = 0;
            iarg = iarg + 1;
        otherwise
            error(['Unknown option: ' varargin{iarg}]) ;
    end
end

%if the type is not set the determine it.
if strcmpi(type, 'unknown') == 1
    
    allfiles = dir('*sine*');

    for x = 1:length(allfiles)
        file = allfiles(x).name;
        info = FileTitleInformation(file);
        ending{x} = info{7};
    end
        
    endings = unique(ending);
   
    if numel(endings) == 0
        error('There are no recognised ''*sine*'' files in the target directory.')
    elseif numel(endings) > 1
        error('There is more then one type of ''*sine*'' file in the target directory.')
    end

    if ~isempty(strfind(lower(endings), 'ssd'))
        type = 'ssd';
    else
        type = 'disp';
    end
end


% process type to set file strings
if strfind(lower(type), 'disp') == 1 % not SSD therefore normal displacements
    
    %string of files to check.
    if exist('file_str','var') == 0
        file_str = '*sine_fits.txt';
        not_str = 'SSD';
    end
    %string for data files. (ending)
    if exist('to_process','var') == 0
        to_process = 'position_change.txt';
    end
    %string for data files. (ending)
    if exist('funct','var') == 0
        funct = 'Phases';
    end
    
else %default SSD settings
    
    %string of files to check.
    if exist('file_str','var') == 0
        file_str = '*SSD_sine_fits.txt';
        not_str = '';
    end
    %string for data files. (ending)
    if exist('to_process','var') == 0
        to_process = 'SSD.mat';
    end
    %string for data files. (ending)
    if exist('funct','var') == 0
        funct = 'PhasesSSD';
    end
    
end

%% sort and read phase files.

%get the files.
allfiles = dir;
allnames = {allfiles.name}';
file_str = strrep(file_str, '*', '\w*');
not_str  = strrep(not_str, '*', '\w*');
k = cellfun(@(s) ~isempty(regexpi(s, file_str)), allnames);
d = cellfun(@(s) isempty(regexpi(s, not_str)), allnames);

allnames = allnames(d.*k==1);
number_files = length(allnames);
file_names = allnames;

%allfiles = dir(file_str);
%for x = 1:length(allfiles)
%allfiles(allfiles.names
%number_files = length(allfiles);
%file_names = {allfiles.name}';

if number_files == 0
    error('There are no identified files to work with!')
end

% get nominal loads, temperatures and periods from the file names
% types = [0 0 0]; %this is used to give a warning if temp turns into power
% during the experiments
for x = 1 : number_files
    file = file_names{x};
    info = FileTitleInformation(file);
    expt_name = info{1};
    load(x) = info{2};
    temp(x) = info{3};
    period(x) = info{4};
    run_num(x) = info{5};
end

%sort to be in the order the data was acquired in
[run_num order] = sort(run_num); %sort by run number
file_names = file_names(order);
load = load(order);
temp = temp(order);
period = period(order);

%correction for some data set. 
wrong = find(temp == 44374);
temp(wrong) = 280;

%fix incase data is not present in the file titles.
period(isnan(period)) = -1;
load(isnan(load)) = -1;
temp(isnan(temp)) = -1;


%unique values.
period(isnan(period)) = 0;
nominal_periods = unique(period);
nominal_loads = unique(load);
nominal_temps = unique(temp);

% reads in all data files listed and does some initial processing of the data.
for i = 1:number_files
    in = file_names{i};
    all_data(i) = importdata(in,',');
    all_numbers(i,:,:) = all_data(i).data;
end

% gets periods from all_numbers    
periods(:,:) = all_numbers(:,3,1);

%determine what to plot the data against
if sum(isnan(temp)) ~= numel(temp) && numel(unique(temp)) ~= 1 %all the temperatures are present
    plot_against = temp;
    x_lab = 'Temperature (�C)';
elseif sum(isnan(load)) ~= numel(load) && numel(unique(load)) ~= 1 %all the load are present
    plot_against = load;
    x_lab = 'Load';
else %plot against run number
    plot_against = 1:number_files;%num_num;
    x_lab = 'Run Number';
end
    
if plots_on == 1
    %figure
    fig = figure(1);
    clf
    colours = get(gca,'ColorOrder');
    hold on
    data_names = [];
    for x = 1:length(nominal_periods)
        to_plot = find(period == nominal_periods(x));
        h = plot(plot_against(to_plot),periods(to_plot),'o', 'MarkerEdgeColor',colours(x,:),'MarkerFaceColor',colours(x,:));
        data_names = [data_names; file_names(to_plot), repmat({h},numel(file_names(to_plot)),1)];
    end
    legend(num2str(nominal_periods'))
    xlabel(x_lab),
    ylabel('Calculated periods (seconds)')
    title('Periods as calculated')
    if max(nominal_periods) - min(nominal_periods) > 90
        set(gca, 'YScale', 'log')
    end
    %change content of data cursor
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,data_names})
    
    pause;
end

for i = 1:length(nominal_periods)
    recalculate = [];
    where_nominal = find(period == nominal_periods(i));
    where = where_nominal;
%     periods(where), pause
    j = 0;
    while j == 0
        average = mean(periods(where));
        deviation = std(periods(where));
        ratio = deviation/average;
        range = max(periods(where)) - min(periods(where));
        proportional_range = range / average;
        
%         plot(temp(where_nominal),periods(where_nominal),'ro'), pause

        % iterates until the range of the data is small enough that it
        % passes. Then the routine moves on and recalcualtes all necessary
        % parts.
        if proportional_range > 1e-2 %ratio >= 0.0005 | range >= 5e-4
            diff = [(abs(periods(where) - average))'; 1:length(where)'];
            diff_sorted = flipud(sortrows(diff',1))';            
            
            recalculate = [recalculate where(diff_sorted(2,1))];

            where(diff_sorted(2,1)) = []; %removes this line from where and continues
        else
            file_names{recalculate}
%             average
            j = 1;
            for k = 1 : length(recalculate)
                to_do = file_names{recalculate(k)};
                
                %rename SSD file so it can be kept.
                [~,p1,p2] = fileparts(to_do);
                copyfile([p1,p2], [p1(1:end-1),'OLD_', strtrim(sprintf('%8.0f', now)), p2])
                
                %cut = strfind(to_do, file_str(2:end-1));
                cut = strfind(to_do, strrep(file_str, '\w*',''));
                name_to_give =  [directory,filesep,to_do(1:cut-1), to_process];%strcat(to_do(1:cut-1), to_process);
                
                %make function name string into function name 
                Recalc = str2func(funct);
                
                %Recalculate the phases for the offending data.
                Recalc(name_to_give,'period',average)

                reload = importdata(to_do,',');
                reload_numbers = reload.data;
                periods(recalculate(k)) = reload_numbers(3,1);
            end
%             pause
        end
    end
    if plots_on == 1
        figure(2)
        plot(temp(where_nominal),periods(where_nominal),'ro'), xlabel('Temperature (�C)'), ylabel('Calculated periods (seconds)'), % pause
    end
end

fclose all;

if plots_on == 1 || 1
    fig2 = figure;
    colours = get(gca,'ColorOrder');
    hold on
    data_names = [];
    for x = 1:length(nominal_periods)
        to_plot = find(period == nominal_periods(x));
        h = plot(plot_against(to_plot),periods(to_plot),'o', 'MarkerEdgeColor',colours(x,:));
        data_names = [data_names; file_names(to_plot), repmat({h},numel(file_names(to_plot)),1)];
    end
    legend(num2str(nominal_periods'))
    xlabel(x_lab),
    ylabel('Calculated periods (seconds)')
    title('Updated periods')
    if max(nominal_periods) - min(nominal_periods) > 90
        set(gca, 'YScale', 'log')
    end
    %change content of data cursor
    dcm_obj = datacursormode(fig2);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,data_names})
end
% 
% figure(fig)
% plot(temp,periods,'o'), xlabel('Temperature (°C)'), ylabel('Calculated periods (seconds)')
% %change content of data cursor
% dcm_obj = datacursormode(fig);
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,file_names})
% 
% 
% if max(nominal_periods) - min(nominal_periods) > 90;
%     set(gca, 'YScale', 'log')
% end

end %funciton

function txt = myupdatefcn(~,event_obj,t)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
targ = get(event_obj, 'Target');

%obj = cell2mat(t(:,2));
obj = t(:,2);
%loc = find(obj == targ);
for x = 1:numel(obj), loc(x) = (obj{x} == targ); end
loc = find(loc==1);

txt = {[t{loc(I)}],...
       ['Calculated Period: ',num2str(pos(2)), 's']};
 %  disp('Something is wrong here and it does not necessarily present the right file name')
end

% v2.3.2  10th April 2019
%   - incremented number for fix for when tempereratures or loads are not present in the file name.
% v2.3.1 8 june 2018
%   - fixed bug in text of data fit for lab. - names now present properly.
% v2.3 - 7 & 8 July 2017
%   - Added settings so that the data is plotted against changing pressure or temperature, depending which one is present
%   - MAde changes so that can read just sine_fits and not the SSD files as
%   well.
% v2.2 - May 2017
%   - Added data cursor function that gives the name of the file for each
%   data point.
%   - Added titles to figures
% v2.1.1  - 6th November 2016
%   - Minor fix to the data directory option.
% v2.1  - 9th June 2016
%   - Added data directory option.
% v2.0  - April 2016
%   - Made useable for both displacement and SSD data sets. 