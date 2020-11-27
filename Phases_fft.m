function [ input_period varargout ] = Phases_fft( data_array, time_array, varargin)
%Phases_fft. Extract most likely period for ImageAnalysis data by fft. 
%   Given the displacment or SSD array from ImageAnalysis the script
%   calculates the most likely period by finding the maximum power of the
%   fft. 
%   Syntax:
%     [period, bigger, smaller, phases, amplitudes] = Phases_fft(data_array, image_time)
%         data_array:  array of displacements or SSD from ImageAnalysis.m
%         image_time:  array of time stamps associated with the data_array
%
%         period:      Most likelt period derived from the fft.
%         bigger:      Larger likely value for period from fft.
%         smaller:     Smaller likely value for period from fft.
%         phases:      Phase(s) of the most likely period. 
%         amplitudes:  Amplitudes of the most likely period.
%
%  Simon Hunt. 2016.
%  Version 0.1
          

fit_type = 'poly';
fig_plot = 0;
% resid    = 0;
make_plot = 0;
iarg = 1 ;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case {'poly', 'power', 'old'}
            fit_type = varargin{iarg};
            iarg = iarg + 1;
        case 'plotfigs'
            fig_plot = 1;
            iarg = iarg + 1;
%         case 'residuals'
%             resid = 1;
%             iarg = iarg + 1;
        otherwise
            error(['Unknown option: ' varargin{iarg}]) ;
    end
end

number_boxes = size(data_array,3);

%% fourier transform
%frequencies
fs = 1/((time_array(1,end)-time_array(1,1))/size(time_array,2));

for x = 1 : number_boxes % process the data for each box individually
    
    % fft data
    Z = data_array(:,:,x);
    if length(Z)/2^nextpow2(size(Z,2)) > .9
        nfft = 2^nextpow2(size(Z,2)); % Next power of 2 from length of Z
    else
        nfft = length(Z);
    end
    Yn = fft(Z',nfft)/nfft;
    
    % extract data from fourier transform
    f   = (0:size(Yn,1)-1).*(fs/nfft); %frequencies
    Pyy = Yn.* conj(Yn); % power
    ang = angle(Yn);     % phase
    A   = (Yn);       % amplitude in freq domain
    
    %half the length of the data to remove repititions
    f   = f(:,1:round(end/2));
    Pyy = Pyy(1:round(end/2),:);
    ang = ang(1:round(end/2),:);
    A   = A(1:round(end/2),:);
    
    % store the data for analysis when all is processed
    all_Pyy(x,:,:) = Pyy;
    all_ang(x,:,:) = ang;
    all_A(x,:,:)   = A;
    
    if fig_plot == 1
        %plot fft if required
        figure%(1)
        subplot(3,1,1)
        plot(f',Pyy); title('Frequency content of data'); xlabel('Frequency (Hz)'); ylabel('Power'); %keyboard
        
        %figure(2)
        subplot(3,1,2)
        plot(f', A); title('Amplitude content of data'); xlabel('Frequency (Hz)'); ylabel('Amplitude in frequency domain'); %keyboard
        
        %figure(3)
        subplot(3,1,3)
        plot(f', ang); title('Phase of data'); xlabel('Frequency (Hz)'); ylabel('Phase in frequency domain'); %keyboard
        
        % keyboard
    end
    
end

%cut ends from data...
cut_percent = 0;
cut = ceil(size(A,1)*cut_percent/100);
if cut == 0, cut = 2; end

%mean power in all fft.
mean_Pyy = squeeze(nanmean(nanmean(all_Pyy,1),3));

%% find period
options = optimset('MaxFunEval', 1E9, 'MaxIter', 1E9);

% for x = 1 : number_boxes
    
%     Pyy_use = mean(all_Pyy(x,:,:),3);
    Pyy_use = mean_Pyy;

    if strcmpi(fit_type, 'old') == 1
        
        %cut ends from data...
        cut_percent = 3;
        cut = ceil(size(A,1)*cut_percent/100);
        
        %mean power in all fft.
        mean_Pyy = squeeze(mean(mean(all_Pyy(:,cut:end,:),1),3));
        f = f(cut:end);
        
        %gradient of mean power at each frequency.
        dPyydf = mean_Pyy(2:end) - mean_Pyy(1:end-1);
        min_dfdp = find(min(dPyydf) == dPyydf, 5);
        max_dfdp = find(max(dPyydf) == dPyydf, 5);
        
        best_guess = min_dfdp; %use minimum in gradient as best frequency
        %this is because the minimum in the gradient is in the same place
        %as the maximum in the power in the test data set (Stsh_006...)
        
        best_guess = find(max(mean_Pyy)==mean_Pyy)+cut-1;
        
    elseif strcmpi(fit_type, 'power') == 1
        power_coefs = fminsearch(@model_misfit_power, [0 0 0], options, f(cut:end), Pyy_use(cut:end));
        mod = power_model(f(cut:end), power_coefs);
        resid = Pyy_use(cut:end) - mod;
        
        %curvature (use as a miniumum frequency cut off)
        K = curvature(power_coefs, f);
        Kmax = find(max(K) == K);
        
        best_guess = Kmax + find(max(resid(Kmax:end))==resid(Kmax:end))+cut-2;
        
    elseif strcmpi(fit_type,'poly') == 1
       
        pp = polyfit(log10(f(cut:end)), log10(Pyy_use(cut:end)),3);
        mod = 10.^(polyval(pp, log10(f(cut:end))));
        
        resid = log10(Pyy_use(cut:end)) - log10(mod);
        
        best_guess = find(max(resid(cut:end))==resid(cut:end))+cut;
        
    else
        error('Fit_type is unrecognised')
    end

    if best_guess == length(f)
        best_guess = best_guess - 2;
    end
input_period = 1./f(best_guess);
    
%define range of periods based on gaps in fft
if size(data_array,2) > 700
    range = 1; %used to be 2
else
    range = 1;
end

bigger = 1/f(best_guess+range);
smaller = 1/f(best_guess-range);

%collate phases of fft to return 
input_phase = mean(rem(squeeze(all_ang(:,best_guess,:))+pi,pi),2);
        % all_ang(:,best_guess,:) -- gets phases at best frequency (see above)
        % squeeze(...)            -- removes the extra dimensions from the data
        % rem(... +pi,pi)         -- removes pi phase shift in data 
        % mean(...)               -- gets average. 
        
%amplitudes are not exctracted because the amplitudes in the SSD do not
%correspond simply to the displacement of the foils
        
%outputs        
varargout = {bigger smaller input_phase};


%debugging function. plots mean power of FFT.
if fig_plot == 1
    figure, clf
    plot(f(cut:end),(mean_Pyy(cut:end)), 'b.-'), hold on,
    if exist('mod','var') 
        plot(f(cut:end),(mod), 'g-')
    end
    plot(f(best_guess-cut+1), (mean_Pyy(best_guess-cut+1)), 'ro')
    title({'FFT of data'},'Interpreter','none'); xlabel('Frequency (Hz)'); ylabel('log Power');
    set(gca,'Yscale', 'log')
    set(gca,'Xscale', 'log')
    
    % keyboard
    %     figure(2), clf
    %     plot(f,K)
    %     hold on
    %     plot(f(Kmax), K(Kmax), 'or')
    %     title({'Curvature of fitted line'},'Interpreter','none'); xlabel('Frequency (Hz)'); ylabel('Curvature');
    
    pause(1)
end


end

 function [sum_sq residuals] = model_misfit_power(coeff, freq, data)

    model = power_model(freq, coeff);
    residuals = data-model;
    residuals = log10(data)-log10(model);
    
    sum_sq = sum(residuals.^2); 

 end

function vals = power_model(freq, coefficients)
    
    a = coefficients(1);
    b = coefficients(2);
    c = coefficients(3);
%     d = coefficients(4);
    
    vals = a.*freq.^b + c;
    
end


function K = curvature(coeffs, freq)

%curvature is defined as:
%K = (d2y/dx2) / [1 + (dy/dx)^2]^(3/2)

%curvature of power function
numerator = coeffs(1) .* coeffs(2) .* (coeffs(2)-1) .* freq.^(coeffs(2)-2);

demonenator = ( 1 + ( coeffs(1) .* coeffs(2) .* freq.^(coeffs(2)-1) ).^2 ).^(3/2);

K = abs(numerator)./demonenator;

end
