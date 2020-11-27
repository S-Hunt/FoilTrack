% FFTImageFilter.
%  Provides function to clean systematic background from Images for the ImageAnalysis processing pipeline. 
%
% Takes images and FFT-filter. Tranforms the image into it's FFT inverse,
% filters the FFT inverse and then returns a restored image.
% 
%   Syntax:
%   FilteredImage, differenceImage = FFTImageFilter(Image,FFT-filter)
%
%   Arguments:
%   1: Image - Unfiltered Image.
%   2: FFT-filter - Binary image same size as unflitered image for filter.
%
%   See also AnalysisOptions, displacement_analysis
%
%   Simon Hunt 2019
%   $ version: 1.0 $ 21st October 2019 $
%       - for details of changes between versions see end of the file.

function varargout = FFTImageFilter(Image, FFT_filter, varargin)  

if nargin == 3
    filt_out = varargin{1};
end


%make sure the image is a double. 
Image = double(Image);

% Fourier transform the image
Fim = fftshift(fft2(Image));

%filter the FFT-Im using the map
% FIX ME: there should be a switch in here for inside or outside the
% filter.
%if filt_out == 0
%    FFTimage_filtered = fftshift(Fim).*~FFT_filter;
%else
    FFTimage_filtered = (Fim).*FFT_filter;
%end

filtered_im = real(ifft2(ifftshift(FFTimage_filtered)));

varargout{1} = filtered_im;
varargout{2} = Image - filtered_im;