%define paths
realPaths = {
    'C:\Users\user\Desktop\norml image\Imtaiz Sir.jpg', ....
    'C:\Users\user\Desktop\norml image\Susmita Sen.jpg',....
    'C:\Users\user\Desktop\norml image\Imtaiz Sir.jpg',....
};

% Define image paths for AI Generated Images
AiPaths = {
    'C:\Users\user\Desktop\ai image\ai1.jpg', ....
    'C:\Users\user\Desktop\ai image\ami nije.jpg',....
    'C:\Users\user\Desktop\ai image\andre.jpg',
};



% Initialize arrays and cell arrays to store properties for Real Images
meanValuesReal = zeros(1, numel(realPaths));
stdDevValuesReal = zeros(1, numel(realPaths));
entropyValuesReal = zeros(1, numel(realPaths));
rmsValuesReal = zeros(1, numel(realPaths));
varianceValuesReal = zeros(1, numel(realPaths));
smoothnessValuesReal = zeros(1, numel(realPaths));
kurtosisValuesReal = zeros(1, numel(realPaths));
skewnessValuesReal = zeros(1, numel(realPaths));
idmValuesReal = zeros(1, numel(realPaths));
contrastValuesReal = zeros(1, numel(realPaths));
correlationValuesReal = zeros(1, numel(realPaths));
energyValuesReal = zeros(1, numel(realPaths));
homogeneityValuesReal = zeros(1, numel(realPaths));
histValuesReal = cell(1, numel(realPaths));
powerSpectrumReal = cell(1, numel(realPaths));
eccentricityReal = zeros(1, numel(realPaths));
solidityReal = zeros(1, numel(realPaths));
snrValuesReal = zeros(1, numel(realPaths));
hueValuesReal = zeros(1, numel(realPaths));
saturationValuesReal = zeros(1, numel(realPaths));
contrastValues2Real = zeros(1, numel(realPaths));
ambianceValuesReal = zeros(1, numel(realPaths));
brightnessValuesReal = zeros(1, numel(realPaths));
numPixelsReal = zeros(1, numel(realPaths));
pixelRateReal = zeros(1, numel(realPaths));

for i = 1:numel(realPaths)
    % Read the image
    inputImage = imread(realPaths{i});

    % Convert to grayscale if necessary
    if size(inputImage, 3) == 3
        grayImage = rgb2gray(inputImage);
    else
        grayImage = inputImage;
    end

    % Calculate Mean for Real Images
    meanValuesReal(i) = mean(double(grayImage(:)));

    % Calculate Standard Deviation for Real Images
    stdDevValuesReal(i) = std(double(grayImage(:)));

    % Calculate Entropy for Real Images
    entropyValuesReal(i) = entropy(grayImage);

    % Calculate Root Mean Square (RMS) for Real Images
    rmsValuesReal(i) = sqrt(mean(double(grayImage(:)).^2));

    % Calculate Variance for Real Images
    varianceValuesReal(i) = var(double(grayImage(:)));

    % Calculate Smoothness for Real Images
    smoothnessValuesReal(i) = 1 - (1 / (1 + varianceValuesReal(i)));

    % Calculate Kurtosis for Real Images
    kurtosisValuesReal(i) = kurtosis(double(grayImage(:)));

    % Calculate Skewness for Real Images
    skewnessValuesReal(i) = skewness(double(grayImage(:)));

    % Calculate Inverse Difference Moment (IDM) for Real Images
    glcm = graycomatrix(grayImage, 'Offset', [0 1; -1 1; -1 0; -1 -1]);
    idmValuesReal(i) = mean(1 ./ (1 + (glcm(:).^2)));

    % Calculate Contrast, Correlation, Energy, Homogeneity for Real Images
    stats = graycoprops(glcm, {'contrast', 'correlation', 'energy', 'homogeneity'});

    % Calculate Color Histogram for Real Images
    histValuesReal{i} = imhist(inputImage);

    % Calculate Power Spectrum for Real Images
    fftImage = fftshift(fft2(double(grayImage)));
    powerSpectrumReal{i} = abs(fftImage).^2;

    % Calculate Eccentricity and Solidity for Real Images
    binaryImage = imbinarize(grayImage);
    stats = regionprops(binaryImage, 'Eccentricity', 'Solidity');
    eccentricityReal(i) = stats.Eccentricity;
    solidityReal(i) = stats.Solidity;

    % Calculate Signal-to-Noise Ratio (SNR) for Real Images
    signal = mean(double(grayImage(:)));
    noise = std(double(grayImage(:)));
    snrValuesReal(i) = 20 * log10(signal / noise);

    % Calculate Hue, Saturation, Contrast, Ambiance, Brightness for Real Images
    hsvImage = rgb2hsv(inputImage);
    hueValuesReal(i) = mean2(hsvImage(:,:,1));
    saturationValuesReal(i) = mean2(hsvImage(:,:,2));
    value = hsvImage(:,:,3);
    contrastValues2Real(i) = std2(value);
    ambianceValuesReal(i) = mean2(value);
    brightnessValuesReal(i) = mean2(value);

    % Calculate Number of Pixels and Pixel Rate for Real Images
    numPixelsReal(i) = numel(grayImage);
    pixelRateReal(i) = snrValuesReal(i) / numPixelsReal(i);
end

% Initialize arrays and cell arrays to store properties for AI Generated Images
meanValuesAI = zeros(1, numel(AiPaths));
stdDevValuesAI = zeros(1, numel(AiPaths));
entropyValuesAI = zeros(1, numel(AiPaths));
rmsValuesAI = zeros(1, numel(AiPaths));
varianceValuesAI = zeros(1, numel(AiPaths));
smoothnessValuesAI = zeros(1, numel(AiPaths));
kurtosisValuesAI = zeros(1, numel(AiPaths));
skewnessValuesAI = zeros(1, numel(AiPaths));
idmValuesAI = zeros(1, numel(AiPaths));
contrastValuesAI = zeros(1, numel(AiPaths));
correlationValuesAI = zeros(1, numel(AiPaths));
energyValuesAI = zeros(1, numel(AiPaths));
homogeneityValuesAI = zeros(1, numel(AiPaths));
histValuesAI = cell(1, numel(AiPaths));
powerSpectrumAI = cell(1, numel(AiPaths));
eccentricityAI = zeros(1, numel(AiPaths));
solidityAI = zeros(1, numel(AiPaths));
snrValuesAI = zeros(1, numel(AiPaths));
hueValuesAI = zeros(1, numel(AiPaths));
saturationValuesAI = zeros(1, numel(AiPaths));
contrastValues2AI = zeros(1, numel(AiPaths));
ambianceValuesAI = zeros(1, numel(AiPaths));
brightnessValuesAI = zeros(1, numel(AiPaths));
numPixelsAI = zeros(1, numel(AiPaths));
pixelRateAI = zeros(1, numel(AiPaths));

for i = 1:numel(AiPaths)
    % Read the image
    inputImage = imread(AiPaths{i});

    % Convert to grayscale if necessary
    if size(inputImage, 3) == 3
        grayImage = rgb2gray(inputImage);
    else
        grayImage = inputImage;
    end

    % Calculate Mean for AI Generated Images
    meanValuesAI(i) = mean(double(grayImage(:)));

    % Calculate Standard Deviation for AI Generated Images
    stdDevValuesAI(i) = std(double(grayImage(:)));

    % Calculate Entropy for AI Generated Images
    entropyValuesAI(i) = entropy(grayImage);

    % Calculate Root Mean Square (RMS) for AI Generated Images
    rmsValuesAI(i) = sqrt(mean(double(grayImage(:)).^2));

    % Calculate Variance for AI Generated Images
    varianceValuesAI(i) = var(double(grayImage(:)));

    % Calculate Smoothness for AI Generated Images
    smoothnessValuesAI(i) = 1 - (1 / (1 + varianceValuesAI(i)));

    % Calculate Kurtosis for AI Generated Images
    kurtosisValuesAI(i) = kurtosis(double(grayImage(:)));

    % Calculate Skewness for AI Generated Images
    skewnessValuesAI(i) = skewness(double(grayImage(:)));

    % Calculate Inverse Difference Moment (IDM) for AI Generated Images
    glcm = graycomatrix(grayImage, 'Offset', [0 1; -1 1; -1 0; -1 -1]);
    idmValuesAI(i) = mean(1 ./ (1 + (glcm(:).^2)));

    % Calculate Contrast, Correlation, Energy, Homogeneity for AI Generated Images
    stats = graycoprops(glcm, {'contrast', 'correlation', 'energy', 'homogeneity'});

    % Calculate Color Histogram for AI Generated Images
    histValuesAI{i} = imhist(inputImage);

    % Calculate Power Spectrum for AI Generated Images
    fftImage = fftshift(fft2(double(grayImage)));
    powerSpectrumAI{i} = abs(fftImage).^2;

    % Calculate Eccentricity and Solidity for AI Generated Images
    binaryImage = imbinarize(grayImage);
    stats = regionprops(binaryImage, 'Eccentricity', 'Solidity');
    eccentricityAI(i) = stats.Eccentricity;
    solidityAI(i) = stats.Solidity;

    % Calculate Signal-to-Noise Ratio (SNR) for AI Generated Images
    signal = mean(double(grayImage(:)));
    noise = std(double(grayImage(:)));
    snrValuesAI(i) = 20 * log10(signal / noise);

    % Calculate Hue, Saturation, Contrast, Ambiance, Brightness for AI Generated Images
    hsvImage = rgb2hsv(inputImage);
    hueValuesAI(i) = mean2(hsvImage(:,:,1));
    saturationValuesAI(i) = mean2(hsvImage(:,:,2));
    value = hsvImage(:,:,3);
    contrastValues2AI(i) = std2(value);
    ambianceValuesAI(i) = mean2(value);
    brightnessValuesAI(i) = mean2(value);

    % Calculate Number of Pixels and Pixel Rate for AI Generated Images
    numPixelsAI(i) = numel(grayImage);
    pixelRateAI(i) = snrValuesAI(i) / numPixelsAI(i);
end

% Display properties for Real Images
% Display properties for Real Images
disp('Properties for Real Images:');
disp(' ');
for i = 1:numel(realPaths)
    disp(['Properties for Real Image------------>>', num2str(i), ':']);
    fprintf('Mean: %.2f\n', meanValuesReal(i));
    fprintf('Standard Deviation: %.2f\n', stdDevValuesReal(i));
    fprintf('Entropy: %.2f\n', entropyValuesReal(i));
    fprintf('Root Mean Square (RMS): %.2f\n', rmsValuesReal(i));
    fprintf('Variance: %.2f\n', varianceValuesReal(i));
    fprintf('Smoothness: %.2f\n', smoothnessValuesReal(i));
    fprintf('Kurtosis: %.2f\n', kurtosisValuesReal(i));
    fprintf('Skewness: %.2f\n', skewnessValuesReal(i));
    fprintf('Inverse Difference Moment (IDM): %.2f\n', idmValuesReal(i));
    fprintf('Number of Pixels: %.2f\n', numPixelsReal(i));
    fprintf('Pixel Rate: %.2f\n', pixelRateReal(i));
    fprintf('Hue: %.2f\n', hueValuesReal(i));
    fprintf('Saturation: %.2f\n', saturationValuesReal(i));
    fprintf('Contrast (2): %.2f\n', contrastValues2Real(i));
    fprintf('Ambiance: %.2f\n', ambianceValuesReal(i));
    fprintf('Brightness: %.2f\n', brightnessValuesReal(i));
    fprintf('Eccentricity: %.2f\n', eccentricityReal(i));
    fprintf('Solidity: %.2f\n', solidityReal(i));
    fprintf('SNR: %.2f\n', snrValuesReal(i));
    disp(' '); % Add an empty line for better readability
end

% Display properties for AI Generated Images
disp('Properties for AI Generated Images:');
disp(' ');
for i = 1:numel(AiPaths)
    disp(['Properties for AI Generated Image----------->> ', num2str(i), ':']);
    fprintf('Mean: %.2f\n', meanValuesAI(i));
    fprintf('Standard Deviation: %.2f\n', stdDevValuesAI(i));
    fprintf('Entropy: %.2f\n', entropyValuesAI(i));
    fprintf('Root Mean Square (RMS): %.2f\n', rmsValuesAI(i));
    fprintf('Variance: %.2f\n', varianceValuesAI(i));
    fprintf('Smoothness: %.2f\n', smoothnessValuesAI(i));
    fprintf('Kurtosis: %.2f\n', kurtosisValuesAI(i));
    fprintf('Skewness: %.2f\n', skewnessValuesAI(i));
    fprintf('Inverse Difference Moment (IDM): %.2f\n', idmValuesAI(i));
    fprintf('Number of Pixels: %.2f\n', numPixelsAI(i));
    fprintf('Pixel Rate: %.2f\n', pixelRateAI(i));
    fprintf('Hue: %.2f\n', hueValuesAI(i));
    fprintf('Saturation: %.2f\n', saturationValuesAI(i));
    fprintf('Contrast (2): %.2f\n', contrastValues2AI(i));
    fprintf('Ambiance: %.2f\n', ambianceValuesAI(i));
    fprintf('Brightness: %.2f\n', brightnessValuesAI(i));
    fprintf('Eccentricity: %.2f\n', eccentricityAI(i));
    fprintf('Solidity: %.2f\n', solidityAI(i));
    fprintf('SNR: %.2f\n', snrValuesAI(i));
    disp(' '); % Add an empty line for better readability
end


% Combine all calculated values in a table
% Define column names
columnNames = {'ImageType', 'Mean', 'StandardDeviation', 'Entropy', 'RMS', ...
    'Variance', 'Smoothness', 'Kurtosis', 'Skewness', 'IDM', ...
    'Eccentricity', 'Solidity', 'SNR', 'Hue', 'Saturation', 'Contrast2', ...
    'Ambiance', 'Brightness', 'NumPixels', 'PixelRate'};

% Create tables for real and AI data
realImageType = repmat({'Real Image'}, numel(realPaths), 1);
aiImageType = repmat({'Ai Generated Image'}, numel(AiPaths), 1);

realTable = array2table([meanValuesReal', stdDevValuesReal', entropyValuesReal', ...
    rmsValuesReal', varianceValuesReal', smoothnessValuesReal', kurtosisValuesReal', ...
    skewnessValuesReal', idmValuesReal', eccentricityReal', solidityReal', ...
    snrValuesReal', hueValuesReal', saturationValuesReal', contrastValues2Real', ...
    ambianceValuesReal', brightnessValuesReal', numPixelsReal', pixelRateReal'], ...
    'VariableNames', columnNames(2:end));

aiTable = array2table([meanValuesAI', stdDevValuesAI', entropyValuesAI', ...
    rmsValuesAI', varianceValuesAI', smoothnessValuesAI', kurtosisValuesAI', ...
    skewnessValuesAI', idmValuesAI', eccentricityAI', solidityAI', ...
    snrValuesAI', hueValuesAI', saturationValuesAI', contrastValues2AI', ...
    ambianceValuesAI', brightnessValuesAI', numPixelsAI', pixelRateAI'], ...
    'VariableNames', columnNames(2:end'));

% Add the ImageType column to the tables
realTable.ImageType = realImageType;
aiTable.ImageType = aiImageType;

% Reorder the columns to have 'ImageType' as the first column
realTable = realTable(:, ['ImageType', realTable.Properties.VariableNames(1:end-1)]);
aiTable = aiTable(:, ['ImageType', aiTable.Properties.VariableNames(1:end-1)]);

% Concatenate tables
dataTable = [realTable; aiTable];


% Assuming you have already calculated the properties and stored them in variables like meanValuesReal, stdDevValuesReal, etc.

% Define the data structure
data = struct();

% Add properties for Real Images
data.RealImages = struct(...
    'Mean', meanValuesReal, ...
    'StandardDeviation', stdDevValuesReal, ...
    'Entropy', entropyValuesReal, ...
    'RootMeanSquare', rmsValuesReal, ...
    'Variance', varianceValuesReal, ...
    'Smoothness', smoothnessValuesReal, ...
    'Kurtosis', kurtosisValuesReal, ...
    'Skewness', skewnessValuesReal, ...
    'InverseDifferenceMoment', idmValuesReal, ...
    'NumberofPixels', numPixelsReal, ...
    'PixelRate', pixelRateReal, ...
    'Hue', hueValuesReal, ...
    'Saturation', saturationValuesReal, ...
    'Contrast2', contrastValues2Real, ...
    'Ambiance', ambianceValuesReal, ...
    'Brightness', brightnessValuesReal, ...
    'Eccentricity', eccentricityReal, ...
    'Solidity', solidityReal, ...
    'SNR', snrValuesReal ...
);

% Add properties for AI Generated Images
data.AiGeneratedImages = struct(...
    'Mean', meanValuesAI, ...
    'StandardDeviation', stdDevValuesAI, ...
    'Entropy', entropyValuesAI, ...
    'RootMeanSquare', rmsValuesAI, ...
    'Variance', varianceValuesAI, ...
    'Smoothness', smoothnessValuesAI, ...
    'Kurtosis', kurtosisValuesAI, ...
    'Skewness', skewnessValuesAI, ...
    'InverseDifferenceMoment', idmValuesAI, ...
    'NumberofPixels', numPixelsAI, ...
    'PixelRate', pixelRateAI, ...
    'Hue', hueValuesAI, ...
    'Saturation', saturationValuesAI, ...
    'Contrast2', contrastValues2AI, ...
    'Ambiance', ambianceValuesAI, ...
    'Brightness', brightnessValuesAI, ...
    'Eccentricity', eccentricityAI, ...
    'Solidity', solidityAI, ...
    'SNR', snrValuesAI ...
);


% Add Image Type column
realTable.ImageType = repmat({'Real Image'}, height(realTable), 1);
aiTable.ImageType = repmat({'AI Generated Image'}, height(aiTable), 1);

% Combine the tables
combinedTable = vertcat(realTable, aiTable);

% Display the combined table
disp('Combined Table:');
disp(combinedTable);




% Save the data to a .mat file
save('image_properties_data.mat', 'data');

% Save the table as a CSV file
writetable(dataTable, 'image_properties.csv');
