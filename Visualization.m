function imagePropertiesComparisonGUI()

    fig = uifigure('Name', 'Image Properties Comparison');
    
    realImageButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Upload Real Image', ...
        'Position', [50 450 150 30], 'Callback', @uploadRealImage);
    aiGeneratedImageButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Upload AI Generated Image', ...
        'Position', [250 450 200 30], 'Callback', @uploadAIImage);
    
    realImageAxes = uiaxes(fig, 'Position', [50 150 200 200]);
    aiGeneratedImageAxes = uiaxes(fig, 'Position', [250 150 200 200]);
    
    realImagePropertiesLabel = uicontrol(fig, 'Style', 'text', 'String', 'Real Image Properties', ...
        'Position', [50 380 200 30], 'HorizontalAlignment', 'center');
    aiGeneratedImagePropertiesLabel = uicontrol(fig, 'Style', 'text', 'String', 'AI Generated Image Properties', ...
        'Position', [250 380 200 30], 'HorizontalAlignment', 'center');
    
    realImagePropertiesTable = uitable(fig, 'Position', [50 50 200 100], ...
        'ColumnName', {'Property', 'Value'}, 'ColumnWidth', {100, 100});
    
    aiGeneratedImagePropertiesTable = uitable(fig, 'Position', [250 50 200 100], ...
        'ColumnName', {'Property', 'Value'}, 'ColumnWidth', {100, 100});
    
    realImageProperties = [];
    aiGeneratedImageProperties = [];
    
    function uploadRealImage(~, ~)
        [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', 'Image Files (*.jpg, *.png, *.bmp)'}, 'Select Real Image');
        if filename
            realImagePath = fullfile(pathname, filename);
            realImage = imread(realImagePath);
            imshow(realImage, 'Parent', realImageAxes);
            realImageProperties = calculateProperties(realImage);
            displayProperties(realImageProperties, realImagePropertiesTable, 'Real Image');
        end
    end

    function uploadAIImage(~, ~)
        [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', 'Image Files (*.jpg, *.png, *.bmp)'}, 'Select AI Generated Image');
        if filename
            aiGeneratedImagePath = fullfile(pathname, filename);
            aiGeneratedImage = imread(aiGeneratedImagePath);
            imshow(aiGeneratedImage, 'Parent', aiGeneratedImageAxes);
            aiGeneratedImageProperties = calculateProperties(aiGeneratedImage);
            displayProperties(aiGeneratedImageProperties, aiGeneratedImagePropertiesTable, 'AI Generated Image');
        end
    end

    function properties = calculateProperties(image)
        grayImage = rgb2gray(image);

        % Calculate Mean for Real Images
        meanValue = mean(double(grayImage(:)));

        % Calculate Standard Deviation for Real Images
        stdDevValue = std(double(grayImage(:)));

        % Calculate Entropy for Real Images
        entropyValue = entropy(grayImage);

        % Calculate Root Mean Square (RMS) for Real Images
        rmsValue = sqrt(mean(double(grayImage(:)).^2));

        % Calculate Variance for Real Images
        varianceValue = var(double(grayImage(:)));

        % Calculate Smoothness for Real Images
        smoothnessValue = 1 - (1 / (1 + varianceValue));

        % Calculate Kurtosis for Real Images
        kurtosisValue = kurtosis(double(grayImage(:)));

        % Calculate Skewness for Real Images
        skewnessValue = skewness(double(grayImage(:)));

        % Calculate Eccentricity and Solidity for Real Images
        binaryImage = imbinarize(grayImage);
        stats = regionprops(binaryImage, 'Eccentricity', 'Solidity');
        eccentricityValue = stats.Eccentricity;
        solidityValue = stats.Solidity;

        % Calculate Signal-to-Noise Ratio (SNR) for Real Images
        signal = mean(double(grayImage(:)));
        noise = std(double(grayImage(:)));
        snrValue = 20 * log10(signal / noise);

        % Calculate Hue, Saturation, Contrast, Ambiance, Brightness for Real Images
        hsvImage = rgb2hsv(image);
        hueValue = mean2(hsvImage(:,:,1));
        saturationValue = mean2(hsvImage(:,:,2));
        value = hsvImage(:,:,3);
        contrastValue2 = std2(value);
        ambianceValue = mean2(value);
        brightnessValue = mean2(value);

        % Calculate Number of Pixels and Pixel Rate for Real Images
        numPixels = numel(grayImage);
        pixelRate = snrValue / numPixels;

        properties = [
            meanValue, ...
            stdDevValue, ...
            entropyValue, ...
            rmsValue, ...
            varianceValue, ...
            smoothnessValue, ...
            kurtosisValue, ...
            skewnessValue, ...
            eccentricityValue, ...
            solidityValue, ...
            snrValue, ...
            hueValue, ...
            saturationValue, ...
            contrastValue2, ...
            ambianceValue, ...
            brightnessValue, ...
            numPixels, ...
            pixelRate
        ];
    end

  function displayProperties(properties, propertyTable, imageType)
    propertyNames = {
        'Mean';
        'Standard Deviation';
        'Entropy';
        'RMS';
        'Variance';
        'Smoothness';
        'Kurtosis';
        'Skewness';
        'Eccentricity';
        'Solidity';
        'SNR';
        'Hue';
        'Saturation';
        'Contrast2';
        'Ambiance';
        'Brightness';
        'NumPixels';
        'PixelRate'
    };

    propertyData = cell(numel(properties)+1, 2);
    propertyData(2:end, 1) = propertyNames;

    propertyData(2:end, 2) = arrayfun(@num2str, properties, 'UniformOutput', false);

    propertyTable.Data = propertyData;

    uicontrol('Parent', propertyTable.Parent, 'Style', 'text', 'String', [imageType, ' Properties'], ...
        'Position', [propertyTable.Position(1), propertyTable.Position(2) + propertyTable.Position(4) + 10, propertyTable.Position(3), 30], ...
        'HorizontalAlignment', 'center');
end


    function updateGUI()
        realImageAxes.Title.String = 'Real Image';
        aiGeneratedImageAxes.Title.String = 'AI Generated Image';
    end

    updateGUI();

end
