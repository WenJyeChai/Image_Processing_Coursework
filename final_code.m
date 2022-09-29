location = 'CW_images\*.bmp';       %  folder in which your images exists
ds = imageDatastore(location);      %  creates a datastore for all images in your folder
i=1;                                %  initialsed iterator

%% Image processing procedure for all image received
while hasdata(ds) 
    img = read(ds) ;                                                                                                        % read image from datastore
    filter_img = noise_cancelling(img);                                                                                     % perform noise canceling approach on the image
    [identified_nucleus_binary,watershed_ridge_line,region_result,number_of_region] = nucleus_identification(filter_img);   % segmented and detect nuclei based on the image
    display_filtered_result(img,watershed_ridge_line,identified_nucleus_binary,i,number_of_region);                         % diplay all image processing result in individual window
    my_field = strcat('Image_',num2str(i));                                                                                 % save the processed image for analysis
    region_result_group.(my_field) = region_result;                                                                         % save the processed image for analysis
    identified_nucleus_binary_group.(my_field) = identified_nucleus_binary;                                                 % save the processed image for analysis
    i=i+1;
end

%% Size distribution for each image
fns = fieldnames(region_result_group);                                                        % identified saved image
nucleus_size_group=[];                                                                        % initialised variable
Legend = cell(2,1);                                                                           % initialised legend of plot
Legend{1} = "Actual distribution";
Legend{2} = "Normal distribution";
for j=1:i-1
    figure;                                                                                   % initialsed figure
    sgtitle(strcat("Nucleus Size Probability Density Distribution Analysis for Image",num2str(j)));
    region_result = region_result_group.(fns{j});                                             
    nucleus_size_analysis = tabulate(region_result(:));                                       % tabulate the size of nucleus (number of pixel in each detected nucleus)
    nucleus_size_analysis(1,:) = [];                                                          % remove data entry of 0 (representing the backgound of binary image)
    data = nucleus_size_analysis(:,2);
    [min_size,max_size,mean_size,standard_deviation,number_of_region,mode_size] =  data_analysis(data);     % statistical analyze the data 
    
    ha = subplot(2,2,1);                                                                                    
    d = {'Mean Size',mean_size;'Maximum Size',max_size;'Minimum Size',min_size;'Standart Deviation',standard_deviation;'Most Frequent Size',mode_size};
    set_uitable(ha,d);                                                                                      % display the statistical analysis result
    
    subplot(2,2,2);
    colormap_plot(number_of_region,data,region_result);                                      % plot the colormap of the detected nuclei based on their size
    
    subplot(2,2,3);
    histogram_plot_with_normal_distribution(data,'Nucleus Size','Probability Density',Legend)        % plot the histogram of the size distribution of nuclei

    subplot(2,2,4);
    boxplot(data,'Notch','on','Labels','Nucleus');                                           % plot the box plot of the data
    ylabel('The Size of Nucleus');
    
    nucleus_size_group = [nucleus_size_group;data];                                          % append the data for cumulative analysis
end
%% Cumulative size distribution for nucleus in all images 
figure;
sgtitle("Cumulative Analysis for Nucleus Size Distribution of All Image");                   % initialised the figure 
[min_size,max_size,mean_size,standard_deviation,number_of_region,mode_size] =  data_analysis(nucleus_size_group);   % statistical analyze the data

ha = subplot(2,2,1);
d = {'Mean Size',mean_size; 'Maximum Size',max_size; 'Minimum Size',min_size; 'Standart Deviation',standard_deviation; 'Most Frequent Size',mode_size}; 
set_uitable(ha,d);                                                              % display the statistical analysis result of the cumulative data

subplot(2,2,[2 4]);
boxplot(data,'Notch','on','Labels','Nucleus');                                  % plot the box plot of the cumulative data
ylabel('The Size of Nucleus');

subplot(2,2,3);
histogram_plot_with_normal_distribution(nucleus_size_group,'Nucleus Size','Probability Density',Legend);                    % plot the histogram of the size distribution of nuclei in all images

%% Intensity distribution for each image
j=1;
ds = imageDatastore(location);                                                  % retrieve the original image
nucleus_intensity_group = [];
while hasdata(ds) 
    img = read(ds) ; % read image from datastore
    figure;
    sgtitle(strcat("Nucleus Mean Intensity Probability Density Distribution Analysis for Image",num2str(j)));
    region_result = identified_nucleus_binary_group.(fns{j});                   % retrieve the processed image 
    region_image = region_result_group.(fns{j});                                % retrieve the processedd image
    im_grey = (img(:,:,2)) .* cast(region_result,'uint8');                      % filter the grayscale value of each detected nuclei in the image
    CC = bwconncomp(im_grey);                                                   % convert the image to connected components 
    stats = (regionprops(CC,im_grey,'MeanIntensity'));                          % calculate the mean intensity of each detected nuclei in the image
    nucleus_intensity = round(cell2mat(struct2cell(stats)));                    % convert the data into matrix
    [smallest_intensity,largest_intensity,mean_intensity,standard_deviation,number_of_region,mode_intensity] =  data_analysis(nucleus_intensity.'); % statistical analyze the data
    
    ha = subplot(2,2,1);
    d = {'Mean Intensity',mean_intensity;'Highest Intensity',largest_intensity;'Lowest Intensity',smallest_intensity;'Standart Deviation',standard_deviation;'Most Frequent Intensity',mode_intensity};
    set_uitable(ha,d);                                                          % display the statistical analysis result
    
    subplot(2,2,2);
    colormap_plot(number_of_region,nucleus_intensity,region_image);             % plot the colormap of the detected nuclei based on their mean intensity
     
    subplot(2,2,3);
    histogram_plot_with_normal_distribution(nucleus_intensity.','Nucleus Brightness','Probability Density',Legend); % plot the histogram of the mean intensity distribution of nuclei
    
    subplot(2,2,4);
    boxplot(nucleus_intensity,'Notch','on','Labels','Nucleus');
    ylabel('The Intensity of Nucleus');                                         % plot the box plot of the data
   
    j=j+1;
    nucleus_intensity_group = [nucleus_intensity_group nucleus_intensity];      % append the data for cumulative analysis
end
%% Cumulative intensity distribution of nucleus in all images
figure;
sgtitle("Cumulative Analysis for Nucleus Mean Intensity Distribution of All Image");
[smallest_intensity,largest_intensity,mean_intensity,standard_deviation,number_of_region,mode_intensity] =  data_analysis(nucleus_intensity_group.');   % statistical analyze the data

ha = subplot(2,2,1);
d = {'Mean Intensity',mean_intensity;'Highest Intensity',largest_intensity;'Lowest Intensity',smallest_intensity;'Standart Deviation',standard_deviation;'Most Frequent Intensity',mode_intensity};
set_uitable(ha,d);                          % display the statistical analysis result                                     

subplot(2,2,[2,4]);
histogram_plot_with_normal_distribution(nucleus_intensity_group.','Nucleus Brightness','Probability Density',Legend);   % plot the histogram of the mean intensity distribution of nuclei in all images

subplot(2,2,3);
boxplot(nucleus_intensity_group,'Notch','on','Labels','Nucleus');
ylabel('The Intensity of Nucleus');         % plot the box plot of the cumulative data

%% Shape distribution of nucleus in each image
j=1;                                        % initialised iterator
ds = imageDatastore(location);              % retrived original images
Eccentricity_group=[];                      % initalised variable
nucleus_compactness_group =[];

Legend{1} = "Nucleus";
Legend{2} = "Multivariate normal distribution";
while hasdata(ds) 
    img = read(ds) ;                        % read image from datastore
    binary_image = identified_nucleus_binary_group.(fns{j});        % retrived the binary image of the detected nuclei
    CC = bwconncomp(binary_image);                                  % convert the image into connected components
    shape_stats = (regionprops(CC,'Eccentricity','Perimeter','Area','ConvexArea'));     % analyse the shape properties of the detect nuclei
    nucleus_shape=cell2mat(struct2cell(shape_stats)).';             % convert the result to array
    a = (4*pi).*(nucleus_shape(:,1));
    b = (nucleus_shape(:,4).^2);
    nucleus_compactness = (a./b);                                   % calculate the compactness of the nuclei 
    shape_mean = mean(cat(2,nucleus_shape(:,2),nucleus_compactness));       % calculate the mean of the eccentricity and compactness of nuclei
    shape_standart_deviation = cov(cat(2,nucleus_shape(:,2),nucleus_compactness));  % calculate the covariance matrix of the eccentricity and compactness of nuclei
    data_size = size(nucleus_shape);
    indentifier = j*ones(data_size(1),1);                           % set up identidier the differentiate the data of each image
    Eccentricity_group = [Eccentricity_group;[nucleus_shape(:,2) indentifier]]; % append the data for cumulative analysis
    nucleus_compactness_group = [nucleus_compactness_group;[nucleus_compactness indentifier]];
    
    figure;
    sgtitle(strcat("Nucleus Shape Probability Density Distribution Analysis for Image",num2str(j)));
    
    subplot(2,2,1);
    histogram_plot_with_normal_distribution(nucleus_shape(:,2),'Nucleus Eccentricity','Probability Density',Legend); % plot the histogram of the nucleus eccentricity
    
    subplot(2,2,3);
    histogram_plot_with_normal_distribution(nucleus_compactness,'Nucleus Compactness','Probability Density',Legend); % plot the histogram of the nucleus compactness
    
    subplot(2,2,2);
    boxplot([nucleus_shape(:,2),nucleus_compactness],'Notch','on','Labels',{'Nucleus Eccentricity','Nucleus Compactness'}); % plot the boxplot of the nucleus compactness and the nucleus eccentricity
    
    Legend{1} = "Nucleus";
    Legend{2} = "Multivariate normal distribution";
    subplot(2,2,4);
    scatter(nucleus_shape(:,2),nucleus_compactness,'.');            % scatter plot regarding the eccentricity and compactness of nuclei along with the mutivariate distribution 
    hold on;
    G = gmdistribution(shape_mean,shape_standart_deviation);
    F =  @(x,y) arrayfun(@(x0,y0) pdf(G,[x0 y0]),x,y);
    fcontour(F,[0 max(nucleus_shape(:,2)) 0 max(nucleus_compactness)]);
    title(strcat('Image', num2str(j)));
    xlabel('Nucleus Eccentricity')
    ylabel('Nucleus Compactness')
    legend(Legend,'Location','southwest');
    j=j+1;
end

%% Cumulative shape distribution of nucleus in all images
group_shape_mean = mean(cat(2,Eccentricity_group(:,1),nucleus_compactness_group(:,1)));     %calculate the mean of the cumulative data
group_shape_standart_deviation = cov(cat(2,Eccentricity_group(:,1),nucleus_compactness_group(:,1)));    %calculate the covariance matrix of the cumulative data
figure;
sgtitle("Nucleus Shape Probability Density Distribution Analysis of All Image");
hold on;

Legend{1} = "Actual distribution";
Legend{2} = "Normal distribution";

subplot(2,2,1);
histogram_plot_with_normal_distribution(Eccentricity_group(:,1),'Nucleus Eccentricity','Probability Density',Legend);   % plot the histogram of the cumulative nucleus eccentricity 

subplot(2,2,3);
histogram_plot_with_normal_distribution(nucleus_compactness_group(:,1),'Nucleus Compactness','Probability Density',Legend); % plot the histogram of the cumulative nucleus compactness

subplot(2,2,2);
boxplot([Eccentricity_group(:,1),nucleus_compactness_group(:,1)],'Notch','on','Labels',{'Nucleus Eccentricity','Nucleus Compactness'}); % plot the boxplot of the cumulative nucleus compactness and the nucleus eccentricity

subplot(2,2,4);
hold on;
c = round(rand(i,3),2);     % configure the color of for the data from each image
for j=1:i-1
    Legend{j} = strcat("Image",num2str(j));     % plot the scatter plot of data from each image with different color
    ind_1 = Eccentricity_group(:,2) == j;
    ind_2 = nucleus_compactness_group(:,2) == j;
    Eccentricity = Eccentricity_group(ind_1,1);
    nucleus_compactness = nucleus_compactness_group(ind_2,1);
    scatter(Eccentricity,nucleus_compactness,[],c(j,:),'.');
end
Legend{j+1} = "Multivariate normal distribution";
G = gmdistribution(group_shape_mean,group_shape_standart_deviation);    % plot the multivariate distribution of the cumulative data
F =  @(x,y) arrayfun(@(x0,y0) pdf(G,[x0 y0]),x,y);
fcontour(F,[0 max(Eccentricity) 0 max(nucleus_compactness)]);
title(strcat('Cumulative Analysis'));
xlabel('Nucleus Eccentricity')
ylabel('Nucleus Compactness')
legend(Legend,'Location','southwest');

%% user-defined function
function result_image = noise_cancelling(original_image)       % perform noise canceling approach to the 'original image'
    
    im_hsv = rgb2hsv(original_image);                          % convert the image into hsv color space
    noise_filter = medfilt2(cast(original_image(:,:,2),'logical'));     % create a binary mask that highlight the main part the image (nuclei)
    first_filter_result = adapthisteq(im_hsv(:,:,1)) .* noise_filter;   % mask the filter over by performing elementwise multiplication
    level = graythresh(first_filter_result);                            % detect the thresholding value for the image to differentiate fore ground and background
    first_filter_result = first_filter_result.*(first_filter_result>level); % suppressing all the value smaller than the thresholding value to be background (reset to 0)
    
    im_smooth = imgaussfilt(first_filter_result, 0.4);        % apply gaussian filter with the sigma value of 0.4
    se = strel('disk', 3);                                    % create a structuring object
    im_open = imopen(im_smooth, se);                          % morphologically opening the image to eliminate noise
    result_image = imreconstruct(im_open, im_smooth);         % reconstruct the opened image with the original image to prevent overfilter 
    
end

function [identified_nucleus_binary,watershed_ridge_line,region_grow_result,number_of_region] = nucleus_identification(original_image)
    se = strel('disk', 3);                                   % create a structuring object
    im_regionmax = imregionalmax(original_image);            % detect the region maxima on the region in the image to detect nuclei
    D = bwdist(im_regionmax);
    DL = watershed(D);                                       % detect the watershed rigdeline between each region maxima
    watershed_ridge_line = DL == 0;
    identified_nucleus_binary = cast(imdilate(im_regionmax,se),'logical').* imcomplement(watershed_ridge_line);     % dilate the region maxima to detect the full shape of the nuclei
    [region_grow_result, number_of_region] = bwlabel(identified_nucleus_binary); % label the binary image
end

function display_filtered_result(original_image,watershed_ridge_line,label_image,i,number_of_nuclei)
    figure;
    main_title = "Image "+num2str(i)+" with "+num2str(number_of_nuclei)+" detected nuclei ";       % title of the figure with the number of detected nuclei
    sgtitle(main_title);
    subplot(2,2,1);
    imshow(original_image);                   % show the original image                  
    title("original image");
    
    subplot(2,2,2);
    imshow(watershed_ridge_line);             % show the watershedridge line between each detected nuclei
    title("watershed ridge line");
    
    subplot(2,2,3);
    imshow(label_image);                      % show the binary image of the detected nuclei
    title("Identified nucleus");
    
    subplot(2,2,4);
    [B,L] = bwboundaries(imbinarize(label_image));      % highlight the detected nuclei over the original image
    imshow(original_image);
    hold on;
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1)
    end
    title("Labeled nucleus over original image");
    
end

function histogram_plot_with_normal_distribution(data,X_label,Y_label,Legend)

    histogram(data,30,'Normalization','pdf')    % plot the histogram of the input data
    xlabel(X_label)
    ylabel(Y_label)
    hold on;
    t=xlim;
    y = linspace(0,t(2));
    GMModel = fitgmdist(data,1);
    pdf_normal = pdf(GMModel,y.');
    plot(y,pdf_normal,'LineWidth',2);       % plot the normal distribution of the data 
    legend(Legend);

end

function [min_value,max_value,mean_value,sigma_value,number_of_data,mode_value] =  data_analysis(data)

    min_value = min(data);      % analyse minimum value of the data
    max_value = max(data);      % analyse maximum value of the data
    mean_value = mean(data);    % analyse the mean value of the data
    sigma_value = std(data);    % analyse the standart deviation of the data
    data_size = size(data);     % analyse the size of the data
    mode_value = mode(data);    % analyse the mode fo the data
    number_of_data = data_size(1);  % analyse the number of row of data
    
end

function colormap_plot(number_of_region,data,label_image)

    for i=1:number_of_region-1
        label_image(label_image==i)=data(i);
    end
    RGB2 = label2rgb(label_image,'parula','k','noshuffle');     % label the detected nuclei based on the data with different color
    imshow(RGB2);
    colorbar;

end

function set_uitable(ha,d)

    pos = get(ha,'Position');
    ht = uitable('Units','normalized','Position',pos);         % create the table to display the statistic result of the data
    ht.Data = d;
    set(ha, 'Visible', 'Off');
    
end
