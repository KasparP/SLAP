imageSize = [512,1024];

% images as they come in from the microscope (signed int 16)
benchNumImages = 2:20;

%gpu = gpuDevice;
processingtimes = [];
for numImages = benchNumImages;
    %gpu.reset();
    
    fprintf('Benchmarking %d images...',numImages);
    
    srcimages = arrayfun(@(i)randi(2^16-1,imageSize,'int16'),1:numImages,'UniformOutput',false);
    
    start = tic();
    images = registerImagesRecursively(srcimages);
    stop = toc(start);
    processingtimes(end+1) = stop;
    fprintf(' %.2fs\n',stop);
end

plot(benchNumImages',processingtimes','*-');
xlabel('Number of source images');
ylabel('Processing time');