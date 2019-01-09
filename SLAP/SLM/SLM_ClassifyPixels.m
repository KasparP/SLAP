function [status,labelspath]=SLM_ClassifyPixels(IMpath) 
    % fast Tiff reading with Tiff and bin image for fast processing
    InfoImage=imfinfo(IMpath);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    bin=4;
    img=zeros(nImage/bin,mImage/bin,NumberImages,'uint8');
 
    TifLink = Tiff(IMpath, 'r');
    for i=1:NumberImages
        TifLink.setDirectory(i);
        img(:,:,i)=imresize(TifLink.read(),1/bin,'bilinear');
    end
    TifLink.close();

    % write binned refIM
    [~,IMname,ext] = fileparts(IMpath);
    basedir='E:\SLAPMidata\SLM_segmentation\';
    writeIMpath=[basedir, IMname,'_binX',num2str(bin), ext];
    options.overwrite=true;
    saveastiff(img,writeIMpath,options);
    
    % pixel classfication
    ilastikPath = 'C:\Program Files\ilastik-1.3.0';
    projectDr = 'C:\Ilastik\';
    [projectFn, projectDr] = uigetfile([projectDr '*.ilp'], 'Select an ilastik project or classifier output');
    projectPath = [projectDr projectFn];
    labelspath = [IMpath(1:end-4) '_Simple Segmentation.h5'];
    command=['"' ilastikPath '\run-ilastik.bat" --headless --export_source="Simple Segmentation" --project="' projectPath '" "' writeIMpath '"'];
    
    t1=tic;
    fprintf(1,'Running Ilastik... ');
    [status,~]  = dos(command);
    disp(['Done! Duration = ' num2str(toc(t1)) 's.'])
end

