function concensus = makeRigidImage (imagesXfast, imagesYfast)

%Globally register an image acquired with X as the fast axis, and a set with Y
%as the fast axis, to produce an image with reduced spatial deformations and slow intensity fluctuations

motion = 10; %estimated global drift amplitude, in pixels; this determines the strip width

%prior for motion registration; favors [0,0]
prior.width = motion/2;
prior.strength = 0.1;

res = size(imagesXfast); 
nIMs = res(3);
res = res(1:2);

disp('Performing coregistration...')

%COREGISTER FASTX
NstripsX = ceil(res(1)/(2*motion));
xbounds = round(linspace(1, res(1), NstripsX+2));
Xo = nan(NstripsX, 2, nIMs); %X offsets
for stripN = 1:NstripsX
    disp([int2str(stripN) ' of ' int2str(NstripsX)])
    strips = imgaussfilt(imagesXfast(xbounds(stripN):xbounds(stripN+2),:,:).^3, 2); %cube to emphasize signal over low-amplitude periodic noise
    ref{1} = median(strips,3);
    refFFT = fft2(ref{1});
    strips = squeeze(mat2cell(strips, size(strips,1), size(strips, 2), ones(1, size(strips,3))));
    %Align the X strips to create a reference
    %[ref, data] = registerImagesRecursively(strips);
    %refFFT = data.fftimages{1};
    for imN = 1:nIMs %get the rigid offset for each strip
       %embed the strip in a double size image
       embed = zeros(size(ref{1})); 
       corner = round(size(embed)/2 - size(strips{imN})/2);
       embed(corner(1)+1:corner(1)+size(strips{imN},1), corner(2)+1:corner(2)+size(strips{imN},2)) = strips{imN};
       result = dftregistration_wprior(refFFT, fft2(embed), 1, prior);
       
       if any(abs(result)>4)
           result
           keyboard
       end
       
       Xo(stripN,1,imN) = result(4); Xo(stripN,2,imN) = result(3);
    end
    Xo(stripN,1,:) = Xo(stripN,1,:) - round(nanmean(Xo(stripN,1,:)));
    Xo(stripN,2,:) = Xo(stripN,2,:) - round(nanmean(Xo(stripN,2,:)));
end

%COREGISTER FASTY
NstripsY = ceil(res(2)/(2*motion));
ybounds = round(linspace(1, res(2), NstripsY+2));
Yo = nan(NstripsY, 2, nIMs); 
for stripN = 1:NstripsY
    disp([int2str(stripN) ' of ' int2str(NstripsY)])
    strips = imgaussfilt(imagesYfast(:, ybounds(stripN):ybounds(stripN+2),:).^3,2); %cube to emphasize signal over low-amplitude periodic noise
    ref{1} = median(strips,3);
    refFFT = fft2(ref{1});
    strips = squeeze(mat2cell(strips, size(strips,1), size(strips, 2), ones(1, size(strips,3))));
    
    %Align the Y strips to create a reference
%     [ref, data] = registerImagesRecursively(strips);
%     refFFT = data.fftimages{1};
    
    for imN = 1:nIMs %get the rigid offset for each strip
       %embed the strip in a double size image
       embed = zeros(size(ref{1})); 
       corner = round(size(embed)/2 - size(strips{imN})/2);
       embed(corner(1)+1:corner(1)+size(strips{imN},1), corner(2)+1:corner(2)+size(strips{imN},2)) = strips{imN};
       result = dftregistration_wprior(refFFT, fft2(embed),1, prior);
       if any(abs(result)>4)
           result
           keyboard
       end
       Yo(stripN,1,imN) = result(4); Yo(stripN,2,imN) = result(3);
    end
    Yo(stripN,1,:) = Yo(stripN,1,:) - round(mean(Yo(stripN,1,:)));
    Yo(stripN,2,:) = Yo(stripN,2,:) - round(mean(Yo(stripN,2,:)));
end

%build the offsets from registrations
offsetsX = smoothOffsets(Xo, xbounds);
offsetsY = smoothOffsets(Yo, ybounds);

figure, plot(squeeze(offsetsY(:,1,:)))
figure, plot(squeeze(offsetsX(:,1,:)))

%create the concensus images
disp('Making concensus image')
embed = 20;
concensusX = makeConcensusImage(imagesXfast,offsetsX,2, embed);
concensusY = makeConcensusImage(imagesYfast,offsetsY,1, embed);
close all
figure('name', 'X1'), imshow(concensusX, [])
figure('name', 'Y1'), imshow(concensusY, [])

%COREGISTER BY LINE
%align each strip by crosscorrelation over the support region
DY = 2; DX = 2;
disp('Coregistration...')
dX = nan(size(imagesXfast,1), size(imagesXfast,3));
dY = nan(size(dX));
for Xim = 1:size(imagesXfast,3)
    disp(['Processing fastX image: ' int2str(Xim)])
    for strip = 1:size(imagesXfast,1)
        Xc = round(strip + offsetsX(strip,2,Xim)+embed);
        Yc = round(offsetsX(strip,1,Xim)+embed);
        
        [dX(strip,Xim),dY(strip,Xim),E] = alignStrip(imagesXfast(strip,:,Xim), concensusX(Xc-DX:Xc+DX, :), DY, Yc);
    end
end
offsetsX(:,1,:) = offsetsX(:,1,:) + reshape(dY, size(offsetsX) - [0 1 0]);
offsetsX(:,2,:) = offsetsX(:,2,:) + reshape(dX, size(offsetsX) - [0 1 0]);
concensusX = makeConcensusImage(imagesXfast,offsetsX,2, embed);

dY = nan(size(imagesYfast,2), size(imagesYfast,3));
dX = nan(size(dY));
for Yim = 1:size(imagesYfast,3)
    disp(['Processing fastY image: ' int2str(Yim)])
    for strip = 1:size(imagesYfast,2)
        Yc = round(strip + offsetsY(strip,1,Yim)+embed);
        Xc = round(offsetsY(strip,2,Yim)+embed);
        
        [dY(strip,Yim),dX(strip,Yim),E] = alignStrip(imagesYfast(:,strip,Yim)', concensusX(:,Yc-DY:Yc+DY)', DX, Xc);
    end
end
offsetsY(:,1,:) = offsetsY(:,1,:) + reshape(dY, size(offsetsY) - [0 1 0]);
offsetsY(:,2,:) = offsetsY(:,2,:) + reshape(dX, size(offsetsY) - [0 1 0]);
concensusY = makeConcensusImage(imagesYfast,offsetsY,1, embed);

%get global offset for concensus images
[result] = register(concensusX, concensusY, prior);
if any(result(3:4)) %if we need to globally align these image
    offsetsY(:,1,:) = offsetsY(:,1,:) + result(4);
    offsetsY(:,2,:) = offsetsY(:,2,:) + result(3);
    concensusY = makeConcensusImage(imagesYfast,offsetsY,1, embed);
    
    figure('name', 'X2'), imshow(concensusX, [])
    figure('name', 'Y2'), imshow(concensusY, [])
end

%register concensus X and concensus Y on a grid
offsets2D = zeros(NstripsX, NstripsY, 2);
for Yloc = 1:NstripsY
    Yc = ybounds(Yloc):ybounds(Yloc+2);
    for Xloc = 1:NstripsX
        %extract the target square from each image and align
        Xc = xbounds(Xloc):xbounds(Xloc+2);
        result = register(concensusX(Xc,Yc), concensusY(Xc,Yc), prior);
        offsets2D(Xloc,Yloc,1) = result(4); offsets2D(Xloc,Yloc,2) = result(3);
    end
end

figure, imshow(offsets2D(:,:,1),[])
figure, imshow(offsets2D(:,:,2),[])

offsetsY = offsetsY + repmat(smoothOffsets(squeeze(mean(offsets2D,1)), ybounds), 1,1,nIMs);
offsetsX = offsetsX - repmat(smoothOffsets(squeeze(mean(offsets2D,2)), xbounds), 1,1,nIMs);

%remake the concensus images
concensusX = makeConcensusImage(imagesXfast,offsetsX,2, embed);
concensusY = makeConcensusImage(imagesYfast,offsetsY,1, embed);

figure('name', 'X3'), imshow(concensusX, [])
figure('name', 'Y3'), imshow(concensusY, [])

%coordinate descent by strip
disp('Prealign done. Starting coordinate descent')

for i = 1:3 %while not converged
    DX = 3; DY = 3;
    dY = zeros(1,size(concensusX,1));
    dX = zeros(size(dY));
    for strip = embed+1:size(concensusX,1)-embed
        [dY(strip),dX(strip),E] = alignStrip(concensusX(strip,:), concensusY(strip-DX:strip+DX,:), DY, 0);
    end
    for strip_ix = 1:size(offsetsX,1)
        for image_ix = 1:size(offsetsX,3)
            delta = dX(strip_ix + embed + round(offsetsX(strip_ix, 2, image_ix)));
            offsetsX(strip_ix, 1, image_ix) =  offsetsX(strip_ix, 1, image_ix) + delta;
            
            delta = dY(strip_ix + embed + round(offsetsX(strip_ix, 2, image_ix)));
            offsetsX(strip_ix, 2, image_ix) =  offsetsX(strip_ix, 2, image_ix) + delta;
        end
    end
    
%     E1 = squeeze(offsetsX(:,1,:) - motionXfast(:,1,:));
%     E2 = squeeze(offsetsX(:,2,:) - motionXfast(:,2,:));
%     figure('name', ['Xfast post, iter: ' int2str(i)])
%     subplot(2,1,1), plot(median(E2,2))
%     hold on, plot(median(E1,2), 'r')
%     subplot(2,1,2), plot(dY)
%     hold on, plot(dX,'r')
    
    dY = zeros(1,size(concensusY,2));
    dX = zeros(size(dY));
    for strip = embed+1:size(concensusY,2)-embed
        [dY(strip),dX(strip),E] = alignStrip(concensusY(:,strip)', concensusX(:,strip-DX:strip+DX)', DX, 0);
    end
    for strip_ix = 1:size(offsetsY,1)
        for image_ix = 1:size(offsetsY,3)
            delta = dY(strip_ix + embed + round(offsetsY(strip_ix, 1, image_ix)));
            offsetsY(strip_ix, 1, image_ix) =  offsetsY(strip_ix, 1, image_ix) + delta;
            
            delta = dX(strip_ix + embed + round(offsetsY(strip_ix, 1, image_ix)));
            offsetsY(strip_ix, 2, image_ix) =  offsetsY(strip_ix, 2, image_ix) + delta;
        end
    end
    
%     E1 = squeeze(offsetsY(:,1,:) - motionYfast(:,1,:));
%     E2 = squeeze(offsetsY(:,2,:) - motionYfast(:,2,:));
%     figure('name', ['Yfast post, iter: ' int2str(i)]), 
%     subplot(2,1,1), plot(median(E2,2))
%     hold on, plot(median(E1,2), 'r')
%     subplot(2,1,2), plot(dY)
%     hold on, plot(dX,'r')
    
    [concensusX, refsX] = makeConcensusImage(imagesXfast,offsetsX,2, embed);
    [concensusY, refsY] = makeConcensusImage(imagesYfast,offsetsY,1, embed);
end

%make grand concensus image
concensus = nanmedian(cat(3, refsX, refsY),3);
% compare = zeros(size(concensus)); compare(embed+(1:size(sourceCompare,1)), embed+(1:size(sourceCompare,2))) = sourceCompare;

%figure, imshow(concensusX,[]), figure, imshow(concensusY,[]), figure, imshow(concensus,[]), figure, imshow(compare,[]), 
keyboard

end

function [result, moved] = register(REF, MOVE, prior)
    if ~all(size(REF)==size(MOVE))
        error('Wrong image size')
    end
    sz = size(REF);
    imR = zeros(2*sz);
    imM = zeros(2*sz);
    corner = floor(sz/2);
    
    imR(corner(1)+1:corner(1)+sz(1),corner(2)+1:corner(2)+sz(2)) = REF; 
    imM(corner(1)+1:corner(1)+sz(1),corner(2)+1:corner(2)+sz(2)) = MOVE;
    
    imM(isnan(imM)) = 0;
    imR(isnan(imR)) = 0;
    
    [result, Greg] = dftregistration_wprior(fft2(imR), fft2(imM), 1, prior);
    moved = real(ifft2(Greg));
    moved = moved(corner(1)+1:corner(1)+sz(1),corner(2)+1:corner(2)+sz(2));
end

function [refIM, refs] = makeConcensusImage(images, offsets, dim, embed)
    %shift the X and Y positions of each strip in the image (integer shifts?), then nanmedian
    %across them
    %create the concensus images
    res = size(images);
    nIMs = size(images,3);
    
    refs = zeros(res+[2*embed 2*embed 0], 'like', images);
   
    [gridX,gridY] = meshgrid(1:res(2), 1:res(1));
    [sampX,sampY] = meshgrid((1-embed):(res(2)+embed), (1-embed):(res(1)+embed));
    for imN= 1:nIMs
        Xoff = repmat(offsets(:,1,imN), 1, res(dim));
        Yoff = repmat(offsets(:,2,imN), 1, res(dim));
        if dim==1
            Xoff = Xoff'; Yoff = Yoff';
        end
        X = gridX + Xoff;
        Y = gridY + Yoff;
        V = double(images(:,:,imN));
        
        F = scatteredInterpolant(X(:),Y(:),V(:), 'linear', 'none');
        
        refs(:,:,imN) = F(sampX, sampY);
        %figure, imshow(refs(:,:,imN),[])
    end
    
    refIM = nanmedian(refs,3);
end

function offsets = smoothOffsets(vals, bounds)
nIMs = size(vals,3);
offsets = zeros(bounds(end), 2, nIMs);
for imN= 1:nIMs
    offsets(:,1,imN) = interp1(bounds, [vals(1,1,imN) ; vals(:,1,imN) ; vals(end,1,imN)], 1:bounds(end), 'linear', 'extrap');
    offsets(:,2,imN) = interp1(bounds, [vals(1,2,imN) ; vals(:,2,imN) ; vals(end,2,imN)], 1:bounds(end), 'linear', 'extrap');
end
end


function [bestX, bestY, E] = alignStrip (strip, nhood, DY, Ycenter)
DX = floor(size(nhood,1)/2); %we alwayswant an odd DX
xs = linspace(1, 2*DX+1, 4*DX+1); %eps to avoid interp2 putting NaNs at edges
ys = (Ycenter - DY) : 0.5: (Ycenter + DY);
E = nan(length(xs), length(ys));
for x = 1:size(E,1)
    strip_nh = (nhood(floor(xs(x)),:)+nhood(ceil(xs(x)),:))/2;
    for y = 1:size(E,2)
        strip2 = interp1(strip_nh, ys(y)+1:ys(y)+length(strip));
        select = ~isnan(strip) & ~isnan(strip2);
        if sum(select)>10
            E(x,y) = corr(strip(select)', strip2(select)');
        end
    end
end

E = imgaussfilt(E,0.6); % - sqrt(a.^2+b.^2)*nanmax(E(:))/20; %smooth the error and add prior
[minval, minix] = max(E(:));
[bestX,bestY] = ind2sub(size(E),minix);

bestX = (bestX-2*DX-1)/2;
bestY = (bestY-2*DY-1)/2;
end