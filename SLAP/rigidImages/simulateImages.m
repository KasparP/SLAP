function [imagesXfast, imagesYfast, motionXfast, motionYfast] = simulateImages(sourceIm, nIMs, res, center)
    
    noise = std(sourceIm(:)) * 0.4; %gaussian noise, fraction of the signal

    amplitude1 = 2; %fast movement
    amplitude2 = 5; %slow drift
    
    [gridX,gridY] = meshgrid(1:res(2), 1:res(1));
    
    imagesXfast = zeros([res nIMs], 'like', sourceIm);
    imagesYfast = zeros([res nIMs], 'like', sourceIm);
    motionXfast = zeros(res(1), 2, nIMs);
    motionYfast = zeros(res(2), 2, nIMs);
    
    for im_ix = 1:nIMs
        Xmotion = (smooth(smooth(amplitude1*randn(1,res(1)), 3),3) + amplitude2*sin(10*rand + (1:res(1))/(40*randn + 300))');
        Ymotion = (smooth(smooth(amplitude1*randn(1,res(1)), 3),3) + amplitude2*sin(10*rand + (1:res(1))/(40*randn + 300))');
        Xsamp = gridX+repmat(Xmotion, 1, res(2));
        Ysamp = gridY+repmat(Ymotion, 1, res(2));
        imagesXfast(:,:,im_ix) = interp2(sourceIm, Xsamp+center(2)-res(2)/2, Ysamp+center(1)-res(1)/2);
        
        imagesXfast(:,:,im_ix) = imagesXfast(:,:,im_ix) - mean(mean(imagesXfast(:,:,im_ix))) + noise*randn(res);
        motionXfast(:, :, im_ix) = [Xmotion Ymotion];
        
        Xmotion = smooth(smooth(amplitude1*randn(1,res(2)), 3),3) + amplitude2*sin(10*rand + (1:res(2))/(20*randn + 500))';
        Ymotion = (smooth(smooth(amplitude1*randn(1,res(2)), 3),3) + amplitude2*sin(10*rand + (1:res(2))/(20*randn + 500))');
        Xsamp = gridX+repmat(Xmotion', res(1),1);
        Ysamp = gridY+repmat(Ymotion', res(1),1);
        imagesYfast(:,:,im_ix) = interp2(sourceIm, Xsamp+center(2)-res(2)/2, Ysamp+center(1)-res(1)/2);        
        
        imagesYfast(:,:,im_ix) = imagesYfast(:,:,im_ix) - mean(mean(imagesYfast(:,:,im_ix))) + noise*randn(res);
        motionYfast(:, :, im_ix) = [Xmotion Ymotion];
%         figure, imshow(imagesXfast(:,:,im_ix), [])
%         figure, imshow(imagesYfast(:,:,im_ix), [])   
    end
    
  
    
