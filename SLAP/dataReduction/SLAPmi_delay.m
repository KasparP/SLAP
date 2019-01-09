function SLAPMi_delay(P, scandata, refIM)

Pcoords = P.coords;
P = P.P;

plane = ceil(length(Pcoords{3})/2);
inds = false(length(Pcoords{1}), length(Pcoords{2}), length(Pcoords{3}));
inds(:,:,plane) = true;

Pslice = P(:, inds);

%figure, imshow(reshape(full(sum(Pslice,1)), length(Pcoords{1}), length(Pcoords{2})),[])

Rplane = ceil(size(refIM.data,4)/2);
Rslice = refIM.data(:,:,1,Rplane);
Rslice(Rslice<3*prctile(Rslice(:), 20))=0;

yR = Pslice * double(Rslice(:)); yR = yR./max(yR);
y = nanmean([scandata.frames.pmtData],2); y = y./max(y);
figure, plot(y), hold on, plot(yR)

for line = 1:4
    select = scandata.line==line & ~isnan(y);
    [xc, lags] = xcorr(y(select), yR(select));
    %figure('name', ['line' int2str(line)]), plot(xc)
    [~, maxix] = max(xc);
    bestlag(line) = lags(maxix);
end

y(isnan(y)) = 0;

PD = reshape(y'*Pslice, [sqrt(size(Pslice,2)) sqrt(size(Pslice,2))]);
%figure, imagesc(PD)

PDerr = reshape((y-yR)'*Pslice, [sqrt(size(Pslice,2)) sqrt(size(Pslice,2))]);
figure, imagesc(PDerr)

bestlag

C_lag = [0 0 0 0]; %0 7 3 0
y = min(y, max(y(:))/2); y = y-median(y);
%optimize for maximum sum(log(P)) over shifts
Pslice = Pslice./max(Pslice(:));
lags = -2:2; nlags = length(lags);
Q = nan(nlags, nlags,nlags,nlags); Q2 = Q;
for lag1 = 1:nlags;
    Pline1  = y(circshift(scandata.line==1, lags(lag1) + C_lag(1),1))'*Pslice(scandata.line==1,:);
    for lag2 = 1:nlags
        Pline2  = y(circshift(scandata.line==2, lags(lag2)+ C_lag(2),1))'*Pslice(scandata.line==2,:);
        cumprod2 = Pline1.*Pline2;
        for lag3 = 1:nlags
            disp([lag1 lag2 lag3])
            Pline3  = y(circshift(scandata.line==3, lags(lag3) + C_lag(3),1))'*Pslice(scandata.line==3,:);
            cumprod3 = cumprod2.*Pline3;
            for lag4 = 1:nlags
                Pline4  = y(circshift(scandata.line==4, lags(lag4) + C_lag(4),1))'*Pslice(scandata.line==4,:);
                Q(lag1,lag2,lag3,lag4) = sum(cumprod3.*Pline4);
                Q2(lag1,lag2,lag3,lag4) = norm(cumprod3.*Pline4);
            end
        end
    end
end

%visualize
f = find(Q==max(Q(:)));
[k(1),k(2),k(3),k(4)] = ind2sub(size(Q), f(1));

%visualize slices through the best lagset
figure('name', ['lag1=' int2str(lags(k(1)))]), imshow3D(squeeze(Q(k(1),:,:,:)),[])
figure('name', ['lag2=' int2str(lags(k(2)))]), imshow3D(squeeze(Q(:,k(2),:,:)),[])
figure('name', ['lag3=' int2str(lags(k(3)))]), imshow3D(squeeze(Q(:,:,k(3),:)),[])
figure('name', ['lag4=' int2str(lags(k(4)))]), imshow3D(squeeze(Q(:,:,:,k(4))),[])

keyboard

%show the backprojection of the best lag
Pline1  = y(circshift(scandata.line==1, lags(k(1))+ C_lag(1),1))'*Pslice(scandata.line==1,:);
Pline2  = y(circshift(scandata.line==2, lags(k(2))+ C_lag(2),1))'*Pslice(scandata.line==2,:);
Pline3  = y(circshift(scandata.line==3, lags(k(3))+ C_lag(3),1))'*Pslice(scandata.line==3,:);
Pline4  = y(circshift(scandata.line==4, lags(k(4))+ C_lag(4),1))'*Pslice(scandata.line==4,:);

figure, imagesc(reshape(Pline1.*Pline2.*Pline3.*Pline4, [400 400]))
figure, imagesc(reshape(Pline1+Pline2+Pline3+Pline4, [400 400]))
 
%Zs = refIM.metadata.hStackManager.zs;

end