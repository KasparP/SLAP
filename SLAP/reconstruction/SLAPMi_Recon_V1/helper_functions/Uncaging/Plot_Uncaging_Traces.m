%%
close all;

t1 = nan(10,10);
t2 = nan(10,10);
hold on;
for filenum = 1:length(filenames)


        if filenum == 1 || filenames(filenum).FOV - filenames(filenum-1).FOV ~= 0
                index = 1;
        end
        FOV = filenames(filenum).FOV;
        if filenames(filenum).FOV == 0
            FOV = 10;
        end

    Trace1 = UncagingLoc1Traces(filenum,filenames(filenum).StimOnset-50:filenames(filenum).StimOnset+500);   
    Trace2 = UncagingLoc2Traces(filenum,filenames(filenum).StimOnset-50:filenames(filenum).StimOnset+500);
    normalizedTrace1 = (Trace1-Trace1(50))/(max(Trace1)-Trace1(50));
    normalizedTrace2 = (Trace2-Trace2(50))/(max(Trace2)-Trace2(50));
    
    
    minTrace1 = Trace1(50);
    minTrace2 = Trace2(50);
    [maxTrace1,I1] = max(Trace1);
    [maxTrace2,I2] = max(Trace2);
    
    mm1 = maxTrace1-minTrace1;
    mm2 = maxTrace2-minTrace2;
    
    t1(FOV,index) = find(normalizedTrace1(50:end) > 0.25,1,'first');
    t2(FOV,index) = find(normalizedTrace2(50:end) > 0.25,1,'first');

%         
%     plot(normalizedTrace1,'b')
%     plot(normalizedTrace2,'r')

    if FOV == 8 && index == 3
        PaperTrace1 = normalizedTrace1; PaperTrace1(PaperTrace1<0) = 0;
        DFF1 = (Trace1-Trace1(60))/(Trace1(60)+1); DFF1(DFF1<0) = 0; 
        StimOnset = filenames(filenum).StimOnset;
        
        PaperTrace2 = normalizedTrace2; PaperTrace2(PaperTrace2<0) = 0;
        DFF2 = (Trace2-Trace2(50))/(Trace2(50)+1); DFF2(DFF2<0) = 0;
    end
    index = index + 1;

%   [1 6 7 9]

end
figure;
plot(PaperTrace1(1:100),'b');
hold on; plot(PaperTrace2(1:100),'r');
line([50 50],[0 1]);
line([60 60],[0 1]);

figure;
plot(DFF1,'b');
hold on; plot(DFF2,'r');

%%
absTdiff = nanmean(abs(t1-t2),2);
absTdiff([1 7 9]) = nan;
timePts = absTdiff(~isnan(absTdiff));
mean(timePts)
% figure; boxplot(timePts,'OutlierSize',10)
figure('color','w'); scatter(zeros(size(timePts)),timePts,200,'k','LineWidth',2.5)

hold on; plot([-0.5,0,0.5],[10 10 10],'b','LineWidth',2.5)
plot([-0.5,0,0.5],[1 1 1]*mean(timePts),'r','LineWidth',2.5)
ylim([0 15]); xlim([-4 4]);