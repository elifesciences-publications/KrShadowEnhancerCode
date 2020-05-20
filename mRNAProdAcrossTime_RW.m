%% Calculate mean mRNA production across time at different AP bins
%load constructs
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','Kr2xProxEmpty','KrBoth'}%,'KrBothEmpty','KrDist32C','KrProx32C','KrBothSep32C','KrBoth32C','Kr2xProx32C','HbEmpty','KrDistDuplicN','KrDist17C','Kr2xDist32C','KrBoth17C','Kr2xDist17C','KrProx17C','Kr2xDistdBcd','Kr2xProxdHb','Kr2xProx17C','Kr2xDistEmpty','Kr2xDistEmpty32C','KrSEdBcd','KrInvSE','Kr2xDistLessMS2','KrSEdHb','KrEndogDist','KrEndogDist32C','KrEndogDist17C','KrSEChrom3','Kr2xDistChrom3','Kr2xProxChrom3','Kr2xDistChrom3_Empty','Kr4_1xBcd','KrSEChr3_Empty','Kr3_1xBcd'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};

[SourcePath, FISHPath, DropboxFolder, MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2, OutputFolder...
 ]= readMovieDatabase('2017-08-03-mKr1_E1')  %just any random dataset to give us the dropbox folder location

% Load each construct
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    % Go through each embryo
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        nc14 = Data(ee).nc14;
        Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        load(Filename);
        APstuff=[SpotDiff.APBin];
        
        for aa = 1:length(APbinID)
            APTimeTable=[];
            APsubset = SpotDiff(APstuff == APbinID(aa));
            if ~isempty(APsubset)
                for bb =1:length(APsubset)
                    APTimeTableEntry = [APsubset(bb).SmoothSpotOne(nc14:end)]'; %save as a row
                    if length(APTimeTableEntry) > 100
                        APTimeTableEntry = APTimeTableEntry(1:100);
                    elseif length(APTimeTableEntry) < 100
                        % fill end with nan's if shorter than 100 frames
                        % (50 min)
                        APTimeTableEntry = [APTimeTableEntry, nan(1, (100-(length(APTimeTableEntry))))];
                    end
                     APTimeTable = [APTimeTable; APTimeTableEntry];
                end
            else
                %If no spots in AP bin, make row of 100 nan's
                APTimeTable=[nan(1,100)];
                %continue
                % start making structure with data
            end
                TimeProduction(cc).Embryo(ee).APbin(aa).TimeFrames = [APTimeTable];
            
        end
    end
    for aa = 1:length(APbinID)
        AllTimeFrames =[];
        for ee = 1:NEmbryos
            AllTimeFrames = [AllTimeFrames; TimeProduction(cc).Embryo(ee).APbin(aa).TimeFrames];
        end
        TimeProduction(cc).APbin(aa).AllTimeFrames = [AllTimeFrames];
        TimeProduction(cc).APbin(aa).MeanFluo = nanmean(TimeProduction(cc).APbin(aa).AllTimeFrames);
        TimeProduction(cc).APbin(aa).StdFluo = nanstd(TimeProduction(cc).APbin(aa).AllTimeFrames);
        %number of total spots at that time point 
        NumbNuc = sum(~isnan(TimeProduction(cc).APbin(aa).AllTimeFrames));
        TimeProduction(cc).APbin(aa).SEFluo = ((TimeProduction(cc).APbin(aa).StdFluo)./(sqrt(NumbNuc)))
        TimeProduction(cc).APbin(aa).Fluo95CI = (TimeProduction(cc).APbin(aa).SEFluo).*1.95;
    end
end
            
%% Plot things
EggLength = APbinID.*100; 

DistalColor=[1 64 172]./255;
DistalEmptyColor=[8 210 238] ./ 255;
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
Proximal32CColor=[251 150 10] ./ 255;
DoubProxColor=[215 183 58] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
BothEmptyColor=[12 250 100] ./ 255;
Both32CColor=[120 195 82] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;
PDistalColor = [129 161 214]./255;


Colors(1).Color=DistalColor; 
Colors(2).Color=ProxColor; 
Colors(3).Color=BothSepColor; 
Colors(4).Color=DistalColor;
Colors(5).Color=ProxColor;
Colors(6).Color=DoubDistColor; 
Colors(7).Color=DoubProxColor;
Colors(8).Color=DoubProxColor;
Colors(9).Color=BothColor;
Colors(10).Color=BothColor;
Colors(11).Color=DistalColor;
Colors(12).Color=ProxColor;
Colors(13).Color=BothSepColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color='k';
Colors(17).Color=DoubDistColor;
Colors(18).Color=DistalColor;
Colors(19).Color=DoubDistColor;
Colors(20).Color=BothColor;
Colors(21).Color=DoubDistColor;
Colors(22).Color=ProxColor;
Colors(23).Color=DoubDistColor;
Colors(24).Color=DoubProxColor;
Colors(25).Color=DoubProxColor;
Colors(26).Color=DoubDistColor;
Colors(27).Color=DoubDistColor;
Colors(28).Color=BothColor;
Colors(29).Color=BothColor;
Colors(30).Color=DoubDistColor;
Colors(31).Color=BothColor;
Colors(32).Color=DistalColor;
Colors(33).Color=DistalColor;
Colors(34).Color=DistalColor;
Colors(35).Color=BothColor;
Colors(36).Color = DoubDistColor;
Colors(37).Color = DoubProxColor;
%Colors(38).Color = DoubDistColor;
Colors(38).Color = DoubDistColor;
Colors(39).Color = DoubDistColor;
Colors(40).Color = BothColor;
Colors(41).Color=BothColor;

%font/figure sizing info
fontsize=10;
fontname='Arial';
x_width=3; y_width=2.25;
x_widthsplit=1.5; y_widthsplit=1.125;
xSize = 7; ySize = 6; xLeft = 0.5; yTop = 0.5;

FigDirect=[DropboxFolder filesep 'Figures'];

%% Compare singles
figure
errorbar([0.5:0.5:50],TimeProduction(1).APbin(21).MeanFluo, TimeProduction(1).APbin(21).Fluo95CI,'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:50],TimeProduction(2).APbin(21).MeanFluo, TimeProduction(2).APbin(21).Fluo95CI,'Color',Colors(2).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
xlabel('Time into nc14 (min)');
title([num2str(EggLength(21)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SinglesmRNAvTime_50EL','.pdf'],'pdf');

% 35% egg length
figure
errorbar([0.5:0.5:50],TimeProduction(1).APbin(15).MeanFluo, TimeProduction(1).APbin(15).Fluo95CI,'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:50],TimeProduction(2).APbin(15).MeanFluo, TimeProduction(2).APbin(15).Fluo95CI,'Color',Colors(2).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
xlabel('Time into nc14 (min)');
title([num2str(EggLength(15)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SinglesmRNAvTime_35EL','.pdf'],'pdf');

% 60% egg length
figure
errorbar([0.5:0.5:50],TimeProduction(1).APbin(25).MeanFluo, TimeProduction(1).APbin(25).Fluo95CI,'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:50],TimeProduction(2).APbin(25).MeanFluo, TimeProduction(2).APbin(25).Fluo95CI,'Color',Colors(2).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
xlabel('Time into nc14 (min)');
title([num2str(EggLength(25)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SinglesmRNAvTime_60EL','.pdf'],'pdf');
%% Duplicated vs SE pair 
figure
errorbar([0.5:0.5:50],TimeProduction(6).APbin(21).MeanFluo, TimeProduction(6).APbin(21).Fluo95CI,'Color',Colors(6).Color,'LineWidth',4.5);
hold on 
errorbar([0.5:0.5:50],TimeProduction(7).APbin(21).MeanFluo, TimeProduction(7).APbin(21).Fluo95CI,'Color',Colors(7).Color,'LineWidth',4.5);
errorbar([0.5:0.5:50],TimeProduction(9).APbin(21).MeanFluo, TimeProduction(9).APbin(21).Fluo95CI,'Color',Colors(9).Color,'LineWidth',4.5);
ylabel('Mean fluorescence (AU)');
ylim([0 30000]);
xlim([0 40]);
xlabel('Time into nc14 (min)');
title([num2str(EggLength(21)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DoubVSE_mRNAvTime_50EL','.pdf'],'pdf');

% 35% egg length
figure
errorbar([0.5:0.5:50],TimeProduction(6).APbin(15).MeanFluo, TimeProduction(6).APbin(15).Fluo95CI,'Color',Colors(6).Color,'LineWidth',4.5);
hold on 
errorbar([0.5:0.5:50],TimeProduction(7).APbin(15).MeanFluo, TimeProduction(7).APbin(15).Fluo95CI,'Color',Colors(7).Color,'LineWidth',4.5);
errorbar([0.5:0.5:50],TimeProduction(9).APbin(15).MeanFluo, TimeProduction(9).APbin(15).Fluo95CI,'Color',Colors(9).Color,'LineWidth',4.5);
ylabel('Mean fluorescence (AU)');
ylim([0 30000]);
xlim([0 40]);
xlabel('Time into nc14 (min)');
title([num2str(EggLength(15)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DoubvsSE_mRNAvTime_35EL','.pdf'],'pdf');

% 60% egg length
figure
errorbar([0.5:0.5:50],TimeProduction(6).APbin(25).MeanFluo, TimeProduction(6).APbin(25).Fluo95CI,'Color',Colors(6).Color,'LineWidth',4.5);
hold on 
errorbar([0.5:0.5:50],TimeProduction(7).APbin(25).MeanFluo, TimeProduction(7).APbin(25).Fluo95CI,'Color',Colors(7).Color,'LineWidth',4.5);
errorbar([0.5:0.5:50],TimeProduction(9).APbin(25).MeanFluo, TimeProduction(9).APbin(25).Fluo95CI,'Color',Colors(9).Color,'LineWidth',4.5);
ylabel('Mean fluorescence (AU)');
ylim([0 30000]);
xlim([0 40]);
xlabel('Time into nc14 (min)');
title([num2str(EggLength(25)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DoubvSE_mRNAvTime_60EL','.pdf'],'pdf');

%% Limit to 40min into nc14 
figure
errorbar([0.5:0.5:40],TimeProduction(6).APbin(21).MeanFluo(1:80), TimeProduction(6).APbin(21).Fluo95CI(1:80),'Color',Colors(6).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:40],TimeProduction(7).APbin(21).MeanFluo(1:80), TimeProduction(7).APbin(21).Fluo95CI(1:80),'Color',Colors(7).Color,'LineWidth',2.5);
errorbar([0.5:0.5:40],TimeProduction(9).APbin(21).MeanFluo(1:80), TimeProduction(9).APbin(21).Fluo95CI(1:80),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
xlabel('Time into nc14 (min)');
title([num2str(EggLength(21)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DoubVSE_mRNAvTime_50EL_40min','.pdf'],'pdf');

% 35% egg length
figure
errorbar([0.5:0.5:40],TimeProduction(6).APbin(15).MeanFluo(1:80), TimeProduction(6).APbin(15).Fluo95CI(1:80),'Color',Colors(6).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:40],TimeProduction(7).APbin(15).MeanFluo(1:80), TimeProduction(7).APbin(15).Fluo95CI(1:80),'Color',Colors(7).Color,'LineWidth',2.5);
errorbar([0.5:0.5:40],TimeProduction(9).APbin(15).MeanFluo(1:80), TimeProduction(9).APbin(15).Fluo95CI(1:80),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
xlabel('Time into nc14 (min)');
title([num2str(EggLength(15)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DoubvsSE_mRNAvTime_35EL_40min','.pdf'],'pdf');

% 60% egg length
figure
errorbar([0.5:0.5:40],TimeProduction(6).APbin(25).MeanFluo(1:80), TimeProduction(6).APbin(25).Fluo95CI(1:80),'Color',Colors(6).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:40],TimeProduction(7).APbin(25).MeanFluo(1:80), TimeProduction(7).APbin(25).Fluo95CI(1:80),'Color',Colors(7).Color,'LineWidth',2.5);
errorbar([0.5:0.5:40],TimeProduction(9).APbin(25).MeanFluo(1:80), TimeProduction(9).APbin(25).Fluo95CI(1:80),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
xlabel('Time into nc14 (min)');
title([num2str(EggLength(25)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DoubvSE_mRNAvTime_60EL_40min','.pdf'],'pdf');

%% Limit 40min into nc14, Singles vs SE 
figure
errorbar([0.5:0.5:40],TimeProduction(1).APbin(21).MeanFluo(1:80), TimeProduction(1).APbin(21).Fluo95CI(1:80),'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:40],TimeProduction(2).APbin(21).MeanFluo(1:80), TimeProduction(2).APbin(21).Fluo95CI(1:80),'Color',Colors(2).Color,'LineWidth',2.5);
errorbar([0.5:0.5:40],TimeProduction(9).APbin(21).MeanFluo(1:80), TimeProduction(9).APbin(21).Fluo95CI(1:80),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
ylim([0 30000]);
xlabel('Time into nc14 (min)');
title([num2str(EggLength(21)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SinglesVSE_mRNAvTime_50EL_40min','.pdf'],'pdf');

% 35% egg length
figure
errorbar([0.5:0.5:40],TimeProduction(1).APbin(15).MeanFluo(1:80), TimeProduction(1).APbin(15).Fluo95CI(1:80),'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:40],TimeProduction(2).APbin(15).MeanFluo(1:80), TimeProduction(2).APbin(15).Fluo95CI(1:80),'Color',Colors(2).Color,'LineWidth',2.5);
errorbar([0.5:0.5:40],TimeProduction(9).APbin(15).MeanFluo(1:80), TimeProduction(9).APbin(15).Fluo95CI(1:80),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
ylim([0 30000]);
xlabel('Time into nc14 (min)');
title([num2str(EggLength(15)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SinglevsSE_mRNAvTime_35EL_40min','.pdf'],'pdf');

% 60% egg length
figure
errorbar([0.5:0.5:40],TimeProduction(1).APbin(25).MeanFluo(1:80), TimeProduction(1).APbin(25).Fluo95CI(1:80),'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:40],TimeProduction(2).APbin(25).MeanFluo(1:80), TimeProduction(2).APbin(25).Fluo95CI(1:80),'Color',Colors(2).Color,'LineWidth',2.5);
errorbar([0.5:0.5:40],TimeProduction(9).APbin(25).MeanFluo(1:80), TimeProduction(9).APbin(25).Fluo95CI(1:80),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
ylim([0 30000]);
xlabel('Time into nc14 (min)');
title([num2str(EggLength(25)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SinglevSE_mRNAvTime_60EL_40min','.pdf'],'pdf');

% 40% egg length
figure
errorbar([0.5:0.5:40],TimeProduction(1).APbin(17).MeanFluo(1:80), TimeProduction(1).APbin(17).Fluo95CI(1:80),'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar([0.5:0.5:40],TimeProduction(2).APbin(17).MeanFluo(1:80), TimeProduction(2).APbin(17).Fluo95CI(1:80),'Color',Colors(2).Color,'LineWidth',2.5);
errorbar([0.5:0.5:40],TimeProduction(9).APbin(17).MeanFluo(1:80), TimeProduction(9).APbin(17).Fluo95CI(1:80),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('Mean fluorescence (AU)');
ylim([0 30000]);
xlabel('Time into nc14 (min)');
title([num2str(EggLength(17)),'%',' ', 'egg length']);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SinglevSE_mRNAvTime_40EL_40min','.pdf'],'pdf');