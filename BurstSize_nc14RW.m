%% Record the burst size (# transcripts per burst) for each construct

ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'HbEmpty';'Kr2xProxEmpty';'KrDist17C';'Kr2xDist32C';'KrBoth17C';'Kr2xDist17C';'KrProx17C';'Kr2xDistdBcd';'Kr2xProxdHb';'Kr2xProx17C';'Kr2xDistEmpty';'Kr2xDistEmpty32C';'KrSEdBcd';'KrInvSE';'KrSEdHb';'KrEndogDist';'KrEndogDist32C';'KrEndogDist17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
SlopeUse=input('Want to use slope calculations?','s');
Halves=input('Do halves calculations?','s');
%if Halves ~='y'
%Count for each construct
AvgAmplitude=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgAmplitudeCon=[];
    firsttime=1;
    ConAmpAllAP=[];
    ConAmpSE=[];
    ConAmpSD=[];
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        if SlopeUse=='y'
            filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        elseif ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        %seperate out by AP bin
        AmpAllAP=[];
        for aa=1:length(APbinID)
            AmpAP=[];
            AmplitudeAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(AmplitudeAP)
                AmpAP=[AmpAP; nan];
            else
            for bb=1:length(AmplitudeAP)
                if ~isempty(BurstProperties(AmplitudeAP(bb)).BurstSize)
                AmpAP=[AmpAP;[BurstProperties(AmplitudeAP(bb)).BurstSize]'];  
                else
                    AmpAP=[AmpAP;nan];
                end
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(AmpAP)
                AmpAllAP(bb,aa,ee)=AmpAP(bb);
            end
            AmpAllAP(AmpAllAP==0)=nan;
            
            AmplitudeSD(ee,aa,cc)=nanstd(AmpAllAP(:,aa,ee));
        AmplitudeSE(ee,aa,cc)=AmplitudeSD(ee,aa,cc)/sqrt(sum(~isnan(AmpAP)));
            clear AmplitudeAP
        end
        %Compile all data for a construct in one long column
        %to put in structure
        for bb=1:size(AmpAllAP,3)
            ConAmpAllAP=[ConAmpAllAP; AmpAllAP(:,:,bb)];
        end
        for bb=1:size(AmplitudeSD,3)
            ConAmpSD=[ConAmpSD; AmplitudeSD(:,:,bb)];
        end
        for bb=1:size(AmplitudeSE,3)
            ConAmpSE=[ConAmpSE;AmplitudeSE(:,:,bb)];
        end
        AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp=nanmean(AmpAllAP(:,:,ee));
        AvgAmpAllAP(cc).EmbryoAmp(ee).SD=AmplitudeSD(ee,:,cc);
        AvgAmpAllAP(cc).EmbryoAmp(ee).SE=AmplitudeSE(ee,:,cc);
        AvgAmpAllAP(cc).EmbryoAmp(ee).Name=PrefixName;
        AvgAmpAllAP(cc).EmbryoAmp(ee).BurstSizes=[AmpAllAP(:,:,ee)];
    end
        AvgAmpAllAP(cc).AvgAmp=nanmean(ConAmpAllAP);  %Avg burst size of all embryos of a construct
        AvgAmpAllAP(cc).AmpSD=nanmean(ConAmpSD);
        AvgAmpAllAP(cc).AmpSE=nanmean(ConAmpSE);
        AvgAmpAllAP(cc).AllAmps=[ConAmpAllAP];
        AvgAmpAllAP(cc).AllSD=nanstd([AvgAmpAllAP(cc).AllAmps]);
        %calculate SE by AP bin since different number of spots in
        %different AP bins 
        for aa=1:length(APbinID)
        AvgAmpAllAP(cc).AllSE(aa)=(AvgAmpAllAP(cc).AllSD(aa))/(sqrt(sum(~isnan(AvgAmpAllAP(cc).AllAmps(:,aa)))));
        end
        AvgAmpAllAP(cc).All95Conf=(AvgAmpAllAP(cc).AllSE) .*1.95;
        AvgAmpAllAP(cc).ConstructName=ConstructList{cc};

end

% Get rid of points that are only one data point for the whole AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgAmpAllAP(cc).AllAmps(:,aa))) ==1
            AvgAmpAllAP(cc).AllSD(aa)=nan;
            AvgAmpAllAP(cc).AvgAmp(aa)=nan;
            AvgAmpAllAP(cc).AllSE(aa)=nan;
            AvgAmpAllAP(cc).All95Conf(aa)=nan;
        end
    end
end


%% ploting

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
HbColor='k';

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DistalColor;
Colors(5).Color=ProxColor;
Colors(6).Color=DoubDistColor;
Colors(7).Color=DoubProxColor;
Colors(8).Color=BothColor;
Colors(9).Color=BothColor;
Colors(10).Color=DistalColor;
Colors(11).Color=ProxColor;
Colors(12).Color=BothColor;
Colors(13).Color=BothColor;
Colors(14).Color=DoubProxColor;
Colors(15).Color=HbColor;
Colors(16).Color=DoubProxColor;
Colors(17).Color=DistalColor;
Colors(18).Color=DoubDistColor;
Colors(19).Color=BothColor;
Colors(20).Color=DoubDistColor;
Colors(21).Color=ProxColor;
Colors(22).Color=DoubDistColor;
Colors(23).Color=DoubProxColor;
Colors(24).Color=DoubProxColor;
Colors(25).Color=DoubDistColor;
Colors(26).Color=BothColor;
Colors(27).Color=BothColor;
Colors(28).Color=BothColor;
Colors(29).Color=BothColor;
Colors(30).Color=DistalColor;
Colors(31).Color=DistalColor;
Colors(32).Color=DistalColor;

fontsize=10;
fontname='Helvetica';
x_width=3; y_width=2.25;
xSize = 6; ySize = 4.5; xLeft = 0.1; yTop = 0.1;
%Set folder to save figures in 
FigDirect=[DropboxFolder filesep 'Figures'];
% Save structure with data to make figures
save([DropboxFolder filesep 'Constructs' filesep 'BurstSize'],'AvgAmpAllAP');


Egglength=APbinID .*100;
%3/3/2019 updated/corrected Frnap to be 379
% MS2Convers=F1 = Frnap * Time of elongation
MS2Convers=(379*((1.275+4.021)/1.5)); %Frnap * Telongation ((length of MS2 +rest of transcript) / elongation rate)

%% Single enhancers vs SE pair 
figure 
errorbar(Egglength, AvgAmpAllAP(30).AvgAmp, AvgAmpAllAP(30).All95Conf, 'Color', Colors(30).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 350000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(30).AvgAmp./MS2Convers), (AvgAmpAllAP(30).All95Conf./MS2Convers), 'Color', Colors(30).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(2).AvgAmp./MS2Convers), (AvgAmpAllAP(2).All95Conf./MS2Convers), 'Color', Colors(2).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(8).AvgAmp./MS2Convers), (AvgAmpAllAP(8).All95Conf./MS2Convers), 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle','-');
ylabel('transcripts produced','Color','k');
set(gca,'YColor','k');
ylim([0 (350000/MS2Convers)]);
%legend('Distal', 'Proximal','Both');
print( [FigDirect filesep 'SinglesvBothBurstSize_ED'],'-dpdf');

%% Endogenous Distal 
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5, 'LineStyle', ':');
hold on
errorbar(Egglength, AvgAmpAllAP(30).AvgAmp, AvgAmpAllAP(30).All95Conf, 'Color', Colors(30).Color, 'LineWidth', 2.5);
%legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 350000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(1).AvgAmp./MS2Convers), (AvgAmpAllAP(1).All95Conf./MS2Convers), 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle',':');
errorbar(Egglength, (AvgAmpAllAP(30).AvgAmp./MS2Convers), (AvgAmpAllAP(30).All95Conf./MS2Convers), 'Color', Colors(30).Color, 'LineWidth', 2.5);
ylabel('number transcripts produced');
ylim([0 (350000/MS2Convers)]);
set(gca, 'FontName',fontname,'FontSize',fontsize,'YColor','k');
set(gcf,'PaperUnits','inches','PaperPosition',[xLeft yTop xSize ySize]);
print( [FigDirect filesep 'EndogvsOrigDistBurstSize'],'-dpdf');

%% All at once
figure 
errorbar(Egglength, AvgAmpAllAP(30).AvgAmp, AvgAmpAllAP(30).All95Conf, 'Color', Colors(30).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(6).AvgAmp, AvgAmpAllAP(6).All95Conf, 'Color', Colors(6).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 350000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(30).AvgAmp./MS2Convers), (AvgAmpAllAP(30).All95Conf./MS2Convers), 'Color', Colors(30).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(2).AvgAmp./MS2Convers), (AvgAmpAllAP(2).All95Conf./MS2Convers), 'Color', Colors(2).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(8).AvgAmp./MS2Convers), (AvgAmpAllAP(8).All95Conf./MS2Convers), 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(6).AvgAmp./MS2Convers), (AvgAmpAllAP(6).All95Conf./MS2Convers), 'Color', Colors(6).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(7).AvgAmp./MS2Convers), (AvgAmpAllAP(7).All95Conf./MS2Convers), 'Color', Colors(7).Color, 'LineWidth', 2.5,'LineStyle','-');
ylabel('transcripts produced','Color','k');
set(gca,'YColor','k');
ylim([0 (350000/MS2Convers)]);