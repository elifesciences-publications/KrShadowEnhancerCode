%% calculate average burst frequency for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'Kr2xProxEmpty';'KrDist17C';'KrBoth17C';'Kr2xDist32C';'Kr2xDist17C';'KrProx17C';'Kr2xDistdBcd';'Kr2xProxdHb';'Kr2xProx17C';'Kr2xDistEmpty';'KrSEdBcd';'KrInvSE';'KrSEdHb';'KrEndogDist';'KrEndogDist32C';'KrEndogDist17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Count for each construct
AvgFreq=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgFreqCon=[];
    firsttime=1;
    ConFreqAllAP=[];
    ConFreqSE=[];
    ConFreqSD=[];
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        load(filename);
        Frequency(ee,cc)=length([BurstProperties.Frequency]);  
          
        %seperate out by AP bin
        FreqAllAP=[];
        for aa=1:length(APbinID)
            FreqAP=[];
            FrequencyAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(FrequencyAP)
                FreqAP=[FreqAP; nan];
            else
            for bb=1:length(FrequencyAP)

                   FreqAP=[FreqAP;[BurstProperties(FrequencyAP(bb)).Frequency]']; 
            
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(FreqAP)
                FreqAllAP(bb,aa,ee)=FreqAP(bb);
            end
            FreqAllAP(FreqAllAP==0)=nan;
            
            FrequencySD(ee,aa,cc)=nanstd(FreqAllAP(:,aa,ee));
        FrequencySE(ee,aa,cc)=FrequencySD(ee,aa,cc)/sqrt(sum(~isnan(FreqAP)));             
            clear FrequencyAP
        end
        %Compile all data for a construct in one long column
        %to put in structure
        for bb=1:size(FreqAllAP,3)
            ConFreqAllAP=[ConFreqAllAP; FreqAllAP(:,:,bb)];
        end
        for bb=1:size(FrequencySD,3)
            ConFreqSD=[ConFreqSD; FrequencySD(:,:,bb)];
        end
        for bb=1:size(FrequencySE,3)
            ConFreqSE=[ConFreqSE;FrequencySE(:,:,bb)];
        end
        AvgFreqAllAP(cc).EmbryoFreq(ee).AvgFreq=nanmean(FreqAllAP(:,:,ee));
        AvgFreqAllAP(cc).EmbryoFreq(ee).Freq=[FreqAllAP(:,:,ee)];
        AvgFreqAllAP(cc).EmbryoFreq(ee).SE=FrequencySE(ee,:,cc);
        AvgFreqAllAP(cc).EmbryoFreq(ee).SD=FrequencySD(ee,:,cc);
    end
        AvgFreqAllAP(cc).AvgFreq=nanmean(ConFreqAllAP);  %Avg frequency of all embryos of a construct
        AvgFreqAllAP(cc).FreqSD=nanmean(ConFreqSD);
        AvgFreqAllAP(cc).FreqSE=nanmean(ConFreqSE);
        AvgFreqAllAP(cc).AllFreq=[ConFreqAllAP];
        AvgFreqAllAP(cc).AllFreqSD=nanstd(AvgFreqAllAP(cc).AllFreq);
        
        for aa=1:length(APbinID)
        AvgFreqAllAP(cc).AllFreqSE(aa)=([AvgFreqAllAP(cc).AllFreqSD(aa)])/(sqrt(sum(~isnan(AvgFreqAllAP(cc).AllFreq(:,aa)))));
        end
        AvgFreqAllAP(cc).AllFreq95Conf=[AvgFreqAllAP(cc).AllFreqSE] .* 1.95;
        AvgFreqAllAP(cc).ConstructName=ConstructList{cc};

end
% Get rid of single data points for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgFreqAllAP(cc).AllFreq(:,aa)))==1
            AvgFreqAllAP(cc).AllFreqSD(aa)=nan;
            AvgFreqAllAP(cc).AvgFreq(aa)=nan;
            AvgProdAllAP(cc).AllFreqSE(aa)=nan;
            AvgProdAllAP(cc).AllFreq95Conf(aa)=nan;
        end
    end
end

%% Plot mean of each embryo at specific AP position for each construct
Egglength=APbinID.*100;


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


Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DistalColor;
Colors(5).Color=ProxColor;
Colors(6).Color=DoubDistColor;
Colors(7).Color=DoubProxColor;
Colors(8).Color=BothColor;
Colors(9).Color=BothColor;
Colors(10).Color=DistalColor
Colors(11).Color=ProxColor;
Colors(12).Color=BothSepColor;
Colors(13).Color=BothColor;
Colors(14).Color=DoubProxColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color=DistalColor;
Colors(17).Color=BothColor;
Colors(18).Color=DoubDistColor;
Colors(19).Color=DoubDistColor;
Colors(20).Color=ProxColor;
Colors(21).Color=DoubDistColor;
Colors(22).Color=DoubProxColor;
Colors(23).Color=DoubProxColor;
Colors(24).Color=DoubDistColor;
Colors(25).Color=BothColor;
Colors(26).Color=BothColor;
Colors(27).Color=BothColor;
Colors(28).Color=DistalColor;
Colors(29).Color=DistalColor;
Colors(30).Color=DistalColor;

%Set display/font settings 
fontsize=12;
fontname='Helvetica';
x_width=3; y_width=2.25;
xSize = 6; ySize = 4.5; xLeft = 0.1; yTop = 0.1;
%Folder to save in 
FigDirect=[DropboxFolder filesep 'Figures' filesep 'Transcriptional dynamics' filesep 'Frequency'];
%Save structure needed to create all graphs
save([DropboxFolder filesep 'Constructs' filesep 'BurstFrequencySlope'],'AvgFreqAllAP');

%% Single/duplicated enhancers vs SE pair 
figure 
errorbar(Egglength, AvgFreqAllAP(28).AvgFreq, AvgFreqAllAP(28).AllFreq95Conf,'Color', Colors(28).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(2).AvgFreq, AvgFreqAllAP(2).AllFreq95Conf,'Color', Colors(2).Color,'LineWidth',2.5)
errorbar(Egglength, AvgFreqAllAP(8).AvgFreq, AvgFreqAllAP(8).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5)
%legend('Distal', 'Proximal','Both')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print([FigDirect filesep 'SinglesvBothFreq_ED'],'-dsvg');

%Duplicates vs SE
figure 
errorbar(Egglength,AvgFreqAllAP(6).AvgFreq,AvgFreqAllAP(6).AllFreqSE,'Color',Colors(6).Color,'LineWidth',2.5);
hold on
errorbar(Egglength,AvgFreqAllAP(7).AvgFreq,AvgFreqAllAP(7).AllFreqSE,'Color',Colors(7).Color,'LineWidth',2.5);
errorbar(Egglength,AvgFreqAllAP(8).AvgFreq,AvgFreqAllAP(8).AllFreqSE,'Color',Colors(8).Color,'LineWidth',2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per min');
%% Endogenous Distal Comparison
figure 
errorbar(Egglength, AvgFreqAllAP(1).AvgFreq, AvgFreqAllAP(1).AllFreq95Conf,'Color', Colors(1).Color,'LineWidth',2.5, 'LineStyle',':')
hold on 
errorbar(Egglength, AvgFreqAllAP(28).AvgFreq, AvgFreqAllAP(28).AllFreq95Conf,'Color', Colors(28).Color,'LineWidth',2.5)
%legend('Both', '1x Both')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
set(gcf,'PaperUnits','inches','PaperPosition',[xLeft yTop xSize ySize]);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'EndogvsOrigDistFreq'],'-dpdf');

%% Compare all at once
figure 
errorbar(Egglength, AvgFreqAllAP(28).AvgFreq, AvgFreqAllAP(28).AllFreq95Conf,'Color', Colors(28).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(2).AvgFreq, AvgFreqAllAP(2).AllFreq95Conf,'Color', Colors(2).Color,'LineWidth',2.5)
errorbar(Egglength, AvgFreqAllAP(6).AvgFreq, AvgFreqAllAP(6).AllFreq95Conf,'Color', Colors(6).Color,'LineWidth',2.5)
errorbar(Egglength, AvgFreqAllAP(7).AvgFreq, AvgFreqAllAP(7).AllFreq95Conf,'Color', Colors(7).Color,'LineWidth',2.5)
errorbar(Egglength, AvgFreqAllAP(8).AvgFreq, AvgFreqAllAP(8).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5)
%legend('Distal', 'Proximal','Both')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');