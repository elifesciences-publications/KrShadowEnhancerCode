%% Estimate ability to discern 2 transcriptional spots in each nucleus 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'Kr2xProxEmpty';'KrDist17C';'KrBoth17C';'Kr2xDist32C';'Kr2xDist17C';'KrProx17C';'Kr2xDistdBcd';'Kr2xProxdHb';'Kr2xProx17C';'Kr2xDistEmpty';'KrSEdBcd';'KrInvSE';'KrSEdHb';'KrEndogDist'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    NumberBursts=nan(15,41);
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        filename2=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(filename); load(filename2);
        for aa=1:length(APbinID)
            NBursts=[];
            BurstsAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(BurstsAP)
                continue
            end
            %All of the nuclei in an AP bin 
            NucleiToUse=[BurstProperties([BurstsAP]).Nucleus]; 
            for nn=1:length(NucleiToUse)
                NucLook=find([BurstProperties.Nucleus]==[NucleiToUse(nn)]); %BurstProperties is done by spots so find the two rows (ie spots) associated with each nucleus
                NBursts=[NBursts,NucLook(1)]; %Count the spots associated with each nucleus, put into a vector for the whole AP bin
                if length(NucLook) > 1 & (~isempty(BurstProperties(NucLook(2)).SmoothTrace))
                    NBursts=[NBursts,NucLook(2)];
                end
            end
            NumberBursts(ee,aa)=length(NBursts); %Group by embryo so not thrown off by some constructs have more/less embryos imaged
        end
    end
    ConstructData(cc).NumberBursts=NumberBursts;
    ConstructData(cc).SDBurstNumber=nanstd(NumberBursts);
    ConstructData(cc).MeanBurstNumber=nanmean(NumberBursts);
    for aa=1:length(APbinID)
        ConstructData(cc).SEBurstNumber(aa)=ConstructData(cc).SDBurstNumber(aa)/(sqrt(sum(~isnan(ConstructData(cc).NumberBursts(:,aa)))));
    end
    ConstructData(cc).NumBursts95CI=ConstructData(cc).SEBurstNumber.*1.95;
end
%% Plotting 
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

fontsize=10;
fontname='Helvetica';
x_width=3; y_width=2.25;
x_widthsplit=1.5; y_widthsplit=1.125;

EggLength=APbinID .* 100;
% Set where to save figures
FigDirect=[DropboxFolder filesep 'Figures'];
% Save structure of data to create figures 
save([DropboxFolder filesep 'Constructs' filesep 'TransvectionEstimates'],'ConstructData');

%Look up the AP bin of max expression for each construct 
load([DropboxFolder filesep 'Constructs' filesep 'ConstructExpressionBins.mat']);
for cc=1:length(ConstructList)
    ContoUse=find(strcmp([ExpressionBins(2,:)],ConstructList{cc}));
    if length(ContoUse) >1
        ContoUse=ContoUse(1);
    end
    ConsList(cc)=ContoUse;
    AvgBursts(cc)=nanmean(ConstructData(cc).NumberBursts(:,ExpressionBins{1,ContoUse}));
end
%%
figure
plot(1,AvgBursts(1),'o','color',Colors(1).Color,'LineWidth',2.5);
hold on
%bar(2,(AvgBursts(4)*2));
plot(2,(nanmean(ConstructData(4).NumberBursts(:,ExpressionBins{1,ConsList(1)}))*2),'o','color',Colors(4).Color,'LineWidth',2.5);
%bar(3,AvgBursts(4));
plot(3,(nanmean(ConstructData(4).NumberBursts(:,ExpressionBins{1,ConsList(1)}))),'o','color',Colors(4).Color,'LineWidth',2.5);
EstmMultFactor(1)=(AvgBursts(4)*2)/AvgBursts(1);
ExpMultFactor(1)=(AvgBursts(1)/(AvgBursts(4)));
ObsvExpect(1)=(1/EstmMultFactor(1));

%Proximal
figure
bar(1,AvgBursts(2));
hold on
bar(2,(AvgBursts(5)*2));
bar(3,AvgBursts(5));
EstmMultFactor(2)=(AvgBursts(5)*2)/AvgBursts(2);
ExpMultFactor(2)=(AvgBursts(2)/(AvgBursts(5)));
ObsvExpect(2)=(1/EstmMultFactor(2));

%2xDistal
figure
bar(1,AvgBursts(6));
hold on
bar(2,(AvgBursts(24)*2));
bar(3,AvgBursts(24));
EstmMultFactor(3)=(AvgBursts(24)*2)/AvgBursts(6);
ExpMultFactor(3)=(AvgBursts(6)/(AvgBursts(24)));
ObsvExpect(3)=(1/EstmMultFactor(3));

%2xProximal
figure
bar(1,AvgBursts(7));
hold on
bar(2,(AvgBursts(15)*2));
bar(3,AvgBursts(15));
EstmMultFactor(4)=(AvgBursts(15)*2)/AvgBursts(7);
ExpMultFactor(4)=(AvgBursts(7)/(AvgBursts(15)));
ObsvExpect(4)=(1/EstmMultFactor(4));

%SE
EstmMultFactor(5)=(AvgBursts(9)*2)/AvgBursts(8);
ExpMultFactor(5)=(AvgBursts(8)/(AvgBursts(9)));
ObsvExpect(5)=(1/EstmMultFactor(5));


figure
ColorstoUse=[1,6,2,7,8]; 
for cc=1:length(ColorstoUse)
plot((cc),EstmMultFactor(cc),'o','color',Colors(ColorstoUse(cc)).Color,'LineWidth',2.5);
hold on 
end
ylabel('Expected/Obs # bursts');

figure
for cc=1:length(ColorstoUse)
    plot(cc,ExpMultFactor(cc),'o','color',Colors(ColorstoUse(cc)).Color,'LineWidth',2.5);
    hold on 
end
ylabel('Homozygote/Hemi # Bursts');

figure
ColorstoUse=[1,6,2,7,8];%,28];
PlaceHolder=0.5;
for cc=1:length(ColorstoUse)
    PlaceHolder=PlaceHolder+1;
errorbar(cc,ObsvExpect(cc),'o','color',Colors(ColorstoUse(cc)).Color,'LineWidth',2.5);
hold on 
end
ylim([0 3]);
xlim([0 (length(ColorstoUse)+1)]);
ylabel('Obs/Expected # of spots');
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
print( [FigDirect filesep 'TransvectionEstm'],'-dsvg');
