%% Plot fluorescence traces with called ON/OFF points for visual inspection 
% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList= {'KrBoth','KrProxDuplic','KrDistDuplicN','KrProx','KrEndogDist'};%{'KrEndogDist','KrProx','KrDistDuplicN','KrProxDuplic','KrBoth'};
% Set our pathway 
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
% Load the max expression data 
load([DropboxFolder filesep 'Constructs' filesep 'ConstructExpressionBins.mat']);

BothColor=[52 119 71]./255;
fontsize=10;
fontname='Helvetica';

%Load BurstProperties and CompiledParticles
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    NEmbryos = length(Data);
    ContoUse=find(strcmp([ExpressionBins(2,:)],ConstructList{cc}));
    if length(ContoUse) >1
        ContoUse=ContoUse(1);
    end
    ConsList(cc)=ExpressionBins{1,ContoUse};
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        load(filename);
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(CompPars);
        APBinUsed=ConsList(cc);
        APBins=[BurstProperties.APBin];
        BurstProperties_AP=BurstProperties(APBins==APbinID(APBinUsed));
        
        if ~isempty(BurstProperties_AP)
            y=datasample([1:length(BurstProperties_AP)],2);
            for nn=y(1):y(end)
                figure
                if ~isempty(BurstProperties_AP(nn).SmoothTrace)
                if ((ElapsedTime(end))-ElapsedTime(nc14)) > 51 
                plot([0:100],BurstProperties_AP(nn).SmoothTrace,'Color',BothColor,'LineWidth',2.5');
                
                hold on
                H1=area([BurstProperties_AP(nn).ONFrames(1):BurstProperties_AP(nn).ONFrames(2)],[BurstProperties_AP(nn).SmoothTrace([BurstProperties_AP(nn).ONFrames(1):BurstProperties_AP(nn).ONFrames(2)])]);
                H1.FaceColor=BothColor;
                plot([BurstProperties_AP(nn).ONFrames],BurstProperties_AP(nn).SmoothTrace(BurstProperties_AP(nn).ONFrames),'o','Color','k','LineWidth',2);
                plot([BurstProperties_AP(nn).OFFFrames],BurstProperties_AP(nn).SmoothTrace(BurstProperties_AP(nn).OFFFrames),'*','Color','r','LineWidth',2);
                xticks([0:20:100]);
                xticklabels([0:20:100]./2);
                
                else
                    plot([1:length(ElapsedTime(nc14:end))],BurstProperties_AP(nn).SmoothTrace,'Color', BothColor,'LineWidth',2.5);
                    hold on
                    H1=area([BurstProperties_AP(nn).ONFrames(1):BurstProperties_AP(nn).ONFrames(2)],[BurstProperties_AP(nn).SmoothTrace([BurstProperties_AP(nn).ONFrames(1):BurstProperties_AP(nn).ONFrames(2)])]);
                    H1.FaceColor=BothColor;
                    plot([BurstProperties_AP(nn).ONFrames],BurstProperties_AP(nn).SmoothTrace(BurstProperties_AP(nn).ONFrames),'o','Color','k','LineWidth',2);
                    plot([BurstProperties_AP(nn).OFFFrames],BurstProperties_AP(nn).SmoothTrace(BurstProperties_AP(nn).OFFFrames),'*','Color','r','LineWidth',2);
                    xticks([0:20:100]);
                    xticklabels([0:20:100]./2);
                end
                xlabel('time into nc14 (min)');
                ylabel('fluorescence intensity (AU)');
                title([PrefixName,'nucleus',num2str(BurstProperties_AP(nn).Nucleus)]);
                set(gca, 'FontSize', fontsize, 'FontName', fontname);
            end
            end
        
        end
        clear BurstProperties BurstProperties_AP ElapsedTime nc14 CompiledParticles ElapsedTime
    end
end