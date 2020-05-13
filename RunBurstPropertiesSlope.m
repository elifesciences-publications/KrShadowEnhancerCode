%% Run SlopeBurstCalling for all data sets
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'Kr2xProxEmpty';'KrDist17C';'KrBoth17C';'Kr2xDist17C';'Kr2xDistdBcd';'Kr2xProxdHb';'KrProx17C';'Kr2xProx17C';'Kr2xDistEmpty';'KrInvSE';'KrSEdBcd';'KrSEdHb';'KrEndogDist';'KrEndogDist32C';'Kr3_1xBcd';'Kr4_1xBcd';'Kr3_1xHb';'Kr3_1xZld';'Kr3_1xStat92E';'Kr4_1xHb';'Kr4_1xZld';'Kr5_1xHb';'Kr5_1xStat92E'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
%'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
for cc=1:length(ConstructList)
    Data=LoadMS2SetsCS(ConstructList{cc});
    NEmbryos = length(Data);
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat']
        load(CompPars);
        Spotz=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat']
        load(Spotz);
        SlopeBurstCallingRW
        FileName=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        save(FileName, 'BurstProperties');
        clear BurstProperties PrefixName CompiledParticles FileName CompPars
    end
end