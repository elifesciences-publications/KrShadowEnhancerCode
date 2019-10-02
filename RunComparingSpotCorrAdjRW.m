%% Run ComparingSpotsRWAdj for all data sets - generates SpotDiff.m files for all embryos approved in DataStatus sheet
%load constructs
ConstructList={'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDistEmpty';'KrProxEmpty';'KrBothEmpty';'Kr2xDistEmpty';'Kr2xProxEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'KrDist17C';'Kr2xDist32C';'KrBoth17C';'Kr2xDist17C';'KrProx17C';'Kr2xDistdBcd';'Kr2xProxdHb';'Kr2xProx17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
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
        FileName2=[DropboxFolder filesep PrefixName filesep 'APDetection.mat']
        FileName3=[DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat']
        load(CompPars);
        load(FileName2);
        load(FileName3);
        ComparingSpotsRWAdj
        FileName=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        save(FileName, 'SpotDiff');
        clear SpotDiff PrefixName CompiledParticles FileName2 FileName3 CompPars
    end
end