% Runs simulated annealing under different assumptions. 
% --Inputs: 
% model: common
% strat: sifdc
% init: initial temperature for annealing
% random: 1 for random initial conditions, 0 for 
% analysis-derived initial conditions
% NOTE - input all parameters as strings 
% --Outputs:
% storage: the fitted parameter set
% out: the error associated with the fitted parameter set
% energies: the SSE over time
% iters: the total number of iterations 
% times: the amount of time spent per iteration
% NOTE - all the above quantities will be written to 
% text files in the same directory.
% Last updated: 04/20/20
function [storage, out, energies, iters, times] = ...
 mainAND(model,strat,init,random)

    tic
    % true for the production rate of T* to be uniform
    flat = false;

    % so that randomness is generated from the current time
    rng('shuffle')
    
    % number of nuclei per bin
    extra.numSpots = 80;
    % number of tries per temperature
    tries = 30;
    % time limit of 24 hours = 86400
    timeLim = 86400;
    % maximum number off tries in total in the anneal
    maxRejectAnneal = Inf;
    % cooling and stopping temperature for the annealing
    cool = 0.8;
    % fraction of the original temperature at which annealing 
    stoppingFraction  = 10;   

    % storing args in their original string form
    randomOrig = random;
    initOrig = init; 
    % arguments to compiled matlab code are passed as
    % strings
    random = str2double(random);
    init = str2double(init);
    
    % range for simulation
    upperEgg = 80; 
    lowerEgg = 20; 

    % we are dividing raw fluorescence data by F_rnap or F_1 since our model
    % operates in transcript count.
    F1 = 1338;

    files = dir(strcat(pwd,'/allData'));
    for myFiles = 3:length(files)
        myStruct(myFiles - 2).files = load(strcat(pwd,'/allData/',files(myFiles).name));
    end

    DistalCorr = myStruct(1).files.DistalCorr;
    ProximalCorr = myStruct(14).files.ProximalCorr;

    ProxBurstSize = myStruct(9).files.ProxBurstSize;
    ProximalDuration = myStruct(10).files.ProximalDuration;
    ProxFreq = myStruct(11).files.ProxFreq;
    ProxInterBurstTime = myStruct(12).files.ProxInterBurstTime;
    ProxIF = myStruct(13).files.ProxTotalProd;
    % fixing the duration vectors that are not the right dimensions
    ProximalDuration = ProximalDuration';

    files = dir(strcat(pwd,'/allDistalData'));
    for myFiles = 3:length(files)
        myStruct(myFiles - 2).files = load(strcat(pwd,'/allDistalData/',files(myFiles).name));
    end

    DistBurstSize = myStruct(1).files.DistBurstSize;
    DistDuration = myStruct(2).files.DistalDuration;
    DistFreq = myStruct(3).files.DistFreq;
    DistInterBurst = myStruct(4).files.DistInterBurstTime;
    DistIF = myStruct(5).files.DistTotalProd; 

    DistDuration = DistDuration';

    files = dir(strcat(pwd,'/allCorr'));
    for myFiles = 3:length(files)
        myStruct(myFiles - 2).files = load(strcat(pwd,'/allCorr/',files(myFiles).name));
    end

    ShadowCorr = myStruct(1).files.SepCorr;


    proxDurIdeal = nanmean(ProximalDuration);
    proxBetweenBurstsIdeal = nanmean(ProxInterBurstTime);
    distDurIdeal = nanmean(DistDuration);
    distBetweenBurstsIdeal = nanmean(DistInterBurst);


    burstLength = proxDurIdeal;
    betweenBursts = proxBetweenBurstsIdeal;
    expFreq = ProxFreq;
    expDur = ProximalDuration;
    expCorr = ProximalCorr;
    expSize = ProxBurstSize/F1;
    expIF = ProxIF/F1;

    burstLength2 = distDurIdeal;
    betweenBursts2 = distBetweenBurstsIdeal;
    expFreq2 = DistFreq;
    expDur2 = DistDuration; 
    expCorr2 = DistalCorr;
    expSize2 = DistBurstSize/F1;
    expIF2 = DistIF/F1;

    expCorrShadow = ShadowCorr;

    % Note that each of the data vectors have 41 entries in
    % them therefore each bin corresponds to around 2.43% in egg length. 
    % This implies that the first and last 20% of the egg length
    % are described by the first and last 8 bins of those 41 bins.
    % bins. Therefore, we only focus on bins 8 to 33.
    lowestBin = 8;
    highestBin = 33; 

    yFreq = expFreq(lowestBin:highestBin);
    yDur = expDur(lowestBin:highestBin);
    ySize = expSize(lowestBin:highestBin);
    yIF = expIF(lowestBin:highestBin); 

    yFreq2 = expFreq2(lowestBin:highestBin);
    yDur2 = expDur2(lowestBin:highestBin);
    ySize2 = expSize2(lowestBin:highestBin);
    yIF2 = expIF2(lowestBin:highestBin); 

    yCorr = expCorr(lowestBin:highestBin);
    yCorr2 = expCorr2(lowestBin:highestBin);
    yCorrShadow = expCorrShadow(lowestBin:highestBin);

    % holds the maximum values for each burst property for 
    % subsequent 
    yGlobal = [ySize/max(ySize) yIF/max(yIF) yFreq/max(yFreq) yDur/max(yDur) yCorr/max(yCorr) ...
              ySize2/max(ySize2) yIF2/max(yIF2) yFreq2/max(yFreq2) yDur2/max(yDur2) yCorr2/max(yCorr2) ...
              yCorrShadow/max(yCorrShadow)];

    yMaxGlobal =  [max(ySize) max(yIF) max(yFreq) max(yDur) max(yCorr) ...
                    max(ySize2) max(yIF2) max(yFreq2) max(yDur2) max(yCorr2) ...
                    max(yCorrShadow)];


    % so we can index later
    sizeIndex = 1;
    IFIndex = 2;
    freqIndex = 3;
    durIndex = 4;
    corrIndex = 5;

    sizeIndex2 = 6;
    IFIndex2 = 7;
    freqIndex2 = 8;
    durIndex2 = 9;
    corrIndex2 = 10;

    corrIndexShadow = 11;

    % k_onIC1 = 10.^((2--3).*rand + -3);
    % sdIC1 = randi(20);
    % nIC1 = floor(10^(2*rand + 1));

    if random
        nIC = floor(10^(1*rand + 1));
        scalingIC = 10.^((1--3).*rand + -3); 
        sdIC = randi([5 10],1);
        degIC =  10.^((1--3).*rand + -3);  
        k_onIC = 10.^((1--3).*rand + -3);
        k_offIC = 10.^((1--3).*rand + -3); 

        nIC2 = floor(10^(1*rand + 1));
        scalingIC2 = 10.^((1--3).*rand + -3); 
        sdIC2 = randi([5 10],1);
        degIC2 =  10.^((1--3).*rand + -3);  
        k_onIC2 = 10.^((1--3).*rand + -3);
        k_offIC2 = 10.^((1--3).*rand + -3); 

        nStarIC = floor(10^(1*rand + 1));
        scalingStarIC = 10.^((1--3).*rand + -3);
        sdStarIC = randi([5 10],1);
        degStarIC = 10.^((1--3).*rand + -3);
    
    else 

        nIC = 5;
        scalingIC = 40; 
        sdIC = 6;
        degIC =  1;  
        k_onIC = 1/betweenBursts;
        k_offIC = 1/burstLength; 

        nIC2 = 5;
        scalingIC2 = 40; 
        sdIC2 = 6;
        degIC2 =  1;  
        k_onIC2 = 1/betweenBursts2;
        k_offIC2 = 1/burstLength2; 

        nStarIC = 5;
        % we are assuming a flat distribution for T*
        scalingStarIC = 1;
        % the sd in this case is a trimFlat
        sdStarIC = 1;
        degStarIC = 1;

    end


    % annealing will stop after reaching a fraction of the original 
    % temperature
    stop = init/stoppingFraction;

    % setting the conditions for annealing
    def = struct(...    
            'CoolSched',@(T) (cool*T),...
            'Generator',@(x,t) myRand(x,t),...
            'InitTemp',init,...
            'MaxConsRej',maxRejectAnneal,...
            'MaxSuccess',Inf,...
            'MaxTries',tries,...
            'StopTemp',stop,...
            'StopVal',0,...
            'Verbosity',1);

    extra.alpha = 1.9540; 
    % length of simulation
    extra.endSpan = 50;
    % mean for beta_1,gamma_1
    enhancer.mean = 50;
    % number of bins to simulate
    extra.numBinsRWdata = 26;
    extra.eggLength = linspace(lowerEgg, upperEgg, extra.numBinsRWdata);
    extra.ElapsedTime = 0:0.5:extra.endSpan;


    % proximal 
    enhancer.rFull =  [0   0   67.0257   60.0482   57.4959   51.5140   53.9837   53.1280 ...
       60.4217   90.3269  104.0844  114.0616  119.2144  101.4468   87.6735   79.4335 ...
       68.7479   63.0808   86.9607   0   0   0   0   0 ...
       0   0];
    % distal at proximal location
    enhancer.rFull2 = [0   0   0   0   0   0   0   0   58.8609 ...
       74.4993   86.2019  111.3550  141.3504  131.6084  127.4489  115.2295   98.5590 ...
       95.3566   81.3611   71.7599   67.8933   59.8787   60.7270   59.1140   62.3468 ... 
       64.2655 ];


    xIC = [degIC, k_offIC, k_onIC, nIC, scalingIC, sdIC, ... 
        degIC2, k_offIC2, k_onIC2, nIC2, scalingIC2, sdIC2, ...
        degStarIC, ...
        nStarIC, scalingStarIC, sdStarIC];


    % indexes of different params in the vector
    degIndex = 1; 
    koffIndex = 2;
    konIndex = 3;   
    nIndex = 4;
    scalingIndex = 5;
    sdIndex = 6;

    degIndex2 = 7; 
    koffIndex2 = 8;
    konIndex2 = 9;   
    nIndex2 = 10;
    scalingIndex2 = 11;
    sdIndex2 = 12;

    degStarIndex = 13;  
    nStarIndex = 14;
    scalingStarIndex = 15;
    sdStarIndex = 16;

    f = @paramOptimizer;
    g = @(p)f(p);

    [storage, out, energies, iters, times] =  anneal(g, xIC, def);

    % writing to file
    format = '%f, ';

    % predefining names for the files based on parameter, model type, strategy, initial temperature
    % and random ICs
    scalingFile = strcat('scaling',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    degFile = strcat('deg',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    sdFile = strcat('sd',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    konFile = strcat('kon',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    koffFile = strcat('koff',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    nFile = strcat('n',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');

    scalingStarFile = strcat('scalingStar',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    degStarFile = strcat('degStar',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    sdStarFile = strcat('sdStar',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    nStarFile = strcat('nStar',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');

    scalingFile2 = strcat('scaling2',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    degFile2 = strcat('deg2',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    sdFile2 = strcat('sd2',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    konFile2 = strcat('kon2',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    koffFile2 = strcat('koff2',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    nFile2 = strcat('n2',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');

    energiesFile = strcat('energies',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    timesFile = strcat('times',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    itersFile = strcat('iters',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');
    finalEnergyFile = strcat('finalEnergies',model,strat,num2str(initOrig), num2str(randomOrig),'.txt');

    % writing to file
    fid = fopen(nFile,'a');
    fprintf(fid,format,storage(nIndex));
    fclose all;
    fid = fopen(scalingFile,'a');
    fprintf(fid,format,storage(scalingIndex));
    fclose all;
    % writing to file
    fid = fopen(sdFile,'a');
    fprintf(fid,format,storage(sdIndex));
    fclose all;
    % writing to file
    fid = fopen(degFile,'a');
    fprintf(fid,format,storage(degIndex));
    fclose all;
    % writing to file
    fid = fopen(konFile,'a');
    fprintf(fid,format,storage(konIndex));
    fclose all;
    % writing to file
    fid = fopen(koffFile,'a');
    fprintf(fid,format,storage(koffIndex));
    fclose all;


    fid = fopen(nStarFile,'a');
    fprintf(fid,format,storage(nStarIndex));
    fclose all;
    fid = fopen(scalingStarFile,'a');
    fprintf(fid,format,storage(scalingStarIndex));
    fclose all;
    % writing to file
    fid = fopen(sdStarFile,'a');
    fprintf(fid,format,storage(sdStarIndex));
    fclose all;
    % writing to file
    fid = fopen(degStarFile,'a');
    fprintf(fid,format,storage(degStarIndex));
    fclose all;



    fid = fopen(nFile2,'a');
    fprintf(fid,format,storage(nIndex2));
    fclose all;
    fid = fopen(scalingFile2,'a');
    fprintf(fid,format,storage(scalingIndex2));
    fclose all;
    % writing to file
    fid = fopen(sdFile2,'a');
    fprintf(fid,format,storage(sdIndex2));
    fclose all;
    % writing to file
    fid = fopen(degFile2,'a');
    fprintf(fid,format,storage(degIndex2));
    fclose all;
    % writing to file
    fid = fopen(konFile2,'a');
    fprintf(fid,format,storage(konIndex2));
    fclose all;
    % writing to file
    fid = fopen(koffFile2,'a');
    fprintf(fid,format,storage(koffIndex2));
    fclose all;

    fid = fopen(energiesFile,'a');
    fprintf(fid,format,energies);
    fprintf(fid,'\n','');
    fclose all;

    fid = fopen(finalEnergyFile,'a');
    fprintf(fid,format,energies(end));
    fprintf(fid,'\n','');
    fclose all;

    fid = fopen(itersFile,'a');
    fprintf(fid,format,iters);
    fprintf(fid,'\n','');
    fclose all;


    fid = fopen(timesFile,'a');
    fprintf(fid,format,times);
    fprintf(fid,'\n','');
    fclose all;

    toc

    % Prepares the error function for anneal.m
    % --Input: 
    % x: the vector of parameters to be fitted
    % --Output:
    % error: the error of the data given by simulating
    % with x
    % Alvaro Fletcher
    function error = paramOptimizer(x) 



        % note that for fixed TF concentrations
        % we have that enhancer.tf represents the
        enhancer.degTF = x(1);  
        enhancer.k_off =  x(2);
        enhancer.k_on = x(3); 
        enhancer.tf = x(4);
        enhancer.scaling = x(5); 
        enhancer.sd = x(6);

        enhancer.degTF2 = x(7);
        enhancer.k_off2 =  x(8);
        enhancer.k_on2 = x(9);
        enhancer.tf2 = x(10);
        enhancer.scaling2 = x(11); 
        enhancer.sd2 = x(12);
  
        % this can also be trimFlat
        enhancer.degTFstar = x(13);   
        enhancer.tfstar = x(14);
        enhancer.scalingstar = x(15); 
        enhancer.sdstar = x(16);               

        simData = mainfunction(enhancer, extra);

        sizes =  simData.sizes/yMaxGlobal(sizeIndex);
        IF = simData.IF/yMaxGlobal(IFIndex);
        freqs = simData.freqs/yMaxGlobal(freqIndex);
        durs = simData.durs/yMaxGlobal(durIndex);
        corrs = simData.corr/yMaxGlobal(corrIndex);   

        sizes2 =  simData.sizes2/yMaxGlobal(sizeIndex2);
        IF2 = simData.IF2/yMaxGlobal(IFIndex2);
        freqs2 = simData.freqs2/yMaxGlobal(freqIndex2);
        durs2 = simData.durs2/yMaxGlobal(durIndex2);
        corrs2 = simData.corr2/yMaxGlobal(corrIndex2);   

        corrsShadow = simData.corrShadow/yMaxGlobal(corrIndexShadow);  

        simVec = [sizes, IF, freqs, durs, corrs, ...
                  sizes2, IF2, freqs2, durs2, corrs2, ...
                  corrsShadow];

        error =  SSEcalculator(simVec, yGlobal);  

    end

    % ANNEAL  Minimizes a function with the method of simulated annealing
    % (Kirkpatrick et al., 1983)
    %   joachim.vandekerckhove@psy.kuleuven.be
    %   $Revision: v5 $  $Date: 2006/04/26 12:54:04 $
    function [minAnneal,fval,energiesAnneal,total,timesAnneal] = anneal(loss, parent, def)

        options=def;

        % main settings
        newsol = options.Generator;      % neighborhood space function
        Tinit = options.InitTemp;        % initial temp
        minT = options.StopTemp;         % stopping temp
        coolAnneal = options.CoolSched;        % annealing schedule
        minF = options.StopVal;
        max_consec_rejections = options.MaxConsRej;
        max_try = options.MaxTries;
        max_success = options.MaxSuccess;
        boltzmannC = 1;                           % boltzmann constant
        
        % counters etc
        itry = 0;
        success = 0;
        finished = 0;
        consec = 0;
        Tanneal = Tinit;
        annealTic = tic;
        initenergy = loss(parent);
        oldenergy = initenergy;

        energiesAnneal = [];
        timesAnneal = [];
        
        total = 0;
        fprintf(1,'\n  Tanneal = %7.5f, loss = %10.5f\n',Tanneal,oldenergy); 

        while ~finished
            energiesAnneal = [energiesAnneal oldenergy];
            timesAnneal = [timesAnneal toc(annealTic)];
            annealTic = tic;
            itry = itry + 1; % just an iteration counter
            current = parent; 
            % stop the program if it takes more 
            % than 24 hours
            if toc > timeLim
                minAnneal = parent;
                fval = oldenergy;

                fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
                fprintf(1, '  Final temperature:       \t%g\n', Tanneal);
                fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
                fprintf(1, '  Number of function calls:\t%i\n', total);
                fprintf(1, '  Total final loss:        \t%g\n', fval);
               
                return
            end
            
            % % Stop / decrement Tanneal criteria
            if itry >= max_try || success >= max_success
                if Tanneal < minT || consec >= max_consec_rejections
                    finished = 1;
                    total = total + itry;
                    break;
                else
                    Tanneal = coolAnneal(Tanneal); % decrease Tanneal according to cooling schedule
                    fprintf(1,'  Tanneal = %7.5f, loss = %10.5f\n',Tanneal,oldenergy);
                    total = total + itry;
                    itry = 1;
                    success = 1;
                end
            end
            
            newparam = newsol(current, Tanneal);
            newenergy = loss(newparam);
            
            if (newenergy < minF)
                parent = newparam; 
                oldenergy = newenergy;
                break
            end
            
            if (oldenergy-newenergy > 0)
                parent = newparam;
                oldenergy = newenergy;
                success = success + 1;
                consec = 0;
            else
                if (rand < exp( (oldenergy-newenergy)/(boltzmannC*Tanneal) ))
                    parent = newparam;
                    oldenergy = newenergy; 
                    % it's counted as a success but also 
                    % a failure since consec is not updated
                    success = success + 1;
                else
                    consec = consec + 1;
                end
            end
        end

        minAnneal = parent;
        fval = oldenergy;

        fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
        fprintf(1, '  Final temperature:       \t%g\n', Tanneal);
        fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
        fprintf(1, '  Number of function calls:\t%i\n', total);
        fprintf(1, '  Total final loss:        \t%g\n', fval);
        

    end

    % Updates parameters for annealing routine
    % -- Input:
    % myParams: vector of parameters
    % currentTemp: the current temperature
    % -- Output: 
    % out: the updated parameter set
    % Alvaro Fletcher
    function out = myRand(myParams, currentTemp) 

        myBool = true;  
        % lower bound tolerance
        tol = 0;
        % lower bound for TFs 
        tolTf = 1;
        % upperTolTf star how much we can trim from each side in flat
        upperTolTfStar = 10;

        while myBool

            tempOut = myParams;
            % sample from normal dist centered at 1 with 
            % standard dev multiplied by current temperature and multiply to a randomly chosen parameter
            newMean = 1;
            newSigma = currentTemp;
            % choose a random parameter and add sampled value
            holderParams = (randperm(length(myParams)) == length(myParams)) * (newMean + randn * newSigma);
            % finding index of parameter to update and value by which it will be multiplied
            holderVal = nonzeros(holderParams);
            holderIndex = find(holderParams);
            tempOut(holderIndex) = myParams(holderIndex) * holderVal;

            % the 1st and 7th elements are the number of TFs in a burst 
            % and we force this number to always be an integer 
            % and greater than 0.This is also the case when we are 
            % fitting to fixed concentrations of TF.
            tempOut(4) = floor(tempOut(4));
            tempOut(10) = floor(tempOut(10));
            tempOut(14) = floor(tempOut(14));
            % flooring sd star index
            if flat
                tempOut(16) = floor(tempOut(16)); 
            end
            % && tempOut(20) >= tol && tempOut(20) <= upperTolTfStar
            % all(tempOut(1:19) > tol) 
            
            % the amount of bins that we trim for the ends of the range
            % for the T* rate of production needs to be an integer
            % note that tolTF > tol
            if flat
                if tempOut(4) > tolTf && tempOut(10) > tolTf && ...
                     tempOut(14) > tolTf  && ...
                     all(tempOut(1:15) > tol) && tempOut(16) >= tol && tempOut(16) <= upperTolTfStar

                    out = tempOut;
                    myBool = false;
                    
                end

            else
                if tempOut(4) > tolTf && tempOut(10) > tolTf && ...
                     tempOut(14) > tolTf  && ...
                     all(tempOut(1:16) > tol) 

                    out = tempOut;
                    myBool = false;
                    
                end

            end
     
        end
        
    end

    % Calculates the sum of squared errors between 
    % two vectors vec1 and vec2.
    % --Inputs:
    % vec1: the simulated data
    % vec2: the observations
    % --Outputs:
    % error: the sum of squared errors
    % Alvaro Fletcher
    function error = SSEcalculator(vec1, vec2)

        % stores total error
        error = 0;

        % checking to see vectors are the same length
        if length(vec1) ~= length(vec2)
            disp('Your simulated vector and your experimental vector are not the same length');
        end

        % the length of vec1 and vec2 should be the same
        for sseIndex = 1:length(vec2)
            
            if ~(isnan(vec2(sseIndex)) || isnan(vec1(sseIndex)))
                
                tempError = (vec2(sseIndex) - vec1(sseIndex))^2;
                error = error + tempError;
                
            else
                
                % for our purposes, nan means that there was no data recorded
                % therefore we just set the error to the square of whatever
                % quantity is not nan
                if isnan(vec2(sseIndex)) && ~isnan(vec1(sseIndex))
                    
                    tempError = vec1(sseIndex)^2;
                    error = error + tempError;

                elseif isnan(vec1(sseIndex)) && ~isnan(vec2(sseIndex))
                    
                    tempError = vec2(sseIndex)^2;
                    error = error + tempError;
        
                else
                    
                    tempError = 0;
                    error = error + tempError;
            
                end
                
            end
            
        end 

    end 

    % Generates a vector of values that correspond to 
    % the values of mRNA production rates based on the assumption
    % of using a bell curve as production rate across the egg.
    % --Inputs:
    % enhancer: struct containing the mean of the curve
    % extra: struct containing a vector of bin positions
    % --Outputs:
    % prodTF: the production rate of the TF
    % Alvaro Fletcher
    function prodTF = TFproduction(enhancer, extra)

        % range of egg length for beta and gamma
        tfProdLength = extra.eggLength;

        prodTF = enhancer.scaling * mygauss(tfProdLength, enhancer.mean, enhancer.sd);

        function fGauss = mygauss(gaussVar, muGauss, sGauss)

            p1Gauss = -.5 * (( (gaussVar)  - muGauss)/sGauss) .^ 2;
            p2Gauss = (sGauss * sqrt(2 * pi));
            fGauss = exp(p1Gauss) ./ p2Gauss; 

        end 

    end

    function prodTF = TFproduction2(enhancer, extra)

        % range of egg length for beta and gamma
        tfProdLength = extra.eggLength;

        prodTF = enhancer.scaling2 * mygauss(tfProdLength, enhancer.mean, enhancer.sd2);

        function fGauss = mygauss(gaussVar, muGauss, sGauss)

            p1Gauss = -.5 * (( (gaussVar)  - muGauss)/sGauss) .^ 2;
            p2Gauss = (sGauss * sqrt(2 * pi));
            fGauss = exp(p1Gauss) ./ p2Gauss; 

        end 

    end

    % Generates a vector of values that correspond to 
    % the values of mRNA production rates for the TF 
    % that both enhancers have in common with a flat distribution
    % and a given range
    % --Inputs:
    % enhancer: struct containing the mean of the curve
    % extra: struct containing a vector of bin positions
    % --Outputs:
    % prodTFstar: the production rate of the TF
    % Alvaro Fletcher
    function prodTFstar = TFproductionstarFlat(enhancer, extra)

        prodTFstar = enhancer.scalingstar * ones(1,length(extra.eggLength));

        prodTFstar(1:enhancer.sdstar) = 0;
        prodTFstar(end - (enhancer.sdstar - 1):end) = 0;

    end

    function prodTFstar = TFproductionstar(enhancer, extra)

        % range of egg length for beta and gamma
        tfProdLength = extra.eggLength;

        prodTFstar = enhancer.scalingstar * mygauss(tfProdLength, enhancer.mean, enhancer.sdstar);

        function fGauss = mygauss(gaussVar, muGauss, sGauss)

            p1Gauss = -.5 * (( (gaussVar)  - muGauss)/sGauss) .^ 2;
            p2Gauss = (sGauss * sqrt(2 * pi));
            fGauss = exp(p1Gauss) ./ p2Gauss;

        end 

    end

    % Runs simulations of the model.
     % --Inputs:
    % enhancer: struct containing the mean of the curve
    % extra: struct containing a vector of bin positions
    % --Outputs:
    % simData: struct containing burst properties and/or allele
    % correlation of simulated data
    % Alvaro Fletcher
    function simData = mainfunction(enhancer, extra)

        holderEggLength = length(extra.eggLength);
        holderElapsedTime = extra.ElapsedTime;
        holderNumSpots = extra.numSpots;
        
        % we ignore bins with less than one data point
        ignoreThreshold = 1;  

        % function handles for parfor loop
        fcnSingleEnhancer = @singleEnhancer;
        fcnSingleEnhancer2 = @singleEnhancer2;
        fcnDupEnhancer = @duplicatedEnhancer;
        fcnDupEnhancer2 = @duplicatedEnhancer2;
        fcnDiscrete = @crndiscrete;
        fcnShadow = @shadowPair;

        prodTF = TFproduction(enhancer, extra);
        prodTF2 = TFproduction2(enhancer, extra);
        if flat
             prodTFstar = TFproductionstarFlat(enhancer, extra);
        else
             prodTFstar = TFproductionstar(enhancer, extra);
        end

        rateRNA = enhancer.rFull;
        rateRNA2 = enhancer.rFull2;

        numSpotsPar = extra.numSpots;

        % holds the mRNA values for SpotDiff
        tempSimR = zeros(holderEggLength, length(extra.ElapsedTime), holderNumSpots);
        tempCorr = zeros(1,holderEggLength);
        tempIF = zeros(1,holderEggLength);

        tempSimR2 = zeros(holderEggLength, length(extra.ElapsedTime), holderNumSpots);
        tempCorr2 = zeros(1,holderEggLength);
        tempIF2 = zeros(1,holderEggLength);

        % going over different bins in the egg
        parfor i = 1:holderEggLength

            production = prodTF(i);
            productionStar = prodTFstar(i);
            production2 = prodTF2(i);
            myParR = rateRNA(i);
            myParR2 = rateRNA2(i);
      
            % limits the number of iterations for the Gillespie algorithm.
            % However, we want all simulations to be limited by the time span of
            % nc14 and just choose a very large number for maxIterations
            maxIterations = 10000000;

            % the locations of the mRNA and transcription factor for our network
            RNAloc = 3;

            parTempSimR = zeros(length(holderElapsedTime), holderNumSpots);
            parTempIF = zeros(1, holderNumSpots);

            parTempSimR2 = zeros(length(holderElapsedTime), holderNumSpots);
            parTempIF2 = zeros(1, holderNumSpots);

            % going over all the spots in the bin
            for spot = 1:numSpotsPar

                % defining the network for single enhancer with 
                % fluctuating TFs
                CRN = feval(fcnSingleEnhancer, enhancer, extra, production, ... 
                    myParR, productionStar);
                CRN2 = feval(fcnSingleEnhancer2, enhancer, extra, production2, ... 
                    myParR2, productionStar);
               
                % runs the Gillespie algorithm 
                expN = feval(fcnDiscrete, CRN, maxIterations);
                expN2 = feval(fcnDiscrete, CRN2, maxIterations);

                % saving the mRNA concentration over time for this simulation
                R = expN(:,RNAloc); 
                R2 = expN2(:,RNAloc); 

                parTempSimR(:,spot) = R;
                parTempIF(spot) = trapz(holderElapsedTime, R);

                parTempSimR2(:,spot) = R2;
                parTempIF2(spot) = trapz(holderElapsedTime, R2);

            end


            tempSimR(i,:,:) = parTempSimR; 
            tempIF(i) = nanmean(parTempIF);

            tempSimR2(i,:,:) = parTempSimR2; 
            tempIF2(i) = nanmean(parTempIF2);

            % note that specifying "prox" and "dist" does not necessarily mean 
            % the proximal enhancer or distal enhancer but rather whichever 
            % enhancers are found in the proximal and distal locations
            RNAlocProx = 5;
            RNAlocDist = 8;

            parTempCorr = zeros(1, holderNumSpots);
            parTempCorr2 = zeros(1, holderNumSpots);
            parTempCorrShadow = zeros(1, holderNumSpots);

            for spot = 1:numSpotsPar

                % doing the proximal-proximal correlation 
                CRN = feval(fcnDupEnhancer, enhancer, extra, production, ...
                            myParR, myParR, productionStar);
                % doing the distal-distal correlation
                CRN2 = feval(fcnDupEnhancer2, enhancer, extra, production2, ...
                            myParR2, myParR2, productionStar);                
                % doing the shadow correlation
                CRNShadow = feval(fcnShadow, enhancer, extra, production, ...
                            production2, myParR, myParR2, productionStar);

                % runs the Gillespie algorithm 
                expN = feval(fcnDiscrete,CRN,maxIterations);      
                expN2 = feval(fcnDiscrete,CRN2,maxIterations);   
                expNShadow =  feval(fcnDiscrete,CRNShadow, maxIterations); 

                % the RNA locations for the correlation models
                % that include R1 and R2 (RNA produced by enhancer 1
                % and enhancer 2 respectively) are selected per 
                % choice of model at the beginning
                proxR = expN(:,RNAlocProx);
                distR = expN(:,RNAlocDist);

                proxR2 = expN2(:,RNAlocProx);
                distR2 = expN2(:,RNAlocDist);

                proxRShadow = expNShadow(:, RNAlocProx);
                distRShadow = expNShadow(:, RNAlocDist);

                % storing the correlation between R1 and R2 
                corrholder = corrcoef(proxR,distR);        
                corrholder2 = corrcoef(proxR2,distR2);
                corrholderShadow = corrcoef(proxRShadow,distRShadow);
                
                parTempCorr(spot) = corrholder(1,2);
                parTempCorr2(spot) = corrholder2(1,2);
                parTempCorrShadow(spot) = corrholderShadow(1,2);

                
            end

            tempCorr(i) = nanmean(parTempCorr);
            tempCorr2(i) = nanmean(parTempCorr2);
            tempCorrShadow(i) = nanmean(parTempCorrShadow);

        end

        simData.corr = tempCorr;
        simData.IF = tempIF;        

        simData.corr2 = tempCorr2;
        simData.IF2 = tempIF2;        

        simData.corrShadow = tempCorrShadow;

        % building SpotDiff and simData.corr after the parallelizations
        SpotDiff = struct('SpotOne', cell(1, holderEggLength * extra.numSpots), ...
                    'APBin',cell(1, holderEggLength * extra.numSpots));        


        SpotDiff2 = struct('SpotOne', cell(1, holderEggLength * extra.numSpots), ...
                    'APBin',cell(1, holderEggLength * extra.numSpots));

        postCounter = 1;
        for i = 1:length(extra.eggLength)
            for j = 1:extra.numSpots
                SpotDiff(postCounter).SpotOne = tempSimR(i,:,j);
                SpotDiff(postCounter).APBin = extra.eggLength(i);

                SpotDiff2(postCounter).SpotOne = tempSimR2(i,:,j);
                SpotDiff2(postCounter).APBin = extra.eggLength(i);

                postCounter = postCounter + 1;
            end
        end  

        % obtaining data from burst calling algorithm
        BurstProperties = modSlopeBurstCallingRW(SpotDiff, extra.ElapsedTime);
        BurstProperties2 = modSlopeBurstCallingRW(SpotDiff2, extra.ElapsedTime);

        % Initializing...
        simData.sizes = [];
        simData.freqs = [];
        simData.durs = [];
        simData.uniqueBins = [];        

        simData.sizes2 = [];
        simData.freqs2 = [];
        simData.durs2 = [];
        simData.uniqueBins2 = [];

        % storing all unique values in the AP bin field of BurstProperties
        uniqueAPBins = unique([BurstProperties(1:end).APBin]);
        uniqueAPBins2 = unique([BurstProperties2(1:end).APBin]);

        % removing unnecessary NaNs
        uniqueAPBins = uniqueAPBins(~isnan(uniqueAPBins));
        uniqueAPBins2 = uniqueAPBins2(~isnan(uniqueAPBins2));

        % putting all the AP bin data into vector form
        APBinData = [BurstProperties(1:end).APBin];
        APBinData2 = [BurstProperties2(1:end).APBin];

        % sometimes not all AP bins are represented or they can be out of order
        for i = 1:length(uniqueAPBins)
            % finding all indexes in BurstProperties for the current AP bin
            currentAPBinVals =  find(APBinData == uniqueAPBins(i));
            % checking number of empty rows, it suffices to check amplitude
            % since the existence of any burst property implies the 
            % existence of all the others
            % tracks the number of empty rows
            nonEmptyCounter = 0;
            % in case nothing happened for every single bin then no 
            % burst property field in the struct is created (unlikely)
            if isfield(BurstProperties, 'BurstSize')    
                for j = currentAPBinVals       
                    if ~isempty(BurstProperties(j).BurstSize) 
                         nonEmptyCounter = nonEmptyCounter + 1;       
                    end      
                end
            end
            % we ignore AP bins with not enough activity
            if nonEmptyCounter > ignoreThreshold

                simData.sizes(end + 1) = ... 
                      nanmean([BurstProperties(currentAPBinVals).BurstSize]); 
                simData.durs(end + 1) = ...
                      nanmean([BurstProperties(currentAPBinVals).Duration]);   
                simData.freqs(end + 1) = ...
                       nanmean([BurstProperties(currentAPBinVals).Frequency]); 
                simData.uniqueBins(end + 1) = uniqueAPBins(i);
            else 
                simData.sizes(end + 1) = nan;       
                simData.durs(end + 1) = nan;
                simData.freqs(end + 1) = nan; 
                simData.uniqueBins(end + 1) = uniqueAPBins(i);
            end
        end

        % sometimes not all AP bins are represented or they can be out of order
        for i = 1:length(uniqueAPBins2)
            % finding all indexes in BurstProperties for the current AP bin
            currentAPBinVals =  find(APBinData2 == uniqueAPBins2(i));
            % checking number of empty rows, it suffices to check amplitude
            % since the existence of any burst property implies the 
            % existence of all the others
            % tracks the number of empty rows
            nonEmptyCounter = 0;
            % in case nothing happened for every single bin then no 
            % burst property field in the struct is created (unlikely)
            if isfield(BurstProperties2, 'BurstSize')    
                for j = currentAPBinVals       
                    if ~isempty(BurstProperties2(j).BurstSize) 
                         nonEmptyCounter = nonEmptyCounter + 1;       
                    end      
                end
            end
            % we ignore AP bins with not enough activity
            if nonEmptyCounter > ignoreThreshold
                simData.sizes2(end + 1) = ... 
                      nanmean([BurstProperties2(currentAPBinVals).BurstSize]); 
                simData.durs2(end + 1) = ...
                      nanmean([BurstProperties2(currentAPBinVals).Duration]);   
                simData.freqs2(end + 1) = ...
                       nanmean([BurstProperties2(currentAPBinVals).Frequency]); 
                simData.uniqueBins2(end + 1) = uniqueAPBins2(i);
            else 
                simData.sizes2(end + 1) = nan;       
                simData.durs2(end + 1) = nan;
                simData.freqs2(end + 1) = nan; 
                simData.uniqueBins2(end + 1) = uniqueAPBins2(i);
            end
        end

    end

    % Defines a single enhancer network to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF: production rate of TF
    % simR: production rate of mRNA
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN = singleEnhancer(params, extra, prodTF, simR, prodTFstar)

        % number of TFs per burst
        n1 = params.tf;
        n3 = params.tfstar;

        % E = enhancer 1, T = TF, C = complex, R=RNA

        %  E + T  <-> C -> C + R
        %  E + T*  <-> C1* -> C1* + R
        %  R -> 0
        %  n*T <- 0 <- T
        %  n3*T* <- 0 <- T*
        
        reaction=[...
            %   E  C  R  T   T*     E  C  R  T   T*  
                1  0  0  1   1      0  1  0  0   0 
                0  1  0  0   0      1  0  0  1   1 
                0  1  0  0   0      0  1  1  0   0    
                0  0  1  0   0      0  0  0  0   0 
                0  0  0  0   0      0  0  0  n1  0  
                0  0  0  1   0      0  0  0  0   0    
                0  0  0  0   0      0  0  0  0   n3  
                0  0  0  0   1      0  0  0  0   0  
        ];

        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum,1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum,speciesNum+1:2*speciesNum); % target complexes

        % default parameters - any vector with m entries
        CRN.params = [params.k_on params.k_off simR extra.alpha prodTF params.degTF ...
                     prodTFstar params.degTFstar];  

        %           E   C   R   T    T*
        CRN.init = [1   0   0   0    0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

    end


    % Defines a single enhancer network to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF: production rate of TF
    % simR: production rate of mRNA
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN = singleEnhancer2(params, extra, prodTF, simR, prodTFstar)

 
        % number of TFs per burst
        n1 = params.tf2;
        n3 = params.tfstar;

        % E = enhancer 1, T = TF, C = complex, R=RNA

        %  E + T  <-> C -> C + R
        %  E + T*  <-> C1* -> C1* + R
        %  R -> 0
        %  n*T <- 0 <- T
        %  n3*T* <- 0 <- T*
        
        reaction=[...
            %   E  C  R  T   T*     E  C  R  T   T*  
                1  0  0  1   1      0  1  0  0   0 
                0  1  0  0   0      1  0  0  1   1 
                0  1  0  0   0      0  1  1  0   0    
                0  0  1  0   0      0  0  0  0   0 
                0  0  0  0   0      0  0  0  n1  0  
                0  0  0  1   0      0  0  0  0   0    
                0  0  0  0   0      0  0  0  0   n3  
                0  0  0  0   1      0  0  0  0   0  
        ];

        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum,1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum,speciesNum + 1:2*speciesNum); % target complexes


        CRN.params = [params.k_on2 params.k_off2 simR extra.alpha prodTF params.degTF2 ...
                     prodTFstar params.degTFstar];  

        %           E   C   R   T  T*
        CRN.init = [1   0   0   0  0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

    end


    % Defines a double enhancer network 
    % to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF: production rate of TF
    % simR1: production rate of mRNA for complex 1
    % simR2: production rate of mRNA for complex 2
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN = duplicatedEnhancer(params, extra, prodTF, simR1, simR2, prodTFstar)


        % number of TFs per burst
        n1 = params.tf;
        n3 = params.tfstar;

        % E1= enhancer 1, E2 = enhancer 2, C=complex, R=RNA

        %  E1 + T1  -> C1 -> C1 + R
        %  E2 + T1  -> C2 -> C2 + R
        %  R -> 0
        %  n*T1 <- 0 <- T1

        % it is convenient to leave T2 so that locations of R1
        % and R2 in the  and shadow case are the same
        
        reaction=[...
          % E1 E2 C1 C2 R1 T1 T2 R2 T*      E1 E2 C1 C2 R1 T1 T2  R2 T*  
            1  0  0  0  0  1  0  0  1       0  0  1  0  0  0   0  0  0     
            0  0  1  0  0  0  0  0  0       1  0  0  0  0  1   0  0  1     
            0  0  1  0  0  0  0  0  0       0  0  1  0  1  0   0  0  0     
            0  1  0  0  0  1  0  0  1       0  0  0  1  0  0   0  0  0     
            0  0  0  1  0  0  0  0  0       0  1  0  0  0  1   0  0  1            
            0  0  0  1  0  0  0  0  0       0  0  0  1  0  0   0  1  0    
            0  0  0  0  1  0  0  0  0       0  0  0  0  0  0   0  0  0     
            0  0  0  0  0  0  0  0  0       0  0  0  0  0  n1  0  0  0     
            0  0  0  0  0  1  0  0  0       0  0  0  0  0  0   0  0  0     
            0  0  0  0  0  0  0  1  0       0  0  0  0  0  0   0  0  0       
            0  0  0  0  0  0  0  0  0       0  0  0  0  0  0   0  0  n3     
            0  0  0  0  0  0  0  0  1       0  0  0  0  0  0   0  0  0      
            
        ];

        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum,1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum,speciesNum + 1:2 * speciesNum); % target complexes

        % default parameters - any vector with reactionsNum entries

        CRN.params=[params.k_on params.k_off simR1 params.k_on params.k_off ...
          simR1 extra.alpha prodTF params.degTF extra.alpha prodTFstar params.degTFstar];  

        %         E1 E2 C1 C2 R1 T1 T2 R2 T* 
        CRN.init=[1  1  0  0  0  0  0  0  0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span
        
    end

    function CRN = duplicatedEnhancer2(params, extra, prodTF, simR1, simR2, prodTFstar)

        % number of TFs per burst
        n1 = params.tf2;
        n3 = params.tfstar;


        % E1= enhancer 1, E2 = enhancer 2, C=complex, R=RNA

        %  E1 + T1  -> C1 -> C1 + R
        %  E2 + T1  -> C2 -> C2 + R
        %  R -> 0
        %  n*T1 <- 0 <- T1

        % it is convenient to leave T2 so that locations of R1
        % and R2 in the  and shadow case are the same
        
        reaction=[...
          % E1 E2 C1 C2 R1 T1 T2 R2 T*      E1 E2 C1 C2 R1 T1 T2  R2 T*  
            1  0  0  0  0  1  0  0  1       0  0  1  0  0  0   0  0  0     
            0  0  1  0  0  0  0  0  0       1  0  0  0  0  1   0  0  1     
            0  0  1  0  0  0  0  0  0       0  0  1  0  1  0   0  0  0     
            0  1  0  0  0  1  0  0  1       0  0  0  1  0  0   0  0  0     
            0  0  0  1  0  0  0  0  0       0  1  0  0  0  1   0  0  1            
            0  0  0  1  0  0  0  0  0       0  0  0  1  0  0   0  1  0    
            0  0  0  0  1  0  0  0  0       0  0  0  0  0  0   0  0  0     
            0  0  0  0  0  0  0  0  0       0  0  0  0  0  n1  0  0  0     
            0  0  0  0  0  1  0  0  0       0  0  0  0  0  0   0  0  0     
            0  0  0  0  0  0  0  1  0       0  0  0  0  0  0   0  0  0       
            0  0  0  0  0  0  0  0  0       0  0  0  0  0  0   0  0  n3     
            0  0  0  0  0  0  0  0  1       0  0  0  0  0  0   0  0  0      
            
        ];


        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum,1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum,speciesNum + 1:2 * speciesNum); % target complexes

        % default parameters - any vector with reactionsNum entries

        CRN.params=[params.k_on params.k_off simR1 params.k_on params.k_off ...
          simR1 extra.alpha prodTF params.degTF extra.alpha prodTFstar params.degTFstar];  

        %         E1 E2 C1 C2 R1 T1 T2 R2 T* 
        CRN.init=[1  1  0  0  0  0  0  0  0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span
        
    end


    % Defines a double enhancer network 
    % to be simulated.
    % --Inputs:
    % params: network parameters
    % extra: struct containing mRNA and length of simulation
    % prodTF1: production rate of TF 1
    % prodTF2: production rate of TF 2
    % simR1: production rate of mRNA for complex 1
    % simR2: production rate of mRNA for complex 2
    % --Outputs:
    % CRN: a struct containing relevant network information.
    % Alvaro Fletcher
    function CRN  = shadowPair(params, extra, prodTF, prodTF2, simR1, simR2, prodTFstar)

        % number of TFs per burst
        n1 = params.tf;
        n2 = params.tf2;
        n3 = params.tfstar;

        % E1= enhancer 1, E2 = enhancer 2, C=complex, R=RNA

        %  E1 + T1  <> C1 -> C1 + R1
        %  E2 + T1  <> C2 -> C2 + R2
        %  E1 + T*  <> C1* -> C1* + R1
        %  E2 + T*  <> C2* -> C2* + R2
        %  R -> 0
        %  n1*T1 <- 0 <- T1
        %  n3*T* <- 0 <- T*

    
        reaction=[...
          % E1 E2 C1 C2 R1 T1 T2 R2 T*      E1 E2 C1 C2 R1 T1 T2  R2 T*  
            1  0  0  0  0  1  0  0  1       0  0  1  0  0  0   0  0  0     
            0  0  1  0  0  0  0  0  0       1  0  0  0  0  1   0  0  1     
            0  0  1  0  0  0  0  0  0       0  0  1  0  1  0   0  0  0     
            0  1  0  0  0  0  1  0  1       0  0  0  1  0  0   0  0  0     
            0  0  0  1  0  0  0  0  0       0  1  0  0  0  0   1  0  1            
            0  0  0  1  0  0  0  0  0       0  0  0  1  0  0   0  1  0    
            0  0  0  0  1  0  0  0  0       0  0  0  0  0  0   0  0  0     
            0  0  0  0  0  0  0  0  0       0  0  0  0  0  n1  0  0  0     
            0  0  0  0  0  1  0  0  0       0  0  0  0  0  0   0  0  0     
            0  0  0  0  0  0  0  1  0       0  0  0  0  0  0   0  0  0       
            0  0  0  0  0  0  0  0  0       0  0  0  0  0  0   0  0  n3     
            0  0  0  0  0  0  0  0  1       0  0  0  0  0  0   0  0  0     
            0  0  0  0  0  0  0  0  0       0  0  0  0  0  0   n2 0  0     
            0  0  0  0  0  0  1  0  0       0  0  0  0  0  0   0  0  0    
            
        ];



        speciesNum=length(reaction(1,:))/2;  % number of speciesNum
        CRN.speciesNum=speciesNum;

        reactionsNum=length(reaction(:,1));  %number of reactions
        CRN.reactionsNum=reactionsNum;

        CRN.source=reaction(1:reactionsNum, 1:speciesNum);  % source complexes

        CRN.target=reaction(1:reactionsNum, speciesNum+1:2*speciesNum); % target complexes

        CRN.params=[params.k_on params.k_off simR1 params.k_on2 params.k_off2 ...
          simR2 extra.alpha prodTF params.degTF extra.alpha prodTFstar params.degTFstar prodTF2 params.degTF2];  

        %         E1 E2 C1 C2 R1 T1 T2 R2 T*
        CRN.init=[1  1  0  0  0  0  0  0  0];  % default initial conditions

        CRN.tspan=[0,extra.endSpan];  % default time span

    end

   
    % This function is given aGE CRN data structure (see samplereaction files), 
    % and optionally an initial state.  
    % It carries out the Gillespie algorithm for a fixed number of steps and 
    % produces a plot.  
    % This file has been corrected to have the right reaction propensities when
    % more than one identical molecule is involved.  For instance, if 2A->B,
    % and A=10, then the propensity is A*(A-1)=90, not A^2=100.  
    % --Inputs:
    % CRN: the network to be simulated
    % maxiterations: the maximum number of iterations
    % -- Outputs: 
    % expN: a vector containing simulated mRNA counts at 30 second intervals
    %  German Enciso
    % Last updated: 2009-11-14
    %
    % With some modifications by Alvaro Fletcher on 07/23/2018:
    % 1. Fixed a bug with the limiting time span of the simulation. 
    % 2. This version of crndiscrete also generates data as is recorded
    % experimentally, in this case, every 30 second intervals. 
    function expN = crndiscrete(CRN,maxiterations)

        nGE=CRN.speciesNum;   % number of species
        mGE=CRN.reactionsNum;   % number of reactions
        source=CRN.source ;  % source is for source
        target=CRN.target; % target is for target
        paramsGE=CRN.params ; %  all the parameters in the system
        initGE=CRN.init;
        tspan=CRN.tspan;

        %default maximum number of iterations
        if nargin==1    
            maxiterations=500*nGE*mGE; 
        end

        % make sure initGE is aGE column vector
        temp = size(initGE);

        if temp(1)==1 

            initGE=initGE'; 

        end

        N(1,:) = initGE;   % the initial state;
        tGE(1) = tspan(1);  %  the initial time;


        % rachel records concentrations of different species
        % at specific times so we have uniform recordings for all
        % simulations

        expN(1,:) = N(1,:);
        myTime = 2;

        for iGE = 1:maxiterations-1

            xGE = N(iGE,:);  % vector of concentrations
      
            % going over every reaction
            for k = 1:mGE

                pow=source(k,:);   % vector of powers
                reactionrate(k) = paramsGE(k) * prod( xGE.^pow );

                % correct reaction rate in case that any of the powers is larger than one    
                % going over every species
                for el = 1:nGE

                    if xGE(el) ~= 0

                        for alpha = 2:source(k,el)

                              reactionrate = reactionrate * (xGE(el) - alpha + 1)/xGE(el);

                        end

                    end

                end   

            end

            rTGE = sum(reactionrate); % Total rate
            u_time = rand;         % Uniformly distributed number

            % time(iGE) is the lapse of time before the next reaction
            time(iGE) = -log(u_time)/rTGE;

            timeHolderGE = tGE(iGE) + time(iGE);
                
            % preventing time from going above tspan(2)
            if timeHolderGE > tspan(2)  
                

                 for mycount = floor(tGE(iGE)) + 0.5:0.5:tspan(2) 
                     
                    % setting the remainder of expN to the last value recorded 
                    expN(myTime,:) = N(iGE,:);      
                    myTime = myTime + 1;


                 end   

                 break

            else 
                    
                  tGE(iGE+1) = timeHolderGE;

                  % we choose 0.5 to resemble the experimental data from Rachel

                  for mycount = floor(tGE(iGE)) + 0.5:0.5:floor(tGE(iGE+1))
                      
                        expN(myTime,:) = N(iGE,:);      
                        myTime = myTime + 1;

                  end

            end



            % choose a number numGE defining which reaction took place
            uGE=rTGE*rand;
            aGE=0;

            % mGE is the number of reactions
            for kk=1:mGE

                aGE=aGE+reactionrate(kk);

                if uGE <= aGE 
                        
                      % we choose the kk^th reaction to take place
                      numGE=kk; 
                      break  

                end

            end
            
            N(iGE+1,:) = N(iGE,:) + target(numGE,:) - source(numGE,:);

        end  

    end

    % Calulates burst properties (among other things) 
    % from experimental or simulated data.
    % --Inputs:
    % SpotDiff: a struct containing simulated data for all bins
    % ElapsedTime: a vector containing the time at which each point  
    % was recorded in the traces. 
    % --Outputs:
    % BurstProperties: a struct containing burst properties
    % Rachel Waymack. 
    %
    % Modifications done by Alvaro Fletcher on 07/23/2018:
    % 1. Made the script into a function.
    % 2. Deleted calculation for SpotTwo (allele correlation 
    % is calculated directly from the simulated data)
    % 3. RNAScale and nc14 parameters set to 1 (model operates in transcript count). 
    % 
    function [BurstProperties] = modSlopeBurstCallingRW(SpotDiffRW,ElapsedTime)

        % MOD: added RNA scale to ON threshold and OFF threshold
        RNAscale = 1;
        % nc14 length just set it to 1
        nc14 = 1;
        
        %Slope thresholds
        ONThreshold = RNAscale;   %6/20/18 thinking ~377AU=1 polymerase so 1500 ~ 5ish polymerases
        OFFThreshold = -RNAscale; 

        for n=1:length(SpotDiffRW)
            if length(SpotDiffRW(n).SpotOne)==length(ElapsedTime)
               %try limiting frames to nc14 to reduce extreme smoothing
               % MOD: got rid of the nc14 variable here
               % did I put the 100 here?
               if length(ElapsedTime) >= 100
               SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:(nc14+100)),SpotDiffRW(n).SpotOne(nc14:(nc14+100)),0.1,'lowess');
               else
                   ShortTrace='y';
                   SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:end),SpotDiffRW(n).SpotOne(nc14:end),0.1,'lowess');
               end
               %SmoothParticles(n).Smoothed=smooth(ElapsedTime(nc14:end),SpotDiffRW(n).SpotOne(nc14:end),0.1,'lowess');
            %fluorescence should never be negative but sometimes with lots of 0s
            %looks like the smoothing f(x) makes it so
               SmoothParticles(n).Smoothed(SmoothParticles(n).Smoothed<0)=0;
               if exist('ShortTrace') & ShortTrace=='y'
               SmoothParticles(n).RawSlope=diff(SmoothParticles(n).Smoothed)./([diff(ElapsedTime(nc14:end))]');
               else
               SmoothParticles(n).RawSlope=diff(SmoothParticles(n).Smoothed)./([diff(ElapsedTime(nc14:(nc14+100)))]');
               end
        %Slope thresholds
               %OFFThresholdSlope=dx<=Q;
        % AboveLine=[SmoothParticles(ss).Smoothed];
        % AboveLine(~ONThresholdSlope)=nan;
        % AboveLine(OFFThresholdSlope)=nan;
        % BelowLine=[SmoothParticles(ss).Smoothed];
        % BelowLine(ONThresholdSlope)=nan;
        % BelowLine(~OFFThresholdSlope)=nan;
        StartFramePool=[];
        StartSlopePool=[];
        EndSlopePool=[];
        Starts=[find(SmoothParticles(n).RawSlope>=ONThreshold)];
        Starts=min(Starts);
        StartFramePool=[StartFramePool,Starts];
        EndPts=[find(SmoothParticles(n).RawSlope<=OFFThreshold)];
        EndPts=min(EndPts);
        EndFramePool=[1000,EndPts];
        if ~isempty(EndPts) & (~isempty(Starts))
        while EndFramePool(end) ~= EndFramePool(end-1)
            for ss=EndFramePool(end):length(SmoothParticles(n).RawSlope)
                if SmoothParticles(n).RawSlope(ss) >=ONThreshold
                    %if StartFramePool(end)~=ss
                    StartFramePool=[StartFramePool,ss];
                    StartSlopePool=[StartSlopePool, SmoothParticles(n).RawSlope(ss)];
                    break
                    %end
                end
            end
            for ss=StartFramePool(end):(length(SmoothParticles(n).RawSlope))
                if SmoothParticles(n).RawSlope(ss) <=OFFThreshold
                    EndFramePool=[EndFramePool, ss];
                    EndSlopePool=[EndSlopePool,SmoothParticles(n).RawSlope(ss)];
                    break
                end
            end
            if sum(SmoothParticles(n).RawSlope(EndFramePool(end):end)>=ONThreshold)==0
        %         EndFramePool(end+1)=EndFramePool(end);
        %         StartFramePool(end+1)=StartFramePool(end);
                PotentialPeaks(n).Starts=[StartFramePool];
                PotentialPeaks(n).Ends=[EndFramePool(2:end)];
                break
            elseif sum(SmoothParticles(n).RawSlope(StartFramePool(end):end)<=OFFThreshold)==0
        %         EndFramePool(end+1)=EndFramePool(end);
        %         StartFramePool(end+1)=StartFramePool(end);
                PotentialPeaks(n).Starts=[StartFramePool];

                PotentialPeaks(n).Ends=[EndFramePool(2:end)];
                break
            end
        end
        if length(PotentialPeaks(n).Ends)>1 & (PotentialPeaks(n).Ends(end)==PotentialPeaks(n).Ends(end-1))
            PotentialPeaks(n).Ends=[PotentialPeaks(n).Ends(1:end-1)];
        end
        if length(PotentialPeaks(n).Starts)>1 & (PotentialPeaks(n).Starts(end)==PotentialPeaks(n).Starts(end-1))
            PotentialPeaks(n).Starts=[PotentialPeaks(n).Starts(1:end-1)];
        end
        BurstProperties(n).SmoothTrace=SmoothParticles(n).Smoothed;
        %Get rid of repeated 1st and 2nd ON frames in case an OFF came before ON in
        %original trace - added 7/27/18
        if (length(PotentialPeaks(n).Starts)>=2)& (PotentialPeaks(n).Starts(1)==PotentialPeaks(n).Starts(2))
            PotentialPeaks(n).Starts=[PotentialPeaks(n).Starts(2:end)];
        end
        % If a trace starts with an OFF frame before an ON frame, need to ignore it
        % 7/27/18
        if PotentialPeaks(n).Ends(1) < PotentialPeaks(n).Starts(1)
            if length(PotentialPeaks(n).Ends)==1
                PotentialPeaks(n).Ends=nan;
                PotentialPeaks(n).Starts=nan;
            else
            PotentialPeaks(n).Ends=[PotentialPeaks(n).Ends(2:end)];
            end
        end
        if (isnan(PotentialPeaks(n).Ends)) & (isnan(PotentialPeaks(n).Starts))
            continue
        end
        %If start frames > end frames, ignore the last start since we can't say how
        %long the burst lasts
        if length(PotentialPeaks(n).Starts) > (length(PotentialPeaks(n).Ends))
            BurstProperties(n).Duration=ElapsedTime(PotentialPeaks(n).Ends)-ElapsedTime(PotentialPeaks(n).Starts(1:end-1));
            BurstProperties(n).Duration=BurstProperties(n).Duration(BurstProperties(n).Duration>=0);
            BurstProperties(n).InterBurst=ElapsedTime(PotentialPeaks(n).Starts(2:end-1))-ElapsedTime(PotentialPeaks(n).Ends(1:end-1));
            BurstProperties(n).ONFrames=[PotentialPeaks(n).Starts(1:end-1)];
        else
        BurstProperties(n).Duration=ElapsedTime(PotentialPeaks(n).Ends)-ElapsedTime(PotentialPeaks(n).Starts);
        BurstProperties(n).Duration=BurstProperties(n).Duration(BurstProperties(n).Duration>=0);
        BurstProperties(n).InterBurst=ElapsedTime(PotentialPeaks(n).Starts(2:end))-ElapsedTime(PotentialPeaks(n).Ends(1:end-1));
        BurstProperties(n).ONFrames=[PotentialPeaks(n).Starts];
        end

        BurstProperties(n).OFFFrames=[PotentialPeaks(n).Ends];

        if PotentialPeaks(n).Ends(1)==1
            BurstProperties(n).BurstAmplitude(1)=0;
        else
            BurstProperties(n).BurstAmplitude(1)=max(SmoothParticles(n).Smoothed([1:PotentialPeaks(n).Ends(1)]));
        end

            for pp=2:length(PotentialPeaks(n).Ends)
            BurstProperties(n).BurstAmplitude(pp)=max(SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends(pp-1):PotentialPeaks(n).Ends(pp)]));%[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends]-1)]';
            BurstProperties(n).BurstAmplitude(pp)=(BurstProperties(n).BurstAmplitude(pp))-(SmoothParticles(n).Smoothed(PotentialPeaks(n).Starts(pp)));
            end

        BurstProperties(n).BurstAmplitude=BurstProperties(n).BurstAmplitude(BurstProperties(n).BurstAmplitude~=0); %get rid of 1st 0 bc above started at pp=2
        
        for pp=1:(length(BurstProperties(n).ONFrames)-1)
        BurstProperties(n).BurstSize(pp)=(trapz([BurstProperties(n).ONFrames(pp):BurstProperties(n).ONFrames(pp+1)],BurstProperties(n).SmoothTrace(BurstProperties(n).ONFrames(pp):BurstProperties(n).ONFrames(pp+1))));
        end
        if isfield(BurstProperties, 'BurstSize')
        BurstProperties(n).BurstSize(end+1)=(trapz([BurstProperties(n).ONFrames(end):BurstProperties(n).OFFFrames(end)], BurstProperties(n).SmoothTrace(BurstProperties(n).ONFrames(end):BurstProperties(n).OFFFrames(end))));
        end
        % if PotentialPeaks(n).Ends(1)==1
        %     BurstProperties(n).BurstAmplitude=[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends(2:end)]-1)]';
        % else
        % BurstProperties(n).BurstAmplitude=[SmoothParticles(n).Smoothed([PotentialPeaks(n).Ends]-1)]';
        % end
       % temp=find([CompiledParticles.Nucleus]==SpotDiffRW(n).Nucleus);
       % BurstProperties(n).TotalmRNAError=CompiledParticles(temp(1)).TotalmRNAError;

        BurstProperties(n).FirstTimeOn=ElapsedTime(PotentialPeaks(n).Starts(1));
        BurstProperties(n).NBursts=length(BurstProperties(n).BurstAmplitude);
        %Frequency set btwn time of first ON and 50min into nc14 
        if ((length(ElapsedTime))-nc14)<100
            BurstProperties(n).Frequency=(BurstProperties(n).NBursts)/(ElapsedTime(end)-ElapsedTime(PotentialPeaks(n).Starts(1)));
        else
        BurstProperties(n).Frequency=(BurstProperties(n).NBursts)/(ElapsedTime(nc14+100)-ElapsedTime(PotentialPeaks(n).Starts(1)));
        end
        BurstProperties(n).FrequencyAct=(BurstProperties(n).NBursts)/(ElapsedTime(PotentialPeaks(n).Ends(end))-ElapsedTime(PotentialPeaks(n).Starts(1)));
        BurstProperties(n).APBin=SpotDiffRW(n).APBin;
       % BurstProperties(n).Nucleus=SpotDiffRW(n).Nucleus;
        %BurstProperties(n).Interburst=
      %  BurstProperties(n).TotalmRNA=SpotDiffRW(n).TotalmRNAOne;
        BurstProperties(n).FractON=sum(BurstProperties(n).Duration)/(length(ElapsedTime(nc14:end))); 
        BurstProperties(n).TotalElapsedTime=ElapsedTime(end)-ElapsedTime(nc14);

        else 
            
            BurstProperties(n).APBin = SpotDiffRW(n).APBin;
            
        end
            end
        end
        
        for bb=1:length(BurstProperties)
            if isempty(BurstProperties(bb).APBin)
                BurstProperties(bb).APBin=nan;
                % BurstProperties(bb).Nucleus=nan;
                % MOD: added uniqueID
                BurstProperties(bb).uniqueID = nan;
            end

        end
        clear ShortTrace
        %save('BurstPropertiesSlope','BurstProperties');

    end

end