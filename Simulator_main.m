% Makes function calls to read annealing files and run simulations 
% with specified parameter sets. Then, it automatically saves figures 
% to folders with the name of the chosen strategies (these folders
% must already be present in the same directory).
% Last updated: 29/07/19
% Alvaro Fletcher

% reserved figure nums:1,3 (boxplots and SSE scatter),
% 11,12 (?)
% 50,51 (exp CV and corr),
% 100,101,102, (sim corr,CV, and props),
%  201 - 203,(SSE histograms) 
% 300-~314 (bootstrapping plots), 
% 400-~406 (pairwise correlation and pca)
% 500 - 503 (quality of fitting for temperatures)
tic
clear all; close all; %#ok<CLALL>

initStrat = 1;
% the top maxStrat strategies will be stored in directories
% that should be created beforehand
maxStrat = 10;

% run simulations with this set of parameters
simBool = true;
% overlap experimental data
experimental = true;
% when true, readerEnergies plots will be ignored
ignoreSave = false;

for topIndex = initStrat:maxStrat

	stratDesired = topIndex; 

	% we only need to save energy plots for one run, after that we only 
	% care about the parameter values for the 2nd best, 3rd best, etc etc
	if topIndex > 1 
		ignoreSave = true;
    end

	strats = {'without','cold', 'colddiff', ...
   			'fixed', 'hot', 'hotdiff', 'sifd', 'sifdc','sifds', ...
			 'sifdr', 'sifdcr', 'sifdrs','sifdcn', 'withoutn'}; 

	% string containing the word results and the rank of the 
	% current parameter set
	resultsFolder = ['results' num2str(stratDesired)];

	for i = strats
	    
		% readerEnergies indices
		energyFigIndex = 1;
		boxFigIndex = 3;
		histFigIndex1 = 201;
		histFigIndex2 = 202;
		histFigIndex3 = 203;
		pcIndex = 400;
		temperIndex1 = 501;
		temperIndex2 = 502;
		temperIndex3 = 503;

		% plotter indices
		corrFigIndex = 100;
		cvFigIndex = 101;
		propsFigIndex = 102;
		expCVFigIndex = 50;
		expCorrFigIndex = 51;
		% note: the experimental burst properties
		% index is determined by fig counter

		strat = char(i);

		doubBool = false;
		randBool = false;
		competition = false;
		boot = false;
		ignoreSaveSim = false;

		switch strat
	        case 'sifd'
	  			stratr = 'sifdr';
	  			randBool = true;
	  			doubBool = true;
	  		% at the moment we have only done fittings of double and boot to sifdc
	        case 'sifdc'
				stratr = 'sifdcr';
				randBool = true;
				doubBool = true; 
				competition = true;
				boot = true;
				stratDoub = 'double';
			% no competition
			 case 'sifdcn'
			 	strat = 'sifdc';
				randBool = true;
				doubBool = true; 
				% so it doesn't overwrite
				% plots from sifdc
				ignoreSave = true;
				ignoreSaveSim = true;
			case 'withoutn'
			 	strat = 'without';
				doubBool = true; 
				% so it doesn't overwrite
				% plots from sifdc
				ignoreSave = true;
				ignoreSaveSim = true;
	        case 'sifds'
	        	stratr = 'sifdrs';
				randBool = true;
				doubBool = true;
			case 'sifdr'
				doubBool = true;
			case 'sifdcr'
				doubBool = true; 
			case 'sifdrs'
				doubBool = true;	
			case 'without'
				doubBool = true; 
				competition = true;		
				stratDoub = 'doubles';
	    end

		% generating figures with reader energies
		[bestVals,bestParams] = readerEnergies(strat,stratDesired);

		if ~ignoreSave
			% saving boxplot
			saveas(figure(boxFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\box' num2str(topIndex) '.jpg']);
			close(figure(boxFigIndex))

			saveas(figure(temperIndex1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\temper1' num2str(topIndex) '.jpg']);
			close(figure(temperIndex1))

			saveas(figure(temperIndex2),[pwd '\' 'results' '\' strat '\' resultsFolder  '\temper2' num2str(topIndex) '.jpg']);
			close(figure(temperIndex2))

			% saving correlation 
			saveas(figure(pcIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\' 'paircorr1' num2str(topIndex) '.jpg']);
			close(figure(pcIndex))
			pcIndex = pcIndex + 1;
			saveas(figure(pcIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\' 'paircorr2' num2str(topIndex) '.jpg']);
			close(figure(pcIndex))
			% saving pca
			pcIndex = pcIndex + 1;
			saveas(figure(pcIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\' 'pca1' num2str(topIndex) '.jpg']);
			close(figure(pcIndex))
			pcIndex = pcIndex + 1;
			saveas(figure(pcIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\' 'pca2' num2str(topIndex) '.jpg']);
			close(figure(pcIndex))
			% saving tables
			pcIndex = pcIndex + 1;
			saveas(figure(pcIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\' 'table1' num2str(topIndex) '.jpg']);
			close(figure(pcIndex))
			pcIndex = pcIndex + 1;
			saveas(figure(pcIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\' 'table2' num2str(topIndex) '.jpg']);
			close(figure(pcIndex))

		else 
			close all
		end

		if randBool

			readerEnergies(stratr,stratDesired);
			if ~ignoreSave
				% adding legends to energy plots 
				figure(energyFigIndex)
				legend('proximal','distal','proximal random', 'distal random')

				% adding legends to histograms
				figure(histFigIndex1)
				legend('Fitted IC','Random IC')
				title('Histogram of SSEs for proximal')  
				figure(histFigIndex2)
				legend('Fitted IC','Random IC')
				title('Histogram of SSEs for distal') 
				 
				% saving figures
				saveas(figure(histFigIndex1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\hist1' num2str(topIndex) '.jpg']);
				saveas(figure(histFigIndex2),[pwd '\' 'results' '\' strat '\' resultsFolder  '\hist2' num2str(topIndex) '.jpg']);
				saveas(figure(energyFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\energies' num2str(topIndex) '.jpg']);
			else 
				close all
			end


		else
			
			% adding legends to energy plots 
			figure(energyFigIndex)
			legend('proximal','distal')

			% adding legends to histograms
			figure(histFigIndex1)
			legend('Fitted IC')
			title('Histogram of SSEs for proximal') 
			figure(histFigIndex2)
			legend('Fitted IC')
			title('Histogram of SSEs for distal') 

			if ~ignoreSave		
				% saving figures
				saveas(figure(histFigIndex1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\hist1' num2str(topIndex) '.jpg']);
				saveas(figure(histFigIndex2),[pwd '\' 'results' '\' strat '\' resultsFolder  '\hist2' num2str(topIndex) '.jpg']);
				saveas(figure(energyFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\energies' num2str(topIndex) '.jpg']);
			else 
				close all
			end

		end

		close all

		% section for running simulations with best values from a strat
		if simBool

			% placeholder
			model3 = 0;

			without = false;
			temperature = false;
			fixedTF = false;

			switch strat
		        case 'fixed'
		        	fixedTF = true;
		        	model1 = 'proxfixed';
					model2 = 'distfixed';
					model3 = 'shadowfixed';
		        case 'without'
		         	model1 = 'proxnequal1';
					model2 = 'distnequal1';
					model3 = 'shadownequal1';
					without = true;
		        case 'sifd'
		         	model1 = 'proximal';
					model2 = 'edist';
		        case {'sifdc','sifdcn'}
		         	model1 = 'proximal';
					model2 = 'distal';
					model3 = 'shadowcorr';
		        case 'sifdr'
		         	model1 = 'proximal';
					model2 = 'edist';
		        case 'sifdcr'
		         	model1 = 'proximal';
					model2 = 'distal';
					model3 = 'shadowcorr';
		        case 'hot'
		        	temperature = true;
		         	model1 = 'proxhot';
					model2 = 'disthot';
		        case 'hotdiff' 
		        	temperature = true;
		        	model1 = 'proxhotdiff';
		        	model2 = 'disthotdiff';
		        case 'cold'
		        	temperature = true;
		         	model1 = 'proxcold';
					model2 = 'distcold';
		        case 'colddiff'
		        	temperature = true;
		        	model1 = 'proxcolddiff';
		        	model2 = 'distcolddiff';
		        case 'sifds'
		         	model1 = 'proximal';
					model2 = 'distal';
					model3 = 'shadowcorr';
		        case 'sifdrs'
		         	model1 = 'proximal';
					model2 = 'distal';
					model3 = 'shadowcorr';
		    end

			% single cases
			figCounter = 10;

			% for future storing of og values
			koff1 = 2;
			kon1 = 3;
		    koff2 = 8;
			kon2 = 9;

			if without || fixedTF
				proxBest = bestVals(1:5);
				distBest = bestVals(6:10);
				% n = 1 for without bursts
				proxBest = [proxBest(1:3),1,proxBest(4:5)];	
				distBest = [distBest(1:3),1,distBest(4:5)];	
			else
				proxBest = bestVals(1:6);
				distBest = bestVals(7:12);	
			end
			% running simulations for single proximal,distal and computing
			% shadow allele correlations
			figCounter = plotter(model1, experimental,figCounter,proxBest, proxBest);
			if ~ignoreSaveSim
				saveas(figure(figCounter - 1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsExp1' num2str(topIndex) '.jpg']);
			end

			figCounter = plotter(model2, experimental,figCounter,distBest,distBest);
		    
		    if ~ignoreSaveSim
				saveas(figure(figCounter - 1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsExp2' num2str(topIndex) '.jpg']);
		    end
			
			if model3 ~= 0
				figCounter = plotter(model3, experimental,figCounter,proxBest,distBest);
			end

			% saving figures
			if ~ignoreSaveSim
				saveas(figure(propsFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsSim' num2str(topIndex) '.jpg']);
				saveas(figure(cvFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\cvSim' num2str(topIndex) '.jpg']);
			end

			if ~temperature
			    saveas(figure(corrFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\corrSim' num2str(topIndex) '.jpg']);
			    saveas(figure(expCorrFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\corrExp' num2str(topIndex) '.jpg']);
			    saveas(figure(expCVFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\cvExp' num2str(topIndex) '.jpg']);
			end

			% currently this is only done for the best strat of sifdc
			if doubBool && competition 

				pcIndex = 400;

				% storing the koff and kon values for distal and proximal
				% order out of readerEnergies : koff1, koff2, kon1, kon2
				og.ogKoff1 = bestVals(koff1);
				og.ogKoff2 = bestVals(koff2);
				og.ogKon1 = bestVals(kon1); 
				og.ogKon2 = bestVals(kon2); 

				% double cases
				[bestValsDoub,bestParamsDoub] = readerEnergies(stratDoub,stratDesired,og);

				if ~ignoreSave		

					saveas(figure(boxFigIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\box' strat '' num2str(topIndex) '.jpg']);
					close(figure(boxFigIndex))

					saveas(figure(temperIndex1),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\temper1' num2str(topIndex) '.jpg']);
					close(figure(temperIndex1))

					saveas(figure(temperIndex2),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\temper2' num2str(topIndex) '.jpg']);
					close(figure(temperIndex2))

					saveas(figure(temperIndex3),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\temper3' num2str(topIndex) '.jpg']);
					close(figure(temperIndex3))

					% saving correlation 
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'paircorr1' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))
					pcIndex = pcIndex + 1;
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'paircorr2' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))
					pcIndex = pcIndex + 1;
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'paircorr3' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))
					pcIndex = pcIndex + 1;
					% saving pca
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'pca1' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))
					pcIndex = pcIndex + 1;
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'pca2' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))
					pcIndex = pcIndex + 1;
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'pca3' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))
					% saving tables
					pcIndex = pcIndex + 1;
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'table1' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))
					pcIndex = pcIndex + 1;
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'table2' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))
					pcIndex = pcIndex + 1;
					saveas(figure(pcIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\' 'table3' num2str(topIndex) '.jpg']);
					close(figure(pcIndex))

					figure(energyFigIndex)
					legend('proximal','distal','shadow')
			 

					figure(histFigIndex1)
			 		title('Histogram of SSEs for proximal 2x')
			 		figure(histFigIndex2)
			 		title('Histogram of SSEs for distal 2x') 
			 		figure(histFigIndex3)
			 		title('Histogram of SSEs for shadow') 
					% saving figures
					saveas(figure(histFigIndex1),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\hist1' strat '' num2str(topIndex) '.jpg']);
		 
					saveas(figure(histFigIndex2),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\hist2' strat '' num2str(topIndex) '.jpg']);
				
					saveas(figure(histFigIndex3),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\hist3' strat '' num2str(topIndex) '.jpg']);

					saveas(figure(energyFigIndex),[pwd '\' 'results' '\' stratDoub '\' resultsFolder  '\energies' strat '' num2str(topIndex) '.jpg']);

				end
				% splitting bestvalues for updating 
				% values of kon and koff
				proxBest1 = proxBest;
				proxBest2 = proxBest;

				distBest1 = distBest;
				distBest2 = distBest;

				% both use the 1 index since we have split 
				% bestVals into two sets of 6 and TF params
				% remain the same as in the single enhancer cases
				% (double annealing only fits kons and koffs)
				proxBest1(koff1) = bestValsDoub(1);
				proxBest2(koff1) = bestValsDoub(2);

				proxBest1(kon1) = bestValsDoub(3);
				proxBest2(kon1) = bestValsDoub(4);

				distBest1(koff1) = bestValsDoub(5);
				distBest2(koff1) = bestValsDoub(6);

				distBest1(kon1) = bestValsDoub(7);
				distBest2(kon1) = bestValsDoub(8);

				figCounter = plotter('prox2x', experimental,figCounter,proxBest1, proxBest2);
				saveas(figure(figCounter - 1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsDoubExp1' strat '' num2str(topIndex) '.jpg']);

				figCounter = plotter('dist2x', experimental,figCounter,distBest1,distBest2);
				saveas(figure(figCounter - 1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsDoubExp2' strat '' num2str(topIndex) '.jpg']);

				% updating from the shadow annealing

				proxBest(koff1) = bestValsDoub(9);
				distBest(koff2) = bestValsDoub(10);

				proxBest(kon1) = bestValsDoub(11);
				distBest(kon2) = bestValsDoub(12);

				figCounter = plotter('shadow', experimental,figCounter,proxBest,distBest);
				saveas(figure(figCounter - 1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsDoubExp3' strat '' num2str(topIndex) '.jpg']);

				% saving figures
				saveas(figure(cvFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\cvSimDoub' strat '' num2str(topIndex) '.jpg']);
			    saveas(figure(propsFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsSimDoub' strat '' num2str(topIndex) '.jpg']);
			    saveas(figure(expCVFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\cvExpDoub' strat '' num2str(topIndex) '.jpg']);


			elseif doubBool

				figCounter = plotter('prox2x', experimental,figCounter,proxBest, proxBest);
				saveas(figure(figCounter - 1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsDoubExp1Ind' strat '' num2str(topIndex) '.jpg']);

				figCounter = plotter('dist2x', experimental,figCounter,distBest,distBest);
				saveas(figure(figCounter - 1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsDoubExp2Ind' strat '' num2str(topIndex) '.jpg']);
				
				figCounter = plotter('shadow', experimental,figCounter,proxBest,distBest);
				saveas(figure(figCounter - 1),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsDoubExp3Ind' strat '' num2str(topIndex) '.jpg']);

				saveas(figure(cvFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\cvSimDoubInd' strat '' num2str(topIndex) '.jpg']);
			    saveas(figure(propsFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\propsSimDoubInd' strat '' num2str(topIndex) '.jpg']);
			    saveas(figure(expCVFigIndex),[pwd '\' 'results' '\' strat '\' resultsFolder  '\cvExpDoubInd' strat '' num2str(topIndex) '.jpg']);

			end

		end

		close all

		% doing bootstrapping case 
		if boot && ~ignoreSave
				strat = 'boot';
				% in bootstrapping we are analyzing 3 models:
				% 2x prox, 2x dist, and shadow
				models = {'proximal','distal','shadow'};
				numModels = length(models); 
				readerEnergies(strat,stratDesired,og);

				plotTracker = 300;

				% saving figures
				for j = 1:numModels

					for k = 1:4
						% 4 histograms for each param rankings
						saveas(figure(plotTracker),[pwd '\' 'results' '\' strat '\' resultsFolder  '\' models{j} 'hist' num2str(k) '' num2str(topIndex) '.jpg']);
						plotTracker = plotTracker + 1;
					end

					% 1 combined scatter plot of all the SSE
					% for each params
					saveas(figure(plotTracker),[pwd '\' 'results' '\' strat '\' resultsFolder  '\' models{j} 'energyScatter' num2str(topIndex) '.jpg']);
					plotTracker = plotTracker + 1;
				end

				% boxplot for storage
				saveas(figure(plotTracker),[pwd '\' 'results' '\' strat '\' resultsFolder  '\box' num2str(topIndex) '.jpg']);
				close all

		end
	 
	end

end

toc