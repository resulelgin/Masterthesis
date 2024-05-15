%%% Workflow for Flux Balance Analysis (FBA) of Recon3D model with CobraToolbox %%%
%%% Workflow proceeds to knockout Hmgcs2 gene and analyse the effects from the knockout %%%
%%% Relaxed FBA analysis is encouraged %%%

% Read the model available in VHM - Ines Thiele lab.
model1 = readCbModel('Recon3D_301.mat');
disp(fieldnames(model1));

modelData = load('Recon3D_301.mat'); 
modelFields = fieldnames(modelData);

% Check if 'model' is one of the fields
if ismember('model', modelFields)
    model = modelData.model; % Extract the model structure
else
    % If the model structure is stored under a different variable, adjust accordingly
    % For example, if stored under 'Recon3D':
    model = modelData.('Recon3D');
end

% Set solver if not already set
changeCobraSolver('mosek', 'all');

% Perform FBA
solutionOriginal = optimizeCbModel(model);

% Check the solution
if solutionOriginal.stat == 1
    disp('The model is feasible and optimal for FBA.');
else
    disp('The model is not feasible for FBA. Please check the model.');
end

% List all genes in the model
allGenes = model.genes;

% Genes of interest
genesOfInterest = {'3158.1', '38.1', '5019.1', '3155.1', '622.1'}; % Hmgcs2, 

% Check if your genes of interest are present in the model
isPresent = ismember(genesOfInterest, allGenes);

% Display results
for i = 1:length(genesOfInterest)
    if isPresent(i)
        disp(['Gene ', genesOfInterest{i}, ' is present in the model.']);
    else
        disp(['Gene ', genesOfInterest{i}, ' is NOT present in the model.']);
    end
end

% Display the FBA results
disp('Original model FBA results:');
disp(solutionOriginal);

% Gene to be knocked out
knockoutGene = '3158.1'; % Hmgcs2

% Knockout the gene in the model
modelKO = deleteModelGenes(model, knockoutGene);

% Perform FBA on the knockout model
solutionKO = optimizeCbModel(modelKO);

% Check if the solution is feasible
if solutionKO.stat == 1
    disp(['The model with knockout of ', knockoutGene, ' is feasible and optimal for FBA.']);
else
    disp(['The model with knockout of ', knockoutGene, ' is not feasible for FBA. Please check the model.']);
end

% Display the FBA results for the knockout model
disp(['Knockout model FBA results for gene ', knockoutGene, ':']);
disp(solutionKO);

% Compare the fluxes
fluxComparison = table(model.rxns, solutionOriginal.x, solutionKO.x, ...
    'VariableNames', {'Reaction', 'OriginalFlux', 'KnockoutFlux'});

% Display the flux comparison
disp('Comparison of reaction fluxes:');
disp(fluxComparison);

% Define the flux change threshold
fluxChangeThreshold = 1e-2; % Adjust the threshold as seen fit

% Calculate the absolute change in flux for each reaction
fluxComparison.FluxChange = (fluxComparison.KnockoutFlux - fluxComparison.OriginalFlux);

% Subset the reactions based on the flux change threshold
significantFluxChanges = fluxComparison(abs(fluxComparison.FluxChange) > fluxChangeThreshold, :);

% Display the subset of reactions with significant changes
disp('Reactions with significant flux changes:');
disp(significantFluxChanges);

% Define the reactions of interest
reactionsOfInterest = {'OCOAT1m', 'HMGCOASim', 'BDHm', 'ACACT1rm', 'HMGLm', ...
    'ADRNCPT1', 'HMGCOASi'};

% Filter fluxComparison for reactions of interest
interestIdx = ismember(fluxComparison.Reaction, reactionsOfInterest);
fluxComparisonInterest = fluxComparison(interestIdx, :);

% Display the filtered table
disp('Flux comparison for reactions of interest:');
disp(fluxComparisonInterest);

% Filter significant changes for reactions of interest
significantChangesInterest = fluxComparison(significantChanges & interestIdx, :);

% Display the filtered table
disp('Significant changes for reactions of interest:');
disp(significantChangesInterest);

% Initialize a list to store genes corresponding to significant reactions
significantGenes = {};

% Loop over each significant reaction
for i = 1:height(significantFluxChanges)
    % Find the index of the reaction in the model
    rxnIndex = find(strcmp(model.rxns, significantFluxChanges.Reaction{i}));
    
    % Get the gene indices from the rxnGeneMat matrix
    geneIndices = find(model.rxnGeneMat(rxnIndex, :));
    
    % Get the gene IDs and add them to the list of significant genes
    for j = 1:length(geneIndices)
        significantGenes{end+1} = model.genes{geneIndices(j)};
    end
end

% Remove duplicate genes
significantGenes = unique(significantGenes);

% Display the list of significant genes
disp('Genes corresponding to reactions with significant flux changes:');
disp(significantGenes);

% Convert the list of significant genes to a table
significantGenesTable = table(significantGenes', 'VariableNames', {'GeneID'});

% Save the gene list to a file (e.g., for use with GSEA tools)
writetable(significantGenesTable, 'significant_genes.txt', 'WriteVariableNames', false);

% Extract the Reaction column from the significantFluxChanges table
reactionList = significantFluxChanges.Reaction;

% Convert the list of significant reactions to a table
significantReactionsTable = table(reactionList, 'VariableNames', {'Reaction'});

% Save the reaction list to a file
writetable(significantReactionsTable, 'significant_reactions.txt', 'WriteVariableNames', false);

% Initialize lists to store up- and down-regulated genes and reactions
upRegulatedGenes = {};
downRegulatedGenes = {};
upRegulatedReactions = {};
downRegulatedReactions = {};

% Loop over each significant reaction
for i = 1:height(significantFluxChanges)
    % Find the index of the reaction in the model
    rxnIndex = find(strcmp(model.rxns, significantFluxChanges.Reaction{i}));
    
    % Get the gene indices from the rxnGeneMat matrix
    geneIndices = find(model.rxnGeneMat(rxnIndex, :));
    
    % Check if the flux change is positive or negative
    if significantFluxChanges.FluxChange(i) > 0
        % Add to up-regulated reactions
        upRegulatedReactions{end+1} = significantFluxChanges.Reaction{i};
        
        % Get the gene IDs and add them to the list of up-regulated genes
        for j = 1:length(geneIndices)
            upRegulatedGenes{end+1} = model.genes{geneIndices(j)};
        end
    else
        % Add to down-regulated reactions
        downRegulatedReactions{end+1} = significantFluxChanges.Reaction{i};
        
        % Get the gene IDs and add them to the list of down-regulated genes
        for j = 1:length(geneIndices)
            downRegulatedGenes{end+1} = model.genes{geneIndices(j)};
        end
    end
end

% Remove duplicate genes
upRegulatedGenes = unique(upRegulatedGenes);
downRegulatedGenes = unique(downRegulatedGenes);

% Display the lists of significant genes
disp('Up-regulated genes:');
disp(upRegulatedGenes);

disp('Down-regulated genes:');
disp(downRegulatedGenes);

% Convert the lists of significant genes to tables
upRegulatedGenesTable = table(upRegulatedGenes', 'VariableNames', {'GeneID'});
downRegulatedGenesTable = table(downRegulatedGenes', 'VariableNames', {'GeneID'});

% Convert the lists of significant reactions to tables
upRegulatedReactionsTable = table(upRegulatedReactions', 'VariableNames', {'Reaction'});
downRegulatedReactionsTable = table(downRegulatedReactions', 'VariableNames', {'Reaction'});

% Save the gene lists to files (e.g., for use with GSEA tools)
writetable(upRegulatedGenesTable, 'up_regulated_genes.txt', 'WriteVariableNames', false);
writetable(downRegulatedGenesTable, 'down_regulated_genes.txt', 'WriteVariableNames', false);

% Save the reaction lists to files
writetable(upRegulatedReactionsTable, 'up_regulated_reactions.txt', 'WriteVariableNames', false);
writetable(downRegulatedReactionsTable, 'down_regulated_reactions.txt', 'WriteVariableNames', false);

disp('Up-regulated genes have been saved to up_regulated_genes.txt');
disp('Down-regulated genes have been saved to down_regulated_genes.txt');
disp('Up-regulated reactions have been saved to up_regulated_reactions.txt');
disp('Down-regulated reactions have been saved to down_regulated_reactions.txt');
