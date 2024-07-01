function cbv = som_batch(somSize, input_data, batchSize)

% Create the Self-Organizing Map
somNet = selforgmap(somSize, 100, 8,'hextop','linkdist'); %change parameters from here

% Parameters for batch training
numBatches = ceil(size(input_data, 1) / batchSize);
idx = randperm(size(input_data,1)); % Random batches

% Train the Self-Organizing Map in batches
for i = 1:numBatches
    batch_idx = idx((i-1)*batchSize + 1 : i*batchSize);
    somNet = train(somNet, input_data(batch_idx, :)');
    fprintf('Trained on batch %d of %d\n', i, numBatches);
end

% Extract the codebook vectors
if length(somSize) == 1
   cbv = somNet.IW{1};
elseif length(somSize) ==2
   cbv = somNet.IW{1,1};
else 
    disp('Error: check SOM dimensionality')

end