function cbv = som(somSize, input_data)

% Create the Self-Organizing Map
somNet = selforgmap(somSize, 150, 3,'tritop','mandist'); %change parameters from here
somNet = train(somNet, input_data');

% Extract the codebook vectors
if length(somSize) == 1
   cbv = somNet.IW{1};
elseif length(somSize) ==2
   cbv = somNet.IW{1,1};
else 
    disp('Error: check SOM dimensionality')

end