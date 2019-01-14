
function n_xpos = poisson_eu(mean,n)

% Input Variables
% ------------------
% mean: mean of the Poissonian distribution
% n: number of output variables
% ------------------
%
% Extracts n*mean particles integer positions in a 1D box of size n
% Counts how many ended in each position
% Faster than pm = poissrnd(mean,1,n);


event_pos = ceil( rand(1,n*mean) *n);
acc = accumarray(event_pos',1)';
n_xpos = zeros(1,n);
n_xpos(1:numel(acc)) = acc;

end

