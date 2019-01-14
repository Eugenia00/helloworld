box_size = 256;
Ntot = 3000;
Ttot = 1000; %50;

% initialize
Occupation = zeros(box_size,box_size,Ttot);
Occupation_i = zeros(box_size,box_size);

ind = ceil(rand(1,Ntot)*box_size^2);

% count repeated indexes
acc = accumarray( reshape(ind,[],1) ,1);

Occupation_i(1:numel(acc)) = acc;
Occupation(:,:,1) = Occupation_i;

%% evolution
for i = 2:Ttot

    evolution = floor(rand(1,Ntot)*5)-2; % -2:2    
    % +/-1 on x same as +/-boxsize on linear indexes
    evolution(evolution == 2) = box_size; 
    evolution(evolution == -2) = -box_size; 
    
    ind = ind + evolution;
    
    % boundary periodic conditions
    ind(ind>box_size^2) = rem(ind(ind>box_size^2) + evolution(ind>box_size^2),box_size^2);
    ind(ind<1) = box_size^2 + ind(ind<1);
    
    % count repeated indexes
    acc = accumarray(reshape(ind,[],1),1);
    
    Occupation_i = zeros(box_size);
    Occupation_i(1:numel(acc)) = acc;
    Occupation(:,:,i) = Occupation_i;
    fprintf('--- %d --- \n',i)
    %max(Ai(:))
    
    imshow(Occupation_i,[0,3]); shg
    %pause
    
end
