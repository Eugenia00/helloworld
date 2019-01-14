
function NandBsimulation(kon,koff)

% Parameters
nParticles =  5;
w_r = 0.5;      %0.35;% -0.4:0.4
w_z = 50;      %0.15; % -0.2:0.2
Size = w_r*10;      %w_r*10^2; % better integer
%kon = 2; % number of exents per second per molecule
%koff = 3;
D = 100; %micron/second
tStep = 2.5 * 10^-4; %seconds
nSteps = 1000;
mean_particle_brightness = 1;

% Initialize variables
Acquire_evry = 1;
TimePoints = Acquire_evry:Acquire_evry:nSteps;
NAqPoints = numel(TimePoints);

% particles
x = zeros(NAqPoints, nParticles);
y = zeros(NAqPoints, nParticles);
z = zeros(NAqPoints, nParticles);
isbound = zeros(NAqPoints, nParticles);
brightness = zeros(NAqPoints, nParticles);

% final image
[xm,ym] = meshgrid(-Size/2:Size/2);
Xm = repmat(reshape(xm,1,[]),nParticles,1);
Ym = repmat(reshape(ym,1,[]),nParticles,1);

% Calculate probabilities
p_bind = kon * tStep; % probability to bind
p_unbind = koff * tStep; % probability to unbind
% Steady state probability of having a bound molecule
p_bound_0 = kon/(kon + koff);


% Initial condition
Initial(:,1:3) =  Size * rand(nParticles,3) - Size/2;
Evolution = Initial;

Test =  rand(nParticles,1);
Bound = Test < p_bound_0; % Gives 1 to bound molecules,
Free = 1-Bound;

%% Temporal evolution
j = 0; % acquisition counter

for i = 1:nSteps
    %%
    fprintf('--- %d --- \n',i)
    % binding evolution
    Test =  rand(nParticles,1);
    Bound = (Bound.*Test > p_unbind) | (Free.*Test > (1 - p_bind));
    Free = 1-Bound;
    
    % displacement evolution of free particles
    Displacement = randn(nParticles,3);
    Evolution(:,1) =  Evolution(:,1) + Free * sqrt(2*D*tStep).*Displacement(:,1);
    Evolution(:,2) = Evolution(:,2) + Free * sqrt(2*D*tStep).*Displacement(:,2);
    Evolution(:,3) = Evolution(:,3) + Free * sqrt(2*D*tStep).*Displacement(:,3);
    
    % Boundary reflection
    Evolution(Evolution > Size/2) = Size - Evolution(Evolution > Size/2);
    Evolution(Evolution < -Size/2) = - Size - Evolution(Evolution < -Size/2);
    
    %     % display evolution
    %     Im_dyn = zeros(Size);
    %     %x_im_dyn = ceil(Evolution( abs(Evolution(:,3))<5 ,1)+Size/2);
    %     %y_im_dyn = ceil(Evolution( abs(Evolution(:,3))<5 ,2)+Size/2);
    %     x_im_dyn = ceil(Evolution(:,1)+Size/2);
    %     y_im_dyn = ceil(Evolution(:,2)+Size/2);
    %     ind = sub2ind(size(Im_dyn),x_im_dyn,y_im_dyn);
    %     acc = accumarray( reshape(ind,[],1) ,1);
    %     Im_dyn(1:numel(acc)) = acc;
    %
    %     figure(1)
    %     subplot(1,2,1)
    %     imshow(Im_dyn,[0,3],'InitialMagnification',200)
    %     shg
    
    %%
    if rem(i,Acquire_evry) == 0
        
        j = j+1;
        x(j,:) = Evolution(:,1);
        y(j,:) = Evolution(:,2);
        z(j,:) = Evolution(:,3);
        
        isbound(j,:) = Bound;
        brightness(j,:) = poisson_eu(mean_particle_brightness,nParticles);
        %x = RepeatArrayElements(Evolution(:,1),brightness(j,:));
        X = repmat(Evolution(:,1),1,numel(xm));
        Y = repmat(Evolution(:,2),1,numel(xm));
        Z = repmat(Evolution(:,3),1,numel(xm));
        B = repmat(brightness(j,:)',1,numel(xm));
        
        M_mic = B .* exp(-4*((X-Xm).^2+(Y-Ym).^2)./w_r^2) .* exp(-4*(Z.^2)./w_z^2);
        Im_mic = reshape(sum(M_mic,1),size(xm));
        
        % display evolution
        Im_dyn = zeros(Size);
        %x_im_dyn = ceil(Evolution( abs(Evolution(:,3))<5 ,1)+Size/2);
        %y_im_dyn = ceil(Evolution( abs(Evolution(:,3))<5 ,2)+Size/2);
        x_im_dyn = ceil(Evolution(:,1)+Size/2);
        y_im_dyn = ceil(Evolution(:,2)+Size/2);
        ind = sub2ind(size(Im_dyn),x_im_dyn,y_im_dyn);
        acc = accumarray( reshape(ind,[],1) ,1);
        Im_dyn(1:numel(acc)) = acc;
        
        figure(1)
        subplot(1,3,1)
        imshow(Im_dyn,[0,3],'InitialMagnification',200)
        shg
        subplot(1,3,2)
        imshow(Im_mic,[0,3],'InitialMagnification',200)
        shg
        pause
        
    end
    
end