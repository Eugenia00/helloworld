function NandBsimulation2species
%% Parameters
nParticles =  500;
w_r = 0.35;      %0.35;% -0.4:0.4
w_z = 1;      %0.15; % -0.2:0.2
Size = w_r*30;      %w_r*10^2; % better integer
kon = 1; % number of events per second per molecule
koff = 2;
D = 10; %micron^2/second
tStep = 1 * 10^-3; %seconds
nSteps = 10000;
mean_particle_brightness = 800;
pix_size = 0.2;

%dz = sqrt(-(w_z^2)/4*log(0.3)); % dept able to detect (until 10% of intensity)

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
[xm,ym] = meshgrid(-Size/2:pix_size:Size/2);
% size of size(xm,1) has to be smaller than sqrt(3.166e+10/8/nSteps)
Im_mic_tot = zeros(size(xm,1),size(ym,1),numel(TimePoints));
%Im_mic_point = zeros(size(xm,1),size(ym,1),nSteps);

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

% Temporal evolution
j = 0; % acquisition counter
%%
tic;
for i = 1:nSteps
    
    
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
    
    %     % image of evolution in visible z depth
    %     figure(1)
    %     subplot(1,3,1)
    %     FlatHist( Evolution( abs(Evolution(:,3))<dz ,1), Evolution( abs(Evolution(:,3))<dz,2), xm(1,:),2)
    
    
    if rem(i,Acquire_evry) == 0
        
        j = j+1;
        x(j,:) = Evolution(:,1);
        y(j,:) = Evolution(:,2);
        z(j,:) = Evolution(:,3);
        
        isbound(j,:) = Bound;
        
        brightness(j,:) = poisson_eu(mean_particle_brightness,nParticles);
        brightness(j,logical(isbound(j,:))) = 4*brightness(j,logical(isbound(j,:)));
        
        % repetition of particle coordinates (rows) times the number of pixels
        % on the new image (col) to avoid for loops on the number of the pixels of the final image.
        % The pixels on the simulated mic image are called with linear indexing
        X = repmat(reshape(x(j,:),[],1),1,numel(xm));
        Y = repmat(reshape(y(j,:),[],1),1,numel(xm));
        Z = repmat(reshape(z(j,:),[],1),1,numel(xm));
        B = repmat(reshape(brightness(j,:),[],1),1,numel(xm));
        % same dimensions as XYZ for the mic image coordinates
        Xm = repmat(reshape(xm,1,[]),nParticles,1);
        Ym = repmat(reshape(ym,1,[]),nParticles,1);
        
        % convolution on x and y integration on z
        M_mic = B .* exp(-4*((X-Xm).^2+(Y-Ym).^2)./w_r^2) .* exp(-4*(Z.^2)./w_z^2);
        Im_mic = reshape(sum(M_mic,1),size(xm))';
        
        
        Im_mic_tot(:,:,j) = Im_mic;
        
        
        
        
        %         figure(1)
        %         subplot(1,3,2)
        %         % repetition of particles number according to their brightness
        
        %xbri = RepeatArrayElements(Evolution(:,1),brightness(j,:));
        %ybri = RepeatArrayElements(Evolution(:,2),brightness(j,:));
        %zbri = RepeatArrayElements(Evolution(:,3),brightness(j,:));
        
        %edges{1} = xm(1,:);
        %edges{2} = xm(1,:);
        %Im_mic_point(:,:,i) = hist3([xbri(abs(zbri)<dz), ybri(abs(zbri)<dz)], 'Edges',edges);
        
        %         FlatHist(xbri(abs(zbri)<dz) , ybri(abs(zbri)<dz), xm(1,:),...
        %             mean_particle_brightness+2*sqrt(mean_particle_brightness))
        %         subplot(1,3,3)
        %         imshow(Im_mic,[],'InitialMagnification',200)
        %         shg
        %         %pause
        
    end
    
end
save(['../output/Im_' sprintf('L%d',mean_particle_brightness) '_kon1_Koff2_tetr.mat'],'Im_mic_tot')
%save(['output\Im_points_' sprintf('L%d',mean_particle_brightness) '.mat'],'Im_mic_point')
toc

end