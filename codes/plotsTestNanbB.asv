%files = dir('output\*.mat');
%Nf = numel(files);
Nf = 1;

x = zeros(1,Nf);
e_m = zeros(1,Nf);
e_std = zeros(1,Nf);
n_m = zeros(1,Nf);
n_std = zeros(1,Nf);

for i = 1:Nf
    
    %name = files(i).name;
    %x(i) = str2double(name(5:end-4));
    %out = load(['output\' name]);
    out = load('output\Im_L800.mat');
    out = out.Im_mic_point;
        
    v = var(out,0,3);
    m = mean(out,3);
    
    B300 = v./m;
    %B_m = mean(B(:));
    %B_std = std(B(:))/sqrt(numel(B(:)));
    
    e = (v-m)./m;
    e_m(i) = mean(e(:));
    e_std(i) = std(e(:))/sqrt(numel(e(:)));
    
    
    N = (m.^2)./v;  
    %N_m = mean(N(:));
    %N_std = std(N(:))/sqrt(numel(N(:)));
    
    n = (m.^2)./(v-m);
    n_m(i) = median(n(:));
    n_std(i) = std(n(:))/sqrt(numel(n(:)));
    
    
end
%%
figure(1)
errorbar(x,e_m,e_std,'.','MarkerSize',15,'Linewidth',1)
xlab = xlabel('\lambda');
ylab = ylabel('\epsilon');
setplot(xlab,ylab,0,'output\e_L')

%%
figure(2)
errorbar(x,e_m./x,e_std./x,'.','MarkerSize',15,'Linewidth',1)
xlab = xlabel('\lambda');
ylab = ylabel('\epsilon / \lambda');
setplot(xlab,ylab,0,'output\e_variations')
%% 
figure(3)
errorbar(x,n_m,n_std,'.','MarkerSize',15,'Linewidth',1)
xlab = xlabel('\lambda');
ylab = ylabel('n');
setplot(xlab,ylab,0,'output\numbers')

%% T steps nedded to reach brightness value

files = dir('output\Im_L800Ts*.mat');

x = 10^(-2)*[0.1,0.5,1,5];
%x = zeros(1,numel(files));
step = 1;
steps = step :step :1000;%size(out,3);
e_m = zeros(numel(files),numel(steps) );
    
for i = 2%[1,4,2,3]%:numel(files)
    
    name = files(i).name;
    %x(i) = str2double(name(5:end-4));
    out = load(['output\' name]);
    out = out.Im_mic_tot;
        
    
    for j = steps;
        
        fprintf('-- %d -- \n',j)
        
        v = var(out(:,:,1:j),0,3);
        m = mean(out(:,:,1:j),3);

        e = (v-m)./m;
        e_m(i,j/step) = mean(e(:));
        
    end
    name
       
end
%%
cmap = winter(numel(files));
figure(1)
for i = 2%: numel(files)
    plot(steps.*5*10^(-2),e_m(i,:)./800,'.-','Markersize',15,'Linewidth',2,'color',cmap(i,:))
    %plot(e_m(i,:)./800,'.-','Markersize',15,'Linewidth',2,'color',cmap(i,:))
    hold on
end
%%
ylim([0,0.35])
xlab = xlabel('time (seconds)');
ylab = ylabel('\epsilon/\lambda');
%legend({'10','20','50','100','300','600','1200','2400'},'Location','Southeast')
setplot(xlab,ylab,0,'output\Ntstep')