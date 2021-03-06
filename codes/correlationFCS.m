%out = load('../output/Im_L800_kon0_dt10-4_1s.mat');
%out = load('../output/Im_L800_10000.mat');
%out = load('../output/Im_L800_kon5');
%out = load('../output/Im_L800kon0');
%out = load('../output/Im_L800_kon0_D1');
out = load('../output/Im_L800_kon1_Koff2_tetr');

out = out.Im_mic_tot;
S = size(out);
out = reshape(out,S(1)^2,S(3));
auto_out = zeros(S(1)^2,S(3)-1);
%
out_m = mean(out,2);
m = out-repmat(out_m,[1,size(out,2)]);
%
for i = 1:size(out,1)  
    
    disp(i)
    
    c = xcorr(m(i,:))/(out_m(i)^2);
    %c = xcorr(out(i,:));
    auto_out(i,:) = c(S(3)+1:end);
      
end

%%
co = mean(auto_out,1);
x = [1:numel(co)];
t = x.*10^(-3);
d_decorr = sqrt(t*10);

figure(2)
plot(t,co,'.-b','MarkerSize',15)
set(gca,'Xscale','log')
ylab = ylabel('G(\tau)');
xlab = xlabel('time (seconds)');
hold on
%%
setplot(xlab,ylab,1,'../output/corr_L800_D1D10')

