function setplot(xlab,ylab,s,filename_save)

set([xlab,ylab],'fontWeight','bold','FontSize',14)
set(gca,'fontWeight','bold','FontSize',14) % bold ticks
set(gca,'TickDir','out')

set(gcf,'paperunits','centimeters');
set(gcf,'paperposition', [5 5 20 6]); % left bottom width height
set(gca,'ticklength',[0.02 0.01])
%set(gca,'box','on')
set(gca,'LineWidth',2)

if s==1
    saveas(gcf,filename_save, 'fig')
    print([filename_save '.eps'], '-depsc');
    system(['ps2pdf -dEPSCrop ' filename_save '.eps ' filename_save '.pdf']);
end

end
