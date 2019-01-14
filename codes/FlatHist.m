function FlatHist(X,Y,edg,maxint)
%%
edges{1} = edg;
edges{2} = edg;
n = hist3([X,Y],'Edges',edges);
pcolor(edg,edg,n);
%colormap(gray(max(n(:))));
colormap(gray(maxint));
axis ij
axis square

end