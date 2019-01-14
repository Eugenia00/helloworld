function Im_dyn = XYcont2Lin(X,Y,Size)

% X and Y range from -Size/2 to Size/2

x_im_dyn = ceil(X + Size/2);
y_im_dyn = ceil(Y + Size/2);
ind = sub2ind([Size,Size],x_im_dyn,y_im_dyn);
acc = accumarray( reshape(ind,[],1) ,1);
Im_dyn = zeros(Size);
Im_dyn(1:numel(acc)) = acc;

end