function visualize_P (P)
    V = reshape(full(sum(P.P,1)), length(P.coords{1}), length(P.coords{2}), length(P.coords{3}));
    figure, imshow3D(V);
    
%     plane = ceil(length(Pcoords{3})/2);
%     inds = false(length(Pcoords{1}), length(Pcoords{2}), length(Pcoords{3}));
%     inds(:,:,plane) = true;
%     V = reshape(full(P(:,inds(:))), size(P,1), length(Pcoords{1}), length(Pcoords{2}));
%     V = permute(V, [2 3 1]);
%     figure, imshow3D(V)
end