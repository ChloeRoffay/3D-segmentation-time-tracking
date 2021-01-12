function out = imadjust_stack(I_3D)

I_3D_2D = reshape(I_3D,[size(I_3D,1),size(I_3D,2)*size(I_3D,3)]);

I_3D_2D = imadjust(I_3D_2D);

out = reshape(I_3D_2D,[size(I_3D,1),size(I_3D,2),size(I_3D,3)]);