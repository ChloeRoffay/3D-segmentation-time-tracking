function detect_bottom_matrix = Stage1_3_detect_bottom(t,nameMovie,pathMovie,zmin,zmax,tmin,zMilieu,nombreCelluleEtudie,detect_bottom_matrix)

%% Make the stack
nombre_z = zmax-zmin+1;

% Read size of the image
I_size= imread([pathMovie filesep 't' num2str(t,'%04d') filesep nameMovie '_t' num2str(t,'%04d') '_z' num2str(zmin,'%04d') '.tif']);

% Create matrix of the z-stack, fixed t
I_3D = zeros(size(I_size,1),size(I_size,2),nombre_z);
I_3D = uint16(I_3D);

for z = zmin:1:zmax
    I = imread([pathMovie filesep 't' num2str(t,'%04d') filesep nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.tif']);
    I_3D(:,:,z-zmin+1) = I;
end

I_3D = imadjust_stack(I_3D);

%% Detect bottom

sum_I_3D = squeeze(sum(sum(I_3D,2),1));
detect_bottom_matrix(:,t) = sum_I_3D;