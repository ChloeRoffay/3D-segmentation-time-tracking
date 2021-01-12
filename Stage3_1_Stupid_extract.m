function  [MatrixCentroid_x,MatrixCentroid_y,MatrixPerimeter,MatrixArea] = Stage3_1_Stupid_extract(nameMovie,pathMovie,tmin,tmax,zmin,zmax,scale1D,scale2D)

%% Initialization
nombreCelluleMaxAttendues = 1000;

MatrixCentroid_x = zeros(tmax,zmax,nombreCelluleMaxAttendues);
MatrixCentroid_y = zeros(tmax,zmax,nombreCelluleMaxAttendues);
MatrixPerimeter = zeros(tmax,zmax,nombreCelluleMaxAttendues);
MatrixArea = zeros(tmax,zmax,nombreCelluleMaxAttendues);

%% Loop

for t=tmin:tmax
    for z = zmin:zmax
        BW = imread([pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_results' filesep 'Unionseg_' nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.png']);
        
        image_CC = bwconncomp(BW,4);
        Stat_All_Regions = regionprops(image_CC,'Area','Centroid','PixelIdxList','Perimeter');
        
        I = imread([pathMovie filesep 't' num2str(t,'%04d') filesep nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.tif']);
        Stat_All_Regions_pixel = regionprops(image_CC,I,'PixelValues');

        Stat_All_Regions_cell = (struct2cell(Stat_All_Regions))';
        
        % Area => volume
        matArea = cell2mat(Stat_All_Regions_cell(:,1))*scale2D;
        MatrixArea(t,z,1:size(matArea,1)) = matArea;
        
        % Centroid
        matCentroid = cell2mat(Stat_All_Regions_cell(:,2))*scale1D;
        matCentroid_x = matCentroid(:,1);
        matCentroid_y = matCentroid(:,2);
        
        MatrixCentroid_x(t,z,1:size(matCentroid,1)) = matCentroid_x;
        MatrixCentroid_y(t,z,1:size(matCentroid,1)) = matCentroid_y;
        
        % Perimeter => Surface
        matPerimeter = cell2mat(Stat_All_Regions_cell(:,4))*scale1D;
        MatrixPerimeter(t,z,1:size(matPerimeter,1)) = matPerimeter;
        
    end
end

%% Delete useless column

Matrix_Zeros = zeros(tmax,zmax,nombreCelluleMaxAttendues);
tf_zeros_inutil = ismember(MatrixArea,Matrix_Zeros);
n_cut = nombreCelluleMaxAttendues;

for n = n_cut:-1:1
    grille = tf_zeros_inutil(:,:,n);
    if size(find(grille<1),1) == 0
        n_cut = n_cut - 1;
    else
        break
    end
end

MatrixCentroid_x = MatrixCentroid_x(:,:,1:n_cut);
MatrixCentroid_y = MatrixCentroid_y(:,:,1:n_cut);
MatrixPerimeter = MatrixPerimeter(:,:,1:n_cut);
MatrixArea = MatrixArea(:,:,1:n_cut);

%% Save file 

dossierSave = [pathMovie filesep 'Data'];
status = mkdir(dossierSave);
save([pathMovie filesep 'Data' filesep nameMovie '.mat'],'MatrixCentroid_x','MatrixCentroid_y','MatrixPerimeter','MatrixArea','-v7.3');

