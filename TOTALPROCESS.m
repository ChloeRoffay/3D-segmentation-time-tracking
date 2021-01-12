clear
clc
%%%% One-click analysis
%%%% Identification number of the Movie Of Interest (MOI)
%%%% Possibility to analyse several movie at the same time
k_min = 1;
k_max = 2;


for k_movie = k_min:1:k_max
    
    % k_animal
    [nameMovie , info_function] = identifier_film(k_movie);
    [t_choc,zmin,zmax,zMilieu,tmin,tmax,pathMovie,scale1D,scale2D,nombreCelluleEtudie,confiance,step_time] = info_function(nameMovie);
    
    %% Etapes
    
    Dispatch = 0;
    Segmentation_Z_3D = 0;
    
    Detect_bottom = 0;
    
    Centroid_time = 1;
    Centroid_z = 1;
    
    Stupid_extract = 1;
    True_extract = 1;
    
    %% Dispatcher les images
    if Dispatch == 1
        for t = tmin:tmax
            Stage0_1_Dispatch(t,nameMovie,pathMovie,zmin,zmax)
            waitbar(t/tmax)
        end
    end
    
    %% Segmentation 3D
    if Segmentation_Z_3D == 1
        for t = tmin:tmax
            Stage1_1_Seg_directmax_3D(t,nameMovie,pathMovie,zmin,zmax,tmin,zMilieu,nombreCelluleEtudie)
            waitbar(t/tmax);
            Stage1_2_Unionseg(t,nameMovie,pathMovie,zmin,zmax)
        end
    end
    
    %% Detect and delete out of focus segmentation
    % The most intense image is the bottom of the dish
    
    if Detect_bottom == 1
        
        detect_bottom_matrix  = zeros(zmax,tmax);
        for t = tmin:tmax
            detect_bottom_matrix = Stage1_3_detect_bottom(t,nameMovie,pathMovie,zmin,zmax,tmin,zMilieu,nombreCelluleEtudie,detect_bottom_matrix);
        end
        
        intensity_z_time = mean(detect_bottom_matrix(:,tmin:tmax),2);
        
        % Detect maximum intensity and cut +1
        max_intensity = max(intensity_z_time);
        tf_bottom = intensity_z_time == max_intensity;
        z_cut_intensity = find(tf_bottom);
        z_cut_intensity = z_cut_intensity(1,1);
        
        % Delete skeleton out of focus
        I_size= imread([pathMovie filesep 't' num2str(tmin,'%04d') filesep nameMovie '_t' num2str(tmin,'%04d') '_z' num2str(zmin,'%04d') '.tif']);
        blanc = ones(size(I_size,1),size(I_size,2));
        
        for t = tmin:tmax
            for z = zmin:z_cut_intensity
                % Acquire image from the top or from the bottom ?
                % for z = z_cut_intensity:zmax
                % Adjust depending on how your images are acquired
                imwrite(blanc,[pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_results' filesep 'Unionseg_' nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.png'],'png');
            end
        end
        
    end
    
    %% Detect centroid of cells in the z-stack
    % First we fixe z and we follow cell in time
    % Second we find the cells in the z-stack
    
    if Centroid_time == 1
        
        detect_bottom_matrix  = zeros(zmax,tmax);
        for t = tmin:tmax
            detect_bottom_matrix = Stage1_3_detect_bottom(t,nameMovie,pathMovie,zmin,zmax,tmin,zMilieu,nombreCelluleEtudie,detect_bottom_matrix);
        end
        
        intensity_z_time = mean(detect_bottom_matrix(:,tmin:tmax),2);
        
        % Detect maximum intensity and select +4
        max_intensity = max(intensity_z_time);
        tf_bottom = intensity_z_time == max_intensity;
        z_cut_intensity = find(tf_bottom);
        zMilieu = z_cut_intensity(1,1)+4;
        
        % Find all the centroid of all the objects in time at zMilieu
        [MatrixCentroid_t_zMilieu,zMilieu_t] = Stage2_1_follow_cell(pathMovie,nameMovie,tmin,tmax,scale1D,nombreCelluleEtudie,zMilieu,k_movie,zmax);
    end
    
    if Centroid_z == 1
        nombreCelluleEtudie = size(MatrixCentroid_t_zMilieu,1);
        MatrixCentroid_t_z = NaN(nombreCelluleEtudie,2,tmax,zmax);
        for t = tmin:tmax
            [MatrixCentroid_final_t_z] = Stage2_2_follow_cell_z(pathMovie,nameMovie,tmin,tmax,scale1D,nombreCelluleEtudie,zMilieu,k_movie,t,zmin,zmax,zMilieu_t);
            MatrixCentroid_t_z(:,:,t,:) = MatrixCentroid_final_t_z;
            save([pathMovie filesep 'Data' filesep 'centroid_identite_cell_TEMPS_Z' '.mat'],'MatrixCentroid_final_final_t_z_FINAL','-v7.3');
        end
    end
    
    %% Extraction
    if Stupid_extract == 1
        [MatrixCentroid_x,MatrixCentroid_y,MatrixPerimeter,MatrixArea] = Stage3_1_Stupid_extract(nameMovie,pathMovie,tmin,tmax,zmin,zmax,scale1D,scale2D);
        disp('Stupid extract DONE')
    end
    
    %% Extraction with cell tracking
    if True_extract == 1
        %%% Load geometry to know cell number
        load([pathMovie filesep 'Data' filesep 'centroid_identite_cell' '.mat']);
        nombreCelluleEtudie = size(xChoixScale,1);
        
        load([pathMovie filesep 'Data' filesep nameMovie '.mat']);
        [MatrixPerimeter_ranger,MatrixArea_ranger] = Stage3_2_True_extract(nameMovie,pathMovie,tmin,tmax,zmax,nombreCelluleEtudie,MatrixCentroid_x,MatrixCentroid_y,MatrixPerimeter,MatrixArea,zmin);
    end
    
end
