
function  [MatrixPerimeter_ranger,MatrixArea_ranger] = Stage3_2_True_extract(nameMovie,pathMovie,tmin,tmax,zmax,nombreCelluleEtudie,MatrixCentroid_x,MatrixCentroid_y,MatrixPerimeter,MatrixArea,zmin)
%% Load centroid

centroid_identite_cell = load([pathMovie filesep 'Data' filesep 'centroid_identite_cell_TEMPS_Z' '.mat']);
centroid_all = centroid_identite_cell.MatrixCentroid_final_final_t_z_FINAL;

%% Initialization

MatrixPerimeter_ranger = NaN(tmax,zmax,nombreCelluleEtudie);
MatrixArea_ranger = NaN(tmax,zmax,nombreCelluleEtudie);

%% Order matrix

for t = tmin:tmax
    for z = zmin:zmax
        
        KNOW_centroid_t_z = centroid_all(:,:,t,z);
        
        centroid_x_t_z = MatrixCentroid_x(t,z,:);
        centroid_y_t_z = MatrixCentroid_y(t,z,:);
        
        for n = 1:size(KNOW_centroid_t_z,1)
            
            colX = ismember(centroid_x_t_z,KNOW_centroid_t_z(n,1));
            colY = ismember(centroid_y_t_z,KNOW_centroid_t_z(n,2));
            
            index_tf = colX==1 & colY==1;
            
            if sum(index_tf) > 1
                disp('Several centroid found = probably OVERSEGMENTATION')
                t
                z
                n
            end
            
            index = find(index_tf);
            
            if size(index,1)==0
                continue
            end
            
            MatrixPerimeter_ranger(t,z,n) = MatrixPerimeter(t,z,index);
            MatrixArea_ranger(t,z,n) = MatrixArea(t,z,index);
        end
    end
end