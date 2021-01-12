function Stage1_2_Unionseg(t,nameMovie,pathMovie,zmin,zmax)

        for z=zmin:zmax
            
            disp(['Creating Unionseg file for frame #' num2str(t) ' and slice #' num2str(z) '...']);
            
            union = [pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_results' filesep 'directskel_' nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.png'];
            
            destinationUnionSeg = [pathMovie filesep 't' num2str(t,'%04d') filesep 'Output_results' filesep 'Unionseg_' nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.png'];
            copyfile (union , destinationUnionSeg)
        end
end