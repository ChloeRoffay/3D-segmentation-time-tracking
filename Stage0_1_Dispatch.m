function Stage0_1_Dispatch(t,nameMovie,pathMovie,zmin,zmax)

% Copy file to make it confortable later
dossier = [pathMovie filesep 't' num2str(t,'%04d')];
status = mkdir(dossier);

for z = zmin:zmax
    source = [pathMovie filesep nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.tif'];
end

destination = [pathMovie filesep 't' num2str(t,'%04d') filesep nameMovie '_t' num2str(t,'%04d') '_z' num2str(z,'%04d') '.tif'];
copyfile (source , destination)
delete(source)

end
