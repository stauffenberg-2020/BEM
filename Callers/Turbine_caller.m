function [output, output_details, BEM] = Turbine_caller(file)
    
    % Temporary initilization of structures for MATLAB Coder compatibility

    [General, op_pts, Blade] = read_turbine_file(file);


    Blade.preflap = zeros(length(Blade.r),1);

    [output_details, output, BEM] = core_bem(General, op_pts, Blade);
    
    %% Writing the outputs into a text file
    outfile = 'output.txt';
    fileID = fopen(outfile,'w');
    for j=1:size(output,2)
        fprintf(fileID,'%s \t',output_details{2,j});
        if j==size(output,2)
            fprintf(fileID,'\n');
        end
    end
    for j=1:size(output,2)
        fprintf(fileID,'%s \t\t',output_details{3,j});
        if j==size(output,2)
            fprintf(fileID,'\n');
        end
    end
    for i=1:size(output,1)
        for j=1:size(output,2)
            fprintf(fileID,'%8.2f \t',output(i,j));
            if j==size(output,2)
                fprintf(fileID,'\n');
            end
        end
    end
    fclose(fileID);
    
end
