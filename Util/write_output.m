function write_output(output_details, output, file)
    %% Writing the outputs into a text file
    tmp = extractBefore(file,'.');
    outfile = plus(tmp,"_output.txt");
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