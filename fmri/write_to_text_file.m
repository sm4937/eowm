function [success] = writeToTextFile(data,filename)
%writeToTextFile Write matrix data with(out) blank entries to a textfile
%   Written to aid in running the AFNI_proc.py GLM pipeline, devised by
%   Tommy & Masih to process fMRI data

    fid = fopen(filename,'w');
    
    for i = 1:size(data, 1)
        
      line = data(i,data(i,:)~=-1);
      line = [line -1];
      line = num2str(line);
        
      fprintf(fid, line);
      fprintf(fid, '\n');
      
    end
    
    fclose(fid);
    
end

