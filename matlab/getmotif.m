function [mat, names] = getmotif(filename, m)
% filename is the meme file that contains the motifs;
% n is the nth motif in the file
mat = cell(length(m),1);
names = cell(length(m),1);
[n ind] = sort(m);
fid = fopen(filename);
if fid < 0
    fprintf('No good')
else
   i=0; 
    for j = 1:length(n)
        while i~=n(j)
            line = fgetl(fid);
            if length(line) >= 5
                if strcmp(line(1:5), 'MOTIF')
                    i = i+1;
                end
            end
        end
        a = strsplit(strip(line), ' ');
        names{j} = a{end};
        line = fgetl(fid);
        line = fgetl(fid);
        mat{j} = [];
        while ~isempty(line)
            mat{j} = [mat{j}; str2num(line)];
            line = fgetl(fid);
        end
    end
end
fclose(fid);
if length(n)==1
    mat = cell2mat(mat);
else
    ind2 = zeros(length(n),1);
    for i = 1:length(n)
        ind2(i) = find(n==m(i));
    end
    mat = mat(ind2);
    names = names(ind2);
end
