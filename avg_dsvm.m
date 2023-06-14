function avg_dsvm(fn,ofn)
%fn:fasta file
%ofn: output fileheader
seq = importdata(fn);
n = length(seq)/2;
for j = 1:n
    dsvm = dlmread([ofn '_' num2str(j) '_dsvm.out']);
    ss = letterconvert(seq{2*j});
    nvar = (length(ss)-(l-1)*2);
    mdsvm = zeros(nvar,1);
    if nvar == length(dsvm)/3
        for i =1:nvar;
            mdsvm(i) = -1*mean(dsvm(3*(i-1)+1:3*i));
        end
        fid2 = fopen([ofn '_' num2str(j) '_mdsvm_pwm.txt'],'w');
        fprintf(fid2, ['# one hot encode \n']);
        fprintf(fid2, 'pos\tA\tC\tG\tT\n');
        for i = 0:(l-2)
            fprintf(fid2, '%d\t%f\t%f\t%f\t%f\n', i, 0,0,0,0);
        end
        for i = l:length(ss)-l
            a = zeros(1,4);
            a(ss(i+1)) = mdsvm(i-9);
            fprintf(fid2, '%d\t%f\t%f\t%f\t%f\n', i, a);
        end
        for i = length(ss)-l+1:length(ss)-1
            fprintf(fid2, '%d\t%f\t%f\t%f\t%f\n', i, 0,0,0,0);
        end
        fclose(fid2);
    end
end

function en = letterconvert(s)

l = length(s);
en = zeros(1,l);
for i = 1:l
    if strcmp(s(i),'A') || strcmp(s(i), 'a')
        en(i) = 1;
    elseif strcmp(s(i),'C') || strcmp(s(i),'c')
        en(i) = 2;
    elseif strcmp(s(i),'G') || strcmp(s(i),'g')
        en(i) = 3;
    else
        en(i) = 4;
    end
end
