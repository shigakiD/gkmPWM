function seq2prob(sfn, wfn, ofn)
%sfn: fasta file
%wfn: kmer weight file
%ofn: output header
fid = fopen(wfn, 'r');
X = textscan(fid,'%s\t%f\n');
fclose(fid);
l = length(X{1}{1});
w = zeros(4^l,1);
pow = (4.^(0:(l-1)))';
disp('calculating indices')
for i = 1:numel(X{2});
    ss = letterconvert(X{1}{i});
    rs = 3-fliplr(ss);
    w(ss*pow+1) = X{2}(i);
    w(rs*pow+1) = X{2}(i);
end
m = mean(X{2});
s = std(X{2});
w = (1/2)*(1+erf((w-m)/s));

seq = importdata(sfn);
n = length(seq)/2;
disp('converting kmers to  probabilities')
for i = 1:n
    if mod(i,100)==0
        disp([num2str(i) ' sequences converted'])
    end
    L = length(seq{2*i})-l+1;
    ss = letterconvert(seq{2*i});
    p = zeros(L,1);
    I = ss(1:l)*pow;
    p(1) = w(I+1);
    for j = 2:L
        I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1);
        p(j) = w(I+1);
    end
    dlmwrite([ofn '_' num2str(i) '_prob.out'],p);
end
disp('Done')

function en = letterconvert(s)

l = length(s);
en = zeros(1,l);
for i = 1:l
    if strcmp(s(i),'A') || strcmp(s(i), 'a')
        en(i) = 0;
    elseif strcmp(s(i),'C') || strcmp(s(i),'c')
        en(i) = 1;
    elseif strcmp(s(i),'G') || strcmp(s(i),'g')
        en(i) = 2;
    else
        en(i) = 3;
    end
end
