function seq2var(sfn, wfn, ofn)
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

seq = importdata(sfn);
n = length(seq)/2;
disp('running dsvm')
mat = [1 2 3;0 2 3;0 1 3;0 1 2];
O = ones(1,l);
for i = 1:n
    if mod(i,100)==0
        disp([num2str(i) ' sequences converted'])
    end
    L = length(seq{2*i})-2*l+2;
    ss = letterconvert(seq{2*i});
    p = zeros(L,1);
    I = ss(1:l)*pow;
    p(1) = w(I+1);
    for j = 2:L+l-1
        I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1);
        p(j) = w(I+1);
    end
    v = zeros(3*L,1);
    a = 1;
    for j = l:L
        S = ss(j-l+1:j+l-1);
        ref = O*p(j-l+1:j);
        for ii = 1:3
            S(l) = mat(ss(j)+1,ii);
            I = S(1:l)*pow;
            v(a) = w(I+1);
            for jj = 2:l
                I = (I-S(jj-1))/4+4^(l-1)*S(jj+l-1);
                v(a)=v(a)+w(I+1);
            end
            v(a) = v(a)-ref;
            a = a+1;
        end
    end
    dlmwrite([ofn '_' num2str(i) '_dsvm.out'],v);
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
