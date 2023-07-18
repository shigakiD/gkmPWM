function [gkmc, GCpos, GCneg, mat, mat2] = getgkmcounts(filename,l,k, lk,RC)
%gets the gapped kmer counts using the alphas to weight each support vector
% filename is the name of the support vector sequence fa ('_svseq.fa')
% l and k are the parameters for the length of the gapped kmer (l) and the number of ungapped positions (k)
if isfile([filename '_svseq.fa']) && isfile([filename '_svalpha.out'])
    filenameseq = [filename '_svseq.fa'];
    filenamealpha = [filename '_svalpha.out'];
    alpha = dlmread([filename '_svalpha.out'],'\t',0,1);
    sequences = importdata([filename '_svseq.fa']);
    ver = 0;
elseif isfile([filename '.model.txt'])
    fid = fopen([filename '.model.txt'], 'r');
    a = 1;
    while a == 1
        line = fgetl(fid);
        if strcmp('SV', line);
            a = 0;
        end
    end
    X = textscan(fid,'%f %s\n');
    sequences = X{2};
    alpha = X{1};
    ver = 1;
else
    error(['Needs ' filename '_svseq.fa and ' filename '_svalpha.out or ' filename '.model.txt']);
end

if RC
    x = encodekmers(l, k);
else
    x = encodekmers_norc(l, k);
end  

[comb,~,~,~,~,rcnum] = genIndex(l,k);
lcnum = length(comb);
n = length(alpha);
disp(['# of support vectors: ' num2str(n)])
alphasum = sum(alpha(alpha>0));
gkmc = zeros(4^k*lcnum,1);
GCpos = 0;
GCneg = 0;
np = sum(alpha>0);
nn = sum(alpha<0);
len = zeros(n,1);
mat = zeros(4);
mat2 = zeros(4);
pow = 4.^(0:l-1).';
n4 = 4^(l-1);
slen = numel(sequences);
per = 10;
a = 2;
s = '';
l2 = l-1;
lcnum2 = lcnum-rcnum;
Ker = zeros(l+1,1);
for i = 0:(l-k)
    Ker(l-i+1) = nchoosek(l-i,l-k-i);
end
for i = 1:n
    if mod(i, floor(n/10))==0
        fprintf('%d...', per);
        per = per+10;
    end
    if ver == 0
        while ~strcmp('>', sequences{a}(1))
            s = [s sequences{a}];
            a = a+1;
            if a > slen 
                break
            end
        end
    else
        s = sequences{i};
    end
    ss = letterconvert(s);
    if RC
        temp = zeros(4^k*(lcnum*2), 1);
    else
        temp = zeros(4^k*(lcnum), 1);
    end
    vec = ss(1:l)*pow;
    en = x(vec+1, :);
    temp(en) = temp(en) + 1;
    for j = 2:length(ss)-l2;
        vec = (vec-ss(j-1))/4+n4*ss(j+l2);
        en = x(vec+1, :);
        temp(en) = temp(en) + 1;
    end
    if RC
        temp = temp(1:4^k*lcnum)+temp(4^k*lcnum+1:4^k*2*lcnum);
    end
    if rcnum > 0 && RC
        temp(4^k*lcnum2+1:4^k*lcnum)=temp(4^k*lcnum2+1:4^k*lcnum)/sqrt(2);
    end
    if lk == 1
        gkmc = gkmc + alpha(i)/sqrt(temp'*temp)*temp;
        len(i) = length(ss);
    else
        C = cell(4,2);
        len(i) = length(ss);
        for ii = 1:4
            for j = 1:2
                C{ii,j} = zeros(len(i)-l+1, l);
            end
        end
        for ii = 1:len(i)
            for j = 0:l-1
                if ii - j <= len(i)-l+1 && ii - j > 0
                    C{ss(ii)+1,1}(ii-j,j+1) = 1;
                end
            end
        end
        if RC
            for ii = 0:3
                C{ii+1,2} = rot90(C{4-ii,1}, 2);
            end
            M2 = zeros(len(i)-l+1);
        end
        M = zeros(len(i)-l+1);
        for ii = 1:4
            M = M + C{ii,1}*C{ii,1}.';
            if RC
                M2 = M2 + C{ii,1}*C{ii,2}.';
            end
        end
        if RC
            norm = sqrt(sum(sum(Ker(1 + M)+Ker(1 + M2))));
        else
            norm = sqrt(sum(sum(Ker(1 + M))));
        end
        gkmc = gkmc + alpha(i)/norm*temp;
    end
    if alpha(i) > 0
        GCpos = GCpos + alpha(i)*(sum(ss==1)+sum(ss==2));
        ss = ss+1;
        for j = 1:length(ss)-1
            mat2(ss(j),ss(j+1)) = mat2(ss(j),ss(j+1)) + abs(alpha(i));
        end
    else
        GCneg = GCneg + alpha(i)*(sum(ss==1)+sum(ss==2));
        ss = ss+1;
        for j = 1:length(ss)-1
            mat(ss(j),ss(j+1)) = mat(ss(j),ss(j+1)) + abs(alpha(i));
        end
    end
    s = '';
    a = a+1;
end
fprintf('\n')
GCpos = GCpos/alphasum/mean(len(1:np));;
GCneg = -1*GCneg/alphasum/mean(len(np+1:end));
if RC
    mat = (mat+rot90(rot90(mat,3)'))/2;
    mat2 = (mat2+rot90(rot90(mat2,3)'))/2;
end
mat = mat./repmat(sum(mat,2),1,4);
mat2 = mat2./repmat(sum(mat2,2),1,4);


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

function mat = encodekmers(l,k);
c = genIndex(l,k);
lcnum = length(c);
seqvec = zeros(4^l, l);
vec = (1:4^l)'-1;
for i = 1:l
    seqvec(:,i) = mod(floor(vec/4^(i-1)), 4);
end
seqvec2 = 3-fliplr(seqvec);
mat = zeros(4^l, 2*lcnum);
pow = 4.^(0:k-1)';
for i = 1:lcnum
    mat(:,i) = seqvec(:,c(i,:))*pow+4^k*(i-1)+1;
    mat(:,i+lcnum) =  seqvec2(:,c(i,:))*pow+4^k*(lcnum+i-1)+1;
end

function mat = encodekmers_norc(l,k);
c = genIndex(l,k);
lcnum = length(c);
seqvec = zeros(4^l, l);
vec = (1:4^l)'-1;
for i = 1:l
    seqvec(:,i) = mod(floor(vec/4^(i-1)), 4);
end
mat = zeros(4^l, lcnum);
pow = 4.^(0:k-1)';
for i = 1:lcnum
    mat(:,i) = seqvec(:,c(i,:))*pow+4^k*(i-1)+1;
end
