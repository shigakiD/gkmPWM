function mapTF(fn, mfn, rnum, GC, l_svm, k_svm, ofn)
%usage: mapTF(fn, mfn, rnum, GC)
%fn: fasta file used for seq2prob
%mfn: file name for motifs from process motifs
%rnum: sequence number in fasta file
%GC: GC content as a fraction between 0 to 1
%ofn: output file header
fid = fopen(fn, 'r');
for i = 1:rnum
    h = fgetl(fid);
    s = fgetl(fid);
end
fclose(fid);
GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2];
ss = letterconvert(s);
L = length(ss)-l_svm+1;
ACGT = 'ACGT';
b = 1;
[p,names,len] = getMOTIF(mfn);
a = numel(len);
PWM = cell(a,1);
PWM2 = cell(a,1);
pwm = cell(a,1);
lpwm = cell(a,1);
lab = zeros(a,1);
LEN = zeros(a,1);
LEN_2 = zeros(a,1);
shift = zeros(a,1);
for i = 1:a
    shift(i) = max([l_svm-len(i) 4]);
    PWM{i} = [repmat(GCmat,shift(i), 1); p{i} ;repmat(GCmat,shift(i), 1)];
    PWM2{i} = p{i};
    LEN_2(i) = len(i);
    LEN(i) = len(i)+2*shift(i)-l_svm+1;
    for j = 1:LEN(i)
        pwm{b} = PWM{i}(j:j+l_svm-1,:);
        lpwm{b} = log((pwm{b}+10^-10)/(1+4*10^-10));
        lab(b) = i;
        b = b+1;
    end
end
lab = [lab; 0];
c = combnk(1:l_svm,k_svm);
pwm_prob = zeros(b-1,L);
vec = zeros(l_svm,1);
LL = nchoosek(l_svm,k_svm);
seqmat = zeros(LL,l_svm);
for i = 1:LL
    for j = 1:k_svm
        seqmat(i,c(i,j)) = 1;
    end
end
for i =1:L
    seq = ss(i:i+l_svm-1);
    for j = 1:b-1
        for jj = 1:l_svm
            vec(jj) = lpwm{j}(jj,seq(jj));
        end
        pwm_prob(j,i) = sum(exp(seqmat*vec));
    end
end
for j = 1:b-1
    %vec = max(lpwm{j}');
    %vec2 = min(lpwm{j}');
    %maxnorm = sum(exp(seqmat*vec'));
    %minnorm = sum(exp(seqmat*vec2'));
    %pwm_prob(j,:) = (pwm_prob(j,:)-minnorm)/(maxnorm-minnorm);
    vec3 = mean(lpwm{j}');
    midnorm = sum(exp(seqmat*vec3'));
    pwm_prob(j,:) = (pwm_prob(j,:).^2)./(pwm_prob(j,:)+midnorm);
end
pwm_prob = log(pwm_prob);
n = sum(LEN)+1;
mat = zeros(n,L)-Inf;
ind = zeros(n,L);
TRANS = zeros(n);
LEN = [0;LEN];
C = cumsum(LEN);

gkmprob = dlmread([ofn '_' num2str(rnum) '_prob.out']);
gkmprob(gkmprob>0.99) = 0.99;
pos = gkmprob;
neg = log(1-gkmprob);
for i = 1:a
    for j = 1:LEN(i+1)-1
        TRANS(C(i)+j,C(i)+j+1)=1;
    end
    TRANS(C(i+1),C(1:end-1)+1) = 1;
    TRANS(C(i+1),n)=1;
end
TRANS(n,C(1:end-1)+1) = 1;
TRANS(n,n) = 1;
TRANS = log(TRANS);
for i = 1:a
    mat(C(i)+1,1) = pwm_prob(C(i)+1,1)+pos(1);
end
mat(n,1) = neg(1);
for i = 2:L
    for j = 1:a
        [mat(C(j)+1,i),ind(C(j)+1,i)] = max(mat(1:end,i-1)+TRANS(1:end,C(j)+1)+pwm_prob(C(j)+1,i));
        for jj = 1:LEN(j+1)-1
            mat(C(j)+jj+1,i) = mat(C(j)+jj,i-1)+TRANS(C(j)+jj,C(j)+jj+1)+pwm_prob(C(j)+jj+1,i);
            ind(C(j)+jj+1,i) = C(j)+jj;
        end
    end
    [mat(n,i),ind(n,i)] = max(mat(1:end,i-1)+TRANS(1:end,end)+neg(i-1));
end
path = zeros(L,1);
path(end) = n;
for i = fliplr(1:L-1)
    path(i) = ind(path(i+1),i+1);
end
PATH = zeros(length(s),1);
C_2 = cumsum([0;LEN_2]);
for i = 1:a
    f = find(path==C(i)+1);
    if ~isempty(f)
        for j = 1:length(f)
            vec = f(j)+shift(i):f(j)+shift(i)+LEN_2(i)-1;
            PATH(vec) = C_2(i)+1:C_2(i)+LEN_2(i);
        end
    end
end
N =sum(LEN_2)+1;
PATH(PATH==0) = N;
MAT = [];
lab = [];
for i = 1:a
    MAT = [MAT;PWM2{i}];
    lab = [lab;ones(LEN_2(i),1)*i];
end
lab = [lab;0];
MAT = [MAT;ones(1,4)/4];
n = sum(LEN_2)+1;
PWM = ones(L,4)/4;
PWM_alt = ones(L,4)/4;

GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2];

I = lab(PATH);
Cmat = [I PATH ones(length(I),4)/4];
for i = 1:L+l_svm-1
    if PATH(i)~=n
        Cmat(i,3:6) = MAT(PATH(i),:);
    else    
        Cmat(i,2) = 0;
    end
end

LEN(1) = [];
L = [];
for i = 1:length(PWM2)
    f = find(path==C(i)+1);
    F = [];
    a = 0;
    if ~isempty(f)
        for j = 1:length(f)
            if f(j)+shift(i)+LEN_2(i)< length(ss)-l_svm+1 && f(j)+shift(i) > l_svm-1
                a = a+1;
                F = [F; f(j)+shift(i) f(j)+shift(i)+LEN_2(i)];
            end
        end
        if a > 0
            L = [L;F ones(a,1)*i];
        end
    end
end
fid = fopen([ofn  '_' num2str(rnum) '_kmer_PWM_locs.out'],'w');
if ~isempty(L)
    [~,b] = sort(L(:,1));
    L = L(b,:);
    for i = 1:numel(L)/3
        if sum(ss(L(i,1):L(i,2)-1)) > L(i,2)-L(i,1) && sum(ss(L(i,1):L(i,2)-1)) < 4*(L(i,2)-L(i,1)-1)
            fprintf(fid, '%s\t%d\t%d\t%d\t%f\n', names{L(i,3)}, L(i,3), L(i,1), L(i,2)-1,0);
        end
    end
end
fclose(fid);
kmer =  scoreseqkmer(fn, rnum, l_svm,k_svm, ofn,Cmat);
omat = PWM_corr(ofn, rnum, mfn,kmer,Cmat);
kmer = scoreseqkmer(fn, rnum, l_svm, k_svm,ofn,omat);
PWM_corr2(ofn, rnum, kmer)

function omat = scoreseqkmer(fn, NUM, l,k, ofn,mat);
%fn: fasta file
%NUM: sequence number
%l,k: gapped kmer parameters
%ofn: output header
seqs = importdata(fn);
seq = seqs{2*NUM};
ss = letterconvert(seq);
L = length(ss);
dsvm = dlmread([ofn '_' num2str(NUM) '_dsvm.out']);
nvar = numel(dsvm);
loc = zeros(nvar, 1);
varvec = zeros(nvar, 1);
varc = cell(4,1);
varc{1} = [2 3 4];
varc{2} = [1 3 4];
varc{3} = [1 2 4];
varc{4} = [1 2 3];
for i = 1:nvar
    loc(i) = floor((i-1)/3)+l;
    varvec(i) = varc{ss(loc(i))}(mod(i-1,3)+1);
end
lab = mat(:,1);
lab = [lab; zeros(l,1)];
path_ref = mat(:,2);
path_ref = [path_ref; zeros(l,1)];
PWM = mat(:,3:6);
c = combnk(1:l,k);
[c1,~] = size(c);
O = ones(1,k);
evec = zeros(nvar, 1);
f = find(lab==0);
GCvec = PWM(f(1),:);
for i = 1:nvar
    if lab(loc(i)) > 0
        vec = (max([1 loc(i)-l+1]):min([L loc(i)+l-1]));
        LAB = lab(loc(i));
        path_vec = path_ref(vec);
        path_loc = path_ref(loc(i));
        scores = zeros(length(vec),1);
        V = PWM(loc(i),varvec(i))-PWM(loc(i), ss(loc(i)));
        for j = 1:length(vec);
            if lab(vec(j)) == LAB && path_loc-path_vec(j)==loc(i)-vec(j)
                scores(j) = PWM(vec(j),ss(vec(j)));
            else
                scores(j) = GCvec(ss(vec(j)));
            end
        end
        if loc(i)-l+1 < 1
            a = l-loc(i);
            scores(loc(i)) = V;
        else
            a = 0;
            scores(l) = V;
        end

        if loc(i)+l-1 > L
            b = l-1+L-loc(i);
        else
            b = l-1;
        end
        for j = a:b
            for jj = 1:c1
                if ismember(l-j , c(jj,:));
                    evec(i) = evec(i) + prod(scores(c(jj,:)+j-a));
                end
            end
        end
    end
end
u = unique(lab);
for i = 1:numel(u)
    if u(i) ~= 0
        d = diff(lab==u(i));
        f = find(d==1);
        F = find(d==-1);
        if numel(f) ~= numel(F)
            F = [F;numel(dsvm)];
        end
        for j = 1:numel(f)
            a = [];
            for jj = f(j)+1:F(j)
                a = [a;find(loc==jj)];
            end
            evec(a) = evec(a)*(evec(a)'*dsvm(a))/(evec(a)'*evec(a));
        end
    end
end
omat = [loc varvec evec];


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

function omat = PWM_corr(fn,NUM,mfn,kmer,omat)
dsvm = dlmread([fn '_' num2str(NUM) '_dsvm.out'],'\t');
[coords, ind] = sort(kmer(:,1));
kmer = kmer(ind, :);
dsvm = dsvm(ind, :);
fid = fopen([fn '_' num2str(NUM) '_kmer_PWM_locs.out'],'r');
X = textscan(fid, '%s\t%f\t%f\t%f\t%f\n');
fclose(fid);
L = length(X{2});
fid1 = fopen([fn '_' num2str(NUM) '_kmer_PWM_locs.out'],'w');
ind = ones(L,1);
for i = 1:L
    F = find(coords > X{3}(i));
    FF = find(coords <= X{4}(i));
    f = intersect(F,FF);
    a1 = ip(kmer(f,end),dsvm(f,end));
    if a1 > 0.4
        fprintf(fid1, '%s\t%d\t%d\t%d\t%f\n', X{1}{i}, X{2}(i), X{3}(i), X{4}(i),a1);
    else
        ind(i) = 0;
        omat(X{3}(i):X{4}(i),1:2) = 0;
        omat(X{3}(i):X{4}(i),3:6) = 0.25;
    end
end
fclose(fid1);
p = getMOTIF(mfn);
LEN = zeros(length(p),1);
for i = 1:length(p)
    [LEN(i),~] = size(p);
end
C = cumsum([0;LEN]);
mat = [zeros(length(omat),2) ones(length(omat),4)/4];
for i = 1:L
    if ind(i) == 1
        num = X{2}(i);
        a = 1;
        for j = X{3}(i):X{4}(i)
            if mat(j,1) == 0
                mat(j,1) = num;
                mat(j,2) = C(num)+1;
                mat(j,3:6) = p{num}(a,:);
            else
                mat(j,3:6) = (mat(j,3:6)+p{num}(a,:))/2;
            end
            a = 1+a;
        end
    end 
end 
    
fid = fopen([fn '_' num2str(NUM) '_kmer_pwm.txt'],'w');
fprintf(fid, ['# one hot encode \npos\tA\tC\tG\tT\n']);
for j = 1:length(mat)
    fprintf(fid, '%d\t%f\t%f\t%f\t%f\n', j-1, mat(j,3:6));
end
fclose(fid);

function c = ip(x,y);
c = x'*y/sqrt(x'*x)/sqrt(y'*y);

function [mat, names, len] = getMOTIF(fn)
mat = {};
names = {};
len = [];
fid = fopen(fn);
if fid < 0
    fprintf('No good')
else
   i=0;
   while ~feof(fid)
        line = fgetl(fid);
        if length(line) >= 5
            if strcmp(line(1:5), 'MOTIF')
                i = i+1;
                names{i} = fgetl(fid);
                len = [len;str2double(fgetl(fid))];
                mat{i} = zeros(len(i),4);
                for j = 1:len(i)
                    line = fgetl(fid);
                    mat{i}(j,:) = str2num(line);
                end
            end
        end
    end
end
fclose(fid);

function PWM_corr2(fn,NUM,kmer)
dsvm = dlmread([fn '_' num2str(NUM) '_dsvm.out']);
[coords, ind] = sort(kmer(:,1));
kmer = kmer(ind, :);
dsvm = dsvm(ind, :);
fid = fopen([fn '_' num2str(NUM) '_kmer_PWM_locs.out'],'r');
X = textscan(fid, '%s\t%f\t%f\t%f\t%f\n');
fclose(fid);
L = length(X{2});
fid1 = fopen([fn '_' num2str(NUM) '_kmer_PWM_locs.out'],'w');
ind = ones(L,1);
for i = 1:L
    F = find(coords > X{3}(i));
    FF = find(coords <= X{4}(i));
    f = intersect(F,FF);
    a1 = ip(kmer(f,end),dsvm(f,end));
    fprintf(fid1, '%s\t%d\t%d\t%d\t%f\n', X{1}{i}, X{2}(i), X{3}(i), X{4}(i), a1);
end
fclose(fid1);

