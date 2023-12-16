function mapTF2(varargin)
% mapTF maps the TFBS motifs found from gkmPWM and gkmPWMlasso to regions at
%     base-pair resolution
% 
%     mapTF(seqfile, wfile, gkmPWMmemefile, gkmPWMlassofile, memefile, outputprefix...)
% 
%     Requires outputs of gkmPWM and gkmPWMlasso for a given model. Any combination
%     of parameters for gkmPWM and gkmPWMlasso will work.
% 
%     Positional Parameters (Required):
% 
%     seqfile         The set of sequences to which the motifs will be mapped
%                     (fasta format)
%     wfile           Two column tab delimited file of kmer weights.  See the
%                     README.md for detail on how to generate this file.
%     gkmPWMmemefile  The gkmPWM meme output file
%     gkmPWMlassofile The gkmPWMlasso output file
%     memefile        The collection of PWMs in meme format
%     outputprefix    The prefix of the output files.  
% 
%     Name Value Pair Parameters (Optional):
% 
%     'l'             The full length of the gapped k-mer.  This NEEDS to be 
%                     the same as the l in the gkmSVM model (default: 11)
%     'k'             The number of ungapped positions of the gapped k-mer.
%                     This NEEDS to be the same as the k in the gkmSVM model
%                     (default: 7)
%     'KmerFrac'      Set the fraction of the total number of gapped k-mers to
%                     use with mapTF.  This reduces the memory and runtime
%                     needed.  If the total number of gapped k-mers is too high
%                     with the given combination of (l,k,KmerFrac), KmerFrac will
%                     be automatically set to a lower value to create a more 
%                     workable number of gapped k-mers
% 
%     Outputs 2 files named:
%     outputprefix_kmer_PWM_locs.out
%         Contains the location of the motifs from gkmPWM and gkmPWMlasso.  
%         See the README.md for details on the format of the output
%     outputprefix_motifs.out
%         A file containing the PWMs that were mapped.  This is NOT in meme
%         format and is made to work with the python script "mapTF_profile.py"
%
%     Example (files in the example_files directory):
%     mapTF('GM12878.fa', 'GM12878_weights.out','GM12878_10_6_0_15_denovo.meme',...
%         'GM12878_10_6_30_gkmPWMlasso.out','combined_db_v4.meme', 'GM12878',...
%         'l', 11, 'k', 7,'KmerFrac', 1)
%         Outputs GM12878_kmer_PWM_locs.oout and GM12878_motif.out
if nargin < 6
    error('Need at least 6 inputs')
end

fn = varargin{1};
wfn = varargin{2};
mfn1 = varargin{3};
mfn2 = varargin{4};
memefn = varargin{5};
ofn = varargin{6};
l_svm = 11;
k_svm = 7;
nfrac = 1;
if nargin > 6
    if mod(nargin,2) == 1
        error('Incorrect number of inputs')
    end
    vec = 7:2:nargin;
    inputlib = {'l', 'k','KmerFrac'};
    for i = 1:length(vec)
        f = strcmp(varargin{vec(i)},inputlib);
        if sum(f) == 0
            error([varargin{vec(i)} ' is not an input option'])
        end
    end
    f = find(strcmp('l', varargin));
    if ~isempty(f);
        l_svm = varargin{f+1};
    end
    f = find(strcmp('k', varargin));
    if ~isempty(f);
        k_svm = varargin{f+1};
    end
    f = find(strcmp('KmerFrac', varargin));
    if ~isempty(f)
        nfrac = varargin{f+1};        
    end
    if l_svm < k_svm
        error('7th argument must be greater or equal to the 8th argument')
    end
end

disp('Processing Motifs')
process_motifs(mfn1, mfn2, memefn, ofn)
mfn = [ofn '_motifs.out'];
[P,V, seqindmat, ss, seq] = seq2pv(fn, wfn,l_svm);
GC = countGC(ss);
GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2];
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
[c,~,~,~,~,rcnum] = genIndex(l_svm,k_svm,nfrac);
c2 = c(1:numel(c)/k_svm-rcnum,:);
c = [c;l_svm+1-fliplr(c2)];
C = numel(c)/k_svm; 
seqmat = zeros(C,l_svm);
for i = 1:C
    seqmat(i,c(i,:)) = 1;
end
Smat = cell(l_svm,1);
for i = 1:l_svm
    f = find(prod(c-i,2)==0);
    Smat{i} = zeros(length(f), l_svm);
    for j = 1:length(f)
        Smat{i}(j,c(f(j),:)) = 1;
    end
end
L = length(ss);
B = b-1;
maxnorm = zeros(B,1);
minnorm = zeros(B,1);
vec = zeros(l_svm,1);
IND = zeros(4^l_svm,1);
kmat = zeros(B,4^l_svm);
for j = 1:B
    vec = max(lpwm{j}');
    vec2 = min(lpwm{j}');
    maxnorm(j) = sum(exp(seqmat*vec'));
    minnorm(j) = sum(exp(seqmat*vec2'));
end
dnorm = maxnorm-minnorm;
vec = zeros(l_svm,1);
disp('Mapping motifs')
tic
for I = 1:length(ss)
    seq2 = ss{I};
    for i = 1:length(seqindmat{I})
        ind = seqindmat{I}(i);
        if IND(ind) == 0
            IND(ind) = 1;
            SEQ = seq2(i:i+l_svm-1);
            for j = 1:B
                for jj = 1:l_svm
                    vec(jj) = lpwm{j}(jj,SEQ(jj));
                end
                kmat(j,ind) = sum(exp(seqmat*vec));
            end
            kmat(:,ind) = log((kmat(:,ind)-minnorm)./dnorm);
            pwm_prob(:,i) = kmat(:,ind);
        else
            pwm_prob(:,i) = kmat(:,ind);
        end
    end
    MAPTF(fn, mfn, seq{I*2}, ss{I}, GC, pwm_prob, Smat, l_svm, k_svm, ofn, PWM, PWM2, pwm, lpwm, lab, LEN, LEN_2, shift, P{I}, V{I}, names, len, a,b, I)
    if mod(I,100)==0
        fprintf('%d out of %d done...\n', I, length(ss));
        toc
    end
end

function MAPTF(fn, mfn, s,ss, GC, pwm_prob, Smat, l_svm, k_svm, ofn, PWM, PWM2, pwm, lpwm, lab, LEN, LEN_2, shift, gkmprob, dsvm, names, len, a, b,rnum)
GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2];
L = length(ss)-l_svm+1;
ACGT = 'ACGT';
lab = [lab; 0];
n = sum(LEN)+1;
mat = zeros(n,L)-Inf;
ind = zeros(n,L);
TRANS = zeros(n);
LEN = [0;LEN];
C = cumsum(LEN);

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
Cmat(:,3:6) = (Cmat(:,3:6)+0.0001)/(1.0001);
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
                F = [F; f(j)+shift(i) f(j)+shift(i)+LEN_2(i)-1 mean(gkmprob(f(j):f(j)+2*shift(i)+LEN_2(i)-l_svm))];
            end
        end
        if a > 0
            L = [L;F ones(a,1)*i];
        end
    end
end
NAME = cell(numel(L)/4,1);
Lmat = zeros(numel(L)/4,4);
if ~isempty(L)
    [~,b] = sort(L(:,1));
    L = L(b,:);
    for i = 1:numel(L)/4
        NAME{i} = names{L(i,4)};
        Lmat(i,:) = [L(i,4) L(i,1) L(i,2) L(i,3)];
    end
end
kmer = scoreseqkmer(fn, rnum, ss, Smat, l_svm,k_svm, ofn,Cmat, dsvm);
[omat, NAME, Lmat] = PWM_corr(ofn, rnum, mfn,kmer,Cmat, dsvm, NAME, Lmat,PWM2);
kmer = scoreseqkmer(fn, rnum, ss, Smat, l_svm, k_svm,ofn,omat, dsvm);
PWM_corr2(ofn, rnum, kmer, dsvm, NAME, Lmat, s)

function omat = scoreseqkmer(fn, NUM, ss, Smat, l,k, ofn,mat, dsvm);
%fn: fasta file
%NUM: sequence number
%l,k: gapped kmer parameters
%ofn: output prefix
L = length(ss);
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
PWM = log(mat(:,3:6));
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
        V = exp(PWM(loc(i),varvec(i)))-exp(PWM(loc(i), ss(loc(i))));
        for j = 1:length(vec);
            if lab(vec(j)) == LAB && path_loc-path_vec(j)==loc(i)-vec(j)
                scores(j) = PWM(vec(j),ss(vec(j)));
            else
                scores(j) = GCvec(ss(vec(j)));
            end
        end
        if loc(i)-l+1 < 1
            a = l-loc(i);
            scores(loc(i)) = 0;
        else
            a = 0;
            scores(l) = 0;
        end

        if loc(i)+l-1 > L
            b = l-1+L-loc(i);
        else
            b = l-1;
        end
        for j = a:b
            evec(i) = evec(i)+sum(exp(Smat{l-j}*scores(1+j-a:l+j-a)));
            %for jj = 1:c1
            %    if ismember(l-j , c(jj,:));
            %        evec(i) = evec(i) + prod(scores(c(jj,:)+j-a));
            %    end
            %end
        end
        evec(i) = evec(i)*V;
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
        en(i) = 0;
    elseif strcmp(s(i),'C') || strcmp(s(i),'c')
        en(i) = 1;
    elseif strcmp(s(i),'G') || strcmp(s(i),'g')
        en(i) = 2;
    else
        en(i) = 3;
    end
end

function [omat, NAME2, Lmat2] = PWM_corr(fn,NUM,mfn,kmer,omat, dsvm, NAME, Lmat, p)
[coords, ind] = sort(kmer(:,1));
kmer = kmer(ind, :);
dsvm = dsvm(ind, :);
L = length(NAME);
ind = ones(L,1);
n = [];
for i = 1:L
    F = find(coords > Lmat(i,2));
    FF = find(coords <= Lmat(i,3));
    f = intersect(F,FF);
    a1 = ip(kmer(f,end),dsvm(f,end));
    if a1 > 0.6
        n = [n;i];
    else
        ind(i) = 0;
        omat(Lmat(i,2):Lmat(i,3),1:2) = 0;
        omat(Lmat(i,2):Lmat(i,3),3:6) = 0.25;
    end
end
Lmat2 = Lmat(n,:);
NAME2 = NAME(n);
LEN = zeros(length(p),1);
for i = 1:length(p)
    [LEN(i),~] = size(p);
end
C = cumsum([0;LEN]);
mat = [zeros(length(omat),2) ones(length(omat),4)/4];
for i = 1:L-length(n)
    if ind(i) == 1
        num = Lmat(i,1);
        a = 1;
        for j = Lmat(i,2):Lmat(i,3)
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

function PWM_corr2(fn,NUM,kmer, dsvm, NAME, Lmat, s)
[coords, ind] = sort(kmer(:,1));
kmer = kmer(ind, :);
dsvm = dsvm(ind, :);
L = length(NAME);
if NUM ==1
    fid1 = fopen([fn '_kmer_PWM_locs.out'],'w');
else
    fid1 = fopen([fn '_kmer_PWM_locs.out'],'a+');
end
ind = ones(L,1);
for i = 1:L
    F = find(coords > Lmat(i,2));
    FF = find(coords <= Lmat(i,3));
    f = intersect(F,FF);
    a1 = ip(kmer(f,end),dsvm(f,end));
    fprintf(fid1, '%d\t%s\t%d\t%d\t%d\t%f\t%f\t%s\n', NUM, NAME{i}, Lmat(i,1), Lmat(i,2), Lmat(i,3), Lmat(i,4), a1, s(Lmat(i,2):Lmat(i,3)));
end
fclose(fid1);

function [P,V,seqindmat,seqout,seq] = seq2pv(sfn, wfn, l_svm)
%sfn: fasta file
%wfn: kmer weight file
%ofn: output prefix
fid = fopen(wfn, 'r');
X = textscan(fid,'%s\t%f\n');
fclose(fid);
l = length(X{1}{1});
if l_svm ~= 0
    error('l is not the same length as the kmers in the weight file')
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
W = (1/2)*(1+erf((w-m)/s/sqrt(2)));
seq = importdata(sfn);
n = length(seq)/2;
seqout = cell(n,1);
seqindmat = cell(n,1);
disp('converting kmers to  probabilities')
P = cell(n,1);
for i = 1:n
    if mod(i,1000)==0
        disp([num2str(i) ' sequences converted'])
    end
    L = length(seq{2*i})-l+1;
    seqindmat{i} = zeros(L,1);
    ss = letterconvert(seq{2*i});
    seqout{i} = ss+1;
    p = zeros(L,1);
    I = ss(1:l)*pow;
    seqindmat{i}(1) = I+1;
    p(1) = W(I+1);
    for j = 2:L
        I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1);
        seqindmat{i}(j) = I+1;
        p(j) = W(I+1);
    end
    P{i} = p;
end
disp('running dsvm')
mat = [1 2 3;0 2 3;0 1 3;0 1 2];
O = ones(1,l);
V = cell(n,1);
for i = 1:n
    if mod(i,1000)==0
        disp([num2str(i) ' sequences converted'])
    end
    L = length(seq{2*i})-2*l+2;
    ss = letterconvert(seq{2*i});
    p = zeros(L+l-1,1);
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
    V{i} = v;
end

function GC = countGC(s)
GC = 0;
L = 0;
for i = 1:length(s)
    S = s{i};
    N = length(S);
    L = N + L;
    for j = 1:N
        if S(j) == 2 || S(j) == 3
            GC = GC + 1;
        end
    end
end
GC = GC/L;

function process_motifs(dfn, lfn, memefn, ofn)
%dfn: file name for denovo motifs
%lfn: file name for lasso motifs
%memefn: file name for the meme input for gkmPWMlasso
%ofn: output filename
a = 1;
LEN = zeros(1,1);
shift = zeros(1,1);
[p,w] = getdenovomotif(dfn);
N = numel(w);
fid = fopen(lfn,'r');
X = textscan(fid,'%f\t%f\t%s\t%f\t%f\t%f\n','delimiter', '\t', 'headerlines', 4);
fclose(fid);
n = X{1}(end);
vec = zeros(n,1);
vec2 = zeros(n,1);
w = [w; zeros(n,1)];
for i = 1:n
    f = find(X{1}==i);
    vec(i) = X{2}(f(1));
    vec2(i) = X{5}(f(1));
    w(i+N) = X{4}(f(1));
end
P = getmotif(memefn,vec);
[pp, info, len] = trim_pwm([p;P],0.25);
PWM2 = cell(1,1);
LEN_2 = zeros(1,1);
I_2 = zeros(1,1);
hal = true;
for ii = 1:length(w)
    if ii > N
        hal = true;
        [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]);
        if cor > 0.7 || vec2(ii-N) < 1.5
            hal = false;
        end
    elseif ii > 1 && a > 2
        hal = true;
        [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]);
        if cor > 0.7
            hal = false;
        end
    end
    if hal && w(ii) > 0 && len(ii) >= 6
        if len(ii) > 10
            if info(ii)/len(ii) > 0.7
                PWM2{a} = pp{ii};
                LEN_2(a) = len(ii);
                I_2(a) = info(ii);
                a = a+1;
                PWM2{a} = rot90(PWM2{a-1},2);
                LEN_2(a) = LEN_2(a-1);
                I_2(a) = I_2(a-1);
                a = a+1;
            end
        elseif info(ii) > 6 || info(ii)/len(ii) > 1
            PWM2{a} = pp{ii};
            LEN_2(a) = len(ii);
            I_2(a) = info(ii);
            a = a+1;
            PWM2{a} = rot90(PWM2{a-1},2);
            LEN_2(a) = LEN_2(a-1);
            I_2(a) = I_2(a-1);
            a = a+1;
        end
    end
end
num2 = length(strfind(fileread(memefn),'MOTIF'));
[p,names] = getmotif(memefn,1:num2);
[p,info,lenvec] = trim_pwm(p,0.25);
fid = fopen([ofn '_motifs.out'], 'w');
a = 1;
for i = 1:length(PWM2)
    [ind, r] = ppmsim([PWM2{i};p], [LEN_2(i);lenvec]);
    if r > 0.80
        fprintf(fid,'MOTIF %d\n%s\n%d\n', a, names{ind},LEN_2(i));
        a = a+1;
        for j = 1:LEN_2(i)
            fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4));
        end
        fprintf(fid, '\n');
    elseif I_2(i)/LEN_2(i) > 1
        fprintf(fid,'MOTIF %d\n%s\n%d\n', a, consen(PWM2{i}, LEN_2(i)), LEN_2(i));
        a = a+1;
        for j = 1:LEN_2(i)
            fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4));
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);

function [ind, M] = ppmsim(mot,lenvec)
n = length(lenvec)-1;
simmat = ones(n-1,1);
for i = 1:n+1
    mot{i} = mot{i}-1/4;
    mot{i} = mot{i}/sqrt(sum(sum(mot{i}.^2)));
end
M = 0;
ind = 1;
for j = 2:n+1
    mat = mot{1}*mot{j}';
    rmat = rot90(mot{1},2)*mot{j}';
    MM = max([sum(spdiags(mat)) sum(spdiags(rmat))]);
    if MM > M
        M = MM;
        ind = j-1;
    end
end

function [mat,w] = getdenovomotif(filename)
% filename is the meme file that contains the motifs;
% n is the nth motif in the file
mat = {};
w = [];
fid = fopen(filename);
if fid < 0
    fprintf('No good')
else
   i=0;
   while ~feof(fid)
        line = fgetl(fid);
        if length(line) >= 5
            if strcmp(line(1:5), 'MOTIF')
                i = i+1;
                line = fgetl(fid);
                mat{i} = [];
                [~, tmp] = strtok(line);
                [tmp] = strtok(tmp);
                w = [w;str2double(tmp)];
                while ~isempty(line)
                    mat{i} = [mat{i}; str2num(line)];
                    line = fgetl(fid);
                end
            end
        end
    end
end
mat = mat';
fclose(fid);

function [pp, info, len] = trim_pwm(p,cut)
l = length(p);
info = zeros(l, 1);
len = zeros(l,1);
for i = 1:l
    mat = p{i}+(p{i}==0);
    vec = 2+sum(mat.*log(mat)/log(2),2);
    while (vec(1) < cut || mean(vec(1:3)) < cut || mean(vec(2:4)) < cut) && length(vec) > 4
        p{i}(1,:) = [];
        vec(1) = [];
    end
    while (vec(end) < cut || mean(vec(end-2:end)) < cut || mean(vec(end-3:end-1)) < cut) && length(vec) > 4
        vec(end) = [];
        p{i}(end,:) = [];
    end
    info(i) = sum(vec);
    [len(i), ~] = size(p{i});
end
pp = p;

function s = consen(p, l)
s='';
for i = 1:l
    s(i)='A';
    r = norm(p(i,:)'- [1 0 0 0]');
    rr =  norm(p(i,:)'- [0 1 0 0]');
    if rr < r
        s(i) = 'C';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 0 1 0]');
    if rr < r
        s(i) = 'G';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 0 0 1]');
    if rr < r
        s(i) = 'T';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/2 1/2 0 0]');
    if rr < r
        s(i) = 'M';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/2 0 1/2 0]');
    if rr < r
        s(i) = 'R';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/2 0 0 1/2]');
    if rr < r
        s(i) = 'W';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 1/2 1/2 0]');
    if rr < r
        s(i) = 'S';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 1/2 0 1/2]');
    if rr < r
        s(i) = 'Y';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 0 1/2 1/2]');
    if rr < r
        s(i) = 'K';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/3 1/3 1/3 0]');
    if rr < r
        s(i) = 'V';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/3 1/3 0 1/3]');
    if rr < r
        s(i) = 'H';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/3 0 1/3 1/3]');
    if rr < r
        s(i) = 'D';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 1/3 1/3 1/3]');
    if rr < r
        s(i) = 'B';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0.25 0.25 0.25 0.25]');
    if rr < r
        s(i) = 'N';
        r = rr;
    end
end
