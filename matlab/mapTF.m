function mapTF(varargin)
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
%     'Negative'      Call motifs that are predictive of the negative set.
%                     Output files will have the "negative" added to the end 
%                     e.g. "_negative_motifs.out" (default: false)
%     'KmerFrac'      Set the fraction of the total number of gapped k-mers to
%                     use with mapTF.  This reduces the memory and runtime
%                     needed.  If the total number of gapped k-mers is too high
%                     with the given combination of (l,k,KmerFrac), KmerFrac will
%                     be automatically set to a lower value to create a more 
%                     workable number of gapped k-mers
%     'LS'            Speeds up mapTF for large numbers of sequences by saving 
%                     the PWM probabilities for each k-mer.  Needs around 15GB
%                     of RAM.  Only set this if you have enough RAM.  
%                     (default: false)
%     'PWMcorrcut'    The correlation cutoff to remove redundant motifs in the 
%                     gkmPWM and gkmPWMlasso list.  Motif selection is prioritized
%                     by the regression weight.  (default: 0.80)
%     'dSVMcorrcut'   The correlation cutoff to remove motifs calls that do not fit
%                     the deltaSVM scores for all variants in the predicted TFBS.
%                     (default: 0.60)
%
% 
%     Outputs 2 files named:
%     outputprefix_TFBS_locations.out
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
%         Outputs GM12878_TFBS_locations.out and GM12878_motif.out
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
PWMcorrcut = 0.8;
dsvmcut = 0.6;
Negative = false;
LS = false;
if nargin > 6
    if mod(nargin,2) == 1
        error('Incorrect number of inputs')
    end
    vec = 7:2:nargin;
    inputlib = {'l', 'k','KmerFrac','LS', 'PWMcorrcut', 'dSVMcorrcut', 'Negative'};
    for i = 1:length(vec)
        f = strcmp(varargin{vec(i)},inputlib);
        if sum(f) == 0
            error([varargin{vec(i)} ' is not an input option'])
        end
    end
    f = find(strcmp('l', varargin));
    if ~isempty(f);
        l_svm = varargin{f+1};
        if ~isa(l_svm, 'double') || round(l_svm)-l_svm ~= 0 || l_svm <= 0
            error(['l must be a positive integer'])
        end
    end
    f = find(strcmp('k', varargin));
    if ~isempty(f);
        k_svm = varargin{f+1};
        if ~isa(k_svm, 'double') || round(k_svm)-k_svm ~= 0 || k_svm <= 0 || k_svm > l_svm
            error(['k must be a positive integer less than or equal to l'])
        end
    end
    f = find(strcmp('Negative', varargin));
    if ~isempty(f);
        Negative = varargin{f+1};
        if ~isa(Negative, 'logical')
            error(['Negative must be set to true or false'])
        end
    end
    f = find(strcmp('KmerFrac', varargin));
    if ~isempty(f);
        nfrac = varargin{f+1};
        if ~isa(nfrac, 'double') || nfrac <= 0 || nfrac >1
            error(['KmerFrac must be a positive fraction in (0 1]'])
        end
    end
    f = find(strcmp('PWMcorrcut', varargin));
    if ~isempty(f);
        PWMcorrcut = varargin{f+1};
        if ~isa(PWMcorrcut , 'double') || PWMcorrcut  < -1 || PWMcorrcut  >1
            error(['PWMcorrcut must be a fraction in [-1 1]'])
        end
    end
    f = find(strcmp('dSVMcorrcut', varargin));
    if ~isempty(f);
        dsvmcut = varargin{f+1};
        if ~isa(dsvmcut , 'double') || dsvmcut  < -1 || dsvmcut  >1
            error(['dSVMcorrcut must be a fraction in [-1 1]'])
        end
    end
    f = find(strcmp('LS', varargin));
    if ~isempty(f);
        LS = varargin{f+1};
        if ~isa(LS, 'logical')
            error(['LS must be set to true or false'])
        end
    end
end

if Negative
    ofn = [ofn '_negative'];
end

disp('processing motifs')
process_motifs(mfn1, mfn2, memefn, ofn, PWMcorrcut, Negative)
mfn = [ofn '_motifs.out'];
[P,V, seqindmat, ss, seq] = seq2pv(fn, wfn,l_svm, Negative);
GC = countGC(ss);
GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2];
b = 1;
[p,names,len] = getMOTIF(mfn);
a = numel(len);
PWM = cell(a,1);
PWM2 = cell(a,1);
lPWM2 = cell(a,1);
pwm = cell(sum(len),1);
lpwm = cell(sum(len),1);
lab = zeros(sum(len),1);
LEN = zeros(a,1);
LEN_2 = zeros(a,1);
shift = zeros(a,1);
for i = 1:a
    shift(i) = max([l_svm-len(i) 4]);
    PWM{i} = [repmat(GCmat,shift(i), 1); p{i} ;repmat(GCmat,shift(i), 1)];
    PWM2{i} = [ones(l_svm-1,4)/4; p{i} ;ones(l_svm-1,4)/4];
    lPWM2{i} = log((PWM2{i}+10^-10)/(1+4*10^-10));
    LEN_2(i) = len(i);
    LEN(i) = len(i)+2*shift(i)-l_svm+1;
    for j = 1:LEN(i)
        pwm{b} = PWM{i}(j:j+l_svm-1,:);
        lpwm{b} = log((pwm{b}+10^-10)/(1+4*10^-10));
        lab(b) = i;
        b = b+1;
    end
end
pwm = pwm(1:(b-1));
lpwm = lpwm(1:(b-1));
lab = lab(1:(b-1));
lab2 = lab;
f = (1:(b-1))';
for i = 2:2:a
    ff = find(lab2==i);
    f(ff) = flipud(f(ff));
end
p1 = find(mod(lab,2)==1);
ff = find(mod(lab,2)==0);
p2 = f(ff);
pl = length(p1);
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
if LS
    kmat = zeros(B/2,4^l_svm);
end
for j = 1:B
    vec = max(lpwm{j}');
    vec2 = min(lpwm{j}');
    maxnorm(j) = sum(exp(seqmat*vec'));
    minnorm(j) = sum(exp(seqmat*vec2'));
end
dnorm = maxnorm-minnorm;
vec = zeros(l_svm,1);
LL = cell(length(V),1);
VV = cell(length(V),1);
NN = cell(length(V),1);
disp('mapping motifs')
tic
for I = 1:length(ss)
    seq2 = ss{I};
    pwm_prob = zeros(B,numel(seqindmat{I})/2);
    for i = 1:numel(seqindmat{I})/2
        ind = seqindmat{I}(i,:);
        if LS
            if IND(ind(1)) == 0
                IND(ind(1)) = 1;
                SEQ = seq2(i:i+l_svm-1);
                for j = 1:pl
                    for jj = 1:l_svm
                        vec(jj) = lpwm{p1(j)}(jj,SEQ(jj));
                    end
                    kmat(j,ind(1)) = sum(exp(seqmat*vec));
                end
                kmat(:,ind(1)) = log((kmat(:,ind(1))-minnorm(p1))./dnorm(p1));
            end
            pwm_prob(p1,i) = kmat(:,ind(1));
            if IND(ind(2)) == 0
                IND(ind(2)) = 1;
                SEQ = seq2(i:i+l_svm-1);
                for j = 1:pl
                    for jj = 1:l_svm
                        vec(jj) = lpwm{p2(j)}(jj,SEQ(jj));
                    end
                    kmat(j,ind(2)) = sum(exp(seqmat*vec));
                end
                kmat(:,ind(2)) = log((kmat(:,ind(2))-minnorm(p2))./dnorm(p2));
            end
            pwm_prob(p2,i) = kmat(:,ind(2));
        else
            SEQ = seq2(i:i+l_svm-1);
            for j = 1:B
                for jj = 1:l_svm
                    vec(jj) = lpwm{j}(jj,SEQ(jj));
                end
                pwm_prob(j,i) = log((sum(exp(seqmat*vec))-minnorm(j))/dnorm(j));
            end
        end
    end
    [LL{I}, NN{I}] = MAPTF(fn, ss{I}, pwm_prob, l_svm, k_svm, LEN, LEN_2, shift, P{I}, names, a, b);
    if numel(LL{I}) > 0
        VV{I} = scoreseqkmer(PWM2, lPWM2, LL{I}, ss{I}, Smat, l_svm, k_svm, ofn, V{I});
    end
    if mod(I,100)==0
        fprintf('%d out of %d sequences done...\n', I, length(ss));
        toc
    end
end
fprintf('%d out of %d sequences done...\n', length(ss), length(ss));
clear kmat
PWM_corr(ofn, VV, NN, LL, seq, dsvmcut);

function [Lmat, NAME] = MAPTF(fn, ss, pwm_prob, l_svm, k_svm, LEN, LEN_2, shift, gkmprob, names, a, b)
L = length(ss)-l_svm+1;
n = sum(LEN)+1;
mat = zeros(n,L)-Inf;
ind = zeros(n,L);
LEN = [0;LEN];
C = cumsum(LEN);
D = setdiff(1:C(end),C+1)';
C2 = [C(2:end);n];
pos = log(gkmprob);
neg = log(1-gkmprob);
pwm_prob = pwm_prob+repmat(pos',n-1,1);
for i = 1:a
    mat(C(i)+1,1) = pwm_prob(C(i)+1,1);
end
mat(n,1) = neg(1);
ind(D,:) = repmat(D-1,1,L);
for i = 2:L
    for j = 1:a
        [mat(C(j)+1,i),ind(C(j)+1,i)] = max(mat(C2,i-1)+pwm_prob(C(j)+1,i));
        ind(C(j)+1,i) = C2(ind(C(j)+1,i));
    end
    mat(D,i) = mat(D-1,i-1)+pwm_prob(D,i);
    [mat(n,i),ind(n,i)] = max(mat(C2,i-1)+neg(i));
    ind(n,i) = C2(ind(n,i));
end
path = zeros(L,1);
path(end) = n;
for i = fliplr(1:L-1)
    path(i) = ind(path(i+1),i+1);
end
L2 = [];
for i = 1:length(LEN_2)
    f = find(path==C(i)+1);
    if ~isempty(f)
        for j = 1:length(f)
            if f(j)+shift(i)+LEN_2(i)< length(ss)-l_svm+1 && f(j)+shift(i) > l_svm-1
                [~,ind] = max(gkmprob(f(j):f(j)+2*shift(i)+LEN_2(i)-l_svm));
                if f(j)-1+ind < 3
                    R = 1:5;
                elseif f(j)-1+ind > L-2
                    R = L-4:L;
                else
                    R = f(j)+ind-3:f(j)+ind+1;
                end
                L2 = [L2; f(j)+shift(i) f(j)+shift(i)+LEN_2(i)-1 mean(gkmprob(R)) i];
            end
        end
    end
end
NAME = cell(numel(L2)/4,1);
Lmat = zeros(numel(L2)/4,4);
if ~isempty(L2)
    [~,b] = sort(L2(:,1));
    L2 = L2(b,:);
    for i = 1:numel(L2)/4
        NAME{i} = names{L2(i,4)};
        Lmat(i,:) = [L2(i,4) L2(i,1) L2(i,2) L2(i,3)];
    end
end

function varscore = scoreseqkmer(PWM2, lPWM2, Lmat, ss, Smat, l_svm, k_svm, ofn, dsvm);
varc = [2 3 4; 1 3 4;1 2 4;1 2 3];
O = ones(1,l_svm-1);
n = numel(Lmat)/4;
varscore = zeros(n,1);
for i = 1:n
    M = Lmat(i,1);
    ind = [O ss(Lmat(i,2):Lmat(i,3)) O];
    DSVM = dsvm(Lmat(i,2):Lmat(i,3),:); 
    L = Lmat(i,3)-Lmat(i,2)+1;
    matscore = zeros(L,3);
    for ii = 1:L
        Lind = l_svm-1+ii;
        scores = zeros(2*l_svm-1,1);
        for j = 1:2*l_svm-1;
            scores(j) = lPWM2{M}(Lind-l_svm+j,ind(Lind-l_svm+j));
        end
        evec = 0;
        scores(l_svm) = 0;
        for j = 1:l_svm
            evec = evec+sum(exp(Smat{l_svm-j+1}*scores(j:l_svm+j-1)));
        end
        for iii = 1:3
            V = PWM2{M}(Lind,varc(ind(Lind),iii))-PWM2{M}(Lind,ind(Lind));
            matscore(ii,iii) = evec*V;            
        end
    end
    varscore(i) = ip(matscore(:), DSVM(:));
end

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

function PWM_corr(fn, VV, NN, LL, seq, dsvmcut);
n = length(VV);
fid1 = fopen([fn '_TFBS_locations.out'],'w');
for j = 1:n
    NAME = NN{j};
    Lmat = LL{j};
    varscore = VV{j};
    s = seq{j*2};
    L = length(NAME);
    for i = 1:L
        if varscore(i) > dsvmcut
            fprintf(fid1, '%d\t%s\t%d\t%d\t%d\t%f\t%f\t%s\n', j, NAME{i}, Lmat(i,1), Lmat(i,2), Lmat(i,3), Lmat(i,4), varscore(i), s(Lmat(i,2):Lmat(i,3)));
        end
    end
end
fclose(fid1);

function [P,V,seqindmat,seqout,seq] = seq2pv(sfn, wfn, l_svm, Negative)
%sfn: fasta file
%wfn: kmer weight file
%ofn: output prefix
fid = fopen(wfn, 'r');
X = textscan(fid,'%s\t%f\n');
fclose(fid);
l = length(X{1}{1});
if l_svm ~= l
    error('l is not the same length as the kmers in the weight file')
end
w = zeros(4^l,1);
pow = (4.^(0:(l-1)))';
pow2 = flipud(pow);
disp('calculating indices')
for i = 1:numel(X{2});
    ss = letterconvert(X{1}{i});
    rs = 3-fliplr(ss);
    w(ss*pow+1) = X{2}(i);
    w(rs*pow+1) = X{2}(i);
end
if Negative
    w = -1*w;
    m = -1*mean(X{2});
else
    m = mean(X{2});
end
s = std(X{2});
W = (1/2)*(1+erf((w-m)/s/sqrt(2)));
W(W>0.99) = 0.99;
W(W<0.01) = 0.01;
seq = importdata(sfn);
n = length(seq)/2;
seqout = cell(n,1);
seqindmat = cell(n,1);
disp('converting kmers to probabilities')
P = cell(n,1);
for i = 1:n
    if mod(i,1000)==0
        disp([num2str(i) ' sequences converted'])
    end
    L = length(seq{2*i})-l+1;
    seqindmat{i} = zeros(L,2);
    ss = letterconvert(seq{2*i});
    rs = 3-ss;
    seqout{i} = ss+1;
    p = zeros(L,1);
    I = ss(1:l)*pow;
    I2 = rs(1:l)*pow2;
    seqindmat{i}(1,1) = I+1;
    seqindmat{i}(1,2) = I2+1;
    p(1) = W(I+1);
    for j = 2:L
        I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1);
        I2 = (I2-rs(j-1)*4^(l-1))*4+rs(j+l-1);
        seqindmat{i}(j,1) = I+1;
        seqindmat{i}(j,2) = I2+1;
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
    L = length(seq{2*i});
    ss = seqout{i}-1;
    p = zeros(L+l-1,1);
    I = ss(1:l)*pow;
    p(1) = w(I+1);
    for j = 2:L-l+1
        I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1);
        p(j) = w(I+1);
    end
    v = zeros(L,3);
    for j = 1:L
        R = max([1 j-l+1]):min([j+l-1 L]);
        RR = max([1 j-l+1]):min([j L+l-1]);
        S = ss(R);
        ref = sum(p(RR));
        if max([1 j-l+1]) == 1
            cen = j;
            cen2 = j;
        elseif min([j+l-1 L]) == L
            cen = l;
            cen2 = L-j+1;
        else
            cen = l;
            cen2 = l;
        end
        for ii = 1:3
            S(cen) = mat(ss(j)+1,ii);
            I = S(1:l)*pow;
            v(j,ii) = w(I+1);
            if length(RR) > 1
                for jj = 2:cen2
                    I = (I-S(jj-1))/4+4^(l-1)*S(jj+l-1);
                    v(j,ii)=v(j,ii)+w(I+1);
                end
            end
            v(j,ii) = v(j,ii)-ref;
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

function process_motifs(dfn, lfn, memefn, ofn, PWMcorrcut,Negative)
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
if Negative
    w = -1*w;
    vec2 = -1*vec2;
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
        if cor > PWMcorrcut || vec2(ii-N) < 1.5
            hal = false;
        end
    elseif ii > 1 && a > 2
        hal = true;
        [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]);
        if cor > PWMcorrcut
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
