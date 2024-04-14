function gkmPWMlasso(varargin)
% gkmPWMlasso finds predictive motifs from sequence-based models of regulatory
%     elements using a database of PWMs.  
% 
%     gkmPWMlasso(fileprefix, memefile, motifNum, ...)
% 
%     Works conveniently with the output of the gkmSVM R package or the lsgkm
%     package (https://github.com/Dongwon-Lee/lsgkm).  Can be leveraged to 
%     extract motifs from other sequence based models as long as you have scores 
%     paired with sequences.  The sequences should be in fasta format with the 
%     suffix *_svseq.fa.  The scores should be in a 2 column tab delimited 
%     file with the suffix *svalpha.out with the same prefix as the *svseq.fa 
%     file.  The first column containing the labels of each sequence and the 
%     scores in the second.  
% 
%     Uses a modified version of MATLAB's lasso function.  The changes made
%     were to reduce the computation time of working with large matrices.
%
%     Positional Parameters (Required):
% 
%     fileprefix      The prefix of the gkmSVM/lsgkm model (FILEHEADER_svseq.fa)
%                     or (FILEHEADER.model.txt)
%     memefile        The collection of PWMs in meme format
%     motifNum        The number of motifs to learn.  If 0 is specified, 
%                     gkmPWMlasso will try to find the optimal number of motifs
% 
%     Name Value Pair Parameters (Optional):
% 
%     'MinLength'     PWMs with lengths shorter than MinLength will be filtered
%                     out (default: 10)
%     'MinInfo'       PWMs with an average information per position less than 
%                     MinInfo will be filtered out (default: 0.5)
%     'CorrCutoff'    PWMs with a pearson correlation greater than CorrCutoff
%                     will be clustered together to prevent linear dependence
%                     and redundant features (default: 0.83)
%     'l'             The full length of the gapped k-mer.  This DOES NOT need
%                     to be the same as the l in the gkmSVM model (default: 11)
%     'k'             The number of ungapped positions of the gapped k-mer.
%                     This DOES NOT need to be the same as the k in the gkmSVM
%                     model (default: 7)
%     'Mode'          Use both the positive and negative set to get the background
%                     distribution of gapped k-mers.  This is best used when
%                     the model is trained with both the positive and negative
%                     sets containing regulatory elements.  (default: off)
%                     To turn on: gkmPWMlasso(fileprefix, memefile, motifNum, 'Mode', 'Compare') 
%     'RC'            If true, treat reverse complementary gapped k-mers as
%                     the same feature.  Otherwise treat them as separate features
%                     (default:true)
%     'KmerFrac'      Set the fraction of the total number of gapped k-mers to
%                     use with gkmPWMlasso.  This reduces the memory and runtime
%                     needed.  If the total number of gapped k-mers is too high
%                     with the given combination of (l,k,KmerFrac), KmerFrac will
%                     be automatically set to a lower value to create a more 
%                     workable number of gapped k-mers
%     'KmerFracLimit' Automatically lower KmerFrac if the number of gapped kmers
%.                    is too large.  Note that this will require more memory and
%                     increase runtime if set to false. (default: true)
% 
%     Outputs a file named fileprefix_l_k_motifNum_gkmPWMlasso.out that contains
%     the most predictive motifs with the given motif number.  If motifNum = 0,
%     the file will be named fileprefix_l_k_gkmPWMlasso.out. See the README.md 
%     for details on the format of the output
% 
%     Example (files in the example_files directory):
%     gkmPWMlasso('GM12878', 'combined_db_v4.meme', 30, 'MinLength',10,...
%         'MinInfo', 0.5, 'CorrCutoff', 0.83, 'l', 11, 'k', 7,'RC',true)
%         Outputs GM12878_11_7_30_gkmPWMlasso.out
if nargin < 3
    error('Need at least 3 inputs')
elseif nargin > 3
    if mod(nargin,2) == 0
        error('Incorrect number of inputs')
     else
        vec = 4:2:nargin;
        inputlib = {'MinLength', 'MinInfo', 'CorrCutoff', 'l', 'k', 'Mode', 'RC', 'KmerFrac', 'KmerFracLimit'};
        for i = 1:length(vec)
            f = strcmp(varargin{vec(i)},inputlib);
            if sum(f) == 0
                error([varargin{vec(i)} ' is not an input option'])
            end
        end
    end
end

filename = varargin{1};
memefile = varargin{2};
d = varargin{3};
minL = 10;
minInfo = 0.5;
corrCut = 0.83;
l_svm = 11;
k_svm = 7;
BG_GC = 0;
RC = true;
nfrac = 1;
nfracLim = true;
lk = 1;
if nargin > 2
    f = find(strcmp('MinLength', varargin));
    if ~isempty(f);
        minL = varargin{f+1};
        if ~isa(minL, 'double') || round(minL)-minL ~= 0 || minL <= 0
            error(['MinLength must be a positive integer'])
        end
    end
    f = find(strcmp('MinInfo', varargin));
    if ~isempty(f);
        minInfo = varargin{f+1};
        if ~isa(minInfo, 'double') || minInfo <= 0
            error(['MinInfo must be a positive float'])
        end
    end
    f = find(strcmp('CorrCutoff', varargin));
    if ~isempty(f);
        corrCut = varargin{f+1};
        if ~isa(corrCut, 'double') || corrCut <= 0
            error(['CorrCutoff must be a positive float'])
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
    f = find(strcmp('Mode', varargin));
    if ~isempty(f) && strcmp('Compare',varargin{f+1});
        BG_GC = 1;
    end
    f = find(strcmp('RC', varargin));
    if ~isempty(f);
        RC = varargin{f+1};
        if ~islogical(RC)
            error(['RC must be a boolean'])
        end
    end
    f = find(strcmp('KmerFrac', varargin));
    if ~isempty(f);
        nfrac = varargin{f+1};
        lk = [l_svm k_svm];
        if ~isa(nfrac, 'double') || nfrac <= 0 || nfrac >1
            error(['KmerFrac must be a positive float in (0 1]'])
        end
    end
    f = find(strcmp('KmerFracLimit', varargin));
    if ~isempty(f);
        nfracLim = varargin{f+1};
        if ~islogical(nfracLim)
            error(['KmerFracLimit must be a boolean'])
        end
    end
end

[comb,comb2,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm,nfrac);%generate gapped positions, adjusted for reverse complements

if nfracLim && numel(comb)/k_svm*4^k_svm > 5*10^5
    nfrac = round(5*10^7/4^k_svm/numel(comb)*k_svm)/100;
    disp(['Combination of (l,k) yields too many gapped kmers.  Using ' num2str(nfrac) ' of the total gapped kmers'])
    l_svm2 = l_svm;
    k_svm2 = k_svm;
    lk = ([l_svm k_svm]);
    [comb,comb2,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm,nfrac);
else
    l_svm2 = l_svm;
    k_svm2 = k_svm;
end

disp('Counting gapped kmers')

[cfile, GCpos1, GCneg1,mat,mat2] = getgkmcounts(filename, l_svm, k_svm, lk, RC, comb,rcnum);
if BG_GC == 1
    mat = (mat+mat2)/2;
    GCpos1 = (GCpos1+GCneg1)/2;
    GCneg1 = GCpos1;
end
disp('Generating negative set')
negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC);

disp('Filtering motifs')
num = length(strfind(fileread(memefile),'MOTIF'));
p = getmotif(memefile,1:num);
for i = 1:num
    [r c] = size(p{i});
    for j = 1:r
        a = sum(p{i}(j,:));
        if abs(a-1)>0
            [b1 loc] = max(p{i}(j,:));
            p{i}(j,loc) = b1-a+1;
        end
    end
end
[p, info, lenvec] = trim_pwm(p,0);
indvec = intersect(find(info./lenvec>=minInfo),find(lenvec>=minL));
n = length(indvec);
disp('Mapping PWMs to gkm space')
lcnum = numel(comb)/k_svm;
A=zeros(lcnum*4^k_svm,n);
GCmat = repmat([0.5-GCpos1/2 GCpos1/2 GCpos1/2 0.5-GCpos1/2],l_svm-1,1);
per = 10;
normvec = zeros(n,1);
for j = 1:n
    if mod(j, floor(n/10))==0
        fprintf('%d...', per);
        per = per+10;
    end
    loc = zeros(l_svm*2-2+lenvec(indvec(j)), 1);
    loc(l_svm:lenvec(indvec(j))+l_svm-1) = 1;
    if RC
        A(:,j) = PWM2kmers([GCmat;p{indvec(j)};GCmat],mat,comb2,diffc,indc,loc,xc,l_svm,k_svm,rcnum);
    else
        A(:,j) = PWM2kmers_norc([GCmat;p{indvec(j)};GCmat],mat,comb2,diffc,indc,loc,xc,l_svm,k_svm,rcnum);
    end
    A(:,j) = A(:,j) - negvec*(l_svm-1+lenvec(indvec(j)));
    normvec(j) = sqrt(A(:,j)'*A(:,j));
end
fprintf('\n')
disp('Clustering motifs')
simmat = diag(normvec.^-1)*(A'*A)*diag(normvec.^-1);
motclus = clus_simmat_eig(simmat,corrCut);
clear simmat
disp(['Number of motif clusters:' num2str(length(motclus))])

disp('Selecting Motifs')
cfile2 = cfile-negvec/sum(negvec)*sum(cfile);
cfile2 = cfile2/std(cfile2);
B = zeros((4^k_svm)*lcnum, length(motclus));
corrvec = zeros(n,1);
zvec = zeros(n,1);
Z = zeros(length(motclus),1);
for i = 1:n
    [~,I] = sort(A(:,i),'descend');
    zvec(i) = mean(cfile2(I(1:lcnum)));
    if zvec(i) < 0
        corrvec(i) = -1*corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10)));
    else
        corrvec(i) = corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10)));
    end
end
for i = 1:length(motclus)
    [a,b] = sort(zvec(motclus{i}),'descend');
    f = find(a == a(1));
    if length(f) > 1
        [~,bb] = sort(abs(corrvec(motclus{i}(b(f)))), 'descend');
        b(1:length(f)) = b(bb);
    end
    B(:,i) = A(:,motclus{i}(b(1)));
    motclus{i} = motclus{i}(b);
    Z(i) = zvec(motclus{i}(1));
end
f = find(abs(Z)>1);
B = B(:,f);
B = B/std(B(:))';
motclus = motclus(f);
Z = Z(f);
clear A loc mat GMmat

if d == 0
disp('Running LASSO')
[weigmat, FitInfo] = lasso_cvmat(B, cfile2,'DFmax', length(Z),'Standardize', false, 'NumLambda', 20);
MSE = zeros(length(FitInfo.DF),1);
cnorm = cfile2'*cfile2;
F = cell(numel(FitInfo.DF),1);
for i = 1:length(FitInfo.DF)
    f = find(weigmat(:,i)~=0);
    F{i} = f;
    OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2);
    res = B(:,f)*OLS;
    if sum(abs(res)>0);
        MSE(i) = (cfile2'*res)^2/(res'*res)/cnorm;
    end
end
FF = cell(1,1);
MSE2 = [];
count = 1;
for i = 1:length(MSE)-1
    if numel(setdiff(F{i}, F{i+1})) > 0
        FF{count} = setdiff(F{i}, F{i+1});
        MSE2(count) = (MSE(i)-MSE(i+1));
        count = count +1;
    end
end
res = B*((B.'*B)^-1*(B.'*cfile2));
csm = (cfile2'*res)^2/(res'*res)/cnorm;
[a,b] = sort(MSE2,'descend');
cs = cumsum(a);
cs = cs/csm;
f = find(cs>0.9);
if isempty(f)
    f = 1:length(OLS);
else
    FF = FF(b(1:f(1)));
    f = [];
    for i = 1:length(FF)
        f = [f;FF{i}];
    end
    f=unique(f);
end
%f = find(weigmat(:,1)~=0);
F = length(f);
disp('Selecting top motifs')
ind = true;
OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2);
Pweig = Z(f);
while ind
    ff = [];
    for i = 1:length(f)    
        if sign(OLS(i)) ~= sign(Pweig(i))
            ff = [ff;i];
        end
    end
    if length(ff) > 0
        f(ff) = [];
        OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2);
        Pweig = Z(f);
    else
        ind = false;
    end
end
BB = B(:,f);
BX = B(:,f)'*B(:,f);
BY = B(:,f)'*cfile2;
E = zeros(length(f),1);
for i = 1:length(f)
    B = BB;
    BBX = BX;
    BBY = BY;
    B(:,i) = [];
    BBX(:,i) = [];
    BBX(i,:) = [];
    BBY(i) = [];
    res = cfile2-B*(BBX^-1*BBY);
    E(i) = sqrt(res'*res);
end
res = cfile2-BB*OLS;
EE = sqrt(res'*res);
E = (E-EE)/EE;
motclus = motclus(f);
f = find(E/max(E) >= 0.01);
B = BB;
BB = B(:,f);
BX = B(:,f)'*B(:,f);
BY = B(:,f)'*cfile2;
E = zeros(length(f),1);
for i = 1:length(f)
    B = BB;
    BBX = BX;
    BBY = BY;
    B(:,i) = [];
    BBX(:,i) = [];
    BBX(i,:) = [];
    BBY(i) = [];
    res = cfile2-B*(BBX^-1*BBY);
    E(i) = sqrt(res'*res);
end
OLS = (BB.'*BB)^-1*(BB.'*cfile2);
Pweig = Pweig(f);
res = cfile2-BB*OLS;
EE = sqrt(res'*res);
E = (E-EE)/EE;
for i = 1:length(motclus)
    motclus{i} = indvec(motclus{i});
end
motclus = motclus(f);
gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus, [filename '_' num2str(l_svm2) '_' num2str(k_svm2)] ,memefile,num,minL, minInfo, corr(cfile2, BB*OLS))

else

disp('Running LASSO (1)')
weigmat = lasso_cvmat(B, cfile2,'DFmax', d,'Standardize', false);
f = find(weigmat(:,1)~=0);
cfile3 = cfile2 - B(:,f)*((B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2));

disp('Running LASSO (2)')
weigmat2 = lasso_cvmat(B, cfile3,'DFmax', d,'Standardize', false);
weigmat = abs(weigmat(:,1)) + abs(weigmat2(:,1));
f = find(weigmat~=0);
motclus2 = clus_simmat_eig(corrcoef(B(:,f)),0.6);
if length(length(motclus2)) ~= length(f)
    f2 = f;
    f = zeros(length(motclus2),1);
    for i = 1:length(motclus2)
        if length(motclus2{i}) > 1
            [~,f3] = max(abs(Z(f2(motclus2{i}))));
            f(i) = f2(motclus2{i}(f3));
        else
            f(i) = f2(motclus2{i}(1));
        end
    end
end
[a,b] = sort(abs(Z(f)), 'descend');
f = f(b);
F = length(f);
disp('Selecting top motifs')
ind = true;
while ind && F >= d
    OLS = (B(:,f(1:d)).'*B(:,f(1:d)))^-1*(B(:,f(1:d)).'*cfile2);
    Pweig = Z(f(1:d));
    ff = [];
    for i = 1:d
        if sign(OLS(i)) ~= sign(Pweig(i))
            ff = [ff;i];
        end
    end
    if length(ff) > 0 && F - length(ff) >= d
        f(ff) = [];
        F = F - length(ff);
    else
        ind = false;
        f = f(1:d);
    end
end
if F < d
    OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2);
    Pweig = Z(f);
end
BB = B(:,f);
BX = B(:,f)'*B(:,f);
BY = B(:,f)'*cfile2;
E = zeros(length(f),1);
for i = 1:length(f)
    B = BB;
    BBX = BX;
    BBY = BY;
    B(:,i) = [];
    BBX(:,i) = [];
    BBX(i,:) = [];
    BBY(i) = [];
    res = cfile2-B*(BBX^-1*BBY);
    E(i) = sqrt(res'*res);
end
res = cfile2-BB*OLS;
EE = sqrt(res'*res);
E = (E-EE)/EE;
for i = 1:length(motclus)
    motclus{i} = indvec(motclus{i});
end
motclus = motclus(f);
gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus, [filename '_' num2str(l_svm2) '_' num2str(k_svm2) '_' num2str(d)] ,memefile,num,minL, minInfo, corr(cfile2, BB*OLS))
end

disp('Done')

function gettopmotifs(weigvec, pvec, E, motclus, filename, memefile, n, minL , minInfo, C)
[a, b] = sort(pvec, 'descend');
c = [weigvec(b) a E(b)];
motclus = motclus(b);
[rr,cc] = size(c);
names = cell(n,1);
i = 0;
fid = fopen(memefile, 'r');
while ~feof(fid)
    line = fgetl(fid);
    if length(line) >= 5
        if strcmp(line(1:5), 'MOTIF')
            i = i+1;
            [~,name] = strtok(line);
            names{i} = strtrim(name);
        end
    end
end
fclose(fid);
fidw = fopen([filename '_gkmPWMlasso.out'], 'w');
fprintf(fidw, 'Minimum PWM Length:\t%d\n', minL);
fprintf(fidw, 'Minimum PWM Information:\t%f\n', minInfo);
fprintf(fidw, 'Correlation with SVM weights:\t%f\n', C);
fprintf(fidw, 'Cluster ID\tMotif ID\tMotif Name\tRegression Weight\tZ-score\tImportance\n');
for j = 1:length(E)
    for l = 1:length(motclus{j})
        fprintf(fidw, '%d\t%d\t%s\t%0.3f\t%0.3f\t%0.3f\n', j, motclus{j}(l),names{motclus{j}(l)}, c(j,1), c(j,2), c(j,3));
    end
end
fclose(fidw);

function motclus = clus_simmat_eig(simmat,r)
bin = conncomp(graph(simmat>r));
n = max(bin);
motclus = cell(n,1);
for i = 1:n
    motclus{i} = find(bin == i);
end

function [pp, info, len] = trim_pwm(p,cut)
l = length(p);
info = zeros(l, 1);
len = zeros(l,1);
for i = 1:l
    mat = p{i}+(p{i}==0);
    vec = 2+sum(mat.*log(mat)/log(2),2);
    mvec = movmean(vec,5);
    while min(vec(1),mvec(1)) < cut && length(vec) > 1
        p{i}(1,:) = [];
        vec(1) = [];
        mvec(1)=[];
        mat(1,:) = [];
    end
    while min(mvec(end),vec(end)) < cut && length(vec) > 1
        vec(end) = [];
        mvec(end)=[];
        mat(end,:) = [];
        p{i}(end,:) = [];
    end
    info(i) = sum(vec);
    [len(i), ~] = size(mat);
end
pp = p;

