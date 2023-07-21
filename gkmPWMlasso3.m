function gkmPWMlasso(varargin)

if nargin < 2
    error('Need at least 2 inputs')
end

filename = varargin{1};
memefile = varargin{2};
minL = 10;
minInfo = 0.5;
corrCut = 0.86;
l_svm = 10;
k_svm = 6;
BG_GC = 0;
RC = true;

if nargin > 2
    f = find(strcmp('MinLength', varargin));
    if ~isempty(f);
        minL = varargin{f+1};
    end
    f = find(strcmp('MinInfo', varargin));
    if ~isempty(f);
        minInfo = varargin{f+1};
    end
    f = find(strcmp('CorrCutoff', varargin));
    if ~isempty(f);
        corrCut = varargin{f+1};
    end
    f = find(strcmp('l', varargin));
    if ~isempty(f);
        l_svm = varargin{f+1};
    end
    f = find(strcmp('k', varargin));
    if ~isempty(f);
        k_svm = varargin{f+1};
    end
    f = find(strcmp('Mode', varargin));
    if ~isempty(f) && strcmp('Compare',varargin{f+1});
        BG_GC = 1;
    end
    f = find(strcmp('RC', varargin));
    if ~isempty(f);
        RC = varargin{f+1};
    end
end

[comb,comb2,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm);%generate gapped positions, adjusted for reverse complements

if length(comb)*4^k_svm > 10^6
    %error([num2str(length(comb)*4^k_svm) ' exceeds the maximum number of gapped kmers allowed (10^6)'])
    disp(['l = 10, k = 6'])
    lk = 0;
    l_svm2 = l_svm;
    k_svm2 = k_svm;
    l_svm = 10;
    k_svm = 6;
    [comb,comb2,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm);
else
    lk = 1;
    l_svm2 = l_svm;
    k_svm2 = k_svm;
end

disp('Counting gapped kmers')

[cfile, GCpos1, GCneg1,mat,mat2] = getgkmcounts(filename, l_svm2, k_svm2, lk, RC);
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
[p, info, lenvec] = trim_pwm(p);
indvec = intersect(find(info./lenvec>=minInfo),find(lenvec>=minL));
n = length(indvec);
disp('Mapping PWMs to gkm space')
lcnum = length(comb);
A=zeros(lcnum*4^k_svm,n);
AA=zeros(lcnum*4^k_svm,n);
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
    normvec(j) = (A(:,j)'*A(:,j))^0.5;
    AA(:,j) = A(:,j)/normvec(j);
    %[~, b] = sort(A(:,j), 'descend');
    %AA(b(1:10*lcnum),j)=1;
end
fprintf('\n')
disp('Clustering motifs')
simmat = AA'*AA;
clear AA
motclus = clus_simmat_eig(simmat,corrCut);
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
clear A AA loc mat GMmat


disp('Running LASSO')
[weigmat, FitInfo] = lasso_cvmat(B, cfile2,'DFmax', length(Z),'Standardize', false);
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
FF = FF(b(1:f(1)));
f = [];
for i = 1:length(FF)
    f = [f;FF{i}];
end
f=unique(f);
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
fidw = fopen([filename '_gkmPWMlasso3.out'], 'w');
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

function [pp, info, len] = trim_pwm(p)
l = length(p);
info = zeros(l, 1);
len = zeros(l,1);
cut = 0;
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

