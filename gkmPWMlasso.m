function gkmPWMlasso(varargin)

if nargin < 3
    error('Need at least 3 inputs')
end

filename = varargin{1};
d = varargin{2};
memefile = varargin{3};
minL = 10;
minInfo = 0.5;
corrCut = 0.86;
l_svm = 10;
k_svm = 6;
BG_GC = 0;

if nargin > 3
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
    if strcmp('Compare',varargin{f+1});
        BG_GC = 1;
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
[cfile, GCpos1, GCneg1,mat,mat2] = getgkmcounts(filename, l_svm2, k_svm2, lk);
if BG_GC == 1
    mat = (mat+mat2)/2;
    GCpos1 = (GCpos1+GCneg1)/2;
    GCneg1 = GCpos1;
end
disp('Generating negative set')
negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm);

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
GCmat = repmat([0.5-GCpos1/2 GCpos1/2 GCpos1/2 0.5-GCpos1/2],l_svm-1,1);
per = 10;
normvec = zeros(n,1);
rev = zeros(n,1);
for j = 1:n
    if mod(j, floor(n/10))==0
        fprintf('%d...', per);
        per = per+10;
    end
    loc = zeros(l_svm*2-2+lenvec(indvec(j)), 1);
    loc(l_svm:lenvec(indvec(j))+l_svm-1) = 1;
    [A(:,j),rev(j)] = PWM2kmers([GCmat;p{indvec(j)};GCmat],mat,comb2,diffc,indc,loc,xc,l_svm,k_svm,rcnum);
    A(:,j) = A(:,j) - negvec*(l_svm-1+lenvec(indvec(j)));
    normvec(j) = (A(:,j)'*A(:,j))^0.5;
    A(:,j) = A(:,j)/normvec(j);
end
fprintf('\n')

disp('Clustering motifs')
simmat = A'*A;
motclus = clus_simmat_eig(simmat,corrCut);
disp(['Number of motif clusters:' num2str(length(motclus))])

disp('Selecting Motifs')
m = mean(cfile);
s = std(cfile);
%cfile2 = (1/2)*(1+erf((cfile-m)/s)).*cfile;
cfile2 = cfile-negvec/sum(negvec)*sum(cfile);
%cfile2 = cfile2/max(abs(cfile2));
cfile2 = cfile2/std(cfile2);
corrvec = zeros(n,1);
B = zeros((4^k_svm)*lcnum, length(motclus));
for i = 1:n
    [~,I] = sort(A(:,i),'descend');
    corrvec(i) = abs(corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10))));
    %corrvec(i) = abs(corr(A(:,i), cfile2));
end
for i = 1:length(motclus)
    [~,b] = sort(corrvec(motclus{i}),'descend');
    B(:,i) = A(:,motclus{i}(b(1)))*normvec(motclus{i}(b(1)));
    motclus{i} = motclus{i}(b);
end
clear A loc mat GMmat
stdB = sqrt(sum(sum(B.^2))/4^k_svm/lcnum/length(motclus));
B = B/stdB;

disp('Running LASSO (1)')
weigmat = lasso(B, cfile2,'DFmax', d,'Standardize', false);
f = find(weigmat(:,1)~=0);
cfile3 = cfile2 - B(:,f)*((B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2));

disp('Running LASSO (2)')
weigmat2 = lasso(B, cfile3,'DFmax', d,'Standardize', false);
ff = find(weigmat2(:,1)~=0);
C = corrcoef(B)-eye(length(motclus));
weigmat = abs(weigmat(:,1)) + abs(weigmat2(:,1));
f = find(weigmat~=0);
SEmat =(B(:,f).'*B(:,f))^-1;
OLS = SEmat*(B(:,f).'*cfile2);
res = cfile2-B(:,f)*OLS;
MSE = sqrt(res'*res)/length(cfile2);
SE = sqrt(diag(SEmat)*MSE);
OLSweig = zeros(length(weigmat(:,1)),1);
OLSweig(f) = OLS./SE;
[~,f] = sort(abs(OLSweig), 'descend');

disp('Selecting top motifs')
SEmat =(B(:,f(1:d)).'*B(:,f(1:d)))^-1;
OLS = SEmat*(B(:,f(1:d)).'*cfile2);
Pweig = zeros(d,1);
S = std(cfile2);
BB = B(:,f(1:d));
E = zeros(d,1);
for i = 1:d
    B = BB;
    B(:,i) = [];
    SEmat =(B'*B)^-1;
    OLS2 = SEmat*(B'*cfile2);
    res = cfile2-B*OLS2;
    E(i) = sqrt(res'*res);
    [~,a] = sort(BB(:,i),'descend');
    Pweig(i) = mean(cfile2(a(1:lcnum)))/S;
end
res = cfile2-BB*OLS;
EE = sqrt(res'*res);
E = (E-EE)/EE;
for i = 1:length(motclus)
    motclus{i} = indvec(motclus{i});
end
motclus = motclus(f);
gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus, [filename '_' num2str(l_svm2) '_' num2str(k_svm2) '_' num2str(d)] ,memefile,num,minL, minInfo, corr(cfile2, BB*OLS))

disp('Done')

function gettopmotifs(weigvec, pvec, E, motclus, filename, memefile, n, minL , minInfo, C)
[a, b] = sort(weigvec, 'descend');
c = [a pvec(b) E(b)];
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
fprintf(fid, 'Cluster ID\tMotif ID\tMotif Name\tRegression Weight\tZ-score\tImportance\n');
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

