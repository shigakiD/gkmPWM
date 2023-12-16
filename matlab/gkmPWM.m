function gkmPWM(varargin)
% gkmPWM finds predictive motifs from sequence-based models of regulatory
%     elements de novo  
% 
%     gkmPWM(fileprefix, wfile, memefile, motifNum, ...)
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
%     Positional Parameters (Required):
% 
%     fileprefix      The prefix of the gkmSVM/lsgkm model (FILEHEADER_svseq.fa)
%                     or (FILEHEADER.model.txt)
%     wfile           Two column tab delimited file of kmer weights.  See the
%                     README.md for detail on how to generate this file.
%     memefile        The collection of PWMs in meme format
%     motifNum        The number of positive set motifs to learn.  
% 
%     Name Value Pair Parameters (Optional):
% 
%     'MaxCorr'       If a PWM has a pearson correlation greater than this value
%                     with another PWM, it will be reseeded.  This is to prevent
%                     linear dependence and redundancy. (default: 0.9)
%     'MaxIter'       Maximum number of iterations gkmPWM will run before exiting.
%                     If the loss converges, gkmPWM will finish prior to this
%                     value.  (default: 200)
%     'PNratio'       The ratio of positive PWMs to negative PWMs that will
%                     be seeded.  e.g., PNratio = 2 means twice as many positive
%                     set PWMs will be seeded as negative set PWMs.  If not set
%                     gkmPWM will automatically set this value based on the model.
%     'RegFrac'       This parameter will push the PWMs to have higher information.
%                     Must be a value in [0,1).  Keeping this value low (<0.05)
%                     works best.  (default: 0).
%     'l'             The full length of the gapped k-mer.  This DOES NOT need
%                     to be the same as the l in the gkmSVM model (default: 10)
%     'k'             The number of ungapped positions of the gapped k-mer.
%                     This DOES NOT need to be the same as the k in the gkmSVM
%                     model (default: 6)
%     'Mode'          Use both the positive and negative set to get the background
%                     distribution of gapped k-mers.  This is best used when
%                     the model is trained with both the positive and negative
%                     sets containing regulatory elements.  (default: off)
%                     To turn on: gkmPWM(fileprefix, wfile,memefile, motifNum, 'Mode', 'Compare') 
%     'RC'            If true, treat reverse complementary gapped k-mers as
%                     the same feature.  Otherwise treat them as separate features
%                     (default:true)
%     'KmerFrac'      Set the fraction of the total number of gapped k-mers to
%                     use with gkmPWM.  This reduces the memory and runtime
%                     needed.  If the total number of gapped k-mers is too high
%                     with the given combination of (l,k,KmerFrac), KmerFrac will
%                     be automatically set to a lower value to create a more 
%                     workable number of gapped k-mers
% 
%     Outputs 3 files named:
%     fileprefix_l_k_RegFrac_motifNum_gkmPWM.out
%         Contains the most predictive motifs with the given motif number.  
%         See the README.md for details on the format of the output
%     fileprefix_l_k_RegFrac_motifNum_denovo.meme
%         A meme file containing the PWMs of each motif in gkmPWM
%     fileprefix_l_k_RegFrac_motifNum_error.out
%         A file containing the loss after each iteration.  Allows you to
%         see how well the model converged over time.
%
%     Example (files in the example_files directory):
%     gkmPWM('GM12878', 'GM12878_weights.out','combined_db_v4.meme', 15,...
%         'MaxCorr',0.9, 'MaxIter', 200, 'RegFrac', 0, 'l', 10,'k', 6,'RC',true)
%         Outputs GM12878_10_6_0_15_gkmPWM.out, GM12878_10_6_0_15_denovo.meme
%         and GM12878_10_6_0_15_error.out
if nargin < 4
    error('Need at least 4 inputs')
elseif nargin > 4
    if mod(nargin,2) == 1
        error('Incorrect number of inputs')
     else
        vec = 5:2:nargin;
        inputlib = {'MaxCorr', 'MaxIter', 'PNratio', 'RegFrac', 'l', 'k', 'Mode', 'RC', 'KmerFrac'};
        for i = 1:length(vec)
            f = strcmp(varargin{vec(i)},inputlib);
            if sum(f) == 0
                error([varargin{vec(i)} ' is not an input option'])
            end
        end
    end
end

fileprefix = varargin{1};
wfile = varargin{2};
memefile = varargin{3};
mnum = varargin{4};
num = 200;
rcorr = 0.9;
reg = 0;
l_svm = 10;
k_svm = 6;
BG_GC = 0;
RC = true;
ipnr = true;
nfrac = 1;
lk = 1;
if nargin > 4
    f = find(strcmp('MaxCorr', varargin));
    if ~isempty(f);
        rcorr = varargin{f+1};
    end
    f = find(strcmp('MaxIter', varargin));
    if ~isempty(f);
        rcorr = varargin{f+1};
    end
    f = find(strcmp('PNratio', varargin));
    if ~isempty(f);
        pnr = varargin{f+1};
        ipnr=false;
    end
    f = find(strcmp('RegFrac', varargin));
    if ~isempty(f);
        reg = varargin{f+1};
        if reg < 0 || reg >= 1
            error('RegFrac must be in [0,1)')
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
    f = find(strcmp('Mode', varargin));
    if ~isempty(f) && strcmp('Compare',varargin{f+1});
        BG_GC = 1;
    end
    f = find(strcmp('RC', varargin));
    if ~isempty(f);
        RC = varargin{f+1};
    end
    f = find(strcmp('KmerFrac', varargin));
    if ~isempty(f);
        nfrac = varargin{f+1};
        lk = [l_svm k_svm];
    end
end

[comb,rc,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm,nfrac);%generate gapped positions, adjusted for reverse complements
if length(comb)*4^k_svm > 6*10^5
    nfrac = round(5*10^7/4^k_svm/length(comb))/100;
    disp(['Combination of (l,k) yields too many gapped kmers.  Using ' num2str(nfrac) ' of the total gapped kmers'])
    lk = [l_svm k_svm];
    [comb,rc,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm,nfrac);
end

disp(['Running gkmPWM on ' fileprefix ' for ' num2str(mnum) ' motifs and ' num2str(num) ' iterations'])
disp('Counting gapped k-mers')

[A, GCpos1, GCneg1,mat,mat2] = getgkmcounts(fileprefix,l_svm,k_svm,lk,RC,comb,rcnum); %count gapped k-mers, get GC content, and dinucleotide distribution
if BG_GC == 1
    mat = (mat+mat2)/2;
    GCpos1 = (GCpos1+GCneg1)/2;
    GCneg1 = GCpos1;
end
negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC);%generate expected gapped kmer distribution of background
GC=[0.5-GCneg1/2 GCneg1/2 GCneg1/2 0.5-GCneg1/2];%GC content vector
if ipnr
    pnr = abs(max(A)/min(A));
end
disp('Finding PWM seeds')
%get PWM seeds using the kmer weight vectors
[kmers, seed,p,c] = seed_kmers(wfile, mnum,'descend', {});
if pnr ~= 0
    [kmers2, seed2,pp, c2] = seed_kmers(wfile, max([floor(mnum/pnr) 2]),'ascend',seed);
    kmers = [kmers;kmers2];
    p = [p;pp];
    tot = c2;
else
    tot = c;
end
disp('Seeding PWMs at the following kmers')
for i = 1:tot
    disp(kmers{i})
end

for i = 1:length(p)
    p{i} = extendPWM(p{i},l_svm+1,GC);
end

disp('Running de novo motif discovery')
m = mean(A);
s = std(A);
cfile = A-negvec/sum(negvec)*sum(A);
cfile = cfile/max(abs(cfile));% normalize to speed up computation
clear A
[pp scorevec C r R E Rd] = gkmPWM_lagrange(cfile,mat,p,negvec,num,rcorr,reg,l_svm,k_svm,RC,rc,diffc,indc,xc,rcnum);

createMEME([fileprefix '_' num2str(l_svm) '_' num2str(k_svm) '_' num2str(reg) '_' num2str(mnum)], pp, memefile,GCneg1, C, r, R, rcorr, E, Rd);
dlmwrite([fileprefix '_' num2str(l_svm) '_' num2str(k_svm) '_' num2str(reg) '_' num2str(mnum) '_error.out'],scorevec');

function createMEME(fileh,PWM, memefile, GC, C, r, R, rcorr, E, Rd)

num = numel(C);
GC = round(GC*100)/100;
GCvec = [0.5-GC/2 GC/2 GC/2 0.5-GC/2];
f = strfind(fileread(memefile),'MOTIF');
num2 = length(strfind(fileread(memefile),'MOTIF'));
[p,names] = getmotif(memefile,1:num2);
fid = fopen([fileh '_denovo.meme'], 'w');
fid2 = fopen([fileh '_gkmPWM.out'], 'w');
lenvec = zeros(num2,1);
for i = 1:num2
    [lenvec(i),~] = size(p{i});
end
fprintf(fid, 'MEME\n\n');
fprintf(fid, 'ALPHABET= ACGT\n\n');
fprintf(fid, 'Correlation with SVM weight vector: %0.3f\n\n', r);
fprintf(fid, 'Max PWM Correlation: %0.3f\n\n', rcorr);
fprintf(fid, 'Background letter frequencies (from negative set)\n');
fprintf(fid, 'A %0.2f C %0.2f G %0.2f T %0.2f\n\n', GCvec(1), GCvec(2), GCvec(3), GCvec(4));
fprintf(fid2, 'Correlation with SVM weight vector:\t%0.3f\n', r);
fprintf(fid2, 'Max PWM Correlation:\t%0.3f\n', rcorr);
fprintf(fid2, 'MOTIF\tID\tSimilarity\tRedundancy\tWeight\tZ\tError\n');
a = 1;
b = 1;
for i = 1:num
    [len,~] = size(PWM{i});
    [M, MM] = matchMotif([PWM{i}; p], [len;lenvec]);% matches motifs to the best motif
    fprintf(fid, 'MOTIF %d %s\n', a, names{MM});
    fprintf(fid2,'%s\t%d\t%0.3f\t%0.3f\t%0.2f\t%0.3f\t%0.3f\n',names{MM},MM,M,Rd(i),C(i),R(i),E(i));
    fprintf(fid, 'weight= %0.3f l= 4 w= %d z-score= %0.2f motifsim= %0.3f\n', C(i), len, R(i), M);
    for j = 1:len
        fprintf(fid, '%0.3f %0.3f %0.3f %0.3f\n',PWM{i}(j,1),PWM{i}(j,2),PWM{i}(j,3),PWM{i}(j,4));
    end
    fprintf(fid, '\n');
    a = a+1;
end
fclose(fid);
fclose(fid2);

function [p, mat,pwms, c] = seed_kmers(fn, num, pn, ik);
fid = fopen(fn, 'r');
a = textscan(fid, '%s\t%f\n');
fclose(fid);
[w, ind] = sort(a{2}, pn);
s = a{1}(ind(1:min([100000 length(a{1})])));
l = length(s{1});
k = round(l/2)+1;
ikl = length(ik);
p = cell(num,1);
c = ikl+1;
p{1} = s{1};
mat = cell(ikl+num,1);
mat(1:ikl) = ik;
mat{c} = letterconvert(s{1});
pwms = cell(num,1);
for i = 1:num
    pwms{i} = zeros(l,4);
end
for i = 1:l
    pwms{1}(i,mat{c}(i)+1) = pwms{1}(i,mat{c}(i)+1)+w(i);
end
B = zeros(9,1);
BB = zeros(9,1);
B(1:5) = (0:4)';
B(6:9) = 0;
BB(1:5) = 0;
BB(6:9) = (1:4)';
CC = [l l-1 l-2 l-3 l-4 l-1 l-2 l-3 l-4];
%this process picks kmers to seed the PWMs.  kmers that match one of the seeds by round(l/2)+1 or more are added that particular seed.  Otherwise, it becomes another seed.
for i = 2:100000
    ss = letterconvert(s{i});
    rs = 3-fliplr(ss);
    M = zeros(c,1);
    D = zeros(c,1);
    DD = zeros(c,1);
    for j = 1:c
        [m,d] = max([sum(mat{j}==ss) sum(mat{j}(2:end)==ss(1:l-1)) sum(mat{j}(3:end)==ss(1:l-2)) sum(mat{j}(4:end)==ss(1:l-3)) sum(mat{j}(5:end)==ss(1:l-4)) sum(mat{j}(1:l-1)==ss(2:end)) sum(mat{j}(1:l-2)==ss(3:end)) sum(mat{j}(1:l-3)==ss(4:end)) sum(mat{j}(1:l-4)==ss(5:end))]);
        [mm,dd] = max([sum(mat{j}==rs) sum(mat{j}(2:end)==rs(1:l-1)) sum(mat{j}(3:end)==rs(1:l-2)) sum(mat{j}(4:end)==rs(1:l-3)) sum(mat{j}(5:end)==rs(1:l-4)) sum(mat{j}(1:l-1)==rs(2:end)) sum(mat{j}(1:l-2)==rs(3:end)) sum(mat{j}(1:l-3)==rs(4:end)) sum(mat{j}(1:l-4)==rs(5:end))]);
        [M(j),ddd] = max([m mm]);
        if ddd == 1
            D(j) = d;
            DD(j) = 1;
        else
            D(j) = dd;
            DD(j) = 2;
        end
    end
    if max(M) < k
        c = c+1;
        p{c-ikl} = s{i};
        mat{c} = ss;
        ss = ss+1;
        for j = 1:l
            pwms{c-ikl}(j,ss(j)) = pwms{c-ikl}(j,ss(j))+w(i);
        end
    else
        [~,d] = max(M);
        if DD(d) == 1 && d > ikl
            ss = ss+1;
            d = d-ikl;
            for j = 1:CC(D(d))
                pwms{d}(j+B(D(d)),ss(j+BB(D(d)))) = pwms{d}(j+B(D(d)),ss(j+BB(D(d))))+w(i);
            end
        elseif DD(d) == 2 && d > ikl
            rs = rs+1;
            d = d-ikl;
            for j = 1:CC(D(d))
                pwms{d}(j+B(D(d)),rs(j+BB(D(d)))) = pwms{d}(j+B(D(d)),rs(j+BB(D(d))))+w(i);
            end
        end
    end
    if c == num+ikl
        break
    end
end
mat = mat(1:c);
p = p(1:c-ikl);
pwms = pwms(1:c-ikl);
for i = 1:c-ikl
    for j = 1:l
        pwms{i}(j,:) = pwms{i}(j,:)/sum(pwms{i}(j,:));
    end
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

function [M, ind] = matchMotif(mot,lenvec)
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

function [PWM, scorevec, C, r, R, E, Rd] = gkmPWM_lagrange(kweig,negmat,PWM,negvec,n,rcorr,reg,l_svm,k_svm,RC,rc,diffc,indc,xc,rcnum)
%Note: This code is rather messy.  I block commmented to the best of my ability, so hopefully this makes enough sense.  If something seems non-trivial, then I probably found a mathematical trick to speed up computation time (in particular dynamic programming).

GC = PWM{1}(1,:);
lcomb = length(diffc);
diffC = zeros(lcomb,l_svm);
for i = 1:l_svm
    ct = rc+i-1;
    f = find(sum(ct==l_svm,2));
    ct = ct(f,:);
    CT = zeros(length(ct),k_svm-1);
    for j = 1:length(ct)
        a = 1;
        for jj = 1:k_svm
            if ct(j,jj) ~= l_svm
                CT(j,a) = ct(j,jj);
                a = a+1;
            end
        end
    end
    for j = 2:length(f)
        a = 1;
        while CT(j,a)==CT(j-1,a)
            a = a+1;
        end
        if a < 2
            a = 2;
        end
        diffC(f(j),i)=a;
    end
    diffC(f(1),i)=2;
end
m = length(PWM);
scorevec = zeros(1,n);
lenvec = zeros(m,1);
loc = cell(m, 1);
for i = 1:m
    lenvec(i) = length(PWM{i})-l_svm*2+2;
    loc{i} = zeros(length(PWM{i}), 1);
    loc{i}(l_svm:lenvec(i)+l_svm-1) = 1;
end

kmat = zeros(lcomb*4^k_svm, m);
KMAT = zeros(lcomb*4^k_svm, m);

disp('Mapping PWMs to gkm space')
if RC
    for i = 1:m
        kmat(:,i) = PWM2kmers(PWM{i},negmat,rc,diffc,indc,loc{i},xc,l_svm,k_svm,rcnum)-negvec*(lenvec(i)+l_svm-1);%map PWMs to gapped kmers
    end
else
    for i = 1:m
        kmat(:,i) = PWM2kmers_norc(PWM{i},negmat,rc,diffc,indc,loc{i},xc,l_svm,k_svm,rcnum)-negvec*(lenvec(i)+l_svm-1);%map PWMs to gapped kmers
    end
end

%the following loop creates indices for the PWM column optimize to utilize dynamic programming.
poscell = cell(k_svm,1);
for i = 1:k_svm
    temp = zeros(4^(k_svm-1),1)';
    vec = 1:4^(i):4^(k_svm);
    for ii = 1:4^(i-1)
        t = length(vec);
        temp(1+(ii-1)*t:t+(ii-1)*t) = vec+ii-1;
    end
    poscell{i} = sort(temp)';
end
disp('Running Recursion')
acount = 0;
i = 0;
scount = 0;
tic
while i < n
    i = i+1;
    if mod(i,10) == 0
        toc
        tic
        fprintf('%d iterations done...\n',i);
    end
    scount = scount + 1;
    if  i >= 10 && scount >= 5 && max(-1*diff(scorevec(i-5:i-1))./scorevec(i-4:i-1)) < 0.001 && acount < 5 && i ~= n
        acount = acount + 1;
        fprintf('adjusting PWMs after %d iterations (%d)\n',i,acount);
        scount = 0;
        for ii = 1:m
            if i/n <= 0.8 && C(ord(ii)) > 0
                [PWM{ord(ii)}, lenvec(ord(ii))] = adjust_PWM(PWM{ord(ii)}(l_svm:(length(PWM{ord(ii)})-l_svm+1),:),GC);
                PWM{ord(ii)} = extendPWM(PWM{ord(ii)}, l_svm-1, GC);
                loc{ord(ii)} = zeros(lenvec(ord(ii))+2*l_svm-2, 1);
                loc{ord(ii)}(l_svm:lenvec(ord(ii))+l_svm-1) = 1;
            end
        end
    end
    C = (kmat'*kmat)^(-1)*(kmat'*kweig);
    res = kweig-kmat*C;
    corrvec = zeros(m,1);
    for ii = 1:m
        [~,ind] = sort(kmat(:,ii), 'descend');
        corrvec(ii) = sum(res(ind(1:lcomb)).^2);
        if mod(i,20) == 0
            if RC
                kmat(:,ii) = PWM2kmers(PWM{ii},negmat,rc,diffc,indc,loc{ii},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(ii));
            else
                kmat(:,ii) = PWM2kmers_norc(PWM{ii},negmat,rc,diffc,indc,loc{ii},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(ii));
            end 
        end
    end
    %The order of PWM optimization is determined by the correlation of its top 110 kmers with the gapped kmer weight vector
    [~,ord] = sort(corrvec, 'descend');
    for ii = 1:m
        %The order of the column optimization is determined by the max probability in each column
        v = max(PWM{ord(ii)}(l_svm:lenvec(ord(ii))+l_svm-1,:)');
        [~,c] = sort(v, 'ascend');
        %C = (kmat'*kmat)^(-1)*(kmat'*kweig);
        %res = kweig-kmat*C;
        for iii = 1:length(c)
            PWMtemp = PWM{ord(ii)}(c(iii):c(iii)+l_svm*2-2,:);
            [kweigdiff,PWM{ord(ii)}(c(iii)+l_svm-1,:)] = getEMprob_v3(PWMtemp,res/C(ord(ii)),negmat,poscell,rc,diffC,indc,loc{ord(ii)}(c(iii):c(iii)+2*l_svm-2),xc,reg,l_svm,k_svm,rcnum,RC);
            kmat(:,ord(ii)) = kmat(:,ord(ii)) + kweigdiff;
            res = res-kweigdiff*C(ord(ii));
        end
    end
    %Reseed PWMs if two or more of them are too highly correlated
    if i/n <= 0.80
        info = avg_info(PWM,l_svm);
        [~,ord] = sort(info,'descend');
        kmat = kmat(:,ord);
        lenvec = lenvec(ord);
        loc = loc(ord);
        PWM = PWM(ord);
        C = C(ord);
        for j = 1:m
            KMAT(:,j) = kmat(:,j)/sqrt(kmat(:,j)'*kmat(:,j));
        end
        MAT = KMAT'*KMAT;
        for j = 1:m-1
            vec = MAT(j+1:end,j);
            [a b] = max(vec);
            if a > rcorr && C(j) > 0
                scount = 0;
                disp('reseeding')
                f = j+find(vec > rcorr);
                for jj = 1:length(f)
                    for jjj = 1:min([lenvec(f(jj)) 12])
                        PWM{f(jj)}(jj+9,:) =  PWM{f(jj)}(jjj+l_svm-1,randperm(4));
                    end
                    if lenvec(f(jj)) >= 12
                        PWM{f(jj)} = PWM{f(jj)}(l_svm:l_svm+11,:);
                        PWM{f(jj)} = extendPWM(PWM{f(jj)}, l_svm-1, GC);
                        lenvec(f(jj)) = 12;
                        loc{f(jj)} = zeros(lenvec(f(jj))+l_svm*2-2, 1);
                        loc{f(jj)}(l_svm:lenvec(f(jj))+l_svm-1) = 1;
                    end
                    if RC
                        kmat(:,f(jj)) = PWM2kmers(PWM{f(jj)},negmat,rc,diffc,indc,loc{f(jj)},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(f(jj)));
                    else
                        kmat(:,f(jj)) = PWM2kmers_norc(PWM{f(jj)},negmat,rc,diffc,indc,loc{f(jj)},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(f(jj)));
                    end
                end
            end
        end
    end
    %Breaks the loop if it looks like it converged
    C = (kmat'*kmat)^(-1)*(kmat'*kweig);
    scorevec(i) = sqrt(res'*res);
    if i >= 10 && acount == 5 && scount >= 10 && max(abs(diff(scorevec(i-9:i))./scorevec(i-8:i))) < 0.0001 
        scorevec = scorevec(1:i);
        break
    end
    if i > 10 && acount == 5 && scount >= 10 && sum(diff(scorevec(i-9:i))>0) > 7
        scorevec = scorevec(1:i);
        break
    end
end
toc
disp(['gkmPWM completed after ' num2str(i) ' iterations'])

for i = 1:length(PWM)
    PWM{i} = PWM{i}(l_svm:(length(PWM{i})-l_svm+1),:);
end
%the following just calculates a few interesting quantities
r = corr(kweig, kmat*C);
M = mean(kweig);
S = std(kweig);
R = zeros(m,1);
E = zeros(m,1);
CM = corrcoef(kmat)-eye(m);
for i = 1:m
    [~,a] = sort(kmat(:,i),'descend');
    R(i) = (mean(kweig(a(1:lcomb)))-M)/S;
    Kmat = kmat;
    Kmat(:,i) = [];
    c = (Kmat'*Kmat)^(-1)*(Kmat'*kweig);
    res = kweig-Kmat*c;
    E(i) = (sqrt(res'*res)-scorevec(end))/scorevec(end);
    if C(i) < 0
        CM(i,:) = 0;
        CM(:,i) = 0;
    end
end
[R,a] = sort(R, 'descend');
PWM = PWM(a);
C = C(a);
E = E(a);
Rd = max(CM);
Rd = Rd(a);

function [kweig,P] = getEMprob_v3(PWM,res,negmat,poscell,rc,diffc,indc,indloc,xc,reg,l_svm,k_svm,rcnum,RC)
%Lagrange optimization (see paper for derivation)
a = true;
posvec = 1:4;
if RC
    A =  ls_kweigtree(PWM,negmat,poscell,rc,diffc,indc,indloc,xc,l_svm,k_svm,rcnum);
else
    A =  ls_kweigtree_norc(PWM,negmat,poscell,rc,diffc,indc,indloc,xc,l_svm,k_svm,rcnum);
end
b = res+A*PWM(l_svm,:)';
mat = A'*A;
y = A'*b;
if reg > 0
    M = min(eig(mat));
    B=(mat-reg*M*eye(4))^-1;
else
    M = 0; 
    B = mat^-1;
end
ps = B*y;
p = ps+(1-sum(ps))*B/sum(sum(B))*ones(4,1);
%The solution to the lagrange optimization problem
%The following deals with the case if the optimal solution has a negative probability.  I cheat a little by only considering the cases where the base with the maximum solution is non-zero.  This works just fine in practice since the sum(p) = 1 constraint forces one of the bases to be positive.  It speeds up computation almost two fold.
if min(p) < 0
    I = 0;
    [~,a] = max(p);
    nvec = posvec;
    nvec(a) = [];
    %Check cases where one of the probabilities is zero
    for i = 1:3
        Posvec = posvec;
        Posvec(nvec(i)) = [];
        MAT = mat;
        MAT(:,nvec(i)) =[];
        MAT(nvec(i),:) =[];
        Y = y(Posvec);
        if reg > 0
            B=(MAT-reg*M*eye(3))^-1;
        else
            B=(MAT)^-1;
        end
        ps = B*Y;
        p = ps+(1-sum(ps))*B/sum(sum(B))*ones(3,1);%solution
        %Checks if solution is permitted.  If so, makes sure that it creates a smaller error than other cases
        if min(p) >= 0
            vec = b-A(:,Posvec)*p;
            if I == 0
                I = 1;
                P = zeros(4,1);
                P(Posvec) = p;
                e = vec'*vec-reg*M*p'*p;
            else
                E = vec'*vec-reg*M*p'*p;
                if E < e
                    e = E;
                    P = zeros(4,1);
                    P(Posvec) = p;
                end
            end
        end
    end
    %Check cases where two of the probabilities are zero
    ind = [nvec(1) nvec(2);nvec(1) nvec(3);nvec(2) nvec(3)];
    for i = 1:3
        Posvec = posvec;
        Posvec(ind(i,:)) = [];
        MAT = mat;
        MAT(:,ind(i,:)) =[];
        MAT(ind(i,:),:) =[];
        Y = y(Posvec);
        if reg > 0
            B=(MAT-reg*M*eye(2))^-1;
        else
            B=MAT^-1;
        end
        ps = B*Y;
        p = ps+(1-sum(ps))*B/sum(sum(B))*ones(2,1);%solution
        %Checks if solution is permitted.  If so, makes sure that it creates a smaller error than other cases
        if min(p) >= 0
            vec = b-A(:,Posvec)*p;
            if I == 0
                I = 1;
                P = zeros(4,1);
                P(Posvec) = p;
                e = vec'*vec-reg*M*p'*p;
            else
                E = vec'*vec-reg*M*p'*p;
                if E < e
                    e = E;
                    P = zeros(4,1);
                    P(Posvec) = p;
                end
            end
        end 
    end
    %Checks to see if one non-zero case is better than the other cases
    if I == 0
        P = zeros(4,1);
        P(a) = 1;
    else
        vec = b-A(:,a);
        E = vec'*vec-reg*M;
        if E < e
            e = E;
            P = zeros(4,1);
            P(Posvec) = p;
        end
    end
else 
    P = p;
end
kweig = A*(P-PWM(l_svm,:)');


function kweig = ls_kweigtree(mat,negmat,poscell,c,s,ind,indloc,x,l,k,rcnum)
%uses dynamic programming to find the matrix A (kweig in this code) such that A*p = expected gapped kmer vector, where p is the middle column of a 2*l-1 PWM.

%I use a 1st order Markov model to model the flanking bases for a PWM
p = cell(l,1);
p{1} = eye(4);
for i = 1:l-1
    p{i+1} = p{i}*negmat;
end
n = 4^k*max(max(x)); %number of possible k-mers
mat2 = rot90(mat,2);
kweig = zeros(n, 4);
ktree = cell(k-1,1);
ktree2 = cell(k-1,1);
[rx,cx] = size(x);
m=rx;
M = l-1;
X = cx*ones(M+1,1);
for i = 1:cx
    X(i) = i;
end
indloc2 = flipud(indloc);
for i = 2:5
    ktree{i} = zeros(4^i,1);
    ktree2{i} = zeros(4^i,1);
end
for i = 0:M
    if i > M-cx+1
        m = length(c);
    end
    %the following loops is basically dynamic programming for tensor multiplication.  there are multiple cases to consider, hence the if statements.
    for ii = 1:m
        if sum((c(ii,:)+i)==l) > 0 && ~(i == M-1 && ii > rx && ii ~= m)
            indvec = c(ii,:)+i;
            f = find(indvec == l);
            indvec(f) = [];
            loc = indloc(indvec);
            loc2 = indloc2(indvec);
            sPWM = mat(indvec,:).';
            sPWM2 = mat2(indvec,:).';
            ktree{1} = sPWM(:,1);
            ktree2{1} = sPWM2(:,1);
            for iii = s(ii,i+1):k-1
                if loc(iii)==0
                    if loc(iii-1)==1 && indvec(iii-1) < l
                        matt = sPWM(:,iii).'*p{indvec(iii)-l};
                        a = ktree2{iii-1}.*sPWM2(:,iii).';
                        ktree2{iii} = a(:);
                        ktree{iii} = repmat(ktree{iii-1},4,1).*repelem(matt', 4^(iii-1));
                    else
                        matt = p{indvec(iii)-indvec(iii-1)+1};
                        ktree{iii} = repmat(ktree{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2));
                        a = ktree2{iii-1}.*sPWM2(:,iii).';
                        ktree2{iii} = a(:);
                    end
                elseif loc2(iii)==0
                    if loc2(iii-1)==1 && indvec(iii-1) < l
                        matt = sPWM2(:,iii).'*p{indvec(iii)-l};
                        a = ktree{iii-1}.*sPWM(:,iii).';
                        ktree{iii} = a(:);
                        ktree2{iii} = repmat(ktree2{iii-1},4,1).*repelem(matt', 4^(iii-1));
                    else
                        matt = p{indvec(iii)-indvec(iii-1)+1};
                        ktree2{iii} = repmat(ktree2{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2));
                        a = ktree{iii-1}.*sPWM(:,iii).';
                        ktree{iii} = a(:);
                    end
                else
                    a = ktree{iii-1}.*sPWM(:,iii).';
                    ktree{iii} = a(:);
                    a = ktree2{iii-1}.*sPWM2(:,iii).';
                    ktree2{iii} = a(:);
                end
            end
            %the weird indexing that I did early in the code comes to fruition.  It is critical to do so to make this computation as fast as possible.
            if ii <= rx
                for j = 1:X(i+1)
                    if x(ii,j) ~= 0 
                        for iii = 1:2
                            indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(x(ii,j))-1);
                            indvec2 = indvec+(5-2*iii)*4^(f-1);
                            kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1};
                            kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1};
                            kweig(indvec2,iii) = kweig(indvec2,iii) + ktree2{k-1};
                            kweig(indvec,5-iii) = kweig(indvec,5-iii) + ktree2{k-1};
                        end
                    end
                end
            else
                for iii = 1:2
                    indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(ii)-1);
                    indvec2 = indvec+(5-2*iii)*4^(f-1);
                    kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1};
                    kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1};
                    kweig(indvec2,iii) = kweig(indvec2,iii) + ktree2{k-1};
                    kweig(indvec,5-iii) = kweig(indvec,5-iii) + ktree2{k-1};
                end 
            end
        end
    end
end
kweig(4^k*(max(max(x))-rcnum)+1:end,:) = kweig(4^k*(max(max(x))-rcnum)+1:end,:)/sqrt(2);


function kweig = ls_kweigtree_norc(mat,negmat,poscell,c,s,ind,indloc,x,l,k,rcnum)
%uses dynamic programming to find the matrix A (kweig in this code) such that A*p = expected gapped kmer vector, where p is the middle column of a 2*l-1 PWM.

%I use a 1st order Markov model to model the flanking bases for a PWM
p = cell(l,1);
p{1} = eye(4);
for i = 1:l-1
    p{i+1} = p{i}*negmat;
end
n = 4^k*max(max(x)); %number of possible k-mers
mat2 = rot90(mat,2);
kweig = zeros(n, 4);
ktree = cell(k-1,1);
[rx,cx] = size(x);
m=rx;
M = l-1;
X = cx*ones(M+1,1);
for i = 1:cx
    X(i) = i;
end
indloc2 = flipud(indloc);
for i = 2:5
    ktree{i} = zeros(4^i,1);
end
for i = 0:M
    if i > M-cx+1
        m = length(c);
    end
    %the following loops is basically dynamic programming for tensor multiplication.  there are multiple cases to consider, hence the if statements.
    for ii = 1:m
        if sum((c(ii,:)+i)==l) > 0 && ~(i == M-1 && ii > rx && ii ~= m)
            indvec = c(ii,:)+i;
            f = find(indvec == l);
            indvec(f) = [];
            loc = indloc(indvec);
            loc2 = indloc2(indvec);
            sPWM = mat(indvec,:).';
            ktree{1} = sPWM(:,1);
            for iii = s(ii,i+1):k-1
                if loc(iii)==0
                    if loc(iii-1)==1 && indvec(iii-1) < l
                        matt = sPWM(:,iii).'*p{indvec(iii)-l};
                        ktree{iii} = repmat(ktree{iii-1},4,1).*repelem(matt', 4^(iii-1));
                    else
                        matt = p{indvec(iii)-indvec(iii-1)+1};
                        ktree{iii} = repmat(ktree{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2));
                    end
                elseif loc2(iii)==0
                    if loc2(iii-1)==1 && indvec(iii-1) < l
                        a = ktree{iii-1}.*sPWM(:,iii).';
                        ktree{iii} = a(:);
                    else
                        a = ktree{iii-1}.*sPWM(:,iii).';
                        ktree{iii} = a(:);
                    end
                else
                    a = ktree{iii-1}.*sPWM(:,iii).';
                    ktree{iii} = a(:);
                end
            end
            %the weird indexing that I did early in the code comes to fruition.  It is critical to do so to make this computation as fast as possible.
            if ii <= rx
                for j = 1:X(i+1)
                    if x(ii,j) ~= 0 
                        for iii = 1:2
                            indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(x(ii,j))-1);
                            indvec2 = indvec+(5-2*iii)*4^(f-1);
                            kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1};
                            kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1};
                        end
                    end
                end
            else
                for iii = 1:2
                    indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(ii)-1);
                    indvec2 = indvec+(5-2*iii)*4^(f-1);
                    kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1};
                    kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1};
                end 
            end
        end
    end
end

function [pp, len] = adjust_PWM(p,GC)
%extends or truncates the PWM based on information.  I try to do it intelligently, so you may disagree on the condition required for adjustment.
[len,~] = size(p);
info = zeros(len, 1);
cut = 0.2;%maximum information of a column to truncate
ext = 0.7;%minimum information needed to extend
%The rest of the code is easy enough to read through quickly
mat = p+(p==0);
vec = 2+sum(mat.*log(mat)/log(2),2);
b = true;
b2 = true;
while ((vec(1) < cut && max(vec(2:3)) <= ext) || mean(vec(1:3) < cut)) && len > 12
    p(1,:) = [];
    vec(1) = [];
    len = len-1;
    b = false;
    if ((vec(end) < cut && max(vec(end-2:end-1)) <= ext) || mean(vec(end-2:end) < cut)) && len > 12
        vec(end) = [];
        p(end,:) = [];
        len = len-1;
        b2 = false;
    end
end
if b && min(vec(1:2)) > ext && len < 20
   mat = [GC ; mat];
   p = [GC;p];
   len = len+1;
end
while ((vec(end) < cut && max(vec(end-2:end-1)) <= ext) || mean(vec(end-2:end) < cut)) && len > 12
    vec(end) = [];
    p(end,:) = [];
    len = len-1;
    b2 = false;
end
if b2 && min(vec(end-1:end)) > ext && len < 20
    p = [p;GC];
    len = len+1;
end
pp = p;

function info = avg_info(p,l_svm)
info = zeros(length(p),1);
for i = 1:length(p)
    [l,~] = size(p{i});
    mat = p{i}(l_svm:end-l_svm+1,:)+(p{i}(l_svm:end-l_svm+1,:)==0);
    info(i) = sum(2+sum(mat.*log(mat)/log(2),2));
end

function ext_pwm = extendPWM(pwm, n,GCmat)
mat = repmat(GCmat, n,1);
ext_pwm = [mat;pwm;mat];
