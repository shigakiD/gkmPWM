function getgkmweights(varargin)
% getgkmweights saves the gapped kmers weights for a given (l,k)
%
%     getgkmweights(fileprefix, l, k ...)
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
%     'l'             The full length of the gapped k-mer.
%     'k'             The number of ungapped positions of the gapped k-mer.
%
%     Name Value Pair Parameters (Optional):
%
%     'RC'            If true, treat reverse complementary gapped k-mers as
%                     the same feature.  Otherwise treat them as separate features
%                     (default:true)
%     'KmerFrac'      Set the fraction of the total number of gapped k-mers to
%                     use with gkmPWM.  This reduces the memory and runtime
%                     needed.  If the total number of gapped k-mers is too high
%                     with the given combination of (l,k,KmerFrac), KmerFrac will
%                     be automatically set to a lower value to create a more
%                     workable number of gapped k-mers
%     'KmerFracLimit' Automatically lower KmerFrac if the number of gapped kmers
%.                    is too large.  Note that this will require more memory and
%                     increase runtime if set to false. (default: true)

if nargin < 3
    error('Need at least 3 inputs')
elseif nargin > 3
    if mod(nargin,2) == 0
        error('Incorrect number of inputs')
     else
        vec = 4:2:nargin;
        inputlib = {'RC', 'KmerFrac','KmerFracLimit'};
        for i = 1:length(vec)
            f = strcmp(varargin{vec(i)},inputlib);
            if sum(f) == 0
                error([varargin{vec(i)} ' is not an input option'])
            end
        end
    end
end
fn = varargin{1};
l_svm = varargin{2};
k_svm = varargin{3};
RC = true;
nfrac = 1;
nfracLim = true;
lk = 1;
if nargin > 3
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
        if nfrac < 1
            lk = [l_svm k_svm];
        end
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
[comb,comb2,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm,nfrac);
if nfracLim && length(comb)*4^k_svm > 5*10^5
    nfrac = round(5*10^7/4^k_svm/numel(comb)*k_svm)/100;
    lk = ([l_svm k_svm]);
    [comb,comb2,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm,nfrac);
end
disp(['Using ' num2str(numel(comb)/k_svm*4^k_svm) ' gapped kmers'])
disp('Counting gapped kmers')
[cfile, GCpos1, GCneg1,mat,mat2] = getgkmcounts(fn, l_svm, k_svm, lk, RC, comb,rcnum);
negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC);
cfile = cfile-negvec/sum(negvec)*sum(cfile);
cfile = cfile/std(cfile);
s = '';
for i = 1:k_svm
    s = [s 'A'];
end
kmer = repmat(s,4^k_svm,1);
vec = 0:4^k_svm-1;
for i = 1:k_svm
    vec2 = mod(floor(vec/4^(i-1)),4);
    f = find(vec2==1);
    kmer(f,i) = 'C';
    f = find(vec2==2);
    kmer(f,i) = 'G';
    f = find(vec2==3);
    kmer(f,i) = 'T';
end
s = '';
for i = 1:l_svm
    s = [s '-'];
end
a = 1;
fid = fopen([fn '_' num2str(l_svm) '_' num2str(k_svm) '_gkmweights.out'], 'w');
for i = 1:numel(comb)/k_svm
    for j = 1:4^k_svm
        S = s;
        S(comb(i,:)) = kmer(j,:);
        fprintf(fid, '%s\t%0.5f\n', S, cfile(a));
        a = a+1;
    end
end
fclose(fid);
dlmwrite([fn '_negmat.out'], [mat; GCpos1 GCneg1 0 0], 'precision',10);
