function [negweights] = BGkmer(mat,GC,c,rcnum,l,k)
len = length(c);
alen = len-rcnum;
negweights = zeros(len*4^k,1);
GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2];
p = cell(l,1);
p{1} = eye(4);
for i = 1:l-1
    p{i+1} = p{i}*mat;
end
seqvec = zeros(4^k, k);
for i = 1:k
    for j = 1:4^k
        seqvec(j,i) = mod(floor((j-1)/4^(i-1)), 4)+1;
    end
end
seqvec2 = 5-fliplr(seqvec);
a = 1;
c2 = l+1-fliplr(c);
for i = 1:len
    dc = diff(c(i,:));
    dc2 = diff(c2(i,:));
    startvec = GCmat*p{c(i,1)};
    startvec2 = GCmat*p{c2(i,1)};
    for ii = 1:4^k
        negweights(a) = startvec(seqvec(ii,1));
        tmp = startvec(seqvec2(ii,1));
        for iii = 1:k-1
            matt = p{dc(iii)+1};
            matt2 = p{dc2(iii)+1};
            negweights(a) = negweights(a)*matt(seqvec(ii,iii), seqvec(ii,iii+1));
            tmp = tmp*matt2(seqvec2(ii,iii), seqvec2(ii,iii+1));
        end
        negweights(a) = (negweights(a)+tmp);
        a = a+1;
    end
end
if rcnum > 0
    negweights(4^k*alen+1:end) = negweights(4^k*alen+1:end)/sqrt(2);
end
