function [negweights] = BGkmer(mat,GC,c,rcnum,l,k,RC)
len = numel(c)/k;
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
if RC
    seqvec2 = 5-fliplr(seqvec);
    c2 = l+1-fliplr(c);
end
a = 1;
for i = 1:len
    dc = diff(c(i,:));
    startvec = GCmat*p{c(i,1)};
    if RC
        dc2 = diff(c2(i,:));
        startvec2 = GCmat*p{c2(i,1)};
    end
    for ii = 1:4^k
        negweights(a) = startvec(seqvec(ii,1));
        if RC
            tmp = startvec(seqvec2(ii,1));
        end
        for iii = 1:k-1
            matt = p{dc(iii)+1};
            negweights(a) = negweights(a)*matt(seqvec(ii,iii), seqvec(ii,iii+1));
            if RC
                matt2 = p{dc2(iii)+1};
                tmp = tmp*matt2(seqvec2(ii,iii), seqvec2(ii,iii+1));
            end
        end
        if RC
            negweights(a) = (negweights(a)+tmp);
        end
        a = a+1;
    end
end
if rcnum > 0 && RC
    negweights(4^k*alen+1:end) = negweights(4^k*alen+1:end)/sqrt(2);
end
