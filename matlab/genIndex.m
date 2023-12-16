function [c,C,I,ind,mat,rcnum] = genIndex(l,k,n_frac)
%l = k-mer length
%k = # of ungapped positions
if n_frac > 1  || n_frac < 0
    error('n_frac must be a fraction in [0,1]')
end
if n_frac ~= 1
    rng(1)
    if mod(l,2) == 0
        c = combnk(1:l,k);
        L = numel(c)/k;
        [c,C,I,ind,mat,rcnum] = sort_comb(c,l,k,n_frac);
        L = numel(c)/k;
        c = genIndex_frac(l,k,c,rcnum,L*n_frac);
        [c,C,I,ind,mat,rcnum] = sort_comb(c,l,k,n_frac);
    else
        c = combnk(1:l,k);
        [c,~,~,~,~,rcnum1] = sort_comb(c,l,k,n_frac);
        c1 = combnk(1:(l-1),k-1);
        [c1,~,~,~,~,rcnum1] = sort_comb(c1,l-1,k-1,n_frac);
        L1 = numel(c1)/(k-1);
        c1 = genIndex_frac(l-1,k-1,c1,rcnum1,L1*n_frac);
        L1 = numel(c1)/(k-1);
        x = zeros(L1, l-1);
        for i = 1:L1
            x(i,c1(i,:)) = x(i,c1(i,:))+1;
            x(i,l-fliplr(c1(i,:))) = x(i,l-fliplr(c1(i,:)))+1;
        end
        L2 = round((L1*(l-1)*2-sum(sum(x)))/k/2);
        c2 = combnk(1:(l-1),k);
        [c2,~,~,~,~,rcnum2] = sort_comb(c2,l-1,k,n_frac);
        c2 = genIndex_frac(l-1,k,c2,rcnum2,L2);
        L2 = numel(c2)/k;
        cc1 = zeros(L1,k);
        M = ceil(l/2);
        for i = 1:L1
            f = find(c1(i,:)>=M);
            c1(i,f) = c1(i,f)+1;
            cc1(i,:) = sort([c1(i,:) M]);
        end
        for i = 1:L2
            f = find(c2(i,:)>=M);
            c2(i,f) = c2(i,f)+1;
        end
        [c,C,I,ind,mat,rcnum] = sort_comb([cc1;c2],l,k,n_frac);
    end
else
    c = combnk(1:l,k);
   [c,C,I,ind,mat,rcnum] = sort_comb(c,l,k,n_frac);
end

function [c,C,I,ind,mat,rcnum] = sort_comb(c,l,k,n_frac)
L = numel(c)/k;
d = ones(L,1);
e = []; 
a = 1;
for i = 1:length(d)-1;
    vec = l+1-fliplr(c(i,:));
    if d(i) ~= 0
        if sum(abs(c(i,:)-vec))==0
            e(a) = i;
            d(i) = 0;
            a = a+1;
        else
            for j = i+1:length(d)
                if sum(abs(c(j,:)-vec))==0
                    d(j) = 0;
                end
            end
        end
    end
end
rcnum = a-1;
f = find(d==1);
c=[c(f,:);c(e,:)];
C = c;
L = numel(c)/k;
for i = 1:L
    c2 = l+1-fliplr(c(i,:));
    if sum(0.5.^c(i,:)) < sum(0.5.^c2)
        C(i,:) = c2;
    end
end
c=C;
[~,ind] = sort(sum(0.5.^C,2),'descend');
C = C(ind,:);
f = find(C(:,1)==1);
S = length(f);
ff = find(C(f,end)~=l);
for i = 1:length(ff)-1
    for j = i+1:length(ff)
        if sum(abs(C(f(ff(i)),:)+1-fliplr(l+1-C(f(ff(j)),:))))==0 
           C(f(ff(j)),:) = fliplr(l+1-C(f(ff(j)),:));
        end
    end
end
[~,ind2] = sort(sum(0.5.^C,2),'descend');
C = C(ind2,:);
ind = ind(ind2);
if n_frac == 1
    S = sum(C(:,1)==1);
    mat = zeros(S, max(c(:,1))-1);
    for i = 1:S
        if C(i,end) ~= l
            for j = 2:max(c(:,1));
                a = S+1;
                while sum(abs(C(i,:)+j-1-C(a,:))) ~= 0
                    a = a + 1;
                    if a > L
                        break
                    end
                end
                if a <= L
                    mat(i,j-1) = a;
                end
            end
        end
    end
    mat = [(1:S)' mat];
else
    mat = (1:L)';
end
I = zeros(L,1);
I(1)=2;
for i=2:length(I)
    a = 1;
    while C(i-1,a)==C(i,a)
        a = a+1;
    end
    I(i) = a;
    if I(i) < 2
        I(i)=2;
    end
end
for i = 1:L
    if sum(abs(c(ind(i),:)-C(i,:)))~=0
        c(ind(i),:) = l+1-fliplr(c(ind(i),:));
    end
end



function mat = genIndex_frac(l,k,c,rc,Lfrac);

L = numel(c)/k;
x = zeros(L, l);
for i = 1:L
    x(i,c(i,:)) = x(i,c(i,:))+1;
    x(i,l+1-fliplr(c(i,:))) = x(i,l+1-fliplr(c(i,:)))+1;
end
M = l/2;
y = sum(x==2,2);
nind = true;
ind = 1:L;
X = zeros(L,1);
count = 0;
i = 0;
while count == 0 || count < Lfrac
    if Lfrac-count < count+M-Lfrac
        break
    end
    count = count + 1;
    r = randi(L-count);
    cvec = c(ind(r),:);
    X(count) = ind(r);
    ind(r) = [];
    num = sum(cvec<M+1);
    for j = 1:M-1
        C = zeros(1,k);
        C(1:num) = sort(mod(cvec(1:num)+(j-1),M)+1);
        C(num+1:end) = sort(l-mod(l+1-cvec(num+1:end)+(j-1),M));
        mat = repmat(C, L-count,1);
        f = find(sum(abs(c(ind,:)-mat),2)==0);
        if isempty(f)
            R = l+1-fliplr(C);
            mat = repmat(R,L-count,1);
            f = find(sum(abs(c(ind,:)-mat),2)==0);
            if isempty(f)
                break
            else
                count = count + 1;
                X(count)=ind(f);
                ind(f)=[];
            end
        else
            count = count + 1;
            X(count)=ind(f);
            ind(f)=[];
        end
    end
end
X = X(1:count);
X = sort(X);
mat = c(X,:);

