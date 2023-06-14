function [c,C,I,ind,mat,rcnum] = genIndex(l,k)
%l = k-mer length
%k = # of ungapped positions

c = combnk(1:l,k);
d = ones(length(c),1);
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
for i = 1:length(c)
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
S = sum(C(:,1)==1);
mat = zeros(S, max(c(:,1))-1);
for i = 1:S
    if C(i,end) ~= l
        for j = 2:max(c(:,1));
            a = S+1;
            while sum(abs(C(i,:)+j-1-C(a,:))) ~= 0
                a = a + 1;
                if a > length(c)
                    break
                end
            end
            if a <= length(c)
                mat(i,j-1) = a;
            end
        end
    end
end
mat = [(1:S)' mat];
I = zeros(length(C),1);
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
for i = 1:length(c)
    if sum(abs(c(ind(i),:)-C(i,:)))~=0
        c(ind(i),:) = l+1-fliplr(c(ind(i),:));
    end
end
