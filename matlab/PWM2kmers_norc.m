function [kweig] = PWM2kmers_norc(mat,negmat,c,s,ind,indloc,x,l,k,rcnum)
%makes gkm-pwm faster
p = cell(l,1);
p{1} = eye(4);
for i = 1:l-1
    p{i+1} = p{i}*negmat;
end
m = length(x);
M = length(mat)-l;
n = 4^k*numel(c)/k; %number of possible k-mers
mat2 = rot90(mat,2);
kweig = zeros(n,1);
ktree = cell(k,1);
[rx,cx] = size(x);
X = cx*ones(length(mat)-l+1,1);
for i = 1:cx
    X(i) = i;
end
for i = 2:k
    ktree{i} = zeros(4^i,1);
end
a = 0;
for i = 0:M
    if i == M-1
        m = numel(c)/k;
    end
    for i2 = 1:m
        if ~(i == M-1 && i2 > rx && i2 ~= m)
            indvec = c(i2,:)+i;
            loc = indloc(indvec);
            sPWM = mat(indvec,:).';
            ktree{1} = sPWM(:,1);
            for i3 = s(i2):k
                if loc(i3)==0
                    if loc(i3-1) == 1
                        matt = mat(1,:)*p{indvec(i3)-indvec(i3-1)+1};
                        for i4 = 1:4
                            ktree{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4)) = ktree{i3-1}*matt(i4);
                        end
                    else
                        matt = p{indvec(i3)-indvec(i3-1)+1};
                        ktree{i3} = repmat(ktree{i3-1}, 4, 1).*repelem(matt(:), 4^(i3-2));
                    end
                else
                    a = ktree{i3-1}.*sPWM(:,i3).';
                    ktree{i3} = a(:);
                end
            end
            if i2 <= rx
                for j = 1:X(i+1)
                    if x(i2,j) ~= 0
                        kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) = kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) + ktree{k};
                    end
                end
            else
                kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) = kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) + ktree{k};
            end
        end
    end
end
%alen = length(c)-rcnum;
%kweig(4^k*alen+1:end) = kweig(4^k*alen+1:end)/sqrt(2);
