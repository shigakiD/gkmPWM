function process_motifs(dfn, lfn, ofn)
%dfn: file name for denovo motifs
%lfn: file name for lasso motifs
%ofn: output filename
a = 1;
LEN = zeros(1,1);
shift = zeros(1,1);
[p,w] = getdenovomotif(dfn);
N = numel(w);
fid = fopen(lfn,'r');
X = textscan(fid,'%f\t%f\t%s\t%f\t%f\t%f\n','delimiter', '\t', 'headerlines', 4);
fclose(fid);
n = X{1}(end);
vec = zeros(n,1);
vec2 = zeros(n,1);
w = [w; zeros(n,1)];
for i = 1:n
    f = find(X{1}==i);
    vec(i) = X{2}(f(1));
    vec2(i) = X{5}(f(1));
    w(i+N) = X{4}(f(1));
end
P = getmotif('combined_db_v4.meme',vec);
[pp, info, len] = trim_pwm([p;P],0.25);
PWM2 = cell(1,1);
LEN_2 = zeros(1,1);
I_2 = zeros(1,1);
hal = true;
for ii = 1:length(w)
    if ii > N
        hal = true;
        [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]);
        if cor > 0.7 || vec2(ii-N) < 1.5
            hal = false;
        end
    elseif ii > 1 && a > 2
        hal = true;
        [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]);
        if cor > 0.7
            hal = false;
        end
    end
    if hal && w(ii) > 0 && len(ii) >= 6
        if len(ii) > 10
            if info(ii)/len(ii) > 0.7
                PWM2{a} = pp{ii};
                LEN_2(a) = len(ii);
                I_2(a) = info(ii);
                a = a+1;
                PWM2{a} = rot90(PWM2{a-1},2);
                LEN_2(a) = LEN_2(a-1);
                I_2(a) = I_2(a-1);
                a = a+1;
            end
        elseif info(ii) > 6 || info(ii)/len(ii) > 1
            PWM2{a} = pp{ii};
            LEN_2(a) = len(ii);
            I_2(a) = info(ii);
            a = a+1;
            PWM2{a} = rot90(PWM2{a-1},2);
            LEN_2(a) = LEN_2(a-1);
            I_2(a) = I_2(a-1);
            a = a+1;
        end
    end
end
p = getmotif('combined_db_v4.meme', 1:1968);
[p,info,lenvec] = trim_pwm(p,0.25);
fid = fopen('motif_logos/motif_names_clus','r');
X = textscan(fid, '%s\n');
fclose(fid);
fid = fopen(ofn, 'w');
a = 1;
for i = 1:length(PWM2)
    [ind, r] = ppmsim([PWM2{i};p], [LEN_2(i);lenvec]);
    if r > 0.80
        fprintf(fid,'MOTIF %d\n%s\n%d\n', a, X{1}{ind},LEN_2(i));
        a = a+1;
        for j = 1:LEN_2(i)
            fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4));
        end
        fprintf(fid, '\n');
    elseif I_2(i)/LEN_2(i) > 1
        fprintf(fid,'MOTIF %d\n%s\n%d\n', a, consen(PWM2{i}, LEN_2(i)), LEN_2(i));
        a = a+1;
        for j = 1:LEN_2(i)
            fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4));
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);

function [ind, M] = ppmsim(mot,lenvec)
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

function [mat,w] = getdenovomotif(filename)
% filename is the meme file that contains the motifs;
% n is the nth motif in the file
mat = {};
w = [];
fid = fopen(filename);
if fid < 0
    fprintf('No good')
else
   i=0;
   while ~feof(fid)
        line = fgetl(fid);
        if length(line) >= 5
            if strcmp(line(1:5), 'MOTIF')
                i = i+1;
                line = fgetl(fid);
                mat{i} = [];
                [~, tmp] = strtok(line);
                [tmp] = strtok(tmp);
                w = [w;str2double(tmp)];
                while ~isempty(line)
                    mat{i} = [mat{i}; str2num(line)];
                    line = fgetl(fid);
                end
            end
        end
    end
end
mat = mat';
fclose(fid);

function [pp, info, len] = trim_pwm(p,cut)
l = length(p);
info = zeros(l, 1);
len = zeros(l,1);
for i = 1:l
    mat = p{i}+(p{i}==0);
    vec = 2+sum(mat.*log(mat)/log(2),2);
    while (vec(1) < cut || mean(vec(1:3)) < cut || mean(vec(2:4)) < cut) && length(vec) > 4
        p{i}(1,:) = [];
        vec(1) = [];
    end
    while (vec(end) < cut || mean(vec(end-2:end)) < cut || mean(vec(end-3:end-1)) < cut) && length(vec) > 4
        vec(end) = [];
        p{i}(end,:) = [];
    end
    info(i) = sum(vec);
    [len(i), ~] = size(p{i});
end
pp = p;

function s = consen(p, l)
s='';
for i = 1:l
    s(i)='A';
    r = norm(p(i,:)'- [1 0 0 0]');
    rr =  norm(p(i,:)'- [0 1 0 0]');
    if rr < r
        s(i) = 'C';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 0 1 0]');
    if rr < r
        s(i) = 'G';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 0 0 1]');
    if rr < r
        s(i) = 'T';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/2 1/2 0 0]');
    if rr < r
        s(i) = 'M';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/2 0 1/2 0]');
    if rr < r
        s(i) = 'R';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/2 0 0 1/2]');
    if rr < r
        s(i) = 'W';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 1/2 1/2 0]');
    if rr < r
        s(i) = 'S';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 1/2 0 1/2]');
    if rr < r
        s(i) = 'Y';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 0 1/2 1/2]');
    if rr < r
        s(i) = 'K';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/3 1/3 1/3 0]');
    if rr < r
        s(i) = 'V';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/3 1/3 0 1/3]');
    if rr < r
        s(i) = 'H';
        r = rr;
    end
    rr =  norm(p(i,:)'- [1/3 0 1/3 1/3]');
    if rr < r
        s(i) = 'D';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0 1/3 1/3 1/3]');
    if rr < r
        s(i) = 'B';
        r = rr;
    end
    rr =  norm(p(i,:)'- [0.25 0.25 0.25 0.25]');
    if rr < r
        s(i) = 'N';
        r = rr;
    end
end
