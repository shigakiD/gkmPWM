function convert2bed(locsheader, bed)
%convert2bed converts the *_kmer_PWM_locs.out files to bed format
%   
%   convert2bed(locsheader, bed)
%   
%   locsheader: prefix of *_kmer_PWM_locs.out file
%   bed: bedfile in the same order as the fasta input to mapTF
%   
%   Example (files in example_files directory):
%   convert2bed('GM12878', 'GM12878.bed')

fid = fopen([locsheader '_kmer_PWM_locs.out'], 'r');
X = textscan(fid, '%d\t%s\t%d\t%d\t%d\t%f\t%f\t%s\n', 'delimiter', '\t');
fclose(fid);
fid = fopen(bed, 'r');
B = textscan(fid, '%s\t%d\t%d%*[^\n]', 'delimiter', '\t');
fclose(fid);
L = length(X{1});
fid = fopen([locsheader '_kmer_PWM_locs.bed'], 'w');
for i = 1:L
    n = X{1}(i);
    s = B{2}(n)+X{4}(i)-1;
    e = B{2}(n)+X{5}(i);
    fprintf(fid, '%s\t%d\t%d\t%s\t%f\t%f\t%s\n', B{1}{n}, s, e, X{2}{i}, X{6}(i), X{7}(i), X{8}{i});
end
fclose(fid);
