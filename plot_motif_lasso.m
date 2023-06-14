function plot_motif_lasso(varargin)

inp=varargin{1};
fid = fopen([varargin{1} '_gkmPWMlasso.out'], 'r');
X = textscan(fid, '%f\t%f\t%s\t%f\t%f\t%f\t%f\n','HeaderLines', 4, 'delimiter', '\t');
fclose(fid);
mf = fopen('motif_logos/short_name', 'r');
names = textscan(mf,'%*d\t%s\t%*s\n');
fclose(mf);
nmot=X{1}(end);
IND = zeros(nmot,1);
for i = 1:nmot
    f = find(X{1}==i);
    IND(i) = f(1);
end
figure('units','inches','position',[8 0 8 7*nmot/30])
hold on
scale=85;
axis([-1200,1500,0,(nmot+1)*scale])
text(-1.8*scale,(nmot+.5)*scale,'I')
text(-4.0*scale,(nmot+.5)*scale,'Z')
text(1.2*scale,(nmot+.5)*scale,'W')
text(2.3*scale,(nmot+.5)*scale,'Z')
text(3.3*scale,(nmot+.5)*scale,'I')
text(-6.2*scale,(nmot+.5)*scale,'W')
text(-12.8*scale,(nmot+.5)*scale,'N')
text(-11.6*scale,(nmot+.5)*scale,'ID')
text(-9.8*scale,(nmot+.5)*scale,'Motif')
zmax = max(abs(X{5}));
emax = max(abs(X{6}));
wmax = max(abs(X{4}));
for i=1:nmot 
    if X{4}(IND(i)) > 0
        rectangle('position',85*[1 nmot-i 1 1],'facecolor',[1-X{4}(IND(i))/wmax 1-X{4}(IND(i))/wmax 1])
    else
        rectangle('position',85*[1 nmot-i 1 1],'facecolor',[1 1+X{4}(IND(i))/wmax 1+X{4}(IND(i))/wmax])
    end
    if X{5}(IND(i)) > 0
        rectangle('position',85*[2 nmot-i 1 1],'facecolor',[1-X{5}(IND(i))/zmax 1-X{5}(IND(i))/zmax 1])
    else
        rectangle('position',85*[2 nmot-i 1 1],'facecolor',[1 1+X{5}(IND(i))/zmax 1+X{5}(IND(i))/zmax])
    end
    rectangle('position',85*[3 nmot-i 1 1],'facecolor',[1-X{6}(IND(i))/emax 1-X{6}(IND(i))/emax 1])
    smot=sprintf('motif_logos/m%d.png',X{2}((IND(i))));
    img=imread(smot);
    p=[2 nmot-i 1 1];
    mlab=names{1}{X{2}((IND(i)))};
    e=min(8,length(mlab));
    text(-9.8*scale,(p(2)+.5)*scale,mlab(1:e))
    text(-12.8*scale,(p(2)+.5)*scale,sprintf('%d',i))
    text(-11.6*scale,(p(2)+.5)*scale,sprintf('%d',X{2}(IND(i))))
    text(-4.0*scale,(p(2)+.5)*scale,sprintf('%0.2f',X{5}(IND(i))))
    text(-1.8*scale,(p(2)+.5)*scale,sprintf('%0.2f',X{6}(IND(i))))
    text(-6.2*scale,(p(2)+.5)*scale,sprintf('%0.2f',X{4}(IND(i))))
    p1=[p(1) p(2)+5]*20;
    p2=p1+[400 100];
    image('XData',(p(1)+2.3)*scale,'YData',(p(2)-0.25)*scale,'CData',flipud(img)) 
end
axis off
hold off
fig=gcf;   
fig.PaperUnits='inches';
fig.PaperSize=[8 7*nmot/30]*1.2;
pfile=[varargin{1} '_gkmPWMlasso.pdf'];
print(pfile,'-dpdf','-fillpage')
close all
