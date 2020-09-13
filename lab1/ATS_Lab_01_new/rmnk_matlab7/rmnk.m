% author - Terentyev Alexander Nikolaevich. 2008. NTUU KPI IASA
%
% команда для запуска программы из командной
% строки Матлаба 7
% rmnk('s_max.txt',3)
% но перед этим нужно в папку ...Matlab 7\work\
% поместить файлы rmnk.m и s_max.txt

function rmnk(file_name_data, n_column);

b=10;

fid_sourse = fopen(file_name_data,'r');
B = fscanf(fid_sourse,'%g',[n_column inf]);
A = B';
y=A(:,1);
x=A(:,2:end);

[x_row, x_column]=size(x);

% чрфрхь эрўры№э√х чэрўхэш  p0
p0=eye(x_column)*b;

% чрфрхь эрўры№э√х чэрўхэш  Є¤Єр (є эрё coef)
for i=1:x_column
    coef(i)=0;    
end;  
coef=coef';

i=1;
for i=1:x_row
    fprintf('╚ЄхЁрЎш  - %i',i);
    p1=p0-(p0*x(i,:)'*x(i,:)*p0)/(1+x(i,:)*p0*x(i,:)'),
    coef=coef+p1*x(i,:)'*(y(i)'-x(i,:)*coef),
    p0=p1;
end;

coef,