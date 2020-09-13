% author - Terentyev Alexander Nikolaevich. 2008. NTUU KPI IASA
%
% ������� ��� ����᪠ �ணࠬ�� �� ���������
% ��ப� ��⫠�� 7
% rmnk('s_max.txt',3)
% �� ��। �⨬ �㦭� � ����� ...Matlab 7\work\
% �������� 䠩�� rmnk.m � s_max.txt

function rmnk(file_name_data, n_column);

b=10;

fid_sourse = fopen(file_name_data,'r');
B = fscanf(fid_sourse,'%g',[n_column inf]);
A = B';
y=A(:,1);
x=A(:,2:end);

[x_row, x_column]=size(x);

% ������ ��������� �������� p0
p0=eye(x_column)*b;

% ������ ��������� �������� ���� (� ��� coef)
for i=1:x_column
    coef(i)=0;    
end;  
coef=coef';

i=1;
for i=1:x_row
    fprintf('�������� - %i',i);
    p1=p0-(p0*x(i,:)'*x(i,:)*p0)/(1+x(i,:)*p0*x(i,:)'),
    coef=coef+p1*x(i,:)'*(y(i)'-x(i,:)*coef),
    p0=p1;
end;

coef,