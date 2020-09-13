{
������ Matrices - ����� �������� ��� ������ � ���������.
��� ������������� ������ ���������� ����� Delphi 6 ��� ����� ������� ������.

������ � �������: ������� �������� � ���������� ���� Matrix.
��� ������ ������ � �������� ���������� ������ �� ������ � �������
������� SetSize (�� ����� � ���������� �������� ��������, ��� ���� ���
���������� ������� ������� ������������ � ��� �������� �����������).
��������� � �������� �������: A.Data[i,j], ��� A - ��� ����������,
i - ����� ������ � �������, j - ����� ������� � �������.
��� ������, ����������� � ������� ������� ������ (��������, ���������),
�������� SetSize �� ���������.
��������������� ����:
TMap - ������������ ������ � ��������� �������� (�����) �������,
Vector - ������������ ������, ���������� � �������� ������ ���� ���������
���������� ����� � �������� � �������, � � ��������� - �������� �������, �
���������� � ������ ��������������� (��� ������ ��� ���������� � ���� ���
�������� �� ����).

������ ������� ��� ������ - ��. ��������� ������

� ������ ������������� ��������� ����������, ��������� � ���������������
������� ������ (��������, ��� ��������� ��� �����������).
}

unit matrices;

interface
Uses SysUtils;

type
  Matrix = record
    M, N   : integer;  { M - ����� �����, N - �������� }
    Data   : array of array of extended;
  end;
{$IFNDEF XP_CMATRIX}
{$DEFINE XP_MATRIX}
  TMap = array of integer;
{$ENDIF}
  Vector = array of extended;

  Procedure SetSize(var A:Matrix; M,N:integer); overload; {������� ������� �������}
  Procedure Zero(var A:Matrix); overload; {���������� ������� ������ }
  Procedure E(var A:Matrix); overload; { ��������� ������� }
  Function Transpose(const A:Matrix):Matrix; overload; { ���������������� ������� }
  Function Add(const A,B:Matrix):Matrix; overload; { �������� ������ }
  Function Sub(const A,B:Matrix):Matrix; overload; { ���������� ������ }
  Function Mul(const A,B:Matrix):Matrix; overload; { ��������� ������ }

  Function AddValue(const A:Matrix; Value:extended):Matrix; overload; { �������� ������� � ������ }
  Function MulValue(const A:Matrix; Value:extended):Matrix; overload; { ���������� ������� �� ����� }
  Function Negate(const A:Matrix):Matrix; overload; { ���������� ����� ��������� �������}

  Procedure DeleteRow(var A:Matrix; Row:integer); overload; { �������� ������ }
  Procedure DeleteCol(var A:Matrix; Col:integer); overload; { �������� ������� }
  Procedure DeleteCross(var A:Matrix; Row:integer); overload; { �������� ������ � ������� � �������� ������� }
  Procedure InsertRow(var A:Matrix; Row:integer); overload; { ���������� ������ ����� � ������� }
  Procedure InsertCol(var A:Matrix; Col:integer); overload; { ���������� ������� ����� � ������� }

  Function MapMatrix(const A:Matrix; const Map:TMap):Matrix; overload; { ��������� ������� ��� �������� ����� � �������� }
  Function UnmapMatrix(const A:Matrix; const Map:TMap):Matrix; overload; { �������������� ������� �� ��������� �������}
  Function MapVector(const A:Matrix; const Map:TMap):Matrix; overload; { ��������� ������� ��� �������� ����� }
  Function UnmapVector(const A:Matrix; const Map:TMap):Matrix; overload; { �������������� ��������� ���������� ����� � �������}

  Function SubMatrix(const A:Matrix; const MapX,MapY:TMap):Matrix; overload; { ��������� ����������, ���������� �������� ������� � ������ }
  Function RangeMatrix(const A:Matrix; StartX,EndX,StartY,EndY:integer):Matrix; overload; { ��������� ����������, ���������� �������� �������� �������� � ����� }

  Function Vectorize(const A:Matrix):Vector; overload; { �������������� ������� � ������, ������������� � ������ ��������������� }
  Function Unvectorize(const V:Vector):Matrix; overload; { �������� ��������������: ������ � ������� }

  Function Inverse(const A:Matrix):Matrix; overload; { ��������� �������� ������� ������� �������-������ }
  Function Trace(const A:Matrix):extended; overload; { ��������� ����� ������� (����� ������������ ���������) }

  Procedure JoinBottom(var A:Matrix; const B:Matrix); overload; { ������������� ������� ����� }
  Procedure JoinRight(var A:Matrix; const B:Matrix); overload; { ������������� ������� ������ }
  Procedure JoinDiag(var A:Matrix; const B:Matrix); overload; { ������������� ������� � ������ ������ ���� }

implementation

  Procedure SetSize(var A:Matrix; M,N:integer);
  var
    i  : integer;
  Begin
    A.M:=M;
    A.N:=N;
    SetLength(A.Data,M);
    for i:=0 to M-1 do begin
      SetLength(A.Data[i],N);
    end;
  End;

  Procedure Zero(var A:Matrix);
  var
    i,j  : integer;
  Begin
   for i:=0 to A.M-1 do for j:=0 to A.N-1 do A.Data[i,j]:=0;
  End;

  Procedure E(var A:Matrix);
  var
    i  : integer;
  Begin
    if (A.M<>A.N) then Raise Exception.Create('������� ������������ ������������ ��������� �������!');
    Zero(A);
    for i:=0 to A.M-1 do A.Data[i,i]:=1;
  End;

  Function Transpose(const A:Matrix):Matrix;
  var
   B      :  Matrix;
   i,j    :  integer;
  Begin
    SetSize(B,A.N,A.M);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do B.Data[j,i]:=A.Data[i,j];
    Result:=B;
  End;

  Function Add(const A,B:Matrix):Matrix;
  var
    C    :  Matrix;
    i,j  :  integer;
  Begin
    if ((A.M<>B.M) or (A.N<>B.N)) then Raise Exception.Create('������� ������ ��� �������� �� ���������!');
    SetSize(C,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do C.Data[i,j]:=A.Data[i,j]+B.Data[i,j];
    Result:=C;
  End;

  Function Sub(const A,B:Matrix):Matrix;
  var
    C    :  Matrix;
    i,j  :  integer;
  Begin
    if ((A.M<>B.M) or (A.N<>B.N)) then Raise Exception.Create('������� ������ ��� �������� �� ���������!');
    SetSize(C,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do C.Data[i,j]:=A.Data[i,j]-B.Data[i,j];
    Result:=C;
  End;

  Function Mul(const A,B:Matrix):Matrix;
  var
    C       : Matrix;
    i,j,k   : integer;
  Begin
    if (A.N<>B.M) then Raise Exception.Create('������� ������ ��� ��������� �� ���������!');
    SetSize(C,A.M,B.N);
    Zero(C);
    for i:=0 to A.M-1 do for j:=0 to B.N-1 do for k:=0 to A.N-1 do C.Data[i,j]:=C.Data[i,j]+A.Data[i,k]*B.Data[k,j];
    Result:=C;
  End;

  Function Inverse(const A:Matrix):Matrix;
  var
    i,j,k    :  integer;
    B        :  Matrix;
    sk, sz   :  extended;
  Begin
    if (A.N<>A.M) then Raise Exception.Create('������� ��������� �������� ������� ��� ������������ �������!');
    SetSize(B,A.M,A.M);
    Zero(B);
    for i:=0 to A.M-1 do B.Data[i,i]:=1;
    for i:=0 to A.M-1 do begin
      sk:=1/A.Data[i,i];
      for j:=0 to A.M-1 do if (i<>j) then begin
        sz:=sk*A.Data[j,i];
        for k:=0 to A.M-1 do A.Data[j,k]:=A.Data[j,k]-sz*A.Data[i,k];
        for k:=0 to A.M-1 do B.Data[j,k]:=B.Data[j,k]-sz*B.Data[i,k];
      end;
      for k:=0 to A.M-1 do A.Data[i,k]:=sk*A.Data[i,k];
      for k:=0 to A.M-1 do B.Data[i,k]:=sk*B.Data[i,k];
    end;
    Result:=B;
  End;

  Function Trace(const A:Matrix):extended;
  var
    i   : integer;
    res : extended;
  Begin
    res:=0;
    if (A.N<>A.M) then Raise Exception.Create('������� ��������� ���� ��� ������������ �������!');
    for i:=0 to A.M-1 do res:=res+A.Data[i,i];
    Result:=res;
  End;

  Procedure JoinBottom(var A:Matrix; const B:Matrix);
  var
    i,j,oldy  : integer;
  Begin
    if (A.N<>B.N) then Raise Exception.Create('���������� ���������� ��� ������� �� ���������. �� ������� �� ���������!');
    oldy:=A.M;
    SetSize(A,A.M+B.M,A.N);
    for i:=0 to B.M-1 do for j:=0 to B.N-1 do A.Data[i+oldY,j]:=B.Data[i,j];
  End;

  Procedure JoinRight(var A:Matrix; const B:Matrix); { ������������� ������� ������ }
  var
   i,j,oldx  : integer;
  Begin
    if (A.M<>B.M) then Raise Exception.Create('���������� ���������� ��� ������� �� �����������. �� ������� �� ���������!');
    oldx:=A.N;
    SetSize(A,A.M,A.N+B.N);
    for i:=0 to B.M-1 do for j:=0 to B.N-1 do A.Data[i,j+oldX]:=B.Data[i,j];
  End;

  Procedure JoinDiag(var A:Matrix; const B:Matrix); { ������������� ������� � ������ ������ ���� }
  var
    i,j,oldx,oldy  : integer;
  Begin
    oldY:=A.M;
    oldX:=A.N;
    SetSize(A,A.M+B.M,A.N+B.N);
    for i:=0 to B.M-1 do for j:=0 to B.N-1 do A.Data[i+oldY,j+oldX]:=B.Data[i,j];
    for i:=0 to oldY-1 do for j:=oldX to A.N-1 do A.Data[i,j]:=0;
    for i:=oldY to A.M-1 do for j:=0 to oldX-1 do A.Data[i,j]:=0;
  End;

  Function Negate(const A:Matrix):Matrix;
  var
    i,j   : integer;
    B     : Matrix;
  Begin
    SetSize(B,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do B.Data[i,j]:=-A.Data[i,j];
    Result:=B;
  End;

  Function AddValue(const A:Matrix; Value:extended):Matrix; { �������� ������� � ������ }
  var
    i,j   : integer;
    B     : Matrix;
  Begin
    SetSize(B,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do B.Data[i,j]:=A.Data[i,j]+Value;
    Result:=B;
  End;

  Function MulValue(const A:Matrix; Value:extended):Matrix; { ���������� ������� �� ����� }
  var
    i,j   : integer;
    B     : Matrix;
  Begin
    SetSize(B,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do B.Data[i,j]:=A.Data[i,j]*Value;
    Result:=B;
  End;

  Procedure DeleteRow(var A:Matrix; Row:integer); { �������� ������ }
  var
    i    :  integer;
  Begin
    if (Row>=A.M) then Raise Exception.Create('������� ������� �� ������� �������������� ������!');
    SetLength(A.Data[Row],0);
    for i:=Row+1 to A.M-1 do A.Data[i-1]:=A.Data[i];
    SetLength(A.Data,A.M-1);
    dec(A.M);
  End;

  Procedure DeleteCol(var A:Matrix; Col:integer); { �������� ������� }
  var
    i,j   : integer;
  Begin
    if (Col>=A.N) then Raise Exception.Create('������� ������� �� ������� �������������� �������!');
    for i:=0 to A.M-1 do begin
      for j:=Col+1 to A.N-1 do A.Data[i,j-1]:=A.Data[i,j];
      SetLength(A.Data[i],A.N-1);
    end;
    dec(A.N);
  End;

  Procedure DeleteCross(var A:Matrix; Row:integer); { �������� ������ � ������� � �������� ������� }
  Begin
    DeleteRow(A,Row);
    DeleteCol(A,Row);
  End;

  Procedure InsertRow(var A:Matrix; Row:integer);
  var
    i  : integer;
  Begin
    SetSize(A,A.M+1,A.N);
    for i:=A.M-1 downto Row do A.Data[i+1]:=A.Data[i];
    SetLength(A.Data[Row],A.N);
    for i:=0 to A.N-1 do A.Data[Row,i]:=0;
  End;

  Procedure InsertCol(var A:Matrix; Col:integer);
  var
    i,j   : integer;
  Begin
    SetSize(A,A.M,A.N+1);
    for i:=0 to A.M-1 do begin
      for j:=A.N-1 downto Col do A.Data[i,j+1]:=A.Data[i,j];
      A.Data[i,Col]:=0;
    end;
  End;

  Function MapMatrix(const A:Matrix; const Map:TMap):Matrix; { �������� �������� ����� � �������� }
  var
    i,j        : integer;
    B          : Matrix;
    CurValue   : integer;
  Begin
    SetSize(B,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do B.Data[i,j]:=A.Data[i,j];
    for i:=0 to length(Map)-2 do begin
     for j:=i to length(Map)-1 do if (Map[j]<Map[i]) then begin
       CurValue:=Map[j];
       Map[j]:=Map[i];
       Map[i]:=CurValue;
     end;
    end;
    for i:=length(Map)-1 downto 0 do DeleteRow(B,Map[i]);
    for i:=length(Map)-1 downto 0 do DeleteCol(B,Map[i]);
    Result:=B;
  End;

  Function MapVector(const A:Matrix; const Map:TMap):Matrix; { �������� �������� �����}
  var
    i,j        : integer;
    B          : Matrix;
    CurValue   : integer;
  Begin
    SetSize(B,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do B.Data[i,j]:=A.Data[i,j];
    for i:=0 to length(Map)-2 do begin
     for j:=i to length(Map)-1 do if (Map[j]<Map[i]) then begin
       CurValue:=Map[j];
       Map[j]:=Map[i];
       Map[i]:=CurValue;
     end;
    end;
    for i:=length(Map)-1 downto 0 do DeleteRow(B,Map[i]);
    Result:=B;
  End;

  Function UnmapMatrix(const A:Matrix; const Map:TMap):Matrix; { �������������� ������� �� ��������� �������}
  var
    i,j        : integer;
    B          : Matrix;
    CurValue   : integer;
  Begin
    SetSize(B,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do B.Data[i,j]:=A.Data[i,j];
    for i:=0 to length(Map)-2 do begin
     for j:=i to length(Map)-1 do if (Map[j]<Map[i]) then begin
       CurValue:=Map[j];
       Map[j]:=Map[i];
       Map[i]:=CurValue;
     end;
    end;
    for i:=0 to length(Map)-1 do InsertRow(B,Map[i]);
    for i:=0 to length(Map)-1 do InsertCol(B,Map[i]);
    Result:=B;
  End;

  Function UnmapVector(const A:Matrix; const Map:TMap):Matrix; { �������������� ���������� ����� ������� �� ��������� �������}
  var
    i,j        : integer;
    B          : Matrix;
    CurValue   : integer;
  Begin
    SetSize(B,A.M,A.N);
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do B.Data[i,j]:=A.Data[i,j];
    for i:=0 to length(Map)-2 do begin
     for j:=i to length(Map)-1 do if (Map[j]<Map[i]) then begin
       CurValue:=Map[j];
       Map[j]:=Map[i];
       Map[i]:=CurValue;
     end;
    end;
    for i:=0 to length(Map)-1 do InsertRow(B,Map[i]);
    Result:=B;
  End;


  Function SubMatrix(const A:Matrix; const MapX,MapY:TMap):Matrix; { ��������� ����������, ���������� �������� ������� � ������ }
  var
    i,j       : integer;
    B         : Matrix;
  Begin
    SetSize(B,Length(MapY),Length(MapX));
    for i:=0 to length(MapY)-1 do begin
      for j:=0 to length(MapX)-1 do B.Data[i,j]:=A.Data[MapY[i],MapX[j]];
    end;
    Result:=B;
  End;

  Function RangeMatrix(const A:Matrix; StartX,EndX,StartY,EndY:integer):Matrix; { ��������� ����������, ���������� �������� �������� �������� � ����� }
  var
    i,j   : integer;
    B     : Matrix;
  Begin
    SetSize(B,EndY-StartY+1,EndX-StartX+1);
    for i:=StartY to EndY do begin
      for j:=StartX to EndX do B.Data[i-StartY,j-StartX]:=A.Data[i,j];
    end;
    Result:=B;
  End;

  Function Vectorize(const A:Matrix):Vector; { �������������� ������� � ������, ������������� � ������ ��������������� }
  var
    i,j,k  : integer;
    V     : Vector;
  Begin
    SetLength(V,A.M*A.N+2);
    V[0]:=A.M;
    V[1]:=A.N;
    k:=2;
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do begin
      V[k]:=A.Data[i,j];
      inc(k);
    end;
    Result:=V;
  End;

  Function Unvectorize(const V:Vector):Matrix; { �������� ��������������: ������ � ������� }
  var
    i,j,k : integer;
    A     : Matrix;
  Begin
    SetSize(A,trunc(V[0]),trunc(V[1]));
    k:=2;
    for i:=0 to A.M-1 do for j:=0 to A.N-1 do begin
      A.Data[i,j]:=V[k];
      inc(k);
    end;
    Result:=A;
  End;

end.
