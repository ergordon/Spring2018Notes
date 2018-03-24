%Simple Square
%{
n=3;
A = zeros(n,n);
%Row = j
%Column = i
for i = 1:n
    for j = 1:n
        %Start at top left, move down, jump back to top one column to the right, repeat
        %A(i,j) = i+(j-1)*n;
        
        %Start at top left, move right, jump back to left one row down,repeat
        A(i,j) = (i-1)*n+j;
    end
end
A
A(2,3) %second row, third column
%}

%Irregular Square
%%{
Nx=5; %x subdivisions
Ny=10; %y subdivisions
A = zeros(2*Ny,2*Nx);
%Row = i
%Column = j

for i = 1:2*Ny
    for j = 1:2*Nx
        %Start at top left, move down, jump back to top one column to the right, repeat
        %A(i,j) = i+(j-1)*ny;
        
        %Start at top left, move right, jump back to left one row down,repeat
        %A(i,j) = (i-1)*nx+j;
        if(i <= Nx +1)
            A(i,j) = (i-1)*(2*Ny+1)+j;
        elseif(j >= Ny + 1)
            A(i,j) = (Nx + 1)*(2*Ny+1)+(i-Nx-2)*(Ny+1)+j-Ny;
        end
      
    end
end

qA=1
qB=2*Ny+1        % N
qC=Nx*(2*Ny+1)+1 % Nx*(qB)+1
qD=qC+Ny
qE=qD+Ny
qF=qE+(Nx-1)*(Ny+1)+1
qG=qF+Ny

A
%}
