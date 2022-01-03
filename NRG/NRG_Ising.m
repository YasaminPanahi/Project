clc 
clear
tic
%=========================
%       Problem 3       
%=========================
%Parameters:
J = 1;
h = 1/4;
bnd = 0; %periodic bc = 1 & %open bc = 0

%pauli matrices:
sigma_x = 1/2*[0 1;1 0];sigma_y = 1/2*[0 -1i;1i 0];sigma_z = 1/2*[1 0;0 -1];I = eye(2);
X = sparse(sigma_x);Y = sparse(sigma_y);Z = sparse(sigma_z);I = sparse(I);
G_Energy = zeros(10,11); %Exact ground state energy 

%% part 1)
for N=10:20
    Sx = cell(1,N);Sz = cell(1,N);
    I_others = 1;
    for i=1:N
        Sx{i} = kron(I_others,X);
        Sz{i} = kron(I_others,Z);
        I_others = kron(I_others,I);
    
        for j=1:i-1
            Sx{j} = kron(Sx{j},I);
            Sz{j} = kron(Sz{j},I);
        end
    
        for j=i+1:N
            Sx{j} = I_others;
            Sz{j} = I_others;
        end
    end
    c=1;
    for h=0.1:0.1:1
        H = -J*bnd*(Sz{N}*Sz{1}) -h*(Sx{N}); %boundary conditions
        for i=1:N-1
            H = H -J*(Sz{i}*Sz{i+1}) -h*(Sx{i});
        end
        H = (H+H')/2;
        [States,Energys] = eigs(H,1,'sa');
        energy = diag(Energys);
        G_Energy(c,N-9) = energy(1);
        c=c+1;
    end
end

%% part 2)
N = 10;
Sx = cell(1,N);Sz = cell(1,N);
Sx_p = cell(1,N);Sz_p = cell(1,N);
I_others = 1;
G_NRG = zeros(1,10); %Ground state energy 
% First step N = 10
for i=1:N %Construction of spin operators
     Sx{i} = kron(I_others,X);
     Sz{i} = kron(I_others,Z);
     I_others = kron(I_others,I);
    
     for j=1:i-1
         Sx{j} = kron(Sx{j},I);
         Sz{j} = kron(Sz{j},I);
     end
    
     for j=i+1:N
         Sx{j} = I_others;
         Sz{j} = I_others;
     end
end
%Hamiltonian
H_10 = -J*bnd*(Sz{N}*Sz{1}) -h*(Sx{N}); %boundary conditions
for i=1:N-1
H_10 = H_10 -J*(Sz{i}*Sz{i+1}) -h*(Sx{i});
end  
H_10 = (H_10+H_10')/2; %making sure that it's hermitian
%Diagonalization
[States,Energys] = eigs(H_10,2^10,'sa');
P_new = States(:,1:2^9);%2^10 *2^9
H_Prime10 = P_new'*H_10*P_new; %2^9 *2^9
%Projection of operators
for i=1:N
    Sz_p{i} = P_new'*Sz{i}*P_new;%2^9 *2^9
    Sx_p{i} = P_new'*Sx{i}*P_new;%2^9 *2^9
end

%Bigger systems (adding more atoms)
for n=1:10 %n=1
    Sx_new = cell(1,N+n);Sz_new = cell(1,N+n);
    for i=1:N+n-1
        Sz_new{i} = kron(Sz_p{i},I);
        Sx_new{i} = kron(Sx_p{i},I);
    end
    I_p =  P_new'*eye(2^N)*P_new;
    Sz_new{N+n}=kron(I_p,Z);
    Sx_new{N+n}=kron(I_p,X);
    
    %Hamiltonian
    H_new = -J*bnd*(Sz_new{N+n}*Sz_new{1}) -h*(Sx_new{N+n}); %boundary conditions
    for i=1:N+n-1
        H_new = H_new -J*(Sz_new{i}*Sz_new{i+1}) -h*(Sx_new{i});
    end
    H_new = (H_new+H_new')/2;
    
    %Diagonalization
    [States,Energys] = eigs(H_new,2^10,'sa');
    energy = diag(Energys); %Energy spectrum
    P_new = States(:,1:2^9);%2^10 *2^9
    G_NRG(:,n) = energy(1);%Ground state energy
    G_State = States(:,1); %Ground state wavefunction
    
    %New projection operators
    Sx_p = cell(1,N+n);Sz_p = cell(1,N+n);
    for i=1:N+n
        Sz_p{i} = P_new'*Sz_new{i}*P_new;%2^9 *2^9
        Sx_p{i} = P_new'*Sx_new{i}*P_new;%2^9 *2^9
    end
end

Error = zeros(1,10);
Abs_Error = zeros(1,10);
exp_x = zeros(1,20);
exp_z = zeros(1,20);
exp_xx = zeros(1,20);
for i=1:10
    Error(i) = abs((G_NRG(i) - G_Energy(i+1))/G_Energy(i+1));
    Abs_Error(i) = abs((G_NRG(i) - G_Energy(i+1)));
end
%% part 4)
for i=1:20
    exp_z(i) = (G_State')*Sz_new{i}*(G_State);
    exp_x(i) = (G_State')*Sx_new{i}*(G_State);
    exp_xx(i) = (G_State')*(Sx_new{10}*Sx_new{i})*(G_State);
end
    

plots:
figure(1)
x_plot=linspace(11,20,10);
plot(x_plot,Error,'-s','Linewidth',2);
t=title('NRG Error');
grid();

figure(2)
x_plot=linspace(11,20,10);
plot(x_plot,Abs_Error,'-s','Linewidth',2);
t=title('NRG Absolute Error');
grid();

x_plot=linspace(1,20,20);
figure(3)
plot(x_plot,exp_x,'-s','Linewidth',2);
xlim([1 20]);
t=title('<S_x>');
grid();

figure(4)
plot(x_plot,exp_z,'-s','Linewidth',2);
xlim([1 20]);
t=title('<S_z>');
grid();

figure(5)
plot(x_plot,exp_xx,'-s','Linewidth',2);
xlim([1 20]);
t=title('<S_x(10)S_x(i)>');
grid();
%% part 3)
N = 10;
Sx = cell(1,N);Sz = cell(1,N);
Sx_p = cell(1,N);Sz_p = cell(1,N);
I_others = 1;
G_NRG = zeros(10,10); %Ground state energy
% First step N = 10
for i=1:N %Construction of spin operators
    Sx{i} = kron(I_others,X);
    Sz{i} = kron(I_others,Z);
    I_others = kron(I_others,I);
    
    for j=1:i-1
        Sx{j} = kron(Sx{j},I);
        Sz{j} = kron(Sz{j},I);
    end
    
    for j=i+1:N
        Sx{j} = I_others;
        Sz{j} = I_others;
    end
end
count = 1;
for h=0.1:0.1:1
    %Hamiltonian
    H_10 = -J*bnd*(Sz{N}*Sz{1}) -h*(Sx{N}); %boundary conditions
    for i=1:N-1
        H_10 = H_10 -J*(Sz{i}*Sz{i+1}) -h*(Sx{i});
    end
    H_10 = (H_10+H_10')/2; %making sure that it's hermitian
    %Diagonalization
    [States,Energys] = eigs(H_10,2^10,'sa');
    P_new = States(:,1:2^9);%2^10 *2^9
    %H_Prime10 = P_new'*H_10*P_new; %2^9 *2^9
    %Projection of operators
    for i=1:N
        Sz_p{i} = P_new'*Sz{i}*P_new;%2^9 *2^9
        Sx_p{i} = P_new'*Sx{i}*P_new;%2^9 *2^9
    end
    %Bigger systems (adding more atoms)
    for n=1:10 %n=1
        Sx_new = cell(1,N+n);Sz_new = cell(1,N+n);
        for i=1:N+n-1
            Sz_new{i} = kron(Sz_p{i},I);
            Sx_new{i} = kron(Sx_p{i},I);
        end
        I_p =  P_new'*eye(2^N)*P_new;
        Sz_new{N+n}=kron(I_p,Z);
        Sx_new{N+n}=kron(I_p,X);
        
        %Hamiltonian
        H_new = -J*bnd*(Sz_new{N+n}*Sz_new{1}) -h*(Sx_new{N+n}); %boundary conditions
        for i=1:N+n-1
            H_new = H_new -J*(Sz_new{i}*Sz_new{i+1}) -h*(Sx_new{i});
        end
        H_new = (H_new+H_new')/2;
        %Diagonalization
        [States,Energys] = eigs(H_new,2^10,'sa');
        energy = diag(Energys); %Energy spectrum
        P_new = States(:,1:2^9);%2^10 *2^9
        G_NRG(:,n) = energy(1);%Ground state energy
        G_State = States(:,1); %Ground state wavefunction
        
        %New projection operators
        Sx_p = cell(1,N+n);Sz_p = cell(1,N+n);
        for i=1:N+n
            Sz_p{i} = P_new'*Sz_new{i}*P_new;%2^9 *2^9
            Sx_p{i} = P_new'*Sx_new{i}*P_new;%2^9 *2^9
        end
    end
    
    Error = zeros(1,10);
    Abs_Error = zeros(1,10);
    for i=1:10
        Error(count,i) = abs((G_NRG(count,i) - G_Energy(count,i+1))/G_Energy(count,i+1));
        Abs_Error(count,i) = abs((G_NRG(count,i) - G_Energy(count,i+1)));
    end
    count = count+1;
end
figure(7)
x_plot=linspace(0.1,1,10);
plot(x_plot,Error(:,10),'-s','Linewidth',2);
t=title('NRG Error');
grid();
    
toc