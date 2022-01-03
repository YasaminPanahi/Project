clear;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preparation
bc = 0; % boundary condition
J = 1;
h = 2;
% Number of states kept
n_kept_states = 40;
% Number of iterations.  Final lattice size is 2*n_iter + 2
n_iter = 50;
% exact energy for the Heisenberg model (from Bethe Ansatz), for comparison
%exact_energy = -log(2) + 0.25;

fprintf('Iter\tEnergy\t\tBondEnergy\tEnergyError\tTrunc\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Intialize local operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I= eye(2);
sz = [1/2 0;0 -1/2]; %Sz
sx = [0 1/2;1/2 0]; %Sx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Initial blocks
%               We assume reflection symmetry so we only need 1 block
%               The operator acts on the inner-most site of the block
%               +---------+    +---------+
%               |        *|    |*        |
%               +---------+    +---------+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


block_sz_last = sz;
block_sx_last = sx;
block_Id  = I;

H_2 = -J*kron(sz,sz) - h*(kron(I,sx) + kron(sx,I));
energy= eigs(H_2,1,'sa');   % energy of 2-site lattice

   
%%%%%%%%%%%%%%%%%%%%
% initializing a few useful quantities
energy_per_bond = zeros(1,n_iter);
truncation_error = zeros(1,n_iter);
system_size = zeros(1,n_iter);
EE = zeros(1,n_iter); % enatnglement entropy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to be used for measureing correlation function
max_length_left_block = n_iter + 1;
x_0 = round(max_length_left_block/2);
Z = cell(1,max_length_left_block);  % sz(i);
ZZ = cell(1,max_length_left_block); % sz(i)sz(x_0) % here we consider x_0 = 12 (26!!)
X = cell(1,max_length_left_block);%@@   
XX = cell(1,max_length_left_block);%@@  
for x = 1 : max_length_left_block
    if x == 1
        Z{x} = sz;
        ZZ{x} = sz;
        X{x} = sx;
        XX{x} = sx;

    else
        Z{x} = I;
        ZZ{x} = I;
        X{x} = I;
        XX{x} = I;
    end
end
%% main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Begin main iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
block_H = -h*sx;
for l = 1:n_iter
    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Get the 2*chi-dimensional operators for the block +
    %           site where chi = n_kept_states
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    block_H = kron(block_H, I) -J*kron(block_sz_last, sz) ...
        -h*( kron(block_sx_last, I));
    
    
    block_sz_last = kron(block_Id, sz);
    block_sx_last = kron(block_Id, sx);
    block_Id  = kron(block_Id, I);
    
    left_block_length = l + 1;
    for x = 1 : max_length_left_block
        if left_block_length == x
            Z{x} = kron(Z{x},sz); %@@: kron(I,sz)
            X{x} = kron(X{x},sx); %@@: kron(I,sx)
        else
            Z{x} = kron(Z{x},I);  %@@:
            X{x} = kron(X{x},I);  %@@:
        end
        
        if left_block_length == x || left_block_length == x_0
            ZZ{x} = kron(ZZ{x},sz);
            XX{x} = kron(XX{x},sx);
        else
            ZZ{x} = kron(ZZ{x},I);
            XX{x} = kron(XX{x},I);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   HAMILTONIAN MATRIX forsuperblock
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    H_superblok = kron(block_H, block_Id) + kron(block_Id, block_H) ...
        - J*kron(block_sz_last, block_sz_last) - J*bc*(kron(Z{1},Z{1}));
    
    H_superblok = 0.5 * (H_superblok + H_superblok');  % ensure H is symmetric
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Diagonalizing the Hamiltonian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    last_energy = energy;
    [Psi,energy] = eigs(H_superblok,1,'sa');
    energy_per_bond(1,l) = (energy - last_energy)/2;
    system_size(1,l) = 2*(1+l);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Form the reduced density matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dim_block = size(block_H,1);
    Psi_matrix = reshape(Psi, dim_block, dim_block);
    block_Rho = Psi_matrix * Psi_matrix';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Diagonalize the reduced density matrix to
    %                   obtain the entangelement entropy an truncation
    %                   operator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_kept_temp = min(dim_block, n_kept_states);
    [V,D] = eigs(block_Rho,dim_block,'la'); % la --> largest eigenvalues
    ES = diag(D); % entanglement spectrum
    ES(ES<10^(-15)) = 10^(-15); % making sure all entanglement spectra are positive
    EE(1,l) = - dot(ES,log(ES)); % entanglement entropy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Construct the truncation operator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    block_T = V(:, 1:n_kept_temp);
    truncation_error(1,l) = 1 - sum(ES(1:n_kept_temp));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Transform the block operators into the truncated basis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if l < n_iter
        block_H  = block_T'*block_H*block_T;
        block_sz_last = block_T'*block_sz_last*block_T;
        block_sx_last = block_T'*block_sx_last*block_T;
        block_Id = block_T'*block_Id*block_T;
        for x = 1 : max_length_left_block
            Z{x} = block_T'*Z{x}*block_T;
            ZZ{x} = block_T'*ZZ{x}*block_T;
            X{x} = block_T'*X{x}*block_T;
            XX{x} = block_T'*XX{x}*block_T;
        end
    end
    time_per_step(1,l) = toc;
    current_step_computation_time = [l , time_per_step(1,l)]
    %current_step_info = [l, energy_per_bond(c,l), exact_energy]
    
end
Z_avg = zeros(1,max_length_left_block);
ZZ_avg = zeros(1,max_length_left_block);
X_avg = zeros(1,max_length_left_block);
XX_avg = zeros(1,max_length_left_block);
ZZ_avg_c = zeros(1,max_length_left_block);
XX_avg_c = zeros(1,max_length_left_block);
for x = 1 : max_length_left_block
    Z_avg(1,x) = trace(block_Rho*Z{x});
    ZZ_avg(1,x) = trace(block_Rho*ZZ{x});
    X_avg(1,x) = trace(block_Rho*X{x});
    XX_avg(1,x) = trace(block_Rho*XX{x});
end
ZZ_avg_c = ZZ_avg - Z_avg*Z_avg(1,x_0); % disconnected correlation function <AB>_c = <AB> - <A><B>
XX_avg_c = XX_avg - X_avg*X_avg(1,x_0);

%% plots
figure(1); %<Sz>
plot(Z_avg,'-s','Linewidth',2);
legend('h=0.5')
%hold on
% plot(Z_avg(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(Z_avg(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
% plot(Z_avg(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(Z_avg(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
% plot(Z_avg(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(Z_avg(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
% plot(Z_avg(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('position');
ylabel('<S_z(i)>');

figure(2); %<Sx>
plot(X_avg,'-s','Linewidth',2);
legend('h=0.5')
xlabel('position');
ylabel('<S_x(i)>');

figure(3); % Correlation z
plot(ZZ_avg_c,'-s','Linewidth',2);
legend('h=0.5')
hold on
% plot(ZZ_avg_c(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(ZZ_avg_c(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
% plot(ZZ_avg_c(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(ZZ_avg_c(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
% plot(ZZ_avg_c(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(ZZ_avg_c(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
% plot(ZZ_avg_c(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('position');
ylabel('<S_z(x)S_z(x_0)>-<S_z(x)><S_z(x_0)>');

figure(4); % Correlation x
plot(XX_avg_c,'-s','Linewidth',2);
legend('h=0.5')
hold on
% plot(XX_avg_c(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(XX_avg_c(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
% plot(XX_avg_c(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(XX_avg_c(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
% plot(XX_avg_c(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(XX_avg_c(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
% plot(XX_avg_c(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('position');
ylabel('<S_x(x)S_x(x_0)>-<S_x(x)><S_x(x_0)>');


figure(5); % truncation error
plot(system_size,truncation_error,'-s','Linewidth',2);
legend('h=0.5')
hold on
% plot(system_size(2,:),truncation_error(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(system_size(3,:),truncation_error(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
% plot(system_size(4,:),truncation_error(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(system_size(5,:),truncation_error(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
% plot(system_size(6,:),truncation_error(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(system_size(7,:),truncation_error(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
% plot(system_size(8,:),truncation_error(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('system size');
ylabel('truncation error');


figure(6); % Entangelement Entropy
plot(log(system_size/2),EE,'-s','Linewidth',2);
legend('h=0.5')
hold on
% plot(log(system_size(2,:)/2),EE(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(log(system_size(3,:)/2),EE(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
% plot(log(system_size(4,:)/2),EE(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(log(system_size(5,:)/2),EE(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
% plot(log(system_size(6,:)/2),EE(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(log(system_size(7,:)/2),EE(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
% plot(log(system_size(8,:)/2),EE(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('log(left_subsystem size)');
ylabel('Entangelement Entropy');
x = log(system_size/2);
y = EE;
p = polyfit(x,y,1);
central_charge = 3*(2-bc)*p(1)

figure(7); % energy per bond accuracy
plot(system_size,energy_per_bond,'-s','Linewidth',2);
%hold on;
%plot(system_size,exact_energy*ones(size(system_size)),'-.','Linewidth',2);
xlabel('system size');
ylabel('energy per bond');

total_computation_time = sum(time_per_step)   