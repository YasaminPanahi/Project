clear;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bc = 0; % boundary condition
J = 1;
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
%block_H  = zeros(2);
energy = zeros(1,8);
last_energy = zeros(1,8);
c=1;
for h=0.25:0.25:2
    H_2 = -J*kron(sz,sz) - h*(kron(I,sx) + kron(sx,I)); 
    energy(1,c) = eigs(H_2,1,'sa');   % energy of 2-site lattice
    c=c+1;
end
   

%%%%%%%%%%%%%%%%%%%%
% initializing a few useful quantities
energy_per_bond = zeros(8,n_iter);
truncation_error = zeros(8,n_iter);
system_size = zeros(8,n_iter);
EE = zeros(8,n_iter); % enatnglement entropy


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Begin main iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=1;
for h=0.25:0.25:2
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
            - J*kron(block_sz_last, block_sz_last) - J*bc*( kron(Z{1},Z{1}));
        
        H_superblok = 0.5 * (H_superblok + H_superblok');  % ensure H is symmetric
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   Diagonalizing the Hamiltonian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        last_energy(1,c) = energy(1,c);
        [Psi,Energy] = eigs(H_superblok,1,'sa');
        energy(1,c) = Energy;
        energy_per_bond(c,l) = (energy(1,c) - last_energy(1,c)) / 2;
        system_size(c,l) = 2*(1+l);
        
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
        EE(c,l) = - dot(ES,log(ES)); % entanglement entropy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   Construct the truncation operator
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        block_T = V(:, 1:n_kept_temp);
        truncation_error(c,l) = 1 - sum(ES(1:n_kept_temp));
        
        
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
    Z_avg = zeros(8,max_length_left_block);
    ZZ_avg = zeros(8,max_length_left_block);
    X_avg = zeros(8,max_length_left_block);
    XX_avg = zeros(8,max_length_left_block);
    ZZ_avg_c = zeros(8,max_length_left_block);
    XX_avg_c = zeros(8,max_length_left_block);
    for x = 1 : max_length_left_block
        Z_avg(8,x) = trace(block_Rho*Z{x});
        ZZ_avg(8,x) = trace(block_Rho*ZZ{x});
        X_avg(8,x) = trace(block_Rho*X{x});
        XX_avg(8,x) = trace(block_Rho*XX{x});
    end
    for x = 1 : max_length_left_block
        ZZ_avg_c(c,x) = ZZ_avg(c,x) - Z_avg(c,x)*Z_avg(c,x_0);
        XX_avg_c(c,x) = XX_avg(c,x) - X_avg(c,x)*X_avg(c,x_0);
    end
    c=c+1;
        
end
   
figure(1); %<Sz>
plot(Z_avg(1,:),'-d','Linewidth',2);
legend('h=0.25')
hold on
plot(Z_avg(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(Z_avg(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
plot(Z_avg(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(Z_avg(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
plot(Z_avg(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(Z_avg(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
plot(Z_avg(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('position');
ylabel('<sz(i)>');

figure(2); % Correlation z
plot(ZZ_avg_c,'-d','Linewidth',2);
legend('h=0.25')
hold on
plot(ZZ_avg_c(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(ZZ_avg_c(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
plot(ZZ_avg_c(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(ZZ_avg_c(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
plot(ZZ_avg_c(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(ZZ_avg_c(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
plot(ZZ_avg_c(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('position');
ylabel('<Sz(x)Sz(x_0)>-<Sz(x)><Sz(x_0)>');

figure(3); % Correlation x
plot(XX_avg_c,'-d','Linewidth',2);
legend('h=0.25')
hold on
plot(XX_avg_c(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(XX_avg_c(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
plot(XX_avg_c(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(XX_avg_c(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
plot(XX_avg_c(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(XX_avg_c(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
plot(XX_avg_c(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('position');
ylabel('<Sx(x)Sx(x_0)>-<Sx(x)><Sx(x_0)>');


figure(4); % truncation error
plot(system_size(1,:),truncation_error(1,:),'-d','Linewidth',2);
legend('h=0.25')
hold on
plot(system_size(2,:),truncation_error(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(system_size(3,:),truncation_error(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
plot(system_size(4,:),truncation_error(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(system_size(5,:),truncation_error(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
plot(system_size(6,:),truncation_error(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(system_size(7,:),truncation_error(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
plot(system_size(8,:),truncation_error(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('system size');
ylabel('truncation error');


figure(5); % Entangelement Entropy
plot(log(system_size(1,:)/2),EE(1,:),'-d','Linewidth',2);
legend('h=0.25')
hold on
plot(log(system_size(2,:)/2),EE(2,:),'-d','Linewidth',2,'DisplayName','h=0.5');plot(log(system_size(3,:)/2),EE(3,:),'-d','Linewidth',2,'DisplayName','h=0.75');
plot(log(system_size(4,:)/2),EE(4,:),'-d','Linewidth',2,'DisplayName','h=1');plot(log(system_size(5,:)/2),EE(5,:),'-d','Linewidth',2,'DisplayName','h=1.25');
plot(log(system_size(6,:)/2),EE(6,:),'-d','Linewidth',2,'DisplayName','h=1.5');plot(log(system_size(7,:)/2),EE(7,:),'-d','Linewidth',2,'DisplayName','h=1.75');
plot(log(system_size(8,:)/2),EE(8,:),'-d','Linewidth',2,'DisplayName','h=2');hold all
xlabel('log(left_subsystem size)');
ylabel('Entangelement Entropy');
central_charge = zeros(1,8);
for i=1:8
    x = log(system_size(i,:)/2);
    y = EE(i,:);
    p = polyfit(x,y,1);
    central_charge(i) = 3*(2-bc)*p(1)
end

% figure(6); % energy per bond accuracy
% plot(system_size(1,,energy_per_bond,'-d','Linewidth',2);
% hold on;
% plot(system_size,exact_energy*ones(size(system_size)),'-.','Linewidth',2);
% xlabel('system size');
% ylabel('energy per bond');


total_computation_time = sum(time_per_step)   