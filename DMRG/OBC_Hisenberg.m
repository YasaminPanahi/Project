clear;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bc = 1; % boundary condition

% Number of states kept
n_kept_states = 40;
% Number of iterations.  Final lattice size is 2*n_iter + 2
n_iter = 50;
% exact energy for the Heisenberg model (from Bethe Ansatz), for comparison
exact_energy = -log(2) + 0.25;

fprintf('Iter\tEnergy\t\tBondEnergy\tEnergyError\tTrunc\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Intialize local operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I= eye(2);
sz = [1/2 0;0 -1/2]; %Sz
sp = [0 0;1 0]; %S+
sm = [0 1 ;0 0];%S-
sx = [0 1/2;1/2 0];
sy = 1i*[0 -1/2; 1/2 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Initial blocks
%               We assume reflection symmetry so we only need 1 block
%               The operator acts on the inner-most site of the block
%               +---------+    +---------+
%               |        *|    |*        |
%               +---------+    +---------+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


block_sz_last = sz;
block_sp_last = sp;
block_sm_last = sm;
block_Id  = I;
block_H  = zeros(2);
H_2 = kron(sz,sz)+0.5*(kron(sp,sm)+kron(sm,sp)); % for 2-sites H = Sz(1)Sz(2) + Sx(1)Sx(2)+Sy(1)Sy(2) = Sz(1)Sz(1)+(Sp(1)Sm(2)+Sm(1)Sp(2))/2;
energy = eigs(H_2,1,'sa');   % energy of 2-site lattice

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
Y = cell(1,max_length_left_block);%@@   
%YY = cell(11,max_length_left_block);%@@  
for x = 1 : max_length_left_block
    if x == 1
        Z{x} = sz;
        ZZ{x} = sz;
        X{x} = sx;
        XX{x} = sx;
        Y{x} = sy;
        %YY{x} = sy;

    else
        Z{x} = I;
        ZZ{x} = I;
        X{x} = I;
        XX{x} = I;
        Y{x} = I;
        %YY{x} = I;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Begin main iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l = 1:n_iter
    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Get the 2*chi-dimensional operators for the block +
    %           site where chi = n_kept_states
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block_H = kron(block_H, I) + kron(block_sz_last, sz) ...
        + 0.5 * ( kron(block_sp_last, sm) + kron(block_sm_last, sp) );
    block_sz_last = kron(block_Id, sz);
    block_sp_last = kron(block_Id, sp);
    block_sm_last = kron(block_Id, sm);
    block_Id  = kron(block_Id, I);

    left_block_length = l + 1;
    for x = 1 : max_length_left_block
        if left_block_length == x
            Z{x} = kron(Z{x},sz); %@@: kron(I,sz)
            X{x} = kron(X{x},sx); %@@: kron(I,sx)
            Y{x} = kron(Y{x},sy); %@@: kron(I,sy)   
        else
            Z{x} = kron(Z{x},I);  %@@: 
            X{x} = kron(X{x},I);  %@@: 
            Y{x} = kron(Y{x},I);  %@@: 
        end

        if left_block_length == x || left_block_length == x_0
            ZZ{x} = kron(ZZ{x},sz);
            XX{x} = kron(XX{x},sx);
            %YY{x} = kron(YY{x},sy);
        else
            ZZ{x} = kron(ZZ{x},I);
            XX{x} = kron(XX{x},I);
            %YY{x} = kron(YY{x},I);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   HAMILTONIAN MATRIX forsuperblock
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    H_superblok = kron(block_H, block_Id) + kron(block_Id, block_H) ...
        + kron(block_sz_last, block_sz_last) ...
        + 0.5 * ( kron(block_sp_last, block_sm_last) + kron(block_sm_last, block_sp_last) )...
        + 0.5 * ( kron(Z{1},Z{1}) + kron(X{1},X{1})+ kron(Y{1},Y{1}) );

    H_superblok = 0.5 * (H_superblok + H_superblok');  % ensure H is symmetric

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Diagonalizing the Hamiltonian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    last_energy = energy;
    [Psi,energy] = eigs(H_superblok,1,'sa');
    energy_per_bond(1,l) = (energy - last_energy) / 2;
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
        block_sp_last = block_T'*block_sp_last*block_T;
        block_sm_last = block_T'*block_sm_last*block_T;
        block_Id = block_T'*block_Id*block_T;
        for x = 1 : max_length_left_block
            Z{x} = block_T'*Z{x}*block_T;
            ZZ{x} = block_T'*ZZ{x}*block_T;
            X{x} = block_T'*X{x}*block_T;
            XX{x} = block_T'*XX{x}*block_T;
            Y{x} = block_T'*Y{x}*block_T;
            %YY{x} = block_T'*YY{x}*block_T;
        end
    end
    time_per_step(1,l) = toc;
    current_step_computation_time = [l , time_per_step(1,l)]
    current_step_info = [l, energy_per_bond(1,l), exact_energy]
   
end
   
Z_avg = zeros(1,max_length_left_block);
ZZ_avg = zeros(1,max_length_left_block);
X_avg = zeros(1,max_length_left_block);
XX_avg = zeros(1,max_length_left_block);
for x = 1 : max_length_left_block
    Z_avg(1,x) = trace(block_Rho*Z{x});
    ZZ_avg(1,x) = trace(block_Rho*ZZ{x});
    X_avg(1,x) = trace(block_Rho*X{x});
    XX_avg(1,x) = trace(block_Rho*XX{x});
end
ZZ_avg_c = ZZ_avg - Z_avg*Z_avg(1,x_0); % disconnected correlation function <AB>_c = <AB> - <A><B>
XX_avg_c = XX_avg - X_avg*X_avg(1,x_0);


figure(1); %<Sz>
plot(Z_avg,'-d','Linewidth',2);
xlabel('position');
ylabel('<sz(i)>');

figure(2); % Correlation z
plot(ZZ_avg_c,'-d','Linewidth',2);
xlabel('position');
ylabel('<Sz(x)Sz(x_0)>-<Sz(x)><Sz(x_0)>');

figure(3); % Correlation x
plot(XX_avg_c,'-d','Linewidth',2);
xlabel('position');
ylabel('<Sx(x)Sx(x_0)>-<Sx(x)><Sx(x_0)>');


figure(4); % truncation error
plot(system_size,truncation_error,'-d','Linewidth',2);
xlabel('system size');
ylabel('truncation error');


figure(5); % Entangelement Entropy
plot(log(system_size/2),EE,'-d','Linewidth',2);
xlabel('log(left_subsystem size)');
ylabel('Entangelement Entropy');
x = log(system_size/2);
y = EE;
p = polyfit(x,y,1);
central_charge = 3*(2-bc)*p(1)

figure(6); % energy per bond accuracy
plot(system_size,energy_per_bond,'-d','Linewidth',2);
hold on;
plot(system_size,exact_energy*ones(size(system_size)),'-.','Linewidth',2);
xlabel('system size');
ylabel('energy per bond');



total_computation_time = sum(time_per_step)   