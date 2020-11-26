function [V,FD_grid,time_steps] = bs_vega_timestepping(sigma,r,...
        forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta)



S_min=0;
% set grid and grid-size.
h = (S_max-S_min)/Nx;
FD_grid = linspace(S_min,S_max,Nx+1); 
inner_grid = FD_grid(1:end-1);

% set number of time steps
deltat=T/Nt;
time_steps = [0 deltat*(1:Nt)];

% init matrix containing the solution at each time step
V=zeros(length(FD_grid),Nt+1);
V(:,1)=initial_cond(FD_grid)';


S = inner_grid';
S_up = inner_grid(1:end-1)'; 
S_down = inner_grid(2:end)'; 

Ssq = inner_grid'.^2;
Ssq_up = inner_grid(1:end-1)'.^2;
Ssq_down = inner_grid(2:end)'.^2;

aa_old_up   = -0.5/h^2*sigma^2.*Ssq_up; 
aa_old_main =  1/h^2*sigma^2.*Ssq;
aa_old_down = -0.5/h^2*sigma^2.*Ssq_down;

bb_up = -r*ones(Nx-1,1)/(2*h).*S_up;
bb_down = r*ones(Nx-1,1)/(2*h).*S_down;

cc_main = r*ones(Nx,1);


f_old = forcing(1:end-1,1);
corr_right_old = 0.5*sigma^2*bc_right(0)/h^2*S_max.^2 +...
    r*bc_right(0)/(2*h)*S_max;
f_old(end) =  f_old(end) + corr_right_old;

for tn = 1:Nt
    
    t_new=tn*deltat;

    % forcing at the new time-step
    f_new = forcing(1:end-1,tn+1);
    corr_right_new = 0.5*sigma^2*bc_right(t_new)/h^2*S(end).^2 +...
        r*bc_right(t_new)/(2*h)*S(end);
    f_new(end) =  f_new(end) + corr_right_new;

    
    A_old_up   = aa_old_up+bb_up;
    A_old_main = aa_old_main+cc_main;
    A_old_down = aa_old_down+bb_down;
    A_old = spdiags([NaN; A_old_up],1,Nx,Nx) +...
        spdiags(A_old_main,0,Nx,Nx) + spdiags([A_old_down; NaN],-1,Nx,Nx);

    f_tot = V(1:end-1,tn) - deltat*(1-theta)*A_old*V(1:end-1,tn) +...
        deltat*theta*f_new + deltat*(1-theta)*f_old;

    
    % now build the diagonals for A new and B
    aa_new_up   = -0.5/h^2*sigma^2.*Ssq_up;
    aa_new_main =  1/h^2*sigma^2.*Ssq;
    aa_new_down = -0.5/h^2*sigma^2.*Ssq_down;
        
    B_up = deltat*theta*(aa_new_up+bb_up);
    B_main = deltat*theta*(aa_new_main+cc_main)+ones(Nx,1);
    B_down = deltat*theta*(aa_new_down+bb_down);
        
         
    % solve the system and store the solution
    sol = thomas(B_down,B_main,B_up,f_tot);
    V(:,tn+1)=[sol; bc_right(t_new)];

    % update the forcing at the previous time step
    f_old = f_new;
    
    % update aa_* vectors
    aa_old_up = aa_new_up;
    aa_old_main = aa_new_main;
    aa_old_down = aa_new_down;

end
