% FLUTTER SPEED
% This function looks for the instability speed from Vwind_initial to
% Vwind_final, up to a tolerance given by 'tol' using BINARY SEARCH.
%
%
% some limitations: 1) it works only for aeroelasticity cases (without
% rigid body degrees of freedom)
% 2) it does not compute the equilibrium condition at each speed.
% For this reason, it only works for wings which the equilibrium condition
% doesn't change with speed! like that one of example 2.
% It is very easy to make it more general, but will make the code much
% slower.
%
% yet to do: save all eigenvectors/eigenvalues during the search;
% this would allow visualizing how the eigenvalues change with speed.

function [unstable_speed, unstable_eig_value, unstable_eig_vec] = ...
           flutter_speed(Vwind_initial, Vwind_final, tol, ap, strain_eq, altitude)
   throttle = 0; deltaflap = 0;
   betaeq = zeros(6,1);
   keq = [0;0;0;altitude];
    [~, Aaeroelast, ~] = linearize(ap, strain_eq, betaeq, keq, throttle, deltaflap, Vwind_initial);
   if max(real(eig(Aaeroelast)))>0
       fprintf('Initial speed is already unstable! \n Cant find instability speed in this interval ');
   else
       [~, Aaeroelast, ~] = linearize(ap, strain_eq, betaeq, keq, throttle, deltaflap, Vwind_final);
       if max(real(eig(Aaeroelast)))<0
           fprintf('Final speed is stable! \n Cant find instability speed in this interval! ');
       else
           diff = Vwind_final - Vwind_initial;
           Vnew = (Vwind_final + Vwind_initial)/2;
           while (diff > tol)
               [~, Aaeroelast, ~] = linearize(ap, strain_eq, betaeq, keq, throttle, deltaflap, Vnew);
               [eig_vec, eig_val] = eig(Aaeroelast);
               eig_val = diag(eig_val);
               [real_max, ind_max] = max(real(eig_val));
               if real_max >0
                   Vwind_final = Vnew;
               else
                   Vwind_initial = Vnew;
               end
               diff = Vwind_final - Vwind_initial;
               Vnew = (Vwind_final + Vwind_initial)/2
           end
           unstable_speed = Vwind_initial;
           unstable_eig_value = eig_val(ind_max);
           unstable_eig_vec = eig_vec(:,ind_max);
       end
   end
end