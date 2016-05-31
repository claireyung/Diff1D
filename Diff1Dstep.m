function [F,FRIC] = Diff1Dstep(Fp,kappa,gam,Hz,Hzw,FFlux,BForce,Nz,dt)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                  DIFF1Dstep                                        %
 % This function advances the field Fp (F previous) by one time step. %
 % This is done using a flux formulation where first the diffusive    %
 % fluxes are calculated using the diffusivity kappa, then forcing    %
 % fluxes are added (FFlux). There is also the option to add a body   %
 % force on the RHS (BForce).                                         %
 %                                                                    %
 % Fp = profile of field at current time step on rho grid.            %
 %                                                                    %
 % kappa = vector of diffusivities of length Nz+1 (on w-points).      %
 %                                                                    %
 % Hz = width of rho blocks.                                          %
 %                                                                    %
 % Hzw = width of w blocks.                                           %
 %                                                                    %
 % FFlux = additional fluxes (such as forcing/T-S fluxes) added.      %
 %                                                                    %
 % BForce = body force (rho-points) added to RHS.                     %
 %                                                                    %
 % Ryan Holmes June 2014                                              %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Centered 2nd-order:%%%
    
    % Calculate fluxes:
   Flux = [0; ...
        -kappa(2:(end-1)).*((Fp(2:end)-Fp(1:(end-1)))./Hzw-gam(2:(end-1))); ...
                       0];

    % Add forcing fluxes:
    Flux = Flux + FFlux;

    % step forward field:
    FRIC = -(Flux(2:end)-Flux(1:(end-1)))./Hz;
    F = Fp + dt*FRIC+dt*BForce;
end