%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% UPWIND_HJ_WENO5_1D() computes the fifth-order upwind HJ WENO 
% approximation to phi_x.
%
% Usage: phi_x = UPWIND_HJ_WENO5_1D(phi, vel_x, ghostcell_width, dx)
%
% Arguments:
% - phi:               function for which to compute upwind 
%                        derivative
% - vel_x:             velocity to use in upwinding
% - ghostcell_width:   number of ghostcells at boundary of 
%                        computational domain
% - dx:                grid cell size
%
% Return values:
% - phi_x:             fifth-order, upwind HJ WENO derivative
%
% NOTES:
% - phi_x has the same ghostcell width as phi.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyrights: (c) 2005 The Trustees of Princeton University and Board of
%                 Regents of the University of Texas.  All rights reserved.
%             (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:   $Revision$
% Modified:   $Date$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
