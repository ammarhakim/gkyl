-- This unit test tests the functionality of sci/diff-recursive.lua, which implements
-- recursive dual numbers for arbitrary order automatic differentiation.

local math = require("sci.math").generic
local root = require("sci.root")
local ffi = require("ffi")
local quad = require("sci.quad")
local xsys = require("xsys")
local diff = require("sci.diff-recursive")
local Unit = require "Unit"
local assert_close = Unit.assert_close
local assert_equal = Unit.assert_equal
local stats = Unit.stats

local deriv = diff.deriv

-- test derivatives up to second order of polynomial function
function test_deriv()
   local k = 2
 
   local function Psi(R, Z)
     return 4*R^2*Z^2/k^2+(R^2-1)^2
   end
 
   -- analytical functions for checking
   local function _dPsidR(R,Z)
     return 4*R*(R^2-1) + 8*R*Z^2/k^2
   end
 
   local function _dPsidZ(R,Z)
     return 8*R^2*Z/k^2
   end
 
   local function _d2PsidR2(R,Z)
     return 8*Z^2/k^2 + 12*R^2 - 4
   end
 
   local function _d2PsidRdZ(R,Z)
     return 16*R*Z/k^2
   end
 
   local function _d2PsidZ2(R,Z)
     return 8*R^2/k^2
   end
 
   -- compute derivatives via automatic differentiation
   local function dPsidR(R,Z)
     -- differentiate Psi wrt arg 1 (R)
     return deriv(Psi,1)(R,Z)
   end
 
   local function dPsidZ(R,Z)
     -- differentiate Psi wrt arg 2 (Z)
     return deriv(Psi,2)(R,Z)
   end
 
   local function d2PsidR2(R,Z)
     return deriv(deriv(Psi,1),1)(R,Z)
   end
 
   local function d2PsidRdZ(R,Z)
     return deriv(deriv(Psi,1),2)(R,Z)
   end
 
   local function d2PsidZ2(R,Z)
     return deriv(deriv(Psi,2),2)(R,Z)
   end
 
   assert_equal(dPsidR(1.1, .9), _dPsidR(1.1, .9), "Checking dPsidR")
   assert_equal(dPsidZ(1.1, .9), _dPsidZ(1.1, .9), "Checking dPsidZ")
   assert_equal(d2PsidR2(1.1, .9), _d2PsidR2(1.1, .9), "Checking d2PsidR2")
   assert_equal(d2PsidRdZ(1.1, .9), _d2PsidRdZ(1.1, .9), "Checking d2PsidRdZ")
   assert_equal(d2PsidZ2(1.1, .9), _d2PsidZ2(1.1, .9), "Checking d2PsidZ2")
end
  
-- start from an analytical Solovev equilibrium and generate a field-aligned coordinate system (Psi, alpha, theta).
function test_solovev()
   -- geometry parameters
   local B0 = 1
   local R0 = 1
   local k = 2
   local Ztop = 2
   
   -- Solovev equilibrium
   -- Psi(R,Z) = poloidal flux function in terms of cylindrical (R,Z) coordiantes
   local function Psi(R, Z)
     return 4*R^2*Z^2/k^2+(R^2-1)^2
   end
   
   -- find separatrix
   local PsiSep = Psi(0,0)
   
   -- analytical expressions, for checking
   local function Z_(R, Psi0)
      return k/(2*R)*math.sqrt(Psi0-(R^2-1)^2)
   end
   
   local function dZdR_(R, Psi0)
      return -k*(Psi0 + R^4 -1)/(2*R^2*math.sqrt(Psi0-(R^2-1)^2))
   end
   
   local function d2ZdR2_(R, Psi0)
      return -k*(-Psi0^2 + (R^2-1)^3 + Psi0*(3*R^4-3*R^2+2))/(R^3*(Psi0 - (R^2-1)^2)^(3/2))
   end
   
   local function Rin_(Psi0)
      if tonumber(tostring(Psi0)) <= tonumber(tostring(PsiSep)) then
         return math.sqrt(1-math.sqrt(Psi0))
      else
         return math.sqrt(1 - 2*Ztop^2/k^2 + math.sqrt(Psi0 - 4*Ztop^2/k^2 + 4*Ztop^4/k^4))
      end
   end
   
   local function Rout_(Psi0)
      return math.sqrt(1+math.sqrt(Psi0))
   end
   
   local function norm_(Psi0)
      local function integrand(t)
         return math.sqrt(1+dZdR_(t, Psi0)^2)
      end
      return quad.dblexp(integrand, Rin_(Psi0), Rout_(Psi0), 1e-10)/math.pi
   end
   
   -- stopping criteria for root-finding
   local function stop(tol)
      return function(x, y, xl, xu, yl, yu)
         if diff.lt(math.abs(y), tol) then return true
         else return false
         end
      end
   end
   
   -- numerically invert Psi to get Z(R,Psi) 
   local tol = 1e-14
   local function Z(R, Psi0)
      -- solve f(Zi) = Psi(R,Zi) - Psi0 = 0 for Zi
      local function f(Zi) 
         return Psi(R,Zi)-Psi0
      end
      local val = root.ridders(f, 0, 10, stop(tol))
      return val
   end
   assert_close(Z(0.1,1.2), Z_(0.1,1.2), tol, "Checking Z(R=0.1,Psi=1.2)")
   assert_close(Z(0.9,0.9), Z_(0.9,0.9), tol, "Checking Z(R=0.9,Psi=0.9)")
   
   local function Rout(Psi0)
      local function f(Ri)
         return Psi(Ri,0)-Psi0
      end
      local val = root.ridders(f, R0, 5*R0, stop(tol))
      return val
   end
   assert_close(Rout(.8), Rout_(.8), tol, "Checking Rout(Psi=0.8)")
   assert_close(Rout(1.1), Rout_(1.1), tol, "Checking Rout(Psi=1.1)")
   assert_close(Rout(2.25), Rout_(2.25), tol, "Checking Rout(Psi=2.25)")
   
   local function Rin(Psi0)
      if tonumber(tostring(Psi0)) <= tonumber(tostring(PsiSep)) then
         local function f(Ri)
            return Psi(Ri,0)-Psi0
         end
         local val = root.ridders(f, 0, 1, stop(tol))
         return val
      else
         local function f(Ri)
            return Psi(Ri,Ztop)-Psi0
         end
         local val = root.ridders(f, 0, Rout(Psi0), stop(tol))
         return val
      end
   end
   assert_close(Rin(.8), Rin_(.8), tol, "Checking Rin(Psi=0.8)")
   assert_close(Rin(1.1), Rin_(1.1), tol, "Checking Rin(Psi=1.1)")
   assert_close(Rin(2.25), Rin_(2.25), tol, "Checking Rin(Psi=2.25)")
   
   local function dRindPsi(Psi)
      return deriv(Rin)(Psi)
   end
   local function dRindPsi_(Psi)
      return deriv(Rin_)(Psi)
   end
   assert_close(dRindPsi(.8), dRindPsi_(.8), 1e-10, "Checking dRindPsi(Psi=0.8)")
   assert_close(dRindPsi(1.1), dRindPsi_(1.1), 1e-10, "Checking dRindPsi(Psi=1.1)")
   assert_close(dRindPsi(2.25), dRindPsi_(2.25), 1e-10, "Checking dRindPsi(Psi=2.25)")
   
   local function dRoutdPsi(Psi)
      return deriv(Rout)(Psi)
   end
   local function dRoutdPsi_(Psi)
      return deriv(Rout_)(Psi)
   end
   assert_close(dRoutdPsi(.8), dRoutdPsi_(.8), 1e-10, "Checking dRoutdPsi(Psi=0.8)")
   assert_close(dRoutdPsi(1.1), dRoutdPsi_(1.1), 1e-10, "Checking dRoutdPsi(Psi=1.1)")
   assert_close(dRoutdPsi(2.25), dRoutdPsi_(2.25), 1e-10, "Checking dRoutdPsi(Psi=2.25)")
   
   local function dZdR(R, Psi)
      return deriv(Z,1)(R,Psi)
   end
   assert_close(dZdR(.9,1.2), dZdR_(.9,1.2), 1e-10, "Checking dZ/dR(R=0.9,Psi=1.2)")
   assert_close(dZdR(.1,1.2), dZdR_(.1,1.2), 1e-10, "Checking dZ/dR(R=0.1,Psi=1.2)")
   assert_close(dZdR(1,2.25), dZdR_(1,2.25), 1e-10, "Checking dZ/dR(R=1,Psi=2.25)")
   
   local function norm(Psi)
      local function integrand(t)
         return math.sqrt(1+dZdR(t, Psi)^2)
      end
      return quad.dblexp(integrand, Rin(Psi), Rout(Psi), 1e-10)/math.pi
   end
   assert_close(norm(0.1), norm_(0.1), 1e-7, "Checking norm(Psi=0.1)")  
   assert_close(norm(1.1), norm_(1.1), 1e-7, "Checking norm(Psi=1.1)")  
   assert_close(norm(2.25), norm_(2.25), 1e-7, "Checking norm(Psi=2.25)") 
   
   local function theta(R,Z)
      local Psi = Psi(R,Z)
      local function integrand(t)
         return math.sqrt(1+dZdR(t, Psi)^2)
      end
      local integral, _ = quad.dblexp(integrand, R, Rout(Psi), 1e-10)
      return math.sign(Z)/norm(Psi)*integral
   end
   assert_close(theta(0.2, 0.0), math.pi, 1e-8, "Checking theta(R=0.2,Z=0.0)")
   assert_close(theta(1.2, 0.0), 0, 1e-8, "Checking theta(R=1.2,Z=0.0)")
   --print(string.format("theta(R=0.2,Z=0.0) =\t %.16f", theta(0.2,0.0)))
   --print(string.format("theta(R=1.2,Z=0.0) =\t %.16f", theta(1.2,0.0)))
   --print(string.format("theta(R=1.2,Z=1e-6) =\t %.16f", theta(1.2,1e-6)))
   --print(string.format("theta(R=0.2,Z=-0.1) =\t %.16f", theta(0.2,-0.1)))
   --print(string.format("theta(R=1.2,Z=0.1) =\t %.16f", theta(1.2,0.1)))
   --print(string.format("theta(R=1.0,Z=1.5) =\t %.16f", theta(1.0,1.5)))
   --print(string.format("theta(R=1.0,Z=1.9) =\t %.16f", theta(1.0,1.9)))
   
   local function dthetadR(R, Z)
      return deriv(theta,1)(R,Z)
   end
   --print(string.format("dtheta/dR(R=1.0,Z=1.5) =\t %.16f", dthetadR(1.0,1.5)))
   
   local function dthetadZ(R, Z)
      return deriv(theta,2)(R,Z)
   end
   --print(string.format("dtheta/dZ(R=1.0,Z=1.5) =\t %.16f", dthetadZ(1.0,1.5)))
   
   local function R(Psi, theta0)
      local function f(R)
         return theta(R,Z(R,Psi))-theta0
      end
      local Rval = root.ridders(f, Rin(Psi)+1e-10, Rout(Psi)-1e-10, stop(1e-10))
      return Rval
   end
   
   local function dRdPsi(Psi, theta)
      return deriv(R,1)(Psi, theta)
   end
   
   local function RZ(Psi, theta0)
      local function f(R)
         return theta(R,Z(R,Psi))-theta0
      end
      local Rval = root.ridders(f, Rin(Psi)+1e-10, Rout(Psi)-1e-10, stop(1e-10))
      local Zval = Z(Rval, Psi)
      return Rval, Zval
   end
   --print(string.format("RZ(Psi=.5,theta=pi/3) =\t %.16f  %.16f", RZ(.5,math.pi/3)))
   
   local function Bphi(R)
      return B0*R0/R
   end
   
   local function grad_Psi(R, Z)
      return deriv(Psi,1)(R,Z), deriv(Psi,2)(R,Z), 0
   end
   
   local function grad_theta(R, Z)
      return deriv(theta,1)(R,Z), deriv(theta,2)(R,Z), 0
   end
   
   local function jacob_PsiThetaPhi(R,Z)
      local gradPsi_R, gradPsi_Z, _ = grad_Psi(R,Z)
      return norm(Psi(R,Z))*R/math.sqrt(gradPsi_R^2 + gradPsi_Z^2)
   end
   
   --jacob_PsiThetaPhi(0.9,0.5)
   
   local function alpha(Ri, Zi, phi)
      local Psi = Psi(Ri,Zi)
      local function integrand(t)
         local Zt = Z(t, Psi) -- need to parametrize Z with R=t for integral
         local gradPsi_R, gradPsi_Z, _ = grad_Psi(t,Zt)
         local gradPsi = math.sqrt(gradPsi_R^2 + gradPsi_Z^2)
         return math.sqrt(1+dZdR(t, Psi)^2)/t/gradPsi
      end
      local integral, _ = quad.dblexp(integrand, Ri, Rout(Psi)-1e-10, 1e-10)
      return phi - Bphi(Ri)*Ri*integral
   end
   
   local function grad_alpha(R, Z, phi)
      return deriv(alpha,1)(R,Z,phi), deriv(alpha,2)(R,Z,phi), deriv(alpha,3)(R,Z,phi)/R
   end

   --print(string.format("grad_Psi(R=0.9,Z=0.5) =\t %.16f %.16f %.16f", grad_Psi(0.9,0.5)))
   --print(string.format("grad_alpha(R=0.9,Z=0.5,phi=Pi) =\t %.16f %.16f %.16f", grad_alpha(0.9,0.5,math.pi)))
   
   -- compute bmag = grad(alpha) X grad(psi)
   local function bmag(R, Z, phi)
      local phi = phi or 0 -- bmag is axisymmetric, so there should be no phi dependence
      local grad_Psi_R, grad_Psi_Z, grad_Psi_phi = grad_Psi(R,Z)
      local grad_alpha_R, grad_alpha_Z, grad_alpha_phi = grad_alpha(R,Z,phi)
   
      local B_R = grad_alpha_Z*grad_Psi_phi - grad_alpha_phi*grad_Psi_Z
      local B_Z = grad_alpha_phi*grad_Psi_R - grad_alpha_R*grad_Psi_phi
      local B_phi = grad_alpha_R*grad_Psi_Z - grad_alpha_Z*grad_Psi_R
      return math.sqrt(B_R^2 + B_Z^2 + B_phi^2), B_R, B_Z, B_phi
   end
   

   local B, B_R, B_Z, B_phi = bmag(0.5, 0.1)
   assert_close(B_phi, Bphi(0.5), 1e-10, "Checking Bphi = grad(alpha)Xgrad(psi).hat(phi) = B0*R0/R")
end

-- Run tests.
test_deriv()
test_solovev()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
