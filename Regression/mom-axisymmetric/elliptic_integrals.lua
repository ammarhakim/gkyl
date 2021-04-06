-- Adapted from https://github.com/ElsevierSoftwareX/SOFTX_2018_246/
-- No sanity or special case chesk are performed.
-- Motivated to implement magnetic field potential in Jackson 3rd ed. eq (3.57).
-- Use like el = require "elliptic_integrals".


local function rf( x, y, z)
--RF evaluates the Carlson's incomplete or complete elliptic
--   integral of the 1st kind. 
--
--                       inf
--                       | |  
--                   1   |            dt
--       rf(x,y,z) = -   |  -----------------------
--                   2   |  sqrt((t+x)*(t+y)*(t+z)) 
--                     | | 
--                      0
--
--   Result:
--       rf(x,y,z) -- Real scalar or NaN if either argument is invalid or
--           convergence failed.
--
--   Arguments:
--       x,y,z  ---  real scalars >0, at most one can be zero
--
--   References:
--   [1] B.C.Carlson,"Numerical computation of real or complex elliptic 
--       integrals", Numerical Algorithms 10, 1995, pp 13-26.
--   [2] B.C.Carlson, "Elliptic Integrals" In: F.W.J.Olver et al., NIST 
--       Handbook of Mathematical Functions, Cambridge, 2010.
--   [3] W.H.Press et al., "Numerical Recipes in Fortran 90", 2nd ed., 
--       Cambridge, 1996

    -- General case (based on rf_s from [3])

    local sqrt = math.sqrt
    local abs = math.abs
    local ERRTOL = 0.001;

    local rx, ry, rz
    local lm, av
    local dx, dy, dz
    while true do
        rx = sqrt(x);
        ry = sqrt(y);
        rz = sqrt(z);
        lm  = rx*(ry + rz) + ry*rz;
        x  = 0.25*(lm + x);
        y  = 0.25*(lm + y);
        z  = 0.25*(lm + z);
        av  = (x + y + z)/3;
        dx = (av - x)/av;
        dy = (av - y)/av;
        dz = (av - z)/av;
        if abs(dx) < ERRTOL and abs(dy) <ERRTOL and abs(dz) < ERRTOL then
            break
        end
    end
    local e2 = dx*dy - dz*dz;
    local e3 = dx*dy*dz;
    
    local C1 = 1/24;
    local C2 = 1/10;
    local C3 = 3/44;
    local C4 = 1/14;
    local s  = e2*(C1*e2 - C2 - C3*e3) + C4*e3;
    
    result = (1 + s)/sqrt(av);

    return result
end


local function rd( x, y, z)
--RD Evaluates the Carlson degenerate case of an elliptic integral of
--   the 2nd kind [1]
--
--                       inf
--                       | |
--                   3   |                  dt
--       rd(x,y,z) = -   |  ------------------------------------
--                   2   |  (t + z)*sqrt((t + x)*(t + y)*(t + z))
--                     | |
--                      0
--
--   Result:
--       rd(x,y,z) -- Real scalar or NaN if either argument is invalid or
--           convergence failed.
--
--   Arguments:
--       x,y    ---  real scalars >0, only one can be zero
--       z      ---  real scalar  >0
--
--   References:
--   [1] B.C.Carlson,"Numerical computation of real or complex elliptic
--       integrals", Numerical Algorithms 10, 1995, pp 13-26.
--   [2] B.C.Carlson, "Elliptic Integrals" In: F.W.J.Olver et al., NIST
--       Handbook of Mathematical Functions, Chap.19, Cambridge, 2010.
--       Handbook of Mathematical Functions, Cambridge, 2010.
--   [3] W.H.Press et al., "Numerical Recipes in Fortran 90", 2nd ed.,
--       Cambridge, 1996
--
    -- General case (based on rd_s from [3])

    local sqrt = math.sqrt
    local abs = math.abs
    local ERRTOL = 0.002; -- By tests

    local sm = 0;
    local fc = 1;
    local rx, ry, rz
    local lm, m
    local dx, dy, dz
    while true do
        rx = sqrt(x);
        ry = sqrt(y);
        rz = sqrt(z);
        lm  = rx*(ry + rz) + ry*rz;
        sm = sm + fc/(rz*(z + lm));
        fc = 0.25*fc;
        x  = 0.25*(lm + x);
        y  = 0.25*(lm + y);
        z  = 0.25*(lm + z);
        m  = (x + y + 3*z)/5;
        dx = (m - x)/m;
        dy = (m - y)/m;
        dz = (m - z)/m;
        if abs(dx) < ERRTOL and abs(dy) < ERRTOL and abs(dz) < ERRTOL then
            break
        end
    end
    local ea = dx*dy;
    local eb = dz*dz;
    local ec = ea - eb;
    local ed = ea - 6*eb;
    local ee = ed + ec + ec;

    local C1  = -3/14;
    local C2  = 1/6;
    local C3  = 9/22;
    local C4  = 3/26;
    local C5  = 9/88;
    local C6  = 9/52;
    local s   = ed*(C1 + C5*ed - C6*dz*ee);
    s   = s + dz*(C2*ee + dz*(dz*C4*ea - C3*ec));

    result  = 3*sm + fc*(1 + s)/(m*sqrt(m));

    return result
end


local function K(k)
--ELK  Evaluates the complete elliptic integral of the 1st kind
--
--                1
--               | |  
--               |              dt
--       K(k) =  |  ------------------------------
--               |  sqrt((1 - t^2)*(1 - k^2*t^2)) 
--             | | 
--              0
--
--   Result:
--       K(k) -- real scalar or NaN if argument is invalid or 
--           convergence failes. 
--
--   Arguments:
--       k  -- real scalar, modulus |k|<=1
   return rf(0, 1-k*k, 1)
end


local function E(k)
--ELE Evaluates the complete elliptic integrals of the 2nd kind
--
--                  1
--                 | |  
--                 |  sqrt(1 - k^2*t^2)
--       E(x,k) =  |  ----------------- dt
--                 |  sqrt(1 -    t^2) 
--               | | 
--                0
--
--   Result:
--       E(k)   -- (complete integral) real scalar or NaN if either 
--           argument is invalid or convergenece failed. 
--
--   Arguments:
--       k    -- real scalar, modulus |k*x| <= 1
   return K(k) - k*k/3 * rd(0, 1-k*k, 1)
end

return {
   rf = rf,
   rd = rd,
   K = K,
   E = E,
}
