;;; -*- Mode: LISP; package:maxima; syntax:common-lisp; -*- 
(in-package :maxima)
(DSKSETQ |$varsC| '((MLIST SIMP) $X $Y $Z)) 
(ADD2LNC '|$varsC| $VALUES) 
(DSKSETQ |$varsP| '((MLIST SIMP) $X $Y $Z $VX $VY)) 
(ADD2LNC '|$varsP| $VALUES) 
(DSKSETQ |$basisC|
         '((MLIST SIMP
            (10.
             "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/basis-precalc/basis-pre-calc.mac"
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $X $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $X $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $Y $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $X $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $X $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $Z 2.)))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))))) 
(ADD2LNC '|$basisC| $VALUES) 
(DSKSETQ |$basisP|
         '((MLIST SIMP
            (10.
             "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/basis-precalc/basis-pre-calc.mac"
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VY)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $X $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $X $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VY $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VY $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VY $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $VY)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $Y $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VY $X $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VY $X $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VY $Y $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $VY $X)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $VY $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $VY $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $X $Y
             $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VY $X $Y
             $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $VY $X
             $Y)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $VY $X
             $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $VY $Y
             $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.)) $VX $VY $X $Y $Z))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VY)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $X $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $X $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VY $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VY $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VY $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $VY)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $Z 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $VX 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $VY 2.)))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $Y $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VY $X $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VY $X $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VY $Y $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $VY $X)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $VY $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $VY $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $VY)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $X)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $X $Y
             $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VY $X $Y
             $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $VY $X
             $Y)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $VY $X
             $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.)) $VX $VY $Y
             $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $X $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $X $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Y)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X)
              ((MTIMES SIMP) $VY $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Z)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Z)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X)
              ((MTIMES SIMP) $VY $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Y)
              ((MTIMES SIMP) $VY $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $VY $X)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $VY $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $VY $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $X $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $X $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $X)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $Z)))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.)) $VX $VY $X $Y $Z)
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Y)
              ((MTIMES SIMP) $VX $X $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $X $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Y $Z)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X $Z)
              ((MTIMES SIMP) $VY $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X $Y)
              ((MTIMES SIMP) $VY $X $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $Y)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $X)
              ((MTIMES SIMP) $VX $VY $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $Z)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $Z)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $X)
              ((MTIMES SIMP) $VX $VY $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $Y)
              ((MTIMES SIMP) $VX $VY $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $VY $X $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $VY $X $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $VY $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $X $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $X $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $X $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $Y $Z)))
            ((MTIMES SIMP) 27. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $Y $Z)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 27. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $X $Z)
              ((MTIMES SIMP) $VX $VY $X
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 27. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY $X $Y)
              ((MTIMES SIMP) $VX $VY $X $Y
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $Z 2.))))
            ((MTIMES SIMP) 27. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $VY $X $Y $Z)))
            ((MTIMES SIMP) 27. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VY 2.)
               $X $Y $Z)))))) 
(ADD2LNC '|$basisP| $VALUES) 