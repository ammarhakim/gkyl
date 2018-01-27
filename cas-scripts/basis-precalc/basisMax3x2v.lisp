;;; -*- Mode: LISP; package:maxima; syntax:common-lisp; -*- 
(in-package :maxima)
(DSKSETQ |$varsC| '((MLIST SIMP) $X $Y $Z)) 
(ADD2LNC '|$varsC| $VALUES) 
(DSKSETQ |$varsP| '((MLIST SIMP) $X $Y $Z $VX $VY)) 
(ADD2LNC '|$varsP| $VALUES) 
(DSKSETQ |$basisC|
         '((MLIST SIMP
            (10.
             "/Users/JunoRavin/gkyl/cas-scripts/basis-precalc/basis-pre-calc.mac"
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
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
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 2.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
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
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 3.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
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
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 3.)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               $Y)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 4.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 4.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 2.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 4.)))))) 
(ADD2LNC '|$basisC| $VALUES) 
(DSKSETQ |$basisP|
         '((MLIST SIMP
            (10.
             "/Users/JunoRavin/gkyl/cas-scripts/basis-precalc/basis-pre-calc.mac"
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
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
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VY))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
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
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VY 2.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
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
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
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
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $VY)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $X)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VY)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VY 3.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
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
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
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
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $VY)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $X)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VY)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VY 3.)))
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
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Y)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X)
              ((MTIMES SIMP) $VY $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Z)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Z)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X)
              ((MTIMES SIMP) $VY $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Y)
              ((MTIMES SIMP) $VY $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $VY)
              ((MTIMES SIMP) $VX $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $VY $X)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $VY $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VY $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $VY $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $X $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $X $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $X)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $Y)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               $Z)))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VY 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VY 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VY 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VY 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               $Y)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $Y)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VY $X)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VY $Y)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VY $Z)
              ((MTIMES SIMP) $VY
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $VY)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $VY)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VY $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 3.)
               $X)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VY $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 3.)
               $Y)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VY $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $VY)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VY 3.))))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -11. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 4.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -11. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Y 4.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -11. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 2.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $Z 4.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -11. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 4.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -11. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VY 2.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VY 4.)))))) 
(ADD2LNC '|$basisP| $VALUES) 