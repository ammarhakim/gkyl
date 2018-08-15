;;; -*- Mode: LISP; package:maxima; syntax:common-lisp; -*- 
(in-package :maxima)
(DSKSETQ |$varsC| '((MLIST SIMP) $X $Y $Z)) 
(ADD2LNC '|$varsC| $VALUES) 
(DSKSETQ |$basisC|
         '((MLIST SIMP
            (9.
             "/Users/JunoRavin/gkyl/cas-scripts/basis-precalc/basis-pre-cdim-calc.mac"
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
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $X $Y)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $X $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)) $Y $Z)
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z))
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
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 27. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 27.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
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
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $Y 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $Z 2.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.)))))))
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
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               $Y $Z)))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 3.)))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 3.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Z)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 3.)))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 3.)))))
            ((MTIMES SIMP) 27. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 27.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
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
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $Y 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $Z 2.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 105. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 105. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 $Z)))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 105. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 105. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 105. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y $Z)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 105. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))))
            ((MTIMES SIMP) 175. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.))))))
            ((MTIMES SIMP) 175. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Z)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))))
            ((MTIMES SIMP) 175. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 $Z)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))))
            ((MTIMES SIMP) 225. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 15.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 2.))))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 2.))))))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 2.))))
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 2.))))))))
            ((MTIMES SIMP) 225. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 15.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Y 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.)
                   $Y)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $Y 3.)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $Y 3.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                  ((MTIMES SIMP) $Y
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 2.))))))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Y)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.)
                   $Y)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                  ((MTIMES SIMP) $Y
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 2.))))))))
            ((MTIMES SIMP) 225. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 15.) $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.)
                   $Z)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 2.)
                   $Z)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $Z 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.)
                   $Z)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $Z 3.)))))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 2.)
                   $Z)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $Z 3.)))))))
            ((MTIMES SIMP) 175. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Y $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 $Z)))))
            ((MTIMES SIMP) 175. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Y $Z)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               $Y
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))))
            ((MTIMES SIMP) 175. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 $Z)))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))))
            ((MTIMES SIMP) 21. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 25.) $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Y)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.))))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Y)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.)
                   $Y)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 3.))))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -9. 25.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Y)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.)
                   $Y)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                  ((MTIMES SIMP) $X $Y
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 2.))))))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Y)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 3.))))
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                  ((MTIMES SIMP) $X $Y
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 2.))))))))
            ((MTIMES SIMP) 21. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 25.) $X $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -9. 25.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 $Z)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.)
                   $Z)))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 2.)
                   $Z)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.)
                   $Z)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 3.))))))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 2.)
                   $Z)))
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 3.))))))))
            ((MTIMES SIMP) 21. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 25.) $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -9. 25.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.)
                   $Y $Z)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 3.)
                   $Z)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $Y $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 3.)
                   $Z)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                  ((MTIMES SIMP) $Y
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 3.))))))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.)
                   $Y $Z)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.)
                 $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
                  ((MTIMES SIMP) $Y
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 3.))))))))
            ((MTIMES SIMP) 125. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -27. 125.) $X $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -9. 25.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Y $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -9. 25.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.)
                   $Y $Z)))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 3.)
                   $Z)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Y 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $Z 3.))
              ((MTIMES SIMP) ((RAT SIMP) -9. 25.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Y $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.)
                   $Y $Z)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.)
                 $Y
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                  ((MTIMES SIMP) $X $Y
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 3.))))))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $X $Y $Z)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                  ((MTIMES SIMP) $X
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Y 3.)
                   $Z)))
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Y 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $Z 3.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
                  ((MTIMES SIMP) $X $Y
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $Z 3.))))))))))) 
(ADD2LNC '|$basisC| $VALUES) 
(DSKSETQ |$basisConstant|
         '((MLIST SIMP
            (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
             |$gsOrthoNorm| 30.))
           ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.)))) 
(ADD2LNC '|$basisConstant| $VALUES) 