;;; -*- Mode: LISP; package:maxima; syntax:common-lisp; -*- 
(in-package :maxima)
(DSKSETQ |$varsC| '((MLIST SIMP) $X)) 
(ADD2LNC '|$varsC| $VALUES) 
(DSKSETQ |$varsP| '((MLIST SIMP) $X $VX)) 
(ADD2LNC '|$varsP| $VALUES) 
(DSKSETQ |$basisC|
         '((MLIST SIMP
            (10.
             "/Users/JunoRavin/gkyl/cas-scripts/basis-precalc/basis-pre-calc.mac"
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
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
               $X 4.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
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
            ((MTIMES SIMP) 63. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 5.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
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
            ((MTIMES SIMP) 63. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) 231. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 6.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
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
            ((MTIMES SIMP) 63. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) 231. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 6.)))
            ((MTIMES SIMP) 429. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 5.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 7.))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
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
            ((MTIMES SIMP) 63. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) 231. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 6.)))
            ((MTIMES SIMP) 429. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 5.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 7.)))
            ((MTIMES SIMP) 6435. ((MEXPT SIMP) 2. ((RAT SIMP) -15. 2.))
             ((MEXPT SIMP) 17. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -20. 33.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -210. 143.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -28. 15.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
                ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 4.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 6.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 8.)))))) 
(ADD2LNC '|$basisC| $VALUES) 
(DSKSETQ |$basisP|
         '((MLIST SIMP
            (10.
             "/Users/JunoRavin/gkyl/cas-scripts/basis-precalc/basis-pre-calc.mac"
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $X))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 45. 8.)
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
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X))))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 45. 8.)
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
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 5.)))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 33. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 5.))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 33. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 5.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X))))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 45. 8.)
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
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 5.)))
            ((MTIMES SIMP) ((RAT SIMP) 175. 8.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
              ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 33. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 5.))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 33. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 5.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 6.)))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 6.)))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 39. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 4.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 6.))))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 39. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 6.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 4.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X))))))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 45. 8.)
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
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 5.)))
            ((MTIMES SIMP) ((RAT SIMP) 175. 8.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
              ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 33. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 5.))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 33. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 5.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 6.)))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 6.)))
            ((MTIMES SIMP) ((RAT SIMP) 75. 32.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 25.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -18. 35.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 3.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 4.))))))
            ((MTIMES SIMP) ((RAT SIMP) 75. 32.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 25.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -18. 35.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 4.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))))))
            ((MTIMES SIMP) ((RAT SIMP) 189. 32.)
             ((MEXPT SIMP) 55. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -10. 27.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 5.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 5.)))))
            ((MTIMES SIMP) ((RAT SIMP) 189. 32.)
             ((MEXPT SIMP) 55. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -10. 27.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 5.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 5.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 3.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))))))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 39. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 4.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 6.))))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 39. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 6.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 4.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))))))
            ((MTIMES SIMP) ((RAT SIMP) 429. 32.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 5.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 7.)))
            ((MTIMES SIMP) ((RAT SIMP) 429. 32.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 5.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 7.)))
            ((MTIMES SIMP) ((RAT SIMP) 1287. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 5.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 7.))))
            ((MTIMES SIMP) ((RAT SIMP) 1287. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 7.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 5.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 3.)
                   $X))))))))
           ((MLIST SIMP
             (32. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
              |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 45. 8.)
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
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
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
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 35. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 5.)))
            ((MTIMES SIMP) ((RAT SIMP) 175. 8.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
              ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 33. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 5.))))
            ((MTIMES SIMP) ((RAT SIMP) 63. 16.)
             ((MEXPT SIMP) 33. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 5.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 6.)))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 6.)))
            ((MTIMES SIMP) ((RAT SIMP) 75. 32.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 25.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -18. 35.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 3.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 4.))))))
            ((MTIMES SIMP) ((RAT SIMP) 75. 32.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 25.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -18. 35.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 4.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))))))
            ((MTIMES SIMP) ((RAT SIMP) 189. 32.)
             ((MEXPT SIMP) 55. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -10. 27.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 5.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 5.)))))
            ((MTIMES SIMP) ((RAT SIMP) 189. 32.)
             ((MEXPT SIMP) 55. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -10. 27.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 5.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 5.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 3.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))))))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 39. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 4.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 6.))))
            ((MTIMES SIMP) ((RAT SIMP) 231. 32.)
             ((MEXPT SIMP) 39. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 6.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 4.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))))))
            ((MTIMES SIMP) ((RAT SIMP) 429. 32.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 5.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 7.)))
            ((MTIMES SIMP) ((RAT SIMP) 429. 32.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 3.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 5.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 7.)))
            ((MTIMES SIMP) ((RAT SIMP) 11025. 128.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 25.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 35.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -6. 35.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -36. 49.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
                ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $VX 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 4.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 4.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                  ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $VX 2.)))
                  ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 4.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 4.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                  ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $VX 2.)))
                  ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 4.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 4.)))))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MEXPT SIMP) 77. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 35.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -3. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -2. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $VX $X)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 3.)
                   $X)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.))))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 3.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 5.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 5.))))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MEXPT SIMP) 77. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 35.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -2. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 5.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 3.)
                   $X)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 5.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 3.))
              ((MTIMES SIMP) ((RAT SIMP) -3. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -9. 25.) $VX $X)
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 3.)
                   $X)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))
                ((MTIMES SIMP) ((RAT SIMP) -3. 5.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.))))))))
            ((MTIMES SIMP) ((RAT SIMP) 693. 64.)
             ((MEXPT SIMP) 65. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 21.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -5. 21.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP) ((RAT SIMP) -5. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                  ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $VX 2.)))
                  ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 4.))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 4.)))))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 6.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
                ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 4.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 6.)))))
            ((MTIMES SIMP) ((RAT SIMP) 693. 64.)
             ((MEXPT SIMP) 65. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 21.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 21.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -5. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
                ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $VX 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 4.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 6.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 6.)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 2.))
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
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
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 15.)
                ((MTIMES SIMP) ((RAT SIMP) -2. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $VX 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 4.)))
                ((MTIMES SIMP) ((RAT SIMP) -1. 5.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 4.)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
                  ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $VX 2.)))
                  ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))))))
            ((MTIMES SIMP) ((RAT SIMP) 1287. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 3.))))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 3.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 5.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 7.))))
            ((MTIMES SIMP) ((RAT SIMP) 1287. 32.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 7.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 3.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $VX $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 5.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 3.)
                   $X)))))))
            ((MTIMES SIMP) ((RAT SIMP) 6435. 256.)
             ((MEXPT SIMP) 17. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -20. 33.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -210. 143.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -28. 15.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
                ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $X 4.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $X 6.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $X 8.)))
            ((MTIMES SIMP) ((RAT SIMP) 6435. 256.)
             ((MEXPT SIMP) 17. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -20. 33.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -210. 143.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -28. 15.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
                ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                       SRC |$calcPowers| 62.))
                     $VX 2.)))
                  ((MEXPT SIMP
                    (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                     SRC |$calcPowers| 62.))
                   $VX 4.)))
                ((MEXPT SIMP
                  (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                   |$calcPowers| 62.))
                 $VX 6.)))
              ((MEXPT SIMP
                (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                 |$calcPowers| 62.))
               $VX 8.)))
            ((MTIMES SIMP) ((RAT SIMP) 6435. 256.)
             ((MEXPT SIMP) 51. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -20. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) ((RAT SIMP) -210. 143.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 4.))))
              ((MTIMES SIMP) ((RAT SIMP) -28. 15.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $VX)
                ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 2.))))
                ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                    ((MTIMES SIMP) $VX
                     ((MEXPT SIMP
                       (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                        SRC |$calcPowers| 62.))
                      $X 2.))))
                  ((MTIMES SIMP) $VX
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $X 4.))))
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $X 6.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $X 8.))))
            ((MTIMES SIMP) ((RAT SIMP) 6435. 256.)
             ((MEXPT SIMP) 51. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                  |$calcPowers| 62.))
                $VX 8.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -20. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP) ((RAT SIMP) -210. 143.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 4.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))))
              ((MTIMES SIMP) ((RAT SIMP) -28. 15.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 7.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac" SRC
                    |$calcPowers| 62.))
                  $VX 6.)
                 $X)
                ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 2.)
                   $X)))
                ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
                  ((MTIMES SIMP)
                   ((MEXPT SIMP
                     (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                      SRC |$calcPowers| 62.))
                    $VX 4.)
                   $X)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                    ((MTIMES SIMP)
                     ((MEXPT SIMP
                       (62. "/Users/JunoRavin/gkyl/cas-scripts/modal-basis.mac"
                        SRC |$calcPowers| 62.))
                      $VX 2.)
                     $X)))))))))))) 
(ADD2LNC '|$basisP| $VALUES) 