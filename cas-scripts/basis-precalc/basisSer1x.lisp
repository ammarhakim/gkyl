;;; -*- Mode: LISP; package:maxima; syntax:common-lisp; -*- 
(in-package :maxima)
(DSKSETQ |$varsC| '((MLIST SIMP) $X)) 
(ADD2LNC '|$varsC| $VALUES) 
(DSKSETQ |$basisC|
         '((MLIST SIMP
            (9.
             #A((90.) BASE-CHAR
                . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/basis-precalc/basis-pre-cdim-calc.mac")
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (32.
              #A((68.) BASE-CHAR
                 . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X))
           ((MLIST SIMP
             (32.
              #A((68.) BASE-CHAR
                 . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.))))
           ((MLIST SIMP
             (32.
              #A((68.) BASE-CHAR
                 . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.))))
           ((MLIST SIMP
             (32.
              #A((68.) BASE-CHAR
                 . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 4.))))
           ((MLIST SIMP
             (32.
              #A((68.) BASE-CHAR
                 . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 4.)))
            ((MTIMES SIMP) 63. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 5.))))
           ((MLIST SIMP
             (32.
              #A((68.) BASE-CHAR
                 . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 4.)))
            ((MTIMES SIMP) 63. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) 231. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62.
                     #A((68.) BASE-CHAR
                        . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 6.))))
           ((MLIST SIMP
             (32.
              #A((68.) BASE-CHAR
                 . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 4.)))
            ((MTIMES SIMP) 63. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) 231. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62.
                     #A((68.) BASE-CHAR
                        . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 6.)))
            ((MTIMES SIMP) 429. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62.
                     #A((68.) BASE-CHAR
                        . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 5.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 7.))))
           ((MLIST SIMP
             (32.
              #A((68.) BASE-CHAR
                 . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 4.)))
            ((MTIMES SIMP) 63. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 11. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 3.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 5.)))
            ((MTIMES SIMP) 231. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 13. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
              ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62.
                     #A((68.) BASE-CHAR
                        . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 4.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 6.)))
            ((MTIMES SIMP) 429. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -35. 33.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 3.)))
              ((MTIMES SIMP) ((RAT SIMP) -21. 13.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 7.) $X)
                ((MTIMES SIMP) ((RAT SIMP) -10. 9.)
                 ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
                  ((MEXPT SIMP
                    (62.
                     #A((68.) BASE-CHAR
                        . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                     SRC |$calcPowers| 62.))
                   $X 3.)))
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 5.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 7.)))
            ((MTIMES SIMP) 6435. ((MEXPT SIMP) 2. ((RAT SIMP) -15. 2.))
             ((MEXPT SIMP) 17. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -20. 33.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -210. 143.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62.
                     #A((68.) BASE-CHAR
                        . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 4.)))
              ((MTIMES SIMP) ((RAT SIMP) -28. 15.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 7.)
                ((MTIMES SIMP) ((RAT SIMP) -5. 7.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                  ((MEXPT SIMP
                    (62.
                     #A((68.) BASE-CHAR
                        . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                     SRC |$calcPowers| 62.))
                   $X 2.)))
                ((MTIMES SIMP) ((RAT SIMP) -15. 11.)
                 ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
                  ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
                   ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                    ((MEXPT SIMP
                      (62.
                       #A((68.) BASE-CHAR
                          . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                       SRC |$calcPowers| 62.))
                     $X 2.)))
                  ((MEXPT SIMP
                    (62.
                     #A((68.) BASE-CHAR
                        . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                     SRC |$calcPowers| 62.))
                   $X 4.)))
                ((MEXPT SIMP
                  (62.
                   #A((68.) BASE-CHAR
                      . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 6.)))
              ((MEXPT SIMP
                (62.
                 #A((68.) BASE-CHAR
                    . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 8.)))))) 
(ADD2LNC '|$basisC| $VALUES) 
(DSKSETQ |$basisConstant|
         '((MLIST SIMP
            (32.
             #A((68.) BASE-CHAR
                . "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac")
             SRC |$gsOrthoNorm| 30.))
           ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.)))) 
(ADD2LNC '|$basisConstant| $VALUES) 