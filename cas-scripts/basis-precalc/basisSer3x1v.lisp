;;; -*- Mode: LISP; package:maxima; syntax:common-lisp; -*- 
(in-package :maxima)
(DSKSETQ |$varsC| '((MLIST SIMP) $X $Y $Z)) 
(ADD2LNC '|$varsC| $VALUES) 
(DSKSETQ |$varsP| '((MLIST SIMP) $X $Y $Z $VX)) 
(ADD2LNC '|$varsP| $VALUES) 
(DSKSETQ |$basisC|
         '((MLIST SIMP
            (10.
             #A((71.) BASE-CHAR
                . "/Users/nmandell/Codes/gkyl/cas-scripts/basis-precalc/basis-pre-calc.mac")
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (32.
              #A((54.) BASE-CHAR
                 . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
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
             (32.
              #A((54.) BASE-CHAR
                 . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
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
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.)))))
           ((MLIST SIMP
             (32.
              #A((54.) BASE-CHAR
                 . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
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
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 3.)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y $Z)))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.)))))
           ((MLIST SIMP
             (32.
              #A((54.) BASE-CHAR
                 . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
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
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 3.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 3.)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) 9. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 45. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 4.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Y 2.)))
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 4.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Z 2.)))
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 4.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y $Z)))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) 15. ((MEXPT SIMP) 2. ((RAT SIMP) -5. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) 35. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y)))))
            ((MTIMES SIMP) 35. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.))))
            ((MTIMES SIMP) 35. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Z)))))
            ((MTIMES SIMP) 35. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) 35. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) 35. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) 315. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y $Z)))))
            ((MTIMES SIMP) 315. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) 315. ((MEXPT SIMP) 2. ((RAT SIMP) -9. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))))) 
(ADD2LNC '|$basisC| $VALUES) 
(DSKSETQ |$basisP|
         '((MLIST SIMP
            (10.
             #A((71.) BASE-CHAR
                . "/Users/nmandell/Codes/gkyl/cas-scripts/basis-precalc/basis-pre-calc.mac")
             SRC |$writeBasisToFile| 7.))
           ((MLIST SIMP
             (32.
              #A((54.) BASE-CHAR
                 . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 4.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $X $Y)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $X $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $Y)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Y)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 9. 4.) $VX $X $Y $Z))
           ((MLIST SIMP
             (32.
              #A((54.) BASE-CHAR
                 . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 4.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $X $Y)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $X $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $Y)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Y)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 4.) $VX $X $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Y)
              ((MTIMES SIMP) $VX $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Y $Z))))
           ((MLIST SIMP
             (32.
              #A((54.) BASE-CHAR
                 . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 4.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $X $Y)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $X $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $Y)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Y)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 4.) $VX $X $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Y)
              ((MTIMES SIMP) $VX $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y $Z)
              ((MTIMES SIMP) $VX $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $X $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $X $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y $Z)
              ((MTIMES SIMP) $VX $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $X $Y $Z))))
           ((MLIST SIMP
             (32.
              #A((54.) BASE-CHAR
                 . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
              SRC |$gsOrthoNorm| 30.))
            ((RAT SIMP) 1. 4.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Y)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $X $Y)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $X $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $X)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $Y)
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.) $VX $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $X $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Y)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $X $Z)
            ((MTIMES SIMP) ((RAT SIMP) 1. 4.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.)) $VX $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Z)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 4.) $VX $X $Y $Z)
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
              ((MTIMES SIMP) $VX $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 45. 16.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 45. 16.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 45. 16.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 45. 16.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 45. 16.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Y 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 45. 16.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Z 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $Y $Z)
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 32.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $X 2.)))
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $X 4.)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 32.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Y 2.)))
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Y 4.)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 32.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $Z 2.)))
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $Z 4.)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 32.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (62.
                   #A((54.) BASE-CHAR
                      . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                   SRC |$calcPowers| 62.))
                 $VX 2.)))
              ((MEXPT SIMP
                (62.
                 #A((54.) BASE-CHAR
                    . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                 SRC |$calcPowers| 62.))
               $VX 4.)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Y)
              ((MTIMES SIMP) $VX $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 9. 8.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $Y)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y)))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $X)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X $Y $Z)
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y $Z)
              ((MTIMES SIMP) $VX $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $X $Y)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $X $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y)))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP) $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 4.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $X)))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 4.)
               $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $Y)))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 135. 16.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 135. 16.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $VX $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y)))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
                ((MTIMES SIMP) $VX $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 135. 16.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
                ((MTIMES SIMP) $VX $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
                ((MTIMES SIMP) $VX $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 135. 16.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 2.)
               $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $Y $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 135. 16.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 2.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $X $Z)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 135. 16.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 9.) $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $X $Y)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 2.)
               $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 2.))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 3.)
               $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 3.)
               $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y $Z)
              ((MTIMES SIMP) $VX $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 15. 8.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 3.)
               $X $Y $Z)))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Z)
              ((MTIMES SIMP) $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                ((MTIMES SIMP) $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                ((MTIMES SIMP) $X $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $Y)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y)))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
                ((MTIMES SIMP) $VX $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.))))
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Z)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X)
                ((MTIMES SIMP) $VX $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y)
                ((MTIMES SIMP) $VX $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $VX $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Y)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 4.)
               $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $X $Y)))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 4.)
               $X $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $X $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 315. 32.)
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 4.)
               $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $Y $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $Y $Z)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $X 4.)
               $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $Y $Z)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $X 2.)
                 $Y $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $X $Z)
              ((MTIMES SIMP) $VX $X
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Y 4.)
               $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Z)
                ((MTIMES SIMP) $VX $X
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Y 2.)
                 $Z)))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX $X $Y)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX $X $Y)
                ((MTIMES SIMP) $VX $X $Y
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $Z 2.))))
              ((MTIMES SIMP) $VX $X $Y
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $Z 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 32.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 5. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X $Y $Z)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (62.
                  #A((54.) BASE-CHAR
                     . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                  SRC |$calcPowers| 62.))
                $VX 4.)
               $X $Y $Z)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X $Y $Z)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (62.
                    #A((54.) BASE-CHAR
                       . "/Users/nmandell/Codes/gkyl/cas-scripts/modal-basis.mac")
                    SRC |$calcPowers| 62.))
                  $VX 2.)
                 $X $Y $Z)))))))) 
(ADD2LNC '|$basisP| $VALUES) 