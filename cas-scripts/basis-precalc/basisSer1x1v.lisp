;;; -*- Mode: LISP; package:maxima; syntax:common-lisp; -*- 
(in-package :maxima)
(DSKSETQ |$varsC| '((MLIST SIMP) $X)) 
(ADD2LNC '|$varsC| $VALUES) 
(DSKSETQ |$varsP| '((MLIST SIMP) $X $VX)) 
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
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 2.))))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 3.))))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
            ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
            ((MTIMES SIMP) ((MEXPT SIMP) 2. ((RAT SIMP) -1. 2.))
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) 3. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 2.)))
            ((MTIMES SIMP) 5. ((MEXPT SIMP) 2. ((RAT SIMP) -3. 2.))
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 3.)))
            ((MTIMES SIMP) 105. ((MEXPT SIMP) 2. ((RAT SIMP) -7. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (61.
                   "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                   SRC |$calcPowers| 61.))
                 $X 2.)))
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 4.)))))) 
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
            ((RAT SIMP) 1. 2.)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $X)
            ((MTIMES SIMP) ((RAT SIMP) 1. 2.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 1. 2.)) $VX)
            ((MTIMES SIMP) ((RAT SIMP) 3. 2.) $VX $X))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
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
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $X))))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
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
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 3.)
               $X))))
           ((MLIST SIMP
             (31.
              "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
              SRC |$gsOrthoNorm| 29.))
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
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 5. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $VX 2.)))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 3. 4.)
             ((MEXPT SIMP) 15. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $X)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 7. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX)
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $VX 3.)))
            ((MTIMES SIMP) ((RAT SIMP) 45. 8.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 9.)
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (61.
                   "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                   SRC |$calcPowers| 61.))
                 $VX 2.)))
              ((MTIMES SIMP) ((RAT SIMP) -1. 3.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (61.
                   "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                   SRC |$calcPowers| 61.))
                 $X 2.)))
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 2.)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 2.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 3.))))
            ((MTIMES SIMP) ((RAT SIMP) 5. 4.)
             ((MEXPT SIMP) 21. ((RAT SIMP) 1. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -3. 5.) $VX $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 3.)
               $X)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (61.
                   "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                   SRC |$calcPowers| 61.))
                 $X 2.)))
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $X 4.)))
            ((MTIMES SIMP) ((RAT SIMP) 105. 16.)
             ((MPLUS SIMP) ((RAT SIMP) -1. 5.)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((RAT SIMP) -1. 3.)
                ((MEXPT SIMP
                  (61.
                   "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                   SRC |$calcPowers| 61.))
                 $VX 2.)))
              ((MEXPT SIMP
                (61.
                 "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                 SRC |$calcPowers| 61.))
               $VX 4.)))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $VX)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $VX)
                ((MTIMES SIMP) $VX
                 ((MEXPT SIMP
                   (61.
                    "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                    SRC |$calcPowers| 61.))
                  $X 2.))))
              ((MTIMES SIMP) $VX
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $X 4.))))
            ((MTIMES SIMP) ((RAT SIMP) 35. 16.)
             ((MEXPT SIMP) 3. ((RAT SIMP) 3. 2.))
             ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 5.) $X)
              ((MTIMES SIMP)
               ((MEXPT SIMP
                 (61.
                  "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                  SRC |$calcPowers| 61.))
                $VX 4.)
               $X)
              ((MTIMES SIMP) ((RAT SIMP) -6. 7.)
               ((MPLUS SIMP) ((MTIMES SIMP) ((RAT SIMP) -1. 3.) $X)
                ((MTIMES SIMP)
                 ((MEXPT SIMP
                   (61.
                    "/Users/ahakim/research/gkyl-project/gkyl/cas-scripts/modal-basis.mac"
                    SRC |$calcPowers| 61.))
                  $VX 2.)
                 $X)))))))) 
(ADD2LNC '|$basisP| $VALUES) 