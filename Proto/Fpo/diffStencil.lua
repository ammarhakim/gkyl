function diffStencilFunc(dt, dx, dy,
			 fPtr, fLPtr, fRPtr, fTPtr, fBPtr,
			 gPtr, gLPtr, gRPtr, gTPtr, gBPtr,
			 isTopEdge, isBotEdge, isLeftEdge, isRightEdge,
			 fOutPtr) 

if isLeftEdge then
    fOutPtr[1] = 0.5*dt*((-(1.40625*fRPtr[4]*gRPtr[4]*dy^2)/dx^3)+(0.46875*fPtr[4]*gRPtr[4]*dy^2)/dx^3+(1.299038105676658*fRPtr[3]*gRPtr[4]*dy^2)/dx^3+(0.3247595264191644*fPtr[3]*gRPtr[4]*dy^2)/dx^3-(0.46875*fRPtr[4]*gPtr[4]*dy^2)/dx^3+(1.40625*fPtr[4]*gPtr[4]*dy^2)/dx^3+(0.3247595264191644*fRPtr[3]*gPtr[4]*dy^2)/dx^3+(1.299038105676658*fPtr[3]*gPtr[4]*dy^2)/dx^3+(0.5412658773652741*gRPtr[3]*fRPtr[4]*dy^2)/dx^3-(0.5412658773652741*gPtr[3]*fRPtr[4]*dy^2)/dx^3-(0.5412658773652741*gRPtr[3]*fPtr[4]*dy^2)/dx^3+(0.5412658773652741*gPtr[3]*fPtr[4]*dy^2)/dx^3-(0.46875*fRPtr[3]*gRPtr[3]*dy^2)/dx^3-(0.46875*fPtr[3]*gRPtr[3]*dy^2)/dx^3+(0.46875*fRPtr[3]*gPtr[3]*dy^2)/dx^3+(0.46875*fPtr[3]*gPtr[3]*dy^2)/dx^3-(5.625*fRPtr[2]*gRPtr[2])/dx^3+(1.875*fPtr[2]*gRPtr[2])/dx^3+(5.196152422706631*fRPtr[1]*gRPtr[2])/dx^3+(1.299038105676658*fPtr[1]*gRPtr[2])/dx^3-(1.875*fRPtr[2]*gPtr[2])/dx^3+(5.625*fPtr[2]*gPtr[2])/dx^3+(1.299038105676658*fRPtr[1]*gPtr[2])/dx^3+(5.196152422706631*fPtr[1]*gPtr[2])/dx^3+(2.165063509461096*gRPtr[1]*fRPtr[2])/dx^3-(2.165063509461096*gPtr[1]*fRPtr[2])/dx^3-(2.165063509461096*gRPtr[1]*fPtr[2])/dx^3+(2.165063509461096*gPtr[1]*fPtr[2])/dx^3-(1.875*fRPtr[1]*gRPtr[1])/dx^3-(1.875*fPtr[1]*gRPtr[1])/dx^3+(1.875*fRPtr[1]*gPtr[1])/dx^3+(1.875*fPtr[1]*gPtr[1])/dx^3)+fOutPtr[1]
    fOutPtr[2] = 0.5*dt*((-(2.002683746251514*fRPtr[4]*gRPtr[4]*dy^2)/dx^3)+(0.3788861141556918*fPtr[4]*gRPtr[4]*dy^2)/dx^3+(1.875*fRPtr[3]*gRPtr[4]*dy^2)/dx^3+(0.1875*fPtr[3]*gRPtr[4]*dy^2)/dx^3-(1.24491151794013*fRPtr[4]*gPtr[4]*dy^2)/dx^3+(2.868709150035952*fPtr[4]*gPtr[4]*dy^2)/dx^3+(0.9375*fRPtr[3]*gPtr[4]*dy^2)/dx^3+(2.625*fPtr[3]*gPtr[4]*dy^2)/dx^3+(0.9375*gRPtr[3]*fRPtr[4]*dy^2)/dx^3-(0.9375*gPtr[3]*fRPtr[4]*dy^2)/dx^3-(0.9375*gRPtr[3]*fPtr[4]*dy^2)/dx^3+(0.9375*gPtr[3]*fPtr[4]*dy^2)/dx^3-(0.8118988160479111*fRPtr[3]*gRPtr[3]*dy^2)/dx^3-(0.8118988160479111*fPtr[3]*gRPtr[3]*dy^2)/dx^3+(0.8118988160479111*fRPtr[3]*gPtr[3]*dy^2)/dx^3+(0.8118988160479111*fPtr[3]*gPtr[3]*dy^2)/dx^3-(8.010734985006055*fRPtr[2]*gRPtr[2])/dx^3+(1.515544456622767*fPtr[2]*gRPtr[2])/dx^3+(7.5*fRPtr[1]*gRPtr[2])/dx^3+(0.75*fPtr[1]*gRPtr[2])/dx^3-(4.979646071760522*fRPtr[2]*gPtr[2])/dx^3+(11.47483660014381*fPtr[2]*gPtr[2])/dx^3+(3.75*fRPtr[1]*gPtr[2])/dx^3+(10.5*fPtr[1]*gPtr[2])/dx^3+(3.75*gRPtr[1]*fRPtr[2])/dx^3-(3.75*gPtr[1]*fRPtr[2])/dx^3-(3.75*gRPtr[1]*fPtr[2])/dx^3+(3.75*gPtr[1]*fPtr[2])/dx^3-(3.247595264191645*fRPtr[1]*gRPtr[1])/dx^3-(3.247595264191645*fPtr[1]*gRPtr[1])/dx^3+(3.247595264191645*fRPtr[1]*gPtr[1])/dx^3+(3.247595264191645*fPtr[1]*gPtr[1])/dx^3)+fOutPtr[2]
    fOutPtr[3] = 0.5*dt*((-(2.8125*fRPtr[2]*gRPtr[4]*dy)/dx^3)+(0.9375*fPtr[2]*gRPtr[4]*dy)/dx^3+(2.598076211353316*fRPtr[1]*gRPtr[4]*dy)/dx^3+(0.6495190528383289*fPtr[1]*gRPtr[4]*dy)/dx^3-(0.9375*fRPtr[2]*gPtr[4]*dy)/dx^3+(2.8125*fPtr[2]*gPtr[4]*dy)/dx^3+(0.6495190528383289*fRPtr[1]*gPtr[4]*dy)/dx^3+(2.598076211353316*fPtr[1]*gPtr[4]*dy)/dx^3-(2.8125*gRPtr[2]*fRPtr[4]*dy)/dx^3-(0.9375*gPtr[2]*fRPtr[4]*dy)/dx^3+(1.082531754730548*gRPtr[1]*fRPtr[4]*dy)/dx^3-(1.082531754730548*gPtr[1]*fRPtr[4]*dy)/dx^3+(0.9375*gRPtr[2]*fPtr[4]*dy)/dx^3+(2.8125*gPtr[2]*fPtr[4]*dy)/dx^3-(1.082531754730548*gRPtr[1]*fPtr[4]*dy)/dx^3+(1.082531754730548*gPtr[1]*fPtr[4]*dy)/dx^3+(1.082531754730548*fRPtr[2]*gRPtr[3]*dy)/dx^3-(1.082531754730548*fPtr[2]*gRPtr[3]*dy)/dx^3-(0.9375*fRPtr[1]*gRPtr[3]*dy)/dx^3-(0.9375*fPtr[1]*gRPtr[3]*dy)/dx^3-(1.082531754730548*fRPtr[2]*gPtr[3]*dy)/dx^3+(1.082531754730548*fPtr[2]*gPtr[3]*dy)/dx^3+(0.9375*fRPtr[1]*gPtr[3]*dy)/dx^3+(0.9375*fPtr[1]*gPtr[3]*dy)/dx^3+(2.598076211353316*gRPtr[2]*fRPtr[3]*dy)/dx^3+(0.6495190528383289*gPtr[2]*fRPtr[3]*dy)/dx^3-(0.9375*gRPtr[1]*fRPtr[3]*dy)/dx^3+(0.9375*gPtr[1]*fRPtr[3]*dy)/dx^3+(0.6495190528383289*gRPtr[2]*fPtr[3]*dy)/dx^3+(2.598076211353316*gPtr[2]*fPtr[3]*dy)/dx^3-(0.9375*gRPtr[1]*fPtr[3]*dy)/dx^3+(0.9375*gPtr[1]*fPtr[3]*dy)/dx^3)+fOutPtr[3]
    fOutPtr[4] = 0.5*dt*((-(4.005367492503027*fRPtr[2]*gRPtr[4]*dy)/dx^3)+(0.7577722283113837*fPtr[2]*gRPtr[4]*dy)/dx^3+(3.75*fRPtr[1]*gRPtr[4]*dy)/dx^3+(0.375*fPtr[1]*gRPtr[4]*dy)/dx^3-(2.489823035880261*fRPtr[2]*gPtr[4]*dy)/dx^3+(5.737418300071904*fPtr[2]*gPtr[4]*dy)/dx^3+(1.875*fRPtr[1]*gPtr[4]*dy)/dx^3+(5.25*fPtr[1]*gPtr[4]*dy)/dx^3-(4.005367492503027*gRPtr[2]*fRPtr[4]*dy)/dx^3-(2.489823035880261*gPtr[2]*fRPtr[4]*dy)/dx^3+(1.875*gRPtr[1]*fRPtr[4]*dy)/dx^3-(1.875*gPtr[1]*fRPtr[4]*dy)/dx^3+(0.7577722283113837*gRPtr[2]*fPtr[4]*dy)/dx^3+(5.737418300071904*gPtr[2]*fPtr[4]*dy)/dx^3-(1.875*gRPtr[1]*fPtr[4]*dy)/dx^3+(1.875*gPtr[1]*fPtr[4]*dy)/dx^3+(1.875*fRPtr[2]*gRPtr[3]*dy)/dx^3-(1.875*fPtr[2]*gRPtr[3]*dy)/dx^3-(1.623797632095822*fRPtr[1]*gRPtr[3]*dy)/dx^3-(1.623797632095822*fPtr[1]*gRPtr[3]*dy)/dx^3-(1.875*fRPtr[2]*gPtr[3]*dy)/dx^3+(1.875*fPtr[2]*gPtr[3]*dy)/dx^3+(1.623797632095822*fRPtr[1]*gPtr[3]*dy)/dx^3+(1.623797632095822*fPtr[1]*gPtr[3]*dy)/dx^3+(3.75*gRPtr[2]*fRPtr[3]*dy)/dx^3+(1.875*gPtr[2]*fRPtr[3]*dy)/dx^3-(1.623797632095822*gRPtr[1]*fRPtr[3]*dy)/dx^3+(1.623797632095822*gPtr[1]*fRPtr[3]*dy)/dx^3+(0.375*gRPtr[2]*fPtr[3]*dy)/dx^3+(5.25*gPtr[2]*fPtr[3]*dy)/dx^3-(1.623797632095822*gRPtr[1]*fPtr[3]*dy)/dx^3+(1.623797632095822*gPtr[1]*fPtr[3]*dy)/dx^3)+fOutPtr[4]
elseif isRightEdge then
    fOutPtr[1] = 0.5*dt*((1.40625*fPtr[4]*gPtr[4]*dy^2)/dx^3-(0.46875*fLPtr[4]*gPtr[4]*dy^2)/dx^3-(1.299038105676658*fPtr[3]*gPtr[4]*dy^2)/dx^3-(0.3247595264191644*fLPtr[3]*gPtr[4]*dy^2)/dx^3+(0.46875*fPtr[4]*gLPtr[4]*dy^2)/dx^3-(1.40625*fLPtr[4]*gLPtr[4]*dy^2)/dx^3-(0.3247595264191644*fPtr[3]*gLPtr[4]*dy^2)/dx^3-(1.299038105676658*fLPtr[3]*gLPtr[4]*dy^2)/dx^3-(0.5412658773652741*gPtr[3]*fPtr[4]*dy^2)/dx^3+(0.5412658773652741*gLPtr[3]*fPtr[4]*dy^2)/dx^3+(0.5412658773652741*gPtr[3]*fLPtr[4]*dy^2)/dx^3-(0.5412658773652741*gLPtr[3]*fLPtr[4]*dy^2)/dx^3+(0.46875*fPtr[3]*gPtr[3]*dy^2)/dx^3+(0.46875*fLPtr[3]*gPtr[3]*dy^2)/dx^3-(0.46875*fPtr[3]*gLPtr[3]*dy^2)/dx^3-(0.46875*fLPtr[3]*gLPtr[3]*dy^2)/dx^3+(5.625*fPtr[2]*gPtr[2])/dx^3-(1.875*fLPtr[2]*gPtr[2])/dx^3-(5.196152422706631*fPtr[1]*gPtr[2])/dx^3-(1.299038105676658*fLPtr[1]*gPtr[2])/dx^3+(1.875*fPtr[2]*gLPtr[2])/dx^3-(5.625*fLPtr[2]*gLPtr[2])/dx^3-(1.299038105676658*fPtr[1]*gLPtr[2])/dx^3-(5.196152422706631*fLPtr[1]*gLPtr[2])/dx^3-(2.165063509461096*gPtr[1]*fPtr[2])/dx^3+(2.165063509461096*gLPtr[1]*fPtr[2])/dx^3+(2.165063509461096*gPtr[1]*fLPtr[2])/dx^3-(2.165063509461096*gLPtr[1]*fLPtr[2])/dx^3+(1.875*fPtr[1]*gPtr[1])/dx^3+(1.875*fLPtr[1]*gPtr[1])/dx^3-(1.875*fPtr[1]*gLPtr[1])/dx^3-(1.875*fLPtr[1]*gLPtr[1])/dx^3)+fOutPtr[1]
    fOutPtr[2] = 0.5*dt*((-(2.868709150035952*fPtr[4]*gPtr[4]*dy^2)/dx^3)+(1.24491151794013*fLPtr[4]*gPtr[4]*dy^2)/dx^3+(2.625*fPtr[3]*gPtr[4]*dy^2)/dx^3+(0.9375*fLPtr[3]*gPtr[4]*dy^2)/dx^3-(0.3788861141556918*fPtr[4]*gLPtr[4]*dy^2)/dx^3+(2.002683746251514*fLPtr[4]*gLPtr[4]*dy^2)/dx^3+(0.1875*fPtr[3]*gLPtr[4]*dy^2)/dx^3+(1.875*fLPtr[3]*gLPtr[4]*dy^2)/dx^3+(0.9375*gPtr[3]*fPtr[4]*dy^2)/dx^3-(0.9375*gLPtr[3]*fPtr[4]*dy^2)/dx^3-(0.9375*gPtr[3]*fLPtr[4]*dy^2)/dx^3+(0.9375*gLPtr[3]*fLPtr[4]*dy^2)/dx^3-(0.8118988160479111*fPtr[3]*gPtr[3]*dy^2)/dx^3-(0.8118988160479111*fLPtr[3]*gPtr[3]*dy^2)/dx^3+(0.8118988160479111*fPtr[3]*gLPtr[3]*dy^2)/dx^3+(0.8118988160479111*fLPtr[3]*gLPtr[3]*dy^2)/dx^3-(11.47483660014381*fPtr[2]*gPtr[2])/dx^3+(4.979646071760522*fLPtr[2]*gPtr[2])/dx^3+(10.5*fPtr[1]*gPtr[2])/dx^3+(3.75*fLPtr[1]*gPtr[2])/dx^3-(1.515544456622767*fPtr[2]*gLPtr[2])/dx^3+(8.010734985006055*fLPtr[2]*gLPtr[2])/dx^3+(0.75*fPtr[1]*gLPtr[2])/dx^3+(7.5*fLPtr[1]*gLPtr[2])/dx^3+(3.75*gPtr[1]*fPtr[2])/dx^3-(3.75*gLPtr[1]*fPtr[2])/dx^3-(3.75*gPtr[1]*fLPtr[2])/dx^3+(3.75*gLPtr[1]*fLPtr[2])/dx^3-(3.247595264191645*fPtr[1]*gPtr[1])/dx^3-(3.247595264191645*fLPtr[1]*gPtr[1])/dx^3+(3.247595264191645*fPtr[1]*gLPtr[1])/dx^3+(3.247595264191645*fLPtr[1]*gLPtr[1])/dx^3)+fOutPtr[2]
    fOutPtr[3] = 0.5*dt*((2.8125*fPtr[2]*gPtr[4]*dy)/dx^3-(0.9375*fLPtr[2]*gPtr[4]*dy)/dx^3-(2.598076211353316*fPtr[1]*gPtr[4]*dy)/dx^3-(0.6495190528383289*fLPtr[1]*gPtr[4]*dy)/dx^3+(0.9375*fPtr[2]*gLPtr[4]*dy)/dx^3-(2.8125*fLPtr[2]*gLPtr[4]*dy)/dx^3-(0.6495190528383289*fPtr[1]*gLPtr[4]*dy)/dx^3-(2.598076211353316*fLPtr[1]*gLPtr[4]*dy)/dx^3+(2.8125*gPtr[2]*fPtr[4]*dy)/dx^3+(0.9375*gLPtr[2]*fPtr[4]*dy)/dx^3-(1.082531754730548*gPtr[1]*fPtr[4]*dy)/dx^3+(1.082531754730548*gLPtr[1]*fPtr[4]*dy)/dx^3-(0.9375*gPtr[2]*fLPtr[4]*dy)/dx^3-(2.8125*gLPtr[2]*fLPtr[4]*dy)/dx^3+(1.082531754730548*gPtr[1]*fLPtr[4]*dy)/dx^3-(1.082531754730548*gLPtr[1]*fLPtr[4]*dy)/dx^3-(1.082531754730548*fPtr[2]*gPtr[3]*dy)/dx^3+(1.082531754730548*fLPtr[2]*gPtr[3]*dy)/dx^3+(0.9375*fPtr[1]*gPtr[3]*dy)/dx^3+(0.9375*fLPtr[1]*gPtr[3]*dy)/dx^3+(1.082531754730548*fPtr[2]*gLPtr[3]*dy)/dx^3-(1.082531754730548*fLPtr[2]*gLPtr[3]*dy)/dx^3-(0.9375*fPtr[1]*gLPtr[3]*dy)/dx^3-(0.9375*fLPtr[1]*gLPtr[3]*dy)/dx^3-(2.598076211353316*gPtr[2]*fPtr[3]*dy)/dx^3-(0.6495190528383289*gLPtr[2]*fPtr[3]*dy)/dx^3+(0.9375*gPtr[1]*fPtr[3]*dy)/dx^3-(0.9375*gLPtr[1]*fPtr[3]*dy)/dx^3-(0.6495190528383289*gPtr[2]*fLPtr[3]*dy)/dx^3-(2.598076211353316*gLPtr[2]*fLPtr[3]*dy)/dx^3+(0.9375*gPtr[1]*fLPtr[3]*dy)/dx^3-(0.9375*gLPtr[1]*fLPtr[3]*dy)/dx^3)+fOutPtr[3]
    fOutPtr[4] = 0.5*dt*((-(5.737418300071904*fPtr[2]*gPtr[4]*dy)/dx^3)+(2.489823035880261*fLPtr[2]*gPtr[4]*dy)/dx^3+(5.25*fPtr[1]*gPtr[4]*dy)/dx^3+(1.875*fLPtr[1]*gPtr[4]*dy)/dx^3-(0.7577722283113837*fPtr[2]*gLPtr[4]*dy)/dx^3+(4.005367492503027*fLPtr[2]*gLPtr[4]*dy)/dx^3+(0.375*fPtr[1]*gLPtr[4]*dy)/dx^3+(3.75*fLPtr[1]*gLPtr[4]*dy)/dx^3-(5.737418300071904*gPtr[2]*fPtr[4]*dy)/dx^3-(0.7577722283113837*gLPtr[2]*fPtr[4]*dy)/dx^3+(1.875*gPtr[1]*fPtr[4]*dy)/dx^3-(1.875*gLPtr[1]*fPtr[4]*dy)/dx^3+(2.489823035880261*gPtr[2]*fLPtr[4]*dy)/dx^3+(4.005367492503027*gLPtr[2]*fLPtr[4]*dy)/dx^3-(1.875*gPtr[1]*fLPtr[4]*dy)/dx^3+(1.875*gLPtr[1]*fLPtr[4]*dy)/dx^3+(1.875*fPtr[2]*gPtr[3]*dy)/dx^3-(1.875*fLPtr[2]*gPtr[3]*dy)/dx^3-(1.623797632095822*fPtr[1]*gPtr[3]*dy)/dx^3-(1.623797632095822*fLPtr[1]*gPtr[3]*dy)/dx^3-(1.875*fPtr[2]*gLPtr[3]*dy)/dx^3+(1.875*fLPtr[2]*gLPtr[3]*dy)/dx^3+(1.623797632095822*fPtr[1]*gLPtr[3]*dy)/dx^3+(1.623797632095822*fLPtr[1]*gLPtr[3]*dy)/dx^3+(5.25*gPtr[2]*fPtr[3]*dy)/dx^3+(0.375*gLPtr[2]*fPtr[3]*dy)/dx^3-(1.623797632095822*gPtr[1]*fPtr[3]*dy)/dx^3+(1.623797632095822*gLPtr[1]*fPtr[3]*dy)/dx^3+(1.875*gPtr[2]*fLPtr[3]*dy)/dx^3+(3.75*gLPtr[2]*fLPtr[3]*dy)/dx^3-(1.623797632095822*gPtr[1]*fLPtr[3]*dy)/dx^3+(1.623797632095822*gLPtr[1]*fLPtr[3]*dy)/dx^3)+fOutPtr[4]
else
    fOutPtr[1] = 0.5*dt*((-(1.40625*fRPtr[4]*gRPtr[4]*dy^2)/dx^3)+(0.46875*fPtr[4]*gRPtr[4]*dy^2)/dx^3+(1.299038105676658*fRPtr[3]*gRPtr[4]*dy^2)/dx^3+(0.3247595264191644*fPtr[3]*gRPtr[4]*dy^2)/dx^3-(0.46875*fRPtr[4]*gPtr[4]*dy^2)/dx^3+(2.8125*fPtr[4]*gPtr[4]*dy^2)/dx^3-(0.46875*fLPtr[4]*gPtr[4]*dy^2)/dx^3+(0.3247595264191644*fRPtr[3]*gPtr[4]*dy^2)/dx^3-(0.3247595264191644*fLPtr[3]*gPtr[4]*dy^2)/dx^3+(0.46875*fPtr[4]*gLPtr[4]*dy^2)/dx^3-(1.40625*fLPtr[4]*gLPtr[4]*dy^2)/dx^3-(0.3247595264191644*fPtr[3]*gLPtr[4]*dy^2)/dx^3-(1.299038105676658*fLPtr[3]*gLPtr[4]*dy^2)/dx^3+(0.5412658773652741*gRPtr[3]*fRPtr[4]*dy^2)/dx^3-(0.5412658773652741*gPtr[3]*fRPtr[4]*dy^2)/dx^3-(0.5412658773652741*gRPtr[3]*fPtr[4]*dy^2)/dx^3+(0.5412658773652741*gLPtr[3]*fPtr[4]*dy^2)/dx^3+(0.5412658773652741*gPtr[3]*fLPtr[4]*dy^2)/dx^3-(0.5412658773652741*gLPtr[3]*fLPtr[4]*dy^2)/dx^3-(0.46875*fRPtr[3]*gRPtr[3]*dy^2)/dx^3-(0.46875*fPtr[3]*gRPtr[3]*dy^2)/dx^3+(0.46875*fRPtr[3]*gPtr[3]*dy^2)/dx^3+(0.9375*fPtr[3]*gPtr[3]*dy^2)/dx^3+(0.46875*fLPtr[3]*gPtr[3]*dy^2)/dx^3-(0.46875*fPtr[3]*gLPtr[3]*dy^2)/dx^3-(0.46875*fLPtr[3]*gLPtr[3]*dy^2)/dx^3-(5.625*fRPtr[2]*gRPtr[2])/dx^3+(1.875*fPtr[2]*gRPtr[2])/dx^3+(5.196152422706631*fRPtr[1]*gRPtr[2])/dx^3+(1.299038105676658*fPtr[1]*gRPtr[2])/dx^3-(1.875*fRPtr[2]*gPtr[2])/dx^3+(11.25*fPtr[2]*gPtr[2])/dx^3-(1.875*fLPtr[2]*gPtr[2])/dx^3+(1.299038105676658*fRPtr[1]*gPtr[2])/dx^3-(1.299038105676658*fLPtr[1]*gPtr[2])/dx^3+(1.875*fPtr[2]*gLPtr[2])/dx^3-(5.625*fLPtr[2]*gLPtr[2])/dx^3-(1.299038105676658*fPtr[1]*gLPtr[2])/dx^3-(5.196152422706631*fLPtr[1]*gLPtr[2])/dx^3+(2.165063509461096*gRPtr[1]*fRPtr[2])/dx^3-(2.165063509461096*gPtr[1]*fRPtr[2])/dx^3-(2.165063509461096*gRPtr[1]*fPtr[2])/dx^3+(2.165063509461096*gLPtr[1]*fPtr[2])/dx^3+(2.165063509461096*gPtr[1]*fLPtr[2])/dx^3-(2.165063509461096*gLPtr[1]*fLPtr[2])/dx^3-(1.875*fRPtr[1]*gRPtr[1])/dx^3-(1.875*fPtr[1]*gRPtr[1])/dx^3+(1.875*fRPtr[1]*gPtr[1])/dx^3+(3.75*fPtr[1]*gPtr[1])/dx^3+(1.875*fLPtr[1]*gPtr[1])/dx^3-(1.875*fPtr[1]*gLPtr[1])/dx^3-(1.875*fLPtr[1]*gLPtr[1])/dx^3)+fOutPtr[1]
    fOutPtr[2] = 0.5*dt*((-(2.002683746251514*fRPtr[4]*gRPtr[4]*dy^2)/dx^3)+(0.3788861141556918*fPtr[4]*gRPtr[4]*dy^2)/dx^3+(1.875*fRPtr[3]*gRPtr[4]*dy^2)/dx^3+(0.1875*fPtr[3]*gRPtr[4]*dy^2)/dx^3-(1.24491151794013*fRPtr[4]*gPtr[4]*dy^2)/dx^3+(1.24491151794013*fLPtr[4]*gPtr[4]*dy^2)/dx^3+(0.9375*fRPtr[3]*gPtr[4]*dy^2)/dx^3+(5.25*fPtr[3]*gPtr[4]*dy^2)/dx^3+(0.9375*fLPtr[3]*gPtr[4]*dy^2)/dx^3-(0.3788861141556918*fPtr[4]*gLPtr[4]*dy^2)/dx^3+(2.002683746251514*fLPtr[4]*gLPtr[4]*dy^2)/dx^3+(0.1875*fPtr[3]*gLPtr[4]*dy^2)/dx^3+(1.875*fLPtr[3]*gLPtr[4]*dy^2)/dx^3+(0.9375*gRPtr[3]*fRPtr[4]*dy^2)/dx^3-(0.9375*gPtr[3]*fRPtr[4]*dy^2)/dx^3-(0.9375*gRPtr[3]*fPtr[4]*dy^2)/dx^3+(1.875*gPtr[3]*fPtr[4]*dy^2)/dx^3-(0.9375*gLPtr[3]*fPtr[4]*dy^2)/dx^3-(0.9375*gPtr[3]*fLPtr[4]*dy^2)/dx^3+(0.9375*gLPtr[3]*fLPtr[4]*dy^2)/dx^3-(0.8118988160479111*fRPtr[3]*gRPtr[3]*dy^2)/dx^3-(0.8118988160479111*fPtr[3]*gRPtr[3]*dy^2)/dx^3+(0.8118988160479111*fRPtr[3]*gPtr[3]*dy^2)/dx^3-(0.8118988160479111*fLPtr[3]*gPtr[3]*dy^2)/dx^3+(0.8118988160479111*fPtr[3]*gLPtr[3]*dy^2)/dx^3+(0.8118988160479111*fLPtr[3]*gLPtr[3]*dy^2)/dx^3-(8.010734985006055*fRPtr[2]*gRPtr[2])/dx^3+(1.515544456622767*fPtr[2]*gRPtr[2])/dx^3+(7.5*fRPtr[1]*gRPtr[2])/dx^3+(0.75*fPtr[1]*gRPtr[2])/dx^3-(4.979646071760522*fRPtr[2]*gPtr[2])/dx^3+(4.979646071760522*fLPtr[2]*gPtr[2])/dx^3+(3.75*fRPtr[1]*gPtr[2])/dx^3+(21.0*fPtr[1]*gPtr[2])/dx^3+(3.75*fLPtr[1]*gPtr[2])/dx^3-(1.515544456622767*fPtr[2]*gLPtr[2])/dx^3+(8.010734985006055*fLPtr[2]*gLPtr[2])/dx^3+(0.75*fPtr[1]*gLPtr[2])/dx^3+(7.5*fLPtr[1]*gLPtr[2])/dx^3+(3.75*gRPtr[1]*fRPtr[2])/dx^3-(3.75*gPtr[1]*fRPtr[2])/dx^3-(3.75*gRPtr[1]*fPtr[2])/dx^3+(7.5*gPtr[1]*fPtr[2])/dx^3-(3.75*gLPtr[1]*fPtr[2])/dx^3-(3.75*gPtr[1]*fLPtr[2])/dx^3+(3.75*gLPtr[1]*fLPtr[2])/dx^3-(3.247595264191645*fRPtr[1]*gRPtr[1])/dx^3-(3.247595264191645*fPtr[1]*gRPtr[1])/dx^3+(3.247595264191645*fRPtr[1]*gPtr[1])/dx^3-(3.247595264191645*fLPtr[1]*gPtr[1])/dx^3+(3.247595264191645*fPtr[1]*gLPtr[1])/dx^3+(3.247595264191645*fLPtr[1]*gLPtr[1])/dx^3)+fOutPtr[2]
    fOutPtr[3] = 0.5*dt*((-(2.8125*fRPtr[2]*gRPtr[4]*dy)/dx^3)+(0.9375*fPtr[2]*gRPtr[4]*dy)/dx^3+(2.598076211353316*fRPtr[1]*gRPtr[4]*dy)/dx^3+(0.6495190528383289*fPtr[1]*gRPtr[4]*dy)/dx^3-(0.9375*fRPtr[2]*gPtr[4]*dy)/dx^3+(5.625*fPtr[2]*gPtr[4]*dy)/dx^3-(0.9375*fLPtr[2]*gPtr[4]*dy)/dx^3+(0.6495190528383289*fRPtr[1]*gPtr[4]*dy)/dx^3-(0.6495190528383289*fLPtr[1]*gPtr[4]*dy)/dx^3+(0.9375*fPtr[2]*gLPtr[4]*dy)/dx^3-(2.8125*fLPtr[2]*gLPtr[4]*dy)/dx^3-(0.6495190528383289*fPtr[1]*gLPtr[4]*dy)/dx^3-(2.598076211353316*fLPtr[1]*gLPtr[4]*dy)/dx^3-(2.8125*gRPtr[2]*fRPtr[4]*dy)/dx^3-(0.9375*gPtr[2]*fRPtr[4]*dy)/dx^3+(1.082531754730548*gRPtr[1]*fRPtr[4]*dy)/dx^3-(1.082531754730548*gPtr[1]*fRPtr[4]*dy)/dx^3+(0.9375*gRPtr[2]*fPtr[4]*dy)/dx^3+(5.625*gPtr[2]*fPtr[4]*dy)/dx^3+(0.9375*gLPtr[2]*fPtr[4]*dy)/dx^3-(1.082531754730548*gRPtr[1]*fPtr[4]*dy)/dx^3+(1.082531754730548*gLPtr[1]*fPtr[4]*dy)/dx^3-(0.9375*gPtr[2]*fLPtr[4]*dy)/dx^3-(2.8125*gLPtr[2]*fLPtr[4]*dy)/dx^3+(1.082531754730548*gPtr[1]*fLPtr[4]*dy)/dx^3-(1.082531754730548*gLPtr[1]*fLPtr[4]*dy)/dx^3+(1.082531754730548*fRPtr[2]*gRPtr[3]*dy)/dx^3-(1.082531754730548*fPtr[2]*gRPtr[3]*dy)/dx^3-(0.9375*fRPtr[1]*gRPtr[3]*dy)/dx^3-(0.9375*fPtr[1]*gRPtr[3]*dy)/dx^3-(1.082531754730548*fRPtr[2]*gPtr[3]*dy)/dx^3+(1.082531754730548*fLPtr[2]*gPtr[3]*dy)/dx^3+(0.9375*fRPtr[1]*gPtr[3]*dy)/dx^3+(1.875*fPtr[1]*gPtr[3]*dy)/dx^3+(0.9375*fLPtr[1]*gPtr[3]*dy)/dx^3+(1.082531754730548*fPtr[2]*gLPtr[3]*dy)/dx^3-(1.082531754730548*fLPtr[2]*gLPtr[3]*dy)/dx^3-(0.9375*fPtr[1]*gLPtr[3]*dy)/dx^3-(0.9375*fLPtr[1]*gLPtr[3]*dy)/dx^3+(2.598076211353316*gRPtr[2]*fRPtr[3]*dy)/dx^3+(0.6495190528383289*gPtr[2]*fRPtr[3]*dy)/dx^3-(0.9375*gRPtr[1]*fRPtr[3]*dy)/dx^3+(0.9375*gPtr[1]*fRPtr[3]*dy)/dx^3+(0.6495190528383289*gRPtr[2]*fPtr[3]*dy)/dx^3-(0.6495190528383289*gLPtr[2]*fPtr[3]*dy)/dx^3-(0.9375*gRPtr[1]*fPtr[3]*dy)/dx^3+(1.875*gPtr[1]*fPtr[3]*dy)/dx^3-(0.9375*gLPtr[1]*fPtr[3]*dy)/dx^3-(0.6495190528383289*gPtr[2]*fLPtr[3]*dy)/dx^3-(2.598076211353316*gLPtr[2]*fLPtr[3]*dy)/dx^3+(0.9375*gPtr[1]*fLPtr[3]*dy)/dx^3-(0.9375*gLPtr[1]*fLPtr[3]*dy)/dx^3)+fOutPtr[3]
    fOutPtr[4] = 0.5*dt*((-(4.005367492503027*fRPtr[2]*gRPtr[4]*dy)/dx^3)+(0.7577722283113837*fPtr[2]*gRPtr[4]*dy)/dx^3+(3.75*fRPtr[1]*gRPtr[4]*dy)/dx^3+(0.375*fPtr[1]*gRPtr[4]*dy)/dx^3-(2.489823035880261*fRPtr[2]*gPtr[4]*dy)/dx^3+(2.489823035880261*fLPtr[2]*gPtr[4]*dy)/dx^3+(1.875*fRPtr[1]*gPtr[4]*dy)/dx^3+(10.5*fPtr[1]*gPtr[4]*dy)/dx^3+(1.875*fLPtr[1]*gPtr[4]*dy)/dx^3-(0.7577722283113837*fPtr[2]*gLPtr[4]*dy)/dx^3+(4.005367492503027*fLPtr[2]*gLPtr[4]*dy)/dx^3+(0.375*fPtr[1]*gLPtr[4]*dy)/dx^3+(3.75*fLPtr[1]*gLPtr[4]*dy)/dx^3-(4.005367492503027*gRPtr[2]*fRPtr[4]*dy)/dx^3-(2.489823035880261*gPtr[2]*fRPtr[4]*dy)/dx^3+(1.875*gRPtr[1]*fRPtr[4]*dy)/dx^3-(1.875*gPtr[1]*fRPtr[4]*dy)/dx^3+(0.7577722283113837*gRPtr[2]*fPtr[4]*dy)/dx^3-(0.7577722283113837*gLPtr[2]*fPtr[4]*dy)/dx^3-(1.875*gRPtr[1]*fPtr[4]*dy)/dx^3+(3.75*gPtr[1]*fPtr[4]*dy)/dx^3-(1.875*gLPtr[1]*fPtr[4]*dy)/dx^3+(2.489823035880261*gPtr[2]*fLPtr[4]*dy)/dx^3+(4.005367492503027*gLPtr[2]*fLPtr[4]*dy)/dx^3-(1.875*gPtr[1]*fLPtr[4]*dy)/dx^3+(1.875*gLPtr[1]*fLPtr[4]*dy)/dx^3+(1.875*fRPtr[2]*gRPtr[3]*dy)/dx^3-(1.875*fPtr[2]*gRPtr[3]*dy)/dx^3-(1.623797632095822*fRPtr[1]*gRPtr[3]*dy)/dx^3-(1.623797632095822*fPtr[1]*gRPtr[3]*dy)/dx^3-(1.875*fRPtr[2]*gPtr[3]*dy)/dx^3+(3.75*fPtr[2]*gPtr[3]*dy)/dx^3-(1.875*fLPtr[2]*gPtr[3]*dy)/dx^3+(1.623797632095822*fRPtr[1]*gPtr[3]*dy)/dx^3-(1.623797632095822*fLPtr[1]*gPtr[3]*dy)/dx^3-(1.875*fPtr[2]*gLPtr[3]*dy)/dx^3+(1.875*fLPtr[2]*gLPtr[3]*dy)/dx^3+(1.623797632095822*fPtr[1]*gLPtr[3]*dy)/dx^3+(1.623797632095822*fLPtr[1]*gLPtr[3]*dy)/dx^3+(3.75*gRPtr[2]*fRPtr[3]*dy)/dx^3+(1.875*gPtr[2]*fRPtr[3]*dy)/dx^3-(1.623797632095822*gRPtr[1]*fRPtr[3]*dy)/dx^3+(1.623797632095822*gPtr[1]*fRPtr[3]*dy)/dx^3+(0.375*gRPtr[2]*fPtr[3]*dy)/dx^3+(10.5*gPtr[2]*fPtr[3]*dy)/dx^3+(0.375*gLPtr[2]*fPtr[3]*dy)/dx^3-(1.623797632095822*gRPtr[1]*fPtr[3]*dy)/dx^3+(1.623797632095822*gLPtr[1]*fPtr[3]*dy)/dx^3+(1.875*gPtr[2]*fLPtr[3]*dy)/dx^3+(3.75*gLPtr[2]*fLPtr[3]*dy)/dx^3-(1.623797632095822*gPtr[1]*fLPtr[3]*dy)/dx^3+(1.623797632095822*gLPtr[1]*fLPtr[3]*dy)/dx^3)+fOutPtr[4]
end

if isTopEdge then
    fOutPtr[1] = 0.5*dt*((1.40625*fPtr[4]*gPtr[4]*dx^2)/dy^3-(0.46875*fBPtr[4]*gPtr[4]*dx^2)/dy^3-(1.299038105676658*fPtr[2]*gPtr[4]*dx^2)/dy^3-(0.3247595264191644*fBPtr[2]*gPtr[4]*dx^2)/dy^3+(0.46875*fPtr[4]*gBPtr[4]*dx^2)/dy^3-(1.40625*fBPtr[4]*gBPtr[4]*dx^2)/dy^3-(0.3247595264191644*fPtr[2]*gBPtr[4]*dx^2)/dy^3-(1.299038105676658*fBPtr[2]*gBPtr[4]*dx^2)/dy^3-(0.5412658773652741*gPtr[2]*fPtr[4]*dx^2)/dy^3+(0.5412658773652741*gBPtr[2]*fPtr[4]*dx^2)/dy^3+(0.5412658773652741*gPtr[2]*fBPtr[4]*dx^2)/dy^3-(0.5412658773652741*gBPtr[2]*fBPtr[4]*dx^2)/dy^3+(0.46875*fPtr[2]*gPtr[2]*dx^2)/dy^3+(0.46875*fBPtr[2]*gPtr[2]*dx^2)/dy^3-(0.46875*fPtr[2]*gBPtr[2]*dx^2)/dy^3-(0.46875*fBPtr[2]*gBPtr[2]*dx^2)/dy^3+(5.625*fPtr[3]*gPtr[3])/dy^3-(1.875*fBPtr[3]*gPtr[3])/dy^3-(5.196152422706631*fPtr[1]*gPtr[3])/dy^3-(1.299038105676658*fBPtr[1]*gPtr[3])/dy^3+(1.875*fPtr[3]*gBPtr[3])/dy^3-(5.625*fBPtr[3]*gBPtr[3])/dy^3-(1.299038105676658*fPtr[1]*gBPtr[3])/dy^3-(5.196152422706631*fBPtr[1]*gBPtr[3])/dy^3-(2.165063509461096*gPtr[1]*fPtr[3])/dy^3+(2.165063509461096*gBPtr[1]*fPtr[3])/dy^3+(2.165063509461096*gPtr[1]*fBPtr[3])/dy^3-(2.165063509461096*gBPtr[1]*fBPtr[3])/dy^3+(1.875*fPtr[1]*gPtr[1])/dy^3+(1.875*fBPtr[1]*gPtr[1])/dy^3-(1.875*fPtr[1]*gBPtr[1])/dy^3-(1.875*fBPtr[1]*gBPtr[1])/dy^3)+fOutPtr[1]
    fOutPtr[2] = 0.5*dt*((2.8125*fPtr[3]*gPtr[4]*dx)/dy^3-(0.9375*fBPtr[3]*gPtr[4]*dx)/dy^3-(2.598076211353316*fPtr[1]*gPtr[4]*dx)/dy^3-(0.6495190528383289*fBPtr[1]*gPtr[4]*dx)/dy^3+(0.9375*fPtr[3]*gBPtr[4]*dx)/dy^3-(2.8125*fBPtr[3]*gBPtr[4]*dx)/dy^3-(0.6495190528383289*fPtr[1]*gBPtr[4]*dx)/dy^3-(2.598076211353316*fBPtr[1]*gBPtr[4]*dx)/dy^3+(2.8125*gPtr[3]*fPtr[4]*dx)/dy^3+(0.9375*gBPtr[3]*fPtr[4]*dx)/dy^3-(1.082531754730548*gPtr[1]*fPtr[4]*dx)/dy^3+(1.082531754730548*gBPtr[1]*fPtr[4]*dx)/dy^3-(0.9375*gPtr[3]*fBPtr[4]*dx)/dy^3-(2.8125*gBPtr[3]*fBPtr[4]*dx)/dy^3+(1.082531754730548*gPtr[1]*fBPtr[4]*dx)/dy^3-(1.082531754730548*gBPtr[1]*fBPtr[4]*dx)/dy^3-(2.598076211353316*fPtr[2]*gPtr[3]*dx)/dy^3-(0.6495190528383289*fBPtr[2]*gPtr[3]*dx)/dy^3-(0.6495190528383289*fPtr[2]*gBPtr[3]*dx)/dy^3-(2.598076211353316*fBPtr[2]*gBPtr[3]*dx)/dy^3-(1.082531754730548*gPtr[2]*fPtr[3]*dx)/dy^3+(1.082531754730548*gBPtr[2]*fPtr[3]*dx)/dy^3+(1.082531754730548*gPtr[2]*fBPtr[3]*dx)/dy^3-(1.082531754730548*gBPtr[2]*fBPtr[3]*dx)/dy^3+(0.9375*fPtr[1]*gPtr[2]*dx)/dy^3+(0.9375*fBPtr[1]*gPtr[2]*dx)/dy^3-(0.9375*fPtr[1]*gBPtr[2]*dx)/dy^3-(0.9375*fBPtr[1]*gBPtr[2]*dx)/dy^3+(0.9375*gPtr[1]*fPtr[2]*dx)/dy^3-(0.9375*gBPtr[1]*fPtr[2]*dx)/dy^3+(0.9375*gPtr[1]*fBPtr[2]*dx)/dy^3-(0.9375*gBPtr[1]*fBPtr[2]*dx)/dy^3)+fOutPtr[2]
    fOutPtr[3] = 0.5*dt*((-(2.868709150035952*fPtr[4]*gPtr[4]*dx^2)/dy^3)+(1.24491151794013*fBPtr[4]*gPtr[4]*dx^2)/dy^3+(2.625*fPtr[2]*gPtr[4]*dx^2)/dy^3+(0.9375*fBPtr[2]*gPtr[4]*dx^2)/dy^3-(0.3788861141556918*fPtr[4]*gBPtr[4]*dx^2)/dy^3+(2.002683746251514*fBPtr[4]*gBPtr[4]*dx^2)/dy^3+(0.1875*fPtr[2]*gBPtr[4]*dx^2)/dy^3+(1.875*fBPtr[2]*gBPtr[4]*dx^2)/dy^3+(0.9375*gPtr[2]*fPtr[4]*dx^2)/dy^3-(0.9375*gBPtr[2]*fPtr[4]*dx^2)/dy^3-(0.9375*gPtr[2]*fBPtr[4]*dx^2)/dy^3+(0.9375*gBPtr[2]*fBPtr[4]*dx^2)/dy^3-(0.8118988160479111*fPtr[2]*gPtr[2]*dx^2)/dy^3-(0.8118988160479111*fBPtr[2]*gPtr[2]*dx^2)/dy^3+(0.8118988160479111*fPtr[2]*gBPtr[2]*dx^2)/dy^3+(0.8118988160479111*fBPtr[2]*gBPtr[2]*dx^2)/dy^3-(11.47483660014381*fPtr[3]*gPtr[3])/dy^3+(4.979646071760522*fBPtr[3]*gPtr[3])/dy^3+(10.5*fPtr[1]*gPtr[3])/dy^3+(3.75*fBPtr[1]*gPtr[3])/dy^3-(1.515544456622767*fPtr[3]*gBPtr[3])/dy^3+(8.010734985006055*fBPtr[3]*gBPtr[3])/dy^3+(0.75*fPtr[1]*gBPtr[3])/dy^3+(7.5*fBPtr[1]*gBPtr[3])/dy^3+(3.75*gPtr[1]*fPtr[3])/dy^3-(3.75*gBPtr[1]*fPtr[3])/dy^3-(3.75*gPtr[1]*fBPtr[3])/dy^3+(3.75*gBPtr[1]*fBPtr[3])/dy^3-(3.247595264191645*fPtr[1]*gPtr[1])/dy^3-(3.247595264191645*fBPtr[1]*gPtr[1])/dy^3+(3.247595264191645*fPtr[1]*gBPtr[1])/dy^3+(3.247595264191645*fBPtr[1]*gBPtr[1])/dy^3)+fOutPtr[3]
    fOutPtr[4] = 0.5*dt*((-(5.737418300071904*fPtr[3]*gPtr[4]*dx)/dy^3)+(2.489823035880261*fBPtr[3]*gPtr[4]*dx)/dy^3+(5.25*fPtr[1]*gPtr[4]*dx)/dy^3+(1.875*fBPtr[1]*gPtr[4]*dx)/dy^3-(0.7577722283113837*fPtr[3]*gBPtr[4]*dx)/dy^3+(4.005367492503027*fBPtr[3]*gBPtr[4]*dx)/dy^3+(0.375*fPtr[1]*gBPtr[4]*dx)/dy^3+(3.75*fBPtr[1]*gBPtr[4]*dx)/dy^3-(5.737418300071904*gPtr[3]*fPtr[4]*dx)/dy^3-(0.7577722283113837*gBPtr[3]*fPtr[4]*dx)/dy^3+(1.875*gPtr[1]*fPtr[4]*dx)/dy^3-(1.875*gBPtr[1]*fPtr[4]*dx)/dy^3+(2.489823035880261*gPtr[3]*fBPtr[4]*dx)/dy^3+(4.005367492503027*gBPtr[3]*fBPtr[4]*dx)/dy^3-(1.875*gPtr[1]*fBPtr[4]*dx)/dy^3+(1.875*gBPtr[1]*fBPtr[4]*dx)/dy^3+(5.25*fPtr[2]*gPtr[3]*dx)/dy^3+(1.875*fBPtr[2]*gPtr[3]*dx)/dy^3+(0.375*fPtr[2]*gBPtr[3]*dx)/dy^3+(3.75*fBPtr[2]*gBPtr[3]*dx)/dy^3+(1.875*gPtr[2]*fPtr[3]*dx)/dy^3-(1.875*gBPtr[2]*fPtr[3]*dx)/dy^3-(1.875*gPtr[2]*fBPtr[3]*dx)/dy^3+(1.875*gBPtr[2]*fBPtr[3]*dx)/dy^3-(1.623797632095822*fPtr[1]*gPtr[2]*dx)/dy^3-(1.623797632095822*fBPtr[1]*gPtr[2]*dx)/dy^3+(1.623797632095822*fPtr[1]*gBPtr[2]*dx)/dy^3+(1.623797632095822*fBPtr[1]*gBPtr[2]*dx)/dy^3-(1.623797632095822*gPtr[1]*fPtr[2]*dx)/dy^3+(1.623797632095822*gBPtr[1]*fPtr[2]*dx)/dy^3-(1.623797632095822*gPtr[1]*fBPtr[2]*dx)/dy^3+(1.623797632095822*gBPtr[1]*fBPtr[2]*dx)/dy^3)+fOutPtr[4]
elseif isBotEdge then
    fOutPtr[1] = 0.5*dt*((-(1.40625*fTPtr[4]*gTPtr[4]*dx^2)/dy^3)+(0.46875*fPtr[4]*gTPtr[4]*dx^2)/dy^3+(1.299038105676658*fTPtr[2]*gTPtr[4]*dx^2)/dy^3+(0.3247595264191644*fPtr[2]*gTPtr[4]*dx^2)/dy^3-(0.46875*fTPtr[4]*gPtr[4]*dx^2)/dy^3+(1.40625*fPtr[4]*gPtr[4]*dx^2)/dy^3+(0.3247595264191644*fTPtr[2]*gPtr[4]*dx^2)/dy^3+(1.299038105676658*fPtr[2]*gPtr[4]*dx^2)/dy^3+(0.5412658773652741*gTPtr[2]*fTPtr[4]*dx^2)/dy^3-(0.5412658773652741*gPtr[2]*fTPtr[4]*dx^2)/dy^3-(0.5412658773652741*gTPtr[2]*fPtr[4]*dx^2)/dy^3+(0.5412658773652741*gPtr[2]*fPtr[4]*dx^2)/dy^3-(0.46875*fTPtr[2]*gTPtr[2]*dx^2)/dy^3-(0.46875*fPtr[2]*gTPtr[2]*dx^2)/dy^3+(0.46875*fTPtr[2]*gPtr[2]*dx^2)/dy^3+(0.46875*fPtr[2]*gPtr[2]*dx^2)/dy^3-(5.625*fTPtr[3]*gTPtr[3])/dy^3+(1.875*fPtr[3]*gTPtr[3])/dy^3+(5.196152422706631*fTPtr[1]*gTPtr[3])/dy^3+(1.299038105676658*fPtr[1]*gTPtr[3])/dy^3-(1.875*fTPtr[3]*gPtr[3])/dy^3+(5.625*fPtr[3]*gPtr[3])/dy^3+(1.299038105676658*fTPtr[1]*gPtr[3])/dy^3+(5.196152422706631*fPtr[1]*gPtr[3])/dy^3+(2.165063509461096*gTPtr[1]*fTPtr[3])/dy^3-(2.165063509461096*gPtr[1]*fTPtr[3])/dy^3-(2.165063509461096*gTPtr[1]*fPtr[3])/dy^3+(2.165063509461096*gPtr[1]*fPtr[3])/dy^3-(1.875*fTPtr[1]*gTPtr[1])/dy^3-(1.875*fPtr[1]*gTPtr[1])/dy^3+(1.875*fTPtr[1]*gPtr[1])/dy^3+(1.875*fPtr[1]*gPtr[1])/dy^3)+fOutPtr[1]
    fOutPtr[2] = 0.5*dt*((-(2.8125*fTPtr[3]*gTPtr[4]*dx)/dy^3)+(0.9375*fPtr[3]*gTPtr[4]*dx)/dy^3+(2.598076211353316*fTPtr[1]*gTPtr[4]*dx)/dy^3+(0.6495190528383289*fPtr[1]*gTPtr[4]*dx)/dy^3-(0.9375*fTPtr[3]*gPtr[4]*dx)/dy^3+(2.8125*fPtr[3]*gPtr[4]*dx)/dy^3+(0.6495190528383289*fTPtr[1]*gPtr[4]*dx)/dy^3+(2.598076211353316*fPtr[1]*gPtr[4]*dx)/dy^3-(2.8125*gTPtr[3]*fTPtr[4]*dx)/dy^3-(0.9375*gPtr[3]*fTPtr[4]*dx)/dy^3+(1.082531754730548*gTPtr[1]*fTPtr[4]*dx)/dy^3-(1.082531754730548*gPtr[1]*fTPtr[4]*dx)/dy^3+(0.9375*gTPtr[3]*fPtr[4]*dx)/dy^3+(2.8125*gPtr[3]*fPtr[4]*dx)/dy^3-(1.082531754730548*gTPtr[1]*fPtr[4]*dx)/dy^3+(1.082531754730548*gPtr[1]*fPtr[4]*dx)/dy^3+(2.598076211353316*fTPtr[2]*gTPtr[3]*dx)/dy^3+(0.6495190528383289*fPtr[2]*gTPtr[3]*dx)/dy^3+(0.6495190528383289*fTPtr[2]*gPtr[3]*dx)/dy^3+(2.598076211353316*fPtr[2]*gPtr[3]*dx)/dy^3+(1.082531754730548*gTPtr[2]*fTPtr[3]*dx)/dy^3-(1.082531754730548*gPtr[2]*fTPtr[3]*dx)/dy^3-(1.082531754730548*gTPtr[2]*fPtr[3]*dx)/dy^3+(1.082531754730548*gPtr[2]*fPtr[3]*dx)/dy^3-(0.9375*fTPtr[1]*gTPtr[2]*dx)/dy^3-(0.9375*fPtr[1]*gTPtr[2]*dx)/dy^3+(0.9375*fTPtr[1]*gPtr[2]*dx)/dy^3+(0.9375*fPtr[1]*gPtr[2]*dx)/dy^3-(0.9375*gTPtr[1]*fTPtr[2]*dx)/dy^3+(0.9375*gPtr[1]*fTPtr[2]*dx)/dy^3-(0.9375*gTPtr[1]*fPtr[2]*dx)/dy^3+(0.9375*gPtr[1]*fPtr[2]*dx)/dy^3)+fOutPtr[2]
    fOutPtr[3] = 0.5*dt*((-(2.002683746251514*fTPtr[4]*gTPtr[4]*dx^2)/dy^3)+(0.3788861141556918*fPtr[4]*gTPtr[4]*dx^2)/dy^3+(1.875*fTPtr[2]*gTPtr[4]*dx^2)/dy^3+(0.1875*fPtr[2]*gTPtr[4]*dx^2)/dy^3-(1.24491151794013*fTPtr[4]*gPtr[4]*dx^2)/dy^3+(2.868709150035952*fPtr[4]*gPtr[4]*dx^2)/dy^3+(0.9375*fTPtr[2]*gPtr[4]*dx^2)/dy^3+(2.625*fPtr[2]*gPtr[4]*dx^2)/dy^3+(0.9375*gTPtr[2]*fTPtr[4]*dx^2)/dy^3-(0.9375*gPtr[2]*fTPtr[4]*dx^2)/dy^3-(0.9375*gTPtr[2]*fPtr[4]*dx^2)/dy^3+(0.9375*gPtr[2]*fPtr[4]*dx^2)/dy^3-(0.8118988160479111*fTPtr[2]*gTPtr[2]*dx^2)/dy^3-(0.8118988160479111*fPtr[2]*gTPtr[2]*dx^2)/dy^3+(0.8118988160479111*fTPtr[2]*gPtr[2]*dx^2)/dy^3+(0.8118988160479111*fPtr[2]*gPtr[2]*dx^2)/dy^3-(8.010734985006055*fTPtr[3]*gTPtr[3])/dy^3+(1.515544456622767*fPtr[3]*gTPtr[3])/dy^3+(7.5*fTPtr[1]*gTPtr[3])/dy^3+(0.75*fPtr[1]*gTPtr[3])/dy^3-(4.979646071760522*fTPtr[3]*gPtr[3])/dy^3+(11.47483660014381*fPtr[3]*gPtr[3])/dy^3+(3.75*fTPtr[1]*gPtr[3])/dy^3+(10.5*fPtr[1]*gPtr[3])/dy^3+(3.75*gTPtr[1]*fTPtr[3])/dy^3-(3.75*gPtr[1]*fTPtr[3])/dy^3-(3.75*gTPtr[1]*fPtr[3])/dy^3+(3.75*gPtr[1]*fPtr[3])/dy^3-(3.247595264191645*fTPtr[1]*gTPtr[1])/dy^3-(3.247595264191645*fPtr[1]*gTPtr[1])/dy^3+(3.247595264191645*fTPtr[1]*gPtr[1])/dy^3+(3.247595264191645*fPtr[1]*gPtr[1])/dy^3)+fOutPtr[3]
    fOutPtr[4] = 0.5*dt*((-(4.005367492503027*fTPtr[3]*gTPtr[4]*dx)/dy^3)+(0.7577722283113837*fPtr[3]*gTPtr[4]*dx)/dy^3+(3.75*fTPtr[1]*gTPtr[4]*dx)/dy^3+(0.375*fPtr[1]*gTPtr[4]*dx)/dy^3-(2.489823035880261*fTPtr[3]*gPtr[4]*dx)/dy^3+(5.737418300071904*fPtr[3]*gPtr[4]*dx)/dy^3+(1.875*fTPtr[1]*gPtr[4]*dx)/dy^3+(5.25*fPtr[1]*gPtr[4]*dx)/dy^3-(4.005367492503027*gTPtr[3]*fTPtr[4]*dx)/dy^3-(2.489823035880261*gPtr[3]*fTPtr[4]*dx)/dy^3+(1.875*gTPtr[1]*fTPtr[4]*dx)/dy^3-(1.875*gPtr[1]*fTPtr[4]*dx)/dy^3+(0.7577722283113837*gTPtr[3]*fPtr[4]*dx)/dy^3+(5.737418300071904*gPtr[3]*fPtr[4]*dx)/dy^3-(1.875*gTPtr[1]*fPtr[4]*dx)/dy^3+(1.875*gPtr[1]*fPtr[4]*dx)/dy^3+(3.75*fTPtr[2]*gTPtr[3]*dx)/dy^3+(0.375*fPtr[2]*gTPtr[3]*dx)/dy^3+(1.875*fTPtr[2]*gPtr[3]*dx)/dy^3+(5.25*fPtr[2]*gPtr[3]*dx)/dy^3+(1.875*gTPtr[2]*fTPtr[3]*dx)/dy^3-(1.875*gPtr[2]*fTPtr[3]*dx)/dy^3-(1.875*gTPtr[2]*fPtr[3]*dx)/dy^3+(1.875*gPtr[2]*fPtr[3]*dx)/dy^3-(1.623797632095822*fTPtr[1]*gTPtr[2]*dx)/dy^3-(1.623797632095822*fPtr[1]*gTPtr[2]*dx)/dy^3+(1.623797632095822*fTPtr[1]*gPtr[2]*dx)/dy^3+(1.623797632095822*fPtr[1]*gPtr[2]*dx)/dy^3-(1.623797632095822*gTPtr[1]*fTPtr[2]*dx)/dy^3+(1.623797632095822*gPtr[1]*fTPtr[2]*dx)/dy^3-(1.623797632095822*gTPtr[1]*fPtr[2]*dx)/dy^3+(1.623797632095822*gPtr[1]*fPtr[2]*dx)/dy^3)+fOutPtr[4]
else
    fOutPtr[1] = 0.5*dt*((-(1.40625*fTPtr[4]*gTPtr[4]*dx^2)/dy^3)+(0.46875*fPtr[4]*gTPtr[4]*dx^2)/dy^3+(1.299038105676658*fTPtr[2]*gTPtr[4]*dx^2)/dy^3+(0.3247595264191644*fPtr[2]*gTPtr[4]*dx^2)/dy^3-(0.46875*fTPtr[4]*gPtr[4]*dx^2)/dy^3+(2.8125*fPtr[4]*gPtr[4]*dx^2)/dy^3-(0.46875*fBPtr[4]*gPtr[4]*dx^2)/dy^3+(0.3247595264191644*fTPtr[2]*gPtr[4]*dx^2)/dy^3-(0.3247595264191644*fBPtr[2]*gPtr[4]*dx^2)/dy^3+(0.46875*fPtr[4]*gBPtr[4]*dx^2)/dy^3-(1.40625*fBPtr[4]*gBPtr[4]*dx^2)/dy^3-(0.3247595264191644*fPtr[2]*gBPtr[4]*dx^2)/dy^3-(1.299038105676658*fBPtr[2]*gBPtr[4]*dx^2)/dy^3+(0.5412658773652741*gTPtr[2]*fTPtr[4]*dx^2)/dy^3-(0.5412658773652741*gPtr[2]*fTPtr[4]*dx^2)/dy^3-(0.5412658773652741*gTPtr[2]*fPtr[4]*dx^2)/dy^3+(0.5412658773652741*gBPtr[2]*fPtr[4]*dx^2)/dy^3+(0.5412658773652741*gPtr[2]*fBPtr[4]*dx^2)/dy^3-(0.5412658773652741*gBPtr[2]*fBPtr[4]*dx^2)/dy^3-(0.46875*fTPtr[2]*gTPtr[2]*dx^2)/dy^3-(0.46875*fPtr[2]*gTPtr[2]*dx^2)/dy^3+(0.46875*fTPtr[2]*gPtr[2]*dx^2)/dy^3+(0.9375*fPtr[2]*gPtr[2]*dx^2)/dy^3+(0.46875*fBPtr[2]*gPtr[2]*dx^2)/dy^3-(0.46875*fPtr[2]*gBPtr[2]*dx^2)/dy^3-(0.46875*fBPtr[2]*gBPtr[2]*dx^2)/dy^3-(5.625*fTPtr[3]*gTPtr[3])/dy^3+(1.875*fPtr[3]*gTPtr[3])/dy^3+(5.196152422706631*fTPtr[1]*gTPtr[3])/dy^3+(1.299038105676658*fPtr[1]*gTPtr[3])/dy^3-(1.875*fTPtr[3]*gPtr[3])/dy^3+(11.25*fPtr[3]*gPtr[3])/dy^3-(1.875*fBPtr[3]*gPtr[3])/dy^3+(1.299038105676658*fTPtr[1]*gPtr[3])/dy^3-(1.299038105676658*fBPtr[1]*gPtr[3])/dy^3+(1.875*fPtr[3]*gBPtr[3])/dy^3-(5.625*fBPtr[3]*gBPtr[3])/dy^3-(1.299038105676658*fPtr[1]*gBPtr[3])/dy^3-(5.196152422706631*fBPtr[1]*gBPtr[3])/dy^3+(2.165063509461096*gTPtr[1]*fTPtr[3])/dy^3-(2.165063509461096*gPtr[1]*fTPtr[3])/dy^3-(2.165063509461096*gTPtr[1]*fPtr[3])/dy^3+(2.165063509461096*gBPtr[1]*fPtr[3])/dy^3+(2.165063509461096*gPtr[1]*fBPtr[3])/dy^3-(2.165063509461096*gBPtr[1]*fBPtr[3])/dy^3-(1.875*fTPtr[1]*gTPtr[1])/dy^3-(1.875*fPtr[1]*gTPtr[1])/dy^3+(1.875*fTPtr[1]*gPtr[1])/dy^3+(3.75*fPtr[1]*gPtr[1])/dy^3+(1.875*fBPtr[1]*gPtr[1])/dy^3-(1.875*fPtr[1]*gBPtr[1])/dy^3-(1.875*fBPtr[1]*gBPtr[1])/dy^3)+fOutPtr[1]
    fOutPtr[2] = 0.5*dt*((-(2.8125*fTPtr[3]*gTPtr[4]*dx)/dy^3)+(0.9375*fPtr[3]*gTPtr[4]*dx)/dy^3+(2.598076211353316*fTPtr[1]*gTPtr[4]*dx)/dy^3+(0.6495190528383289*fPtr[1]*gTPtr[4]*dx)/dy^3-(0.9375*fTPtr[3]*gPtr[4]*dx)/dy^3+(5.625*fPtr[3]*gPtr[4]*dx)/dy^3-(0.9375*fBPtr[3]*gPtr[4]*dx)/dy^3+(0.6495190528383289*fTPtr[1]*gPtr[4]*dx)/dy^3-(0.6495190528383289*fBPtr[1]*gPtr[4]*dx)/dy^3+(0.9375*fPtr[3]*gBPtr[4]*dx)/dy^3-(2.8125*fBPtr[3]*gBPtr[4]*dx)/dy^3-(0.6495190528383289*fPtr[1]*gBPtr[4]*dx)/dy^3-(2.598076211353316*fBPtr[1]*gBPtr[4]*dx)/dy^3-(2.8125*gTPtr[3]*fTPtr[4]*dx)/dy^3-(0.9375*gPtr[3]*fTPtr[4]*dx)/dy^3+(1.082531754730548*gTPtr[1]*fTPtr[4]*dx)/dy^3-(1.082531754730548*gPtr[1]*fTPtr[4]*dx)/dy^3+(0.9375*gTPtr[3]*fPtr[4]*dx)/dy^3+(5.625*gPtr[3]*fPtr[4]*dx)/dy^3+(0.9375*gBPtr[3]*fPtr[4]*dx)/dy^3-(1.082531754730548*gTPtr[1]*fPtr[4]*dx)/dy^3+(1.082531754730548*gBPtr[1]*fPtr[4]*dx)/dy^3-(0.9375*gPtr[3]*fBPtr[4]*dx)/dy^3-(2.8125*gBPtr[3]*fBPtr[4]*dx)/dy^3+(1.082531754730548*gPtr[1]*fBPtr[4]*dx)/dy^3-(1.082531754730548*gBPtr[1]*fBPtr[4]*dx)/dy^3+(2.598076211353316*fTPtr[2]*gTPtr[3]*dx)/dy^3+(0.6495190528383289*fPtr[2]*gTPtr[3]*dx)/dy^3+(0.6495190528383289*fTPtr[2]*gPtr[3]*dx)/dy^3-(0.6495190528383289*fBPtr[2]*gPtr[3]*dx)/dy^3-(0.6495190528383289*fPtr[2]*gBPtr[3]*dx)/dy^3-(2.598076211353316*fBPtr[2]*gBPtr[3]*dx)/dy^3+(1.082531754730548*gTPtr[2]*fTPtr[3]*dx)/dy^3-(1.082531754730548*gPtr[2]*fTPtr[3]*dx)/dy^3-(1.082531754730548*gTPtr[2]*fPtr[3]*dx)/dy^3+(1.082531754730548*gBPtr[2]*fPtr[3]*dx)/dy^3+(1.082531754730548*gPtr[2]*fBPtr[3]*dx)/dy^3-(1.082531754730548*gBPtr[2]*fBPtr[3]*dx)/dy^3-(0.9375*fTPtr[1]*gTPtr[2]*dx)/dy^3-(0.9375*fPtr[1]*gTPtr[2]*dx)/dy^3+(0.9375*fTPtr[1]*gPtr[2]*dx)/dy^3+(1.875*fPtr[1]*gPtr[2]*dx)/dy^3+(0.9375*fBPtr[1]*gPtr[2]*dx)/dy^3-(0.9375*fPtr[1]*gBPtr[2]*dx)/dy^3-(0.9375*fBPtr[1]*gBPtr[2]*dx)/dy^3-(0.9375*gTPtr[1]*fTPtr[2]*dx)/dy^3+(0.9375*gPtr[1]*fTPtr[2]*dx)/dy^3-(0.9375*gTPtr[1]*fPtr[2]*dx)/dy^3+(1.875*gPtr[1]*fPtr[2]*dx)/dy^3-(0.9375*gBPtr[1]*fPtr[2]*dx)/dy^3+(0.9375*gPtr[1]*fBPtr[2]*dx)/dy^3-(0.9375*gBPtr[1]*fBPtr[2]*dx)/dy^3)+fOutPtr[2]
    fOutPtr[3] = 0.5*dt*((-(2.002683746251514*fTPtr[4]*gTPtr[4]*dx^2)/dy^3)+(0.3788861141556918*fPtr[4]*gTPtr[4]*dx^2)/dy^3+(1.875*fTPtr[2]*gTPtr[4]*dx^2)/dy^3+(0.1875*fPtr[2]*gTPtr[4]*dx^2)/dy^3-(1.24491151794013*fTPtr[4]*gPtr[4]*dx^2)/dy^3+(1.24491151794013*fBPtr[4]*gPtr[4]*dx^2)/dy^3+(0.9375*fTPtr[2]*gPtr[4]*dx^2)/dy^3+(5.25*fPtr[2]*gPtr[4]*dx^2)/dy^3+(0.9375*fBPtr[2]*gPtr[4]*dx^2)/dy^3-(0.3788861141556918*fPtr[4]*gBPtr[4]*dx^2)/dy^3+(2.002683746251514*fBPtr[4]*gBPtr[4]*dx^2)/dy^3+(0.1875*fPtr[2]*gBPtr[4]*dx^2)/dy^3+(1.875*fBPtr[2]*gBPtr[4]*dx^2)/dy^3+(0.9375*gTPtr[2]*fTPtr[4]*dx^2)/dy^3-(0.9375*gPtr[2]*fTPtr[4]*dx^2)/dy^3-(0.9375*gTPtr[2]*fPtr[4]*dx^2)/dy^3+(1.875*gPtr[2]*fPtr[4]*dx^2)/dy^3-(0.9375*gBPtr[2]*fPtr[4]*dx^2)/dy^3-(0.9375*gPtr[2]*fBPtr[4]*dx^2)/dy^3+(0.9375*gBPtr[2]*fBPtr[4]*dx^2)/dy^3-(0.8118988160479111*fTPtr[2]*gTPtr[2]*dx^2)/dy^3-(0.8118988160479111*fPtr[2]*gTPtr[2]*dx^2)/dy^3+(0.8118988160479111*fTPtr[2]*gPtr[2]*dx^2)/dy^3-(0.8118988160479111*fBPtr[2]*gPtr[2]*dx^2)/dy^3+(0.8118988160479111*fPtr[2]*gBPtr[2]*dx^2)/dy^3+(0.8118988160479111*fBPtr[2]*gBPtr[2]*dx^2)/dy^3-(8.010734985006055*fTPtr[3]*gTPtr[3])/dy^3+(1.515544456622767*fPtr[3]*gTPtr[3])/dy^3+(7.5*fTPtr[1]*gTPtr[3])/dy^3+(0.75*fPtr[1]*gTPtr[3])/dy^3-(4.979646071760522*fTPtr[3]*gPtr[3])/dy^3+(4.979646071760522*fBPtr[3]*gPtr[3])/dy^3+(3.75*fTPtr[1]*gPtr[3])/dy^3+(21.0*fPtr[1]*gPtr[3])/dy^3+(3.75*fBPtr[1]*gPtr[3])/dy^3-(1.515544456622767*fPtr[3]*gBPtr[3])/dy^3+(8.010734985006055*fBPtr[3]*gBPtr[3])/dy^3+(0.75*fPtr[1]*gBPtr[3])/dy^3+(7.5*fBPtr[1]*gBPtr[3])/dy^3+(3.75*gTPtr[1]*fTPtr[3])/dy^3-(3.75*gPtr[1]*fTPtr[3])/dy^3-(3.75*gTPtr[1]*fPtr[3])/dy^3+(7.5*gPtr[1]*fPtr[3])/dy^3-(3.75*gBPtr[1]*fPtr[3])/dy^3-(3.75*gPtr[1]*fBPtr[3])/dy^3+(3.75*gBPtr[1]*fBPtr[3])/dy^3-(3.247595264191645*fTPtr[1]*gTPtr[1])/dy^3-(3.247595264191645*fPtr[1]*gTPtr[1])/dy^3+(3.247595264191645*fTPtr[1]*gPtr[1])/dy^3-(3.247595264191645*fBPtr[1]*gPtr[1])/dy^3+(3.247595264191645*fPtr[1]*gBPtr[1])/dy^3+(3.247595264191645*fBPtr[1]*gBPtr[1])/dy^3)+fOutPtr[3]
    fOutPtr[4] = 0.5*dt*((-(4.005367492503027*fTPtr[3]*gTPtr[4]*dx)/dy^3)+(0.7577722283113837*fPtr[3]*gTPtr[4]*dx)/dy^3+(3.75*fTPtr[1]*gTPtr[4]*dx)/dy^3+(0.375*fPtr[1]*gTPtr[4]*dx)/dy^3-(2.489823035880261*fTPtr[3]*gPtr[4]*dx)/dy^3+(2.489823035880261*fBPtr[3]*gPtr[4]*dx)/dy^3+(1.875*fTPtr[1]*gPtr[4]*dx)/dy^3+(10.5*fPtr[1]*gPtr[4]*dx)/dy^3+(1.875*fBPtr[1]*gPtr[4]*dx)/dy^3-(0.7577722283113837*fPtr[3]*gBPtr[4]*dx)/dy^3+(4.005367492503027*fBPtr[3]*gBPtr[4]*dx)/dy^3+(0.375*fPtr[1]*gBPtr[4]*dx)/dy^3+(3.75*fBPtr[1]*gBPtr[4]*dx)/dy^3-(4.005367492503027*gTPtr[3]*fTPtr[4]*dx)/dy^3-(2.489823035880261*gPtr[3]*fTPtr[4]*dx)/dy^3+(1.875*gTPtr[1]*fTPtr[4]*dx)/dy^3-(1.875*gPtr[1]*fTPtr[4]*dx)/dy^3+(0.7577722283113837*gTPtr[3]*fPtr[4]*dx)/dy^3-(0.7577722283113837*gBPtr[3]*fPtr[4]*dx)/dy^3-(1.875*gTPtr[1]*fPtr[4]*dx)/dy^3+(3.75*gPtr[1]*fPtr[4]*dx)/dy^3-(1.875*gBPtr[1]*fPtr[4]*dx)/dy^3+(2.489823035880261*gPtr[3]*fBPtr[4]*dx)/dy^3+(4.005367492503027*gBPtr[3]*fBPtr[4]*dx)/dy^3-(1.875*gPtr[1]*fBPtr[4]*dx)/dy^3+(1.875*gBPtr[1]*fBPtr[4]*dx)/dy^3+(3.75*fTPtr[2]*gTPtr[3]*dx)/dy^3+(0.375*fPtr[2]*gTPtr[3]*dx)/dy^3+(1.875*fTPtr[2]*gPtr[3]*dx)/dy^3+(10.5*fPtr[2]*gPtr[3]*dx)/dy^3+(1.875*fBPtr[2]*gPtr[3]*dx)/dy^3+(0.375*fPtr[2]*gBPtr[3]*dx)/dy^3+(3.75*fBPtr[2]*gBPtr[3]*dx)/dy^3+(1.875*gTPtr[2]*fTPtr[3]*dx)/dy^3-(1.875*gPtr[2]*fTPtr[3]*dx)/dy^3-(1.875*gTPtr[2]*fPtr[3]*dx)/dy^3+(3.75*gPtr[2]*fPtr[3]*dx)/dy^3-(1.875*gBPtr[2]*fPtr[3]*dx)/dy^3-(1.875*gPtr[2]*fBPtr[3]*dx)/dy^3+(1.875*gBPtr[2]*fBPtr[3]*dx)/dy^3-(1.623797632095822*fTPtr[1]*gTPtr[2]*dx)/dy^3-(1.623797632095822*fPtr[1]*gTPtr[2]*dx)/dy^3+(1.623797632095822*fTPtr[1]*gPtr[2]*dx)/dy^3-(1.623797632095822*fBPtr[1]*gPtr[2]*dx)/dy^3+(1.623797632095822*fPtr[1]*gBPtr[2]*dx)/dy^3+(1.623797632095822*fBPtr[1]*gBPtr[2]*dx)/dy^3-(1.623797632095822*gTPtr[1]*fTPtr[2]*dx)/dy^3+(1.623797632095822*gPtr[1]*fTPtr[2]*dx)/dy^3-(1.623797632095822*gTPtr[1]*fPtr[2]*dx)/dy^3+(1.623797632095822*gBPtr[1]*fPtr[2]*dx)/dy^3-(1.623797632095822*gPtr[1]*fBPtr[2]*dx)/dy^3+(1.623797632095822*gBPtr[1]*fBPtr[2]*dx)/dy^3)+fOutPtr[4]
end





end

return diffStencilFunc
   
