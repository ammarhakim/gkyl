function dragStencilFunc(dt, dx, dy,
			 fPtr, fLPtr, fRPtr, fTPtr, fBPtr,
			 hPtr, hLPtr, hRPtr, hTPtr, hBPtr,
			 isTopEdge, isBotEdge, isLeftEdge, isRightEdge,
			 fOutPtr) 
local hRecovL = {} 
hRecovL[1] = 0.08333333333333333*((-3.464101615137754*hPtr[2])+3.464101615137754*hLPtr[2]+3.0*hPtr[1]+3.0*hLPtr[1]) 
hRecovL[2] = (0.5*((-2.0*hPtr[4])+2.0*hLPtr[4]+1.732050807568877*(hPtr[3]+hLPtr[3])))/dy 
hRecovL[3] = -(0.125*(8.660254037844386*hPtr[2]+8.660254037844386*hLPtr[2]-9.0*hPtr[1]+9.0*hLPtr[1]))/dx 
hRecovL[4] = -(0.25*(15.0*hPtr[4]+15.0*hLPtr[4]+1.732050807568877*(9.0*hLPtr[3]-9.0*hPtr[3])))/(dx*dy) 
hRecovL[5] = -(0.5*(1.732050807568877*hLPtr[2]-1.732050807568877*hPtr[2]))/dx^2 
hRecovL[6] = -(1.0*(3.0*hLPtr[4]-3.0*hPtr[4]))/(dx^2*dy) 
hRecovL[7] = (0.25*(8.660254037844386*hPtr[2]+8.660254037844386*hLPtr[2]-5.0*hPtr[1]+5.0*hLPtr[1]))/dx^3 
hRecovL[8] = (0.5*(15.0*hPtr[4]+15.0*hLPtr[4]+1.732050807568877*(5.0*hLPtr[3]-5.0*hPtr[3])))/(dx^3*dy) 

local hRecovR = {} 
hRecovR[1] = 0.08333333333333333*((-3.464101615137754*hRPtr[2])+3.464101615137754*hPtr[2]+3.0*hRPtr[1]+3.0*hPtr[1]) 
hRecovR[2] = (0.5*((-2.0*hRPtr[4])+2.0*hPtr[4]+1.732050807568877*(hRPtr[3]+hPtr[3])))/dy 
hRecovR[3] = -(0.125*(8.660254037844386*hRPtr[2]+8.660254037844386*hPtr[2]-9.0*hRPtr[1]+9.0*hPtr[1]))/dx 
hRecovR[4] = -(0.25*(15.0*hRPtr[4]+15.0*hPtr[4]+1.732050807568877*(9.0*hPtr[3]-9.0*hRPtr[3])))/(dx*dy) 
hRecovR[5] = -(0.5*(1.732050807568877*hPtr[2]-1.732050807568877*hRPtr[2]))/dx^2 
hRecovR[6] = -(1.0*(3.0*hPtr[4]-3.0*hRPtr[4]))/(dx^2*dy) 
hRecovR[7] = (0.25*(8.660254037844386*hRPtr[2]+8.660254037844386*hPtr[2]-5.0*hRPtr[1]+5.0*hPtr[1]))/dx^3 
hRecovR[8] = (0.5*(15.0*hRPtr[4]+15.0*hPtr[4]+1.732050807568877*(5.0*hPtr[3]-5.0*hRPtr[3])))/(dx^3*dy) 

local hxfRecovL = {} 
hxfRecovL[1] = 0.001388888888888889*(((27.0*fLPtr[4]-27.0*fPtr[4])*hRecovL[8]*dx^2+((-30.0*fPtr[4])-30.0*fLPtr[4]+17.32050807568877*fPtr[3]-17.32050807568877*fLPtr[3])*hRecovL[6]*dx+((-60.0*fPtr[4])+60.0*fLPtr[4]+51.96152422706631*fPtr[3]+51.96152422706631*fLPtr[3])*hRecovL[4])*dy+(93.53074360871933*fLPtr[2]-93.53074360871933*fPtr[2])*hRecovL[7]*dx^2+((-103.9230484541326*fPtr[2])-103.9230484541326*fLPtr[2]+60.0*fPtr[1]-60.0*fLPtr[1])*hRecovL[5]*dx+((-207.8460969082653*fPtr[2])+207.8460969082653*fLPtr[2]+180.0*fPtr[1]+180.0*fLPtr[1])*hRecovL[3]) 
hxfRecovL[2] = (0.008333333333333333*(((15.58845726811989*fLPtr[2]-15.58845726811989*fPtr[2])*hRecovL[8]*dx^2+((-17.32050807568877*fPtr[2])-17.32050807568877*fLPtr[2]+10.0*fPtr[1]-10.0*fLPtr[1])*hRecovL[6]*dx+((-34.64101615137754*fPtr[2])+34.64101615137754*fLPtr[2]+30.0*fPtr[1]+30.0*fLPtr[1])*hRecovL[4])*dy+(54.0*fLPtr[4]-54.0*fPtr[4])*hRecovL[7]*dx^2+((-60.0*fPtr[4])-60.0*fLPtr[4]+34.64101615137754*fPtr[3]-34.64101615137754*fLPtr[3])*hRecovL[5]*dx+hRecovL[3]*((-120.0*fPtr[4])+103.9230484541326*fPtr[3]+103.9230484541326*fLPtr[3])+120.0*hRecovL[3]*fLPtr[4]))/dy 
hxfRecovL[3] = -(0.01041666666666667*(((9.0*fPtr[4]+9.0*fLPtr[4]-5.196152422706631*fPtr[3]+5.196152422706631*fLPtr[3])*hRecovL[8]*dx^2+(12.0*fPtr[4]-12.0*fLPtr[4]-13.85640646055102*fPtr[3]-13.85640646055102*fLPtr[3])*hRecovL[6]*dx+(30.0*fPtr[4]+30.0*fLPtr[4]-31.17691453623978*fPtr[3]+31.17691453623978*fLPtr[3])*hRecovL[4])*dy+(31.17691453623978*fPtr[2]+31.17691453623978*fLPtr[2]-18.0*fPtr[1]+18.0*fLPtr[1])*hRecovL[7]*dx^2+(41.56921938165305*fPtr[2]-41.56921938165305*fLPtr[2]-48.0*fPtr[1]-48.0*fLPtr[1])*hRecovL[5]*dx+(103.9230484541326*fPtr[2]+103.9230484541326*fLPtr[2]-108.0*fPtr[1]+108.0*fLPtr[1])*hRecovL[3]))/dx 
hxfRecovL[4] = -(0.0625*(((5.196152422706631*fPtr[2]+5.196152422706631*fLPtr[2]-3.0*fPtr[1]+3.0*fLPtr[1])*hRecovL[8]*dx^2+(6.928203230275509*fPtr[2]-6.928203230275509*fLPtr[2]-8.0*fPtr[1]-8.0*fLPtr[1])*hRecovL[6]*dx+(17.32050807568877*fPtr[2]+17.32050807568877*fLPtr[2]-18.0*fPtr[1]+18.0*fLPtr[1])*hRecovL[4])*dy+(18.0*fPtr[4]+18.0*fLPtr[4]-10.39230484541326*fPtr[3]+10.39230484541326*fLPtr[3])*hRecovL[7]*dx^2+(24.0*fPtr[4]-24.0*fLPtr[4]-27.71281292110204*fPtr[3]-27.71281292110204*fLPtr[3])*hRecovL[5]*dx+hRecovL[3]*(60.0*fPtr[4]-62.35382907247956*fPtr[3]+62.35382907247956*fLPtr[3])+60.0*hRecovL[3]*fLPtr[4]))/(dx*dy) 
hxfRecovL[5] = -(0.008333333333333333*((((-36.0*fPtr[4])+36.0*fLPtr[4]-25.98076211353316*fPtr[3]-25.98076211353316*fLPtr[3])*hRecovL[8]*dx^2+((-30.0*fPtr[4])-30.0*fLPtr[4]-17.32050807568877*fPtr[3]+17.32050807568877*fLPtr[3])*hRecovL[6]*dx+(30.0*fLPtr[4]-30.0*fPtr[4])*hRecovL[4])*dy+((-124.7076581449591*fPtr[2])+124.7076581449591*fLPtr[2]-90.0*fPtr[1]-90.0*fLPtr[1])*hRecovL[7]*dx^2+((-103.9230484541326*fPtr[2])-103.9230484541326*fLPtr[2]-60.0*fPtr[1]+60.0*fLPtr[1])*hRecovL[5]*dx+(103.9230484541326*fLPtr[2]-103.9230484541326*fPtr[2])*hRecovL[3]))/dx^2 
hxfRecovL[6] = -(0.05*((((-20.78460969082652*fPtr[2])+20.78460969082652*fLPtr[2]-15.0*fPtr[1]-15.0*fLPtr[1])*hRecovL[8]*dx^2+((-17.32050807568877*fPtr[2])-17.32050807568877*fLPtr[2]-10.0*fPtr[1]+10.0*fLPtr[1])*hRecovL[6]*dx+(17.32050807568877*fLPtr[2]-17.32050807568877*fPtr[2])*hRecovL[4])*dy+((-72.0*fPtr[4])+72.0*fLPtr[4]-51.96152422706631*fPtr[3]-51.96152422706631*fLPtr[3])*hRecovL[7]*dx^2+((-60.0*fPtr[4])-60.0*fLPtr[4]-34.64101615137754*fPtr[3]+34.64101615137754*fLPtr[3])*hRecovL[5]*dx-60.0*hRecovL[3]*fPtr[4]+60.0*hRecovL[3]*fLPtr[4]))/(dx^2*dy) 
hxfRecovL[7] = (0.02083333333333333*(((21.0*fPtr[4]+21.0*fLPtr[4]+8.660254037844386*fPtr[3]-8.660254037844386*fLPtr[3])*hRecovL[8]*dx^2+(20.0*fPtr[4]-20.0*fLPtr[4])*hRecovL[6]*dx+(30.0*fPtr[4]+30.0*fLPtr[4]-17.32050807568877*fPtr[3]+17.32050807568877*fLPtr[3])*hRecovL[4])*dy+(72.74613391789283*fPtr[2]+72.74613391789283*fLPtr[2]+30.0*fPtr[1]-30.0*fLPtr[1])*hRecovL[7]*dx^2+(69.28203230275508*fPtr[2]-69.28203230275508*fLPtr[2])*hRecovL[5]*dx+(103.9230484541326*fPtr[2]+103.9230484541326*fLPtr[2]-60.0*fPtr[1]+60.0*fLPtr[1])*hRecovL[3]))/dx^3 
hxfRecovL[8] = (0.04166666666666666*(((36.37306695894642*fPtr[2]+36.37306695894642*fLPtr[2]+15.0*fPtr[1]-15.0*fLPtr[1])*hRecovL[8]*dx^2+(34.64101615137754*fPtr[2]-34.64101615137754*fLPtr[2])*hRecovL[6]*dx+(51.96152422706631*fPtr[2]+51.96152422706631*fLPtr[2]-30.0*fPtr[1]+30.0*fLPtr[1])*hRecovL[4])*dy+(126.0*fPtr[4]+126.0*fLPtr[4]+51.96152422706631*fPtr[3]-51.96152422706631*fLPtr[3])*hRecovL[7]*dx^2+(120.0*fPtr[4]-120.0*fLPtr[4])*hRecovL[5]*dx+hRecovL[3]*(180.0*fPtr[4]-103.9230484541326*fPtr[3]+103.9230484541326*fLPtr[3])+180.0*hRecovL[3]*fLPtr[4]))/(dx^3*dy) 

local hxfRecovR = {} 
hxfRecovR[1] = 0.001388888888888889*(((27.0*fPtr[4]-27.0*fRPtr[4])*hRecovR[8]*dx^2+((-30.0*fRPtr[4])-30.0*fPtr[4]+17.32050807568877*fRPtr[3]-17.32050807568877*fPtr[3])*hRecovR[6]*dx+((-60.0*fRPtr[4])+60.0*fPtr[4]+51.96152422706631*fRPtr[3]+51.96152422706631*fPtr[3])*hRecovR[4])*dy+(93.53074360871933*fPtr[2]-93.53074360871933*fRPtr[2])*hRecovR[7]*dx^2+((-103.9230484541326*fRPtr[2])-103.9230484541326*fPtr[2]+60.0*fRPtr[1]-60.0*fPtr[1])*hRecovR[5]*dx+((-207.8460969082653*fRPtr[2])+207.8460969082653*fPtr[2]+180.0*fRPtr[1]+180.0*fPtr[1])*hRecovR[3]) 
hxfRecovR[2] = (0.008333333333333333*(((15.58845726811989*fPtr[2]-15.58845726811989*fRPtr[2])*hRecovR[8]*dx^2+((-17.32050807568877*fRPtr[2])-17.32050807568877*fPtr[2]+10.0*fRPtr[1]-10.0*fPtr[1])*hRecovR[6]*dx+((-34.64101615137754*fRPtr[2])+34.64101615137754*fPtr[2]+30.0*fRPtr[1]+30.0*fPtr[1])*hRecovR[4])*dy+(54.0*fPtr[4]-54.0*fRPtr[4])*hRecovR[7]*dx^2+((-60.0*fRPtr[4])-60.0*fPtr[4]+34.64101615137754*fRPtr[3]-34.64101615137754*fPtr[3])*hRecovR[5]*dx+hRecovR[3]*((-120.0*fRPtr[4])+103.9230484541326*fRPtr[3]+103.9230484541326*fPtr[3])+120.0*hRecovR[3]*fPtr[4]))/dy 
hxfRecovR[3] = -(0.01041666666666667*(((9.0*fRPtr[4]+9.0*fPtr[4]-5.196152422706631*fRPtr[3]+5.196152422706631*fPtr[3])*hRecovR[8]*dx^2+(12.0*fRPtr[4]-12.0*fPtr[4]-13.85640646055102*fRPtr[3]-13.85640646055102*fPtr[3])*hRecovR[6]*dx+(30.0*fRPtr[4]+30.0*fPtr[4]-31.17691453623978*fRPtr[3]+31.17691453623978*fPtr[3])*hRecovR[4])*dy+(31.17691453623978*fRPtr[2]+31.17691453623978*fPtr[2]-18.0*fRPtr[1]+18.0*fPtr[1])*hRecovR[7]*dx^2+(41.56921938165305*fRPtr[2]-41.56921938165305*fPtr[2]-48.0*fRPtr[1]-48.0*fPtr[1])*hRecovR[5]*dx+(103.9230484541326*fRPtr[2]+103.9230484541326*fPtr[2]-108.0*fRPtr[1]+108.0*fPtr[1])*hRecovR[3]))/dx 
hxfRecovR[4] = -(0.0625*(((5.196152422706631*fRPtr[2]+5.196152422706631*fPtr[2]-3.0*fRPtr[1]+3.0*fPtr[1])*hRecovR[8]*dx^2+(6.928203230275509*fRPtr[2]-6.928203230275509*fPtr[2]-8.0*fRPtr[1]-8.0*fPtr[1])*hRecovR[6]*dx+(17.32050807568877*fRPtr[2]+17.32050807568877*fPtr[2]-18.0*fRPtr[1]+18.0*fPtr[1])*hRecovR[4])*dy+(18.0*fRPtr[4]+18.0*fPtr[4]-10.39230484541326*fRPtr[3]+10.39230484541326*fPtr[3])*hRecovR[7]*dx^2+(24.0*fRPtr[4]-24.0*fPtr[4]-27.71281292110204*fRPtr[3]-27.71281292110204*fPtr[3])*hRecovR[5]*dx+hRecovR[3]*(60.0*fRPtr[4]-62.35382907247956*fRPtr[3]+62.35382907247956*fPtr[3])+60.0*hRecovR[3]*fPtr[4]))/(dx*dy) 
hxfRecovR[5] = -(0.008333333333333333*((((-36.0*fRPtr[4])+36.0*fPtr[4]-25.98076211353316*fRPtr[3]-25.98076211353316*fPtr[3])*hRecovR[8]*dx^2+((-30.0*fRPtr[4])-30.0*fPtr[4]-17.32050807568877*fRPtr[3]+17.32050807568877*fPtr[3])*hRecovR[6]*dx+(30.0*fPtr[4]-30.0*fRPtr[4])*hRecovR[4])*dy+((-124.7076581449591*fRPtr[2])+124.7076581449591*fPtr[2]-90.0*fRPtr[1]-90.0*fPtr[1])*hRecovR[7]*dx^2+((-103.9230484541326*fRPtr[2])-103.9230484541326*fPtr[2]-60.0*fRPtr[1]+60.0*fPtr[1])*hRecovR[5]*dx+(103.9230484541326*fPtr[2]-103.9230484541326*fRPtr[2])*hRecovR[3]))/dx^2 
hxfRecovR[6] = -(0.05*((((-20.78460969082652*fRPtr[2])+20.78460969082652*fPtr[2]-15.0*fRPtr[1]-15.0*fPtr[1])*hRecovR[8]*dx^2+((-17.32050807568877*fRPtr[2])-17.32050807568877*fPtr[2]-10.0*fRPtr[1]+10.0*fPtr[1])*hRecovR[6]*dx+(17.32050807568877*fPtr[2]-17.32050807568877*fRPtr[2])*hRecovR[4])*dy+((-72.0*fRPtr[4])+72.0*fPtr[4]-51.96152422706631*fRPtr[3]-51.96152422706631*fPtr[3])*hRecovR[7]*dx^2+((-60.0*fRPtr[4])-60.0*fPtr[4]-34.64101615137754*fRPtr[3]+34.64101615137754*fPtr[3])*hRecovR[5]*dx-60.0*hRecovR[3]*fRPtr[4]+60.0*hRecovR[3]*fPtr[4]))/(dx^2*dy) 
hxfRecovR[7] = (0.02083333333333333*(((21.0*fRPtr[4]+21.0*fPtr[4]+8.660254037844386*fRPtr[3]-8.660254037844386*fPtr[3])*hRecovR[8]*dx^2+(20.0*fRPtr[4]-20.0*fPtr[4])*hRecovR[6]*dx+(30.0*fRPtr[4]+30.0*fPtr[4]-17.32050807568877*fRPtr[3]+17.32050807568877*fPtr[3])*hRecovR[4])*dy+(72.74613391789283*fRPtr[2]+72.74613391789283*fPtr[2]+30.0*fRPtr[1]-30.0*fPtr[1])*hRecovR[7]*dx^2+(69.28203230275508*fRPtr[2]-69.28203230275508*fPtr[2])*hRecovR[5]*dx+(103.9230484541326*fRPtr[2]+103.9230484541326*fPtr[2]-60.0*fRPtr[1]+60.0*fPtr[1])*hRecovR[3]))/dx^3 
hxfRecovR[8] = (0.04166666666666666*(((36.37306695894642*fRPtr[2]+36.37306695894642*fPtr[2]+15.0*fRPtr[1]-15.0*fPtr[1])*hRecovR[8]*dx^2+(34.64101615137754*fRPtr[2]-34.64101615137754*fPtr[2])*hRecovR[6]*dx+(51.96152422706631*fRPtr[2]+51.96152422706631*fPtr[2]-30.0*fRPtr[1]+30.0*fPtr[1])*hRecovR[4])*dy+(126.0*fRPtr[4]+126.0*fPtr[4]+51.96152422706631*fRPtr[3]-51.96152422706631*fPtr[3])*hRecovR[7]*dx^2+(120.0*fRPtr[4]-120.0*fPtr[4])*hRecovR[5]*dx+hRecovR[3]*(180.0*fRPtr[4]-103.9230484541326*fRPtr[3]+103.9230484541326*fPtr[3])+180.0*hRecovR[3]*fPtr[4]))/(dx^3*dy) 

if isLeftEdge then
    fOutPtr[1] = fPtr[1]-hxfRecovR[1]*dt 
    fOutPtr[2] = fPtr[2]-1.732050807568877*hxfRecovR[1]*dt 
    fOutPtr[3] = fPtr[3]-0.2886751345948129*hxfRecovR[2]*dt*dy 
    fOutPtr[4] = fPtr[4]-0.5*hxfRecovR[2]*dt*dy 
elseif isRightEdge then
    fOutPtr[1] = 1.0*hxfRecovL[1]*dt+fPtr[1] 
    fOutPtr[2] = fPtr[2]-1.732050807568877*hxfRecovL[1]*dt 
    fOutPtr[3] = 0.2886751345948129*hxfRecovL[2]*dt*dy+fPtr[3] 
    fOutPtr[4] = fPtr[4]-0.5*hxfRecovL[2]*dt*dy 
else
    fOutPtr[1] = fPtr[1]-(hxfRecovR[1]-1.0*hxfRecovL[1])*dt 
    fOutPtr[2] = fPtr[2]-(1.732050807568877*hxfRecovR[1]+1.732050807568877*hxfRecovL[1])*dt 
    fOutPtr[3] = fPtr[3]-dt*(0.2886751345948129*hxfRecovR[2]*dy-0.2886751345948129*hxfRecovL[2]*dy) 
    fOutPtr[4] = fPtr[4]-dt*(0.5*hxfRecovR[2]*dy+0.5*hxfRecovL[2]*dy) 
end

fOutPtr[1] = fOutPtr[1] 
fOutPtr[2] = dt*((1.732050807568877*hRecovR[2]*fPtr[4]*dy)/dx^2+(1.732050807568877*hRecovL[2]*fPtr[4]*dy)/dx^2+(hRecovR[2]*fPtr[3]*dy)/dx^2-(1.0*hRecovL[2]*fPtr[3]*dy)/dx^2-(6.0*hPtr[3]*fPtr[4])/dx^2+(6.0*hRecovR[1]*fPtr[2])/dx^2+(6.0*hRecovL[1]*fPtr[2])/dx^2-(6.0*hPtr[1]*fPtr[2])/dx^2+(3.464101615137754*fPtr[1]*hRecovR[1])/dx^2-(3.464101615137754*fPtr[1]*hRecovL[1])/dx^2)+fOutPtr[2] 
fOutPtr[3] = fOutPtr[3] 
fOutPtr[4] = dt*((1.732050807568877*fPtr[2]*hRecovR[2]*dy)/dx^2+(fPtr[1]*hRecovR[2]*dy)/dx^2+(1.732050807568877*fPtr[2]*hRecovL[2]*dy)/dx^2-(1.0*fPtr[1]*hRecovL[2]*dy)/dx^2+(6.0*hRecovR[1]*fPtr[4])/dx^2+(6.0*hRecovL[1]*fPtr[4])/dx^2-(6.0*hPtr[1]*fPtr[4])/dx^2-(6.0*fPtr[2]*hPtr[3])/dx^2+(3.464101615137754*hRecovR[1]*fPtr[3])/dx^2-(3.464101615137754*hRecovL[1]*fPtr[3])/dx^2)+fOutPtr[4] 

   local hRecovT = {} 
hRecovT[1] = 0.08333333333333333*(1.732050807568877*(2.0*hPtr[3]-2.0*hTPtr[3])+3.0*hTPtr[1]+3.0*hPtr[1]) 
hRecovT[2] = -(0.125*(1.732050807568877*(5.0*hTPtr[3]+5.0*hPtr[3])-9.0*hTPtr[1]+9.0*hPtr[1]))/dy 
hRecovT[3] = -(0.8660254037844386*(hPtr[3]-1.0*hTPtr[3]))/dy^2 
hRecovT[4] = (0.25*(1.732050807568877*(5.0*hTPtr[3]+5.0*hPtr[3])-5.0*hTPtr[1]+5.0*hPtr[1]))/dy^3 
hRecovT[5] = (0.5*((-2.0*hTPtr[4])+2.0*hPtr[4]+1.732050807568877*hTPtr[2]+1.732050807568877*hPtr[2]))/dx 
hRecovT[6] = -(0.25*(15.0*hTPtr[4]+15.0*hPtr[4]-15.58845726811989*hTPtr[2]+15.58845726811989*hPtr[2]))/(dx*dy) 
hRecovT[7] = -(1.0*(3.0*hPtr[4]-3.0*hTPtr[4]))/(dx*dy^2) 
hRecovT[8] = (0.5*(15.0*hTPtr[4]+15.0*hPtr[4]-8.660254037844386*hTPtr[2]+8.660254037844386*hPtr[2]))/(dx*dy^3) 

local hRecovB = {} 
hRecovB[1] = 0.08333333333333333*(1.732050807568877*(2.0*hBPtr[3]-2.0*hPtr[3])+3.0*hPtr[1]+3.0*hBPtr[1]) 
hRecovB[2] = -(0.125*(1.732050807568877*(5.0*hPtr[3]+5.0*hBPtr[3])-9.0*hPtr[1]+9.0*hBPtr[1]))/dy 
hRecovB[3] = -(0.8660254037844386*(hBPtr[3]-1.0*hPtr[3]))/dy^2 
hRecovB[4] = (0.25*(1.732050807568877*(5.0*hPtr[3]+5.0*hBPtr[3])-5.0*hPtr[1]+5.0*hBPtr[1]))/dy^3 
hRecovB[5] = (0.5*((-2.0*hPtr[4])+2.0*hBPtr[4]+1.732050807568877*hPtr[2]+1.732050807568877*hBPtr[2]))/dx 
hRecovB[6] = -(0.25*(15.0*hPtr[4]+15.0*hBPtr[4]-15.58845726811989*hPtr[2]+15.58845726811989*hBPtr[2]))/(dx*dy) 
hRecovB[7] = -(1.0*(3.0*hBPtr[4]-3.0*hPtr[4]))/(dx*dy^2) 
hRecovB[8] = (0.5*(15.0*hPtr[4]+15.0*hBPtr[4]-8.660254037844386*hPtr[2]+8.660254037844386*hBPtr[2]))/(dx*dy^3) 

local hyfRecovT = {} 
hyfRecovT[1] = -0.001388888888888889*(((27.0*fTPtr[4]-27.0*fPtr[4])*hRecovT[8]*dx+(93.53074360871933*fTPtr[3]-93.53074360871933*fPtr[3])*hRecovT[4])*dy^2+((30.0*fTPtr[4]+30.0*fPtr[4]-17.32050807568877*fTPtr[2]+17.32050807568877*fPtr[2])*hRecovT[7]*dx+(103.9230484541326*fTPtr[3]+103.9230484541326*fPtr[3]-60.0*fTPtr[1]+60.0*fPtr[1])*hRecovT[3])*dy+(60.0*fTPtr[4]-60.0*fPtr[4]-51.96152422706631*fTPtr[2]-51.96152422706631*fPtr[2])*hRecovT[6]*dx+207.8460969082653*hRecovT[2]*fTPtr[3]+hRecovT[2]*((-207.8460969082653*fPtr[3])-180.0*fTPtr[1]-180.0*fPtr[1])) 
hyfRecovT[2] = -(0.01041666666666667*(((9.0*fTPtr[4]+9.0*fPtr[4]-5.196152422706631*fTPtr[2]+5.196152422706631*fPtr[2])*hRecovT[8]*dx+(31.17691453623978*fTPtr[3]+31.17691453623978*fPtr[3]-18.0*fTPtr[1]+18.0*fPtr[1])*hRecovT[4])*dy^2+((12.0*fTPtr[4]-12.0*fPtr[4]-13.85640646055102*fTPtr[2]-13.85640646055102*fPtr[2])*hRecovT[7]*dx+(41.56921938165305*fTPtr[3]-41.56921938165305*fPtr[3]-48.0*fTPtr[1]-48.0*fPtr[1])*hRecovT[3])*dy+(30.0*fTPtr[4]+30.0*fPtr[4]-31.17691453623978*fTPtr[2]+31.17691453623978*fPtr[2])*hRecovT[6]*dx+103.9230484541326*hRecovT[2]*fTPtr[3]+hRecovT[2]*(103.9230484541326*fPtr[3]-108.0*fTPtr[1]+108.0*fPtr[1])))/dy 
hyfRecovT[3] = (0.008333333333333333*(((36.0*fTPtr[4]-36.0*fPtr[4]+25.98076211353316*fTPtr[2]+25.98076211353316*fPtr[2])*hRecovT[8]*dx+(124.7076581449591*fTPtr[3]-124.7076581449591*fPtr[3]+90.0*fTPtr[1]+90.0*fPtr[1])*hRecovT[4])*dy^2+((30.0*fTPtr[4]+30.0*fPtr[4]+17.32050807568877*fTPtr[2]-17.32050807568877*fPtr[2])*hRecovT[7]*dx+(103.9230484541326*fTPtr[3]+103.9230484541326*fPtr[3]+60.0*fTPtr[1]-60.0*fPtr[1])*hRecovT[3])*dy+(30.0*fTPtr[4]-30.0*fPtr[4])*hRecovT[6]*dx+103.9230484541326*hRecovT[2]*fTPtr[3]-103.9230484541326*hRecovT[2]*fPtr[3]))/dy^2 
hyfRecovT[4] = (0.02083333333333333*(((21.0*fTPtr[4]+21.0*fPtr[4]+8.660254037844386*fTPtr[2]-8.660254037844386*fPtr[2])*hRecovT[8]*dx+(72.74613391789283*fTPtr[3]+72.74613391789283*fPtr[3]+30.0*fTPtr[1]-30.0*fPtr[1])*hRecovT[4])*dy^2+((20.0*fTPtr[4]-20.0*fPtr[4])*hRecovT[7]*dx+(69.28203230275508*fTPtr[3]-69.28203230275508*fPtr[3])*hRecovT[3])*dy+(30.0*fTPtr[4]+30.0*fPtr[4]-17.32050807568877*fTPtr[2]+17.32050807568877*fPtr[2])*hRecovT[6]*dx+103.9230484541326*hRecovT[2]*fTPtr[3]+hRecovT[2]*(103.9230484541326*fPtr[3]-60.0*fTPtr[1]+60.0*fPtr[1])))/dy^3 
hyfRecovT[5] = -(0.008333333333333333*(((15.58845726811989*fTPtr[3]-15.58845726811989*fPtr[3])*hRecovT[8]*dx+(54.0*fTPtr[4]-54.0*fPtr[4])*hRecovT[4])*dy^2+((17.32050807568877*fTPtr[3]+17.32050807568877*fPtr[3]-10.0*fTPtr[1]+10.0*fPtr[1])*hRecovT[7]*dx+60.0*hRecovT[3]*fTPtr[4]+hRecovT[3]*(60.0*fPtr[4]-34.64101615137754*fTPtr[2]+34.64101615137754*fPtr[2]))*dy+(34.64101615137754*fTPtr[3]-34.64101615137754*fPtr[3]-30.0*fTPtr[1]-30.0*fPtr[1])*hRecovT[6]*dx+120.0*hRecovT[2]*fTPtr[4]+hRecovT[2]*((-120.0*fPtr[4])-103.9230484541326*fTPtr[2]-103.9230484541326*fPtr[2])))/dx 
hyfRecovT[6] = -(0.0625*(((5.196152422706631*fTPtr[3]+5.196152422706631*fPtr[3]-3.0*fTPtr[1]+3.0*fPtr[1])*hRecovT[8]*dx+(18.0*fTPtr[4]+18.0*fPtr[4]-10.39230484541326*fTPtr[2]+10.39230484541326*fPtr[2])*hRecovT[4])*dy^2+((6.928203230275509*fTPtr[3]-6.928203230275509*fPtr[3]-8.0*fTPtr[1]-8.0*fPtr[1])*hRecovT[7]*dx+24.0*hRecovT[3]*fTPtr[4]+hRecovT[3]*((-24.0*fPtr[4])-27.71281292110204*fTPtr[2]-27.71281292110204*fPtr[2]))*dy+(17.32050807568877*fTPtr[3]+17.32050807568877*fPtr[3]-18.0*fTPtr[1]+18.0*fPtr[1])*hRecovT[6]*dx+60.0*hRecovT[2]*fTPtr[4]+hRecovT[2]*(60.0*fPtr[4]-62.35382907247956*fTPtr[2]+62.35382907247956*fPtr[2])))/(dx*dy) 
hyfRecovT[7] = (0.05*(((20.78460969082652*fTPtr[3]-20.78460969082652*fPtr[3]+15.0*fTPtr[1]+15.0*fPtr[1])*hRecovT[8]*dx+(72.0*fTPtr[4]-72.0*fPtr[4]+51.96152422706631*fTPtr[2]+51.96152422706631*fPtr[2])*hRecovT[4])*dy^2+((17.32050807568877*fTPtr[3]+17.32050807568877*fPtr[3]+10.0*fTPtr[1]-10.0*fPtr[1])*hRecovT[7]*dx+60.0*hRecovT[3]*fTPtr[4]+hRecovT[3]*(60.0*fPtr[4]+34.64101615137754*fTPtr[2]-34.64101615137754*fPtr[2]))*dy+(17.32050807568877*fTPtr[3]-17.32050807568877*fPtr[3])*hRecovT[6]*dx+60.0*hRecovT[2]*fTPtr[4]-60.0*hRecovT[2]*fPtr[4]))/(dx*dy^2) 
hyfRecovT[8] = (0.04166666666666666*(((36.37306695894642*fTPtr[3]+36.37306695894642*fPtr[3]+15.0*fTPtr[1]-15.0*fPtr[1])*hRecovT[8]*dx+(126.0*fTPtr[4]+126.0*fPtr[4]+51.96152422706631*fTPtr[2]-51.96152422706631*fPtr[2])*hRecovT[4])*dy^2+((34.64101615137754*fTPtr[3]-34.64101615137754*fPtr[3])*hRecovT[7]*dx+120.0*hRecovT[3]*fTPtr[4]-120.0*hRecovT[3]*fPtr[4])*dy+(51.96152422706631*fTPtr[3]+51.96152422706631*fPtr[3]-30.0*fTPtr[1]+30.0*fPtr[1])*hRecovT[6]*dx+180.0*hRecovT[2]*fTPtr[4]+hRecovT[2]*(180.0*fPtr[4]-103.9230484541326*fTPtr[2]+103.9230484541326*fPtr[2])))/(dx*dy^3) 

local hyfRecovB = {} 
hyfRecovB[1] = -0.001388888888888889*(((27.0*fPtr[4]-27.0*fBPtr[4])*hRecovB[8]*dx+(93.53074360871933*fPtr[3]-93.53074360871933*fBPtr[3])*hRecovB[4])*dy^2+((30.0*fPtr[4]+30.0*fBPtr[4]-17.32050807568877*fPtr[2]+17.32050807568877*fBPtr[2])*hRecovB[7]*dx+(103.9230484541326*fPtr[3]+103.9230484541326*fBPtr[3]-60.0*fPtr[1]+60.0*fBPtr[1])*hRecovB[3])*dy+(60.0*fPtr[4]-60.0*fBPtr[4]-51.96152422706631*fPtr[2]-51.96152422706631*fBPtr[2])*hRecovB[6]*dx+207.8460969082653*hRecovB[2]*fPtr[3]+hRecovB[2]*((-207.8460969082653*fBPtr[3])-180.0*fPtr[1]-180.0*fBPtr[1])) 
hyfRecovB[2] = -(0.01041666666666667*(((9.0*fPtr[4]+9.0*fBPtr[4]-5.196152422706631*fPtr[2]+5.196152422706631*fBPtr[2])*hRecovB[8]*dx+(31.17691453623978*fPtr[3]+31.17691453623978*fBPtr[3]-18.0*fPtr[1]+18.0*fBPtr[1])*hRecovB[4])*dy^2+((12.0*fPtr[4]-12.0*fBPtr[4]-13.85640646055102*fPtr[2]-13.85640646055102*fBPtr[2])*hRecovB[7]*dx+(41.56921938165305*fPtr[3]-41.56921938165305*fBPtr[3]-48.0*fPtr[1]-48.0*fBPtr[1])*hRecovB[3])*dy+(30.0*fPtr[4]+30.0*fBPtr[4]-31.17691453623978*fPtr[2]+31.17691453623978*fBPtr[2])*hRecovB[6]*dx+103.9230484541326*hRecovB[2]*fPtr[3]+hRecovB[2]*(103.9230484541326*fBPtr[3]-108.0*fPtr[1]+108.0*fBPtr[1])))/dy 
hyfRecovB[3] = (0.008333333333333333*(((36.0*fPtr[4]-36.0*fBPtr[4]+25.98076211353316*fPtr[2]+25.98076211353316*fBPtr[2])*hRecovB[8]*dx+(124.7076581449591*fPtr[3]-124.7076581449591*fBPtr[3]+90.0*fPtr[1]+90.0*fBPtr[1])*hRecovB[4])*dy^2+((30.0*fPtr[4]+30.0*fBPtr[4]+17.32050807568877*fPtr[2]-17.32050807568877*fBPtr[2])*hRecovB[7]*dx+(103.9230484541326*fPtr[3]+103.9230484541326*fBPtr[3]+60.0*fPtr[1]-60.0*fBPtr[1])*hRecovB[3])*dy+(30.0*fPtr[4]-30.0*fBPtr[4])*hRecovB[6]*dx+103.9230484541326*hRecovB[2]*fPtr[3]-103.9230484541326*hRecovB[2]*fBPtr[3]))/dy^2 
hyfRecovB[4] = (0.02083333333333333*(((21.0*fPtr[4]+21.0*fBPtr[4]+8.660254037844386*fPtr[2]-8.660254037844386*fBPtr[2])*hRecovB[8]*dx+(72.74613391789283*fPtr[3]+72.74613391789283*fBPtr[3]+30.0*fPtr[1]-30.0*fBPtr[1])*hRecovB[4])*dy^2+((20.0*fPtr[4]-20.0*fBPtr[4])*hRecovB[7]*dx+(69.28203230275508*fPtr[3]-69.28203230275508*fBPtr[3])*hRecovB[3])*dy+(30.0*fPtr[4]+30.0*fBPtr[4]-17.32050807568877*fPtr[2]+17.32050807568877*fBPtr[2])*hRecovB[6]*dx+103.9230484541326*hRecovB[2]*fPtr[3]+hRecovB[2]*(103.9230484541326*fBPtr[3]-60.0*fPtr[1]+60.0*fBPtr[1])))/dy^3 
hyfRecovB[5] = -(0.008333333333333333*(((15.58845726811989*fPtr[3]-15.58845726811989*fBPtr[3])*hRecovB[8]*dx+(54.0*fPtr[4]-54.0*fBPtr[4])*hRecovB[4])*dy^2+((17.32050807568877*fPtr[3]+17.32050807568877*fBPtr[3]-10.0*fPtr[1]+10.0*fBPtr[1])*hRecovB[7]*dx+60.0*hRecovB[3]*fPtr[4]+hRecovB[3]*(60.0*fBPtr[4]-34.64101615137754*fPtr[2]+34.64101615137754*fBPtr[2]))*dy+(34.64101615137754*fPtr[3]-34.64101615137754*fBPtr[3]-30.0*fPtr[1]-30.0*fBPtr[1])*hRecovB[6]*dx+120.0*hRecovB[2]*fPtr[4]+hRecovB[2]*((-120.0*fBPtr[4])-103.9230484541326*fPtr[2]-103.9230484541326*fBPtr[2])))/dx 
hyfRecovB[6] = -(0.0625*(((5.196152422706631*fPtr[3]+5.196152422706631*fBPtr[3]-3.0*fPtr[1]+3.0*fBPtr[1])*hRecovB[8]*dx+(18.0*fPtr[4]+18.0*fBPtr[4]-10.39230484541326*fPtr[2]+10.39230484541326*fBPtr[2])*hRecovB[4])*dy^2+((6.928203230275509*fPtr[3]-6.928203230275509*fBPtr[3]-8.0*fPtr[1]-8.0*fBPtr[1])*hRecovB[7]*dx+24.0*hRecovB[3]*fPtr[4]+hRecovB[3]*((-24.0*fBPtr[4])-27.71281292110204*fPtr[2]-27.71281292110204*fBPtr[2]))*dy+(17.32050807568877*fPtr[3]+17.32050807568877*fBPtr[3]-18.0*fPtr[1]+18.0*fBPtr[1])*hRecovB[6]*dx+60.0*hRecovB[2]*fPtr[4]+hRecovB[2]*(60.0*fBPtr[4]-62.35382907247956*fPtr[2]+62.35382907247956*fBPtr[2])))/(dx*dy) 
hyfRecovB[7] = (0.05*(((20.78460969082652*fPtr[3]-20.78460969082652*fBPtr[3]+15.0*fPtr[1]+15.0*fBPtr[1])*hRecovB[8]*dx+(72.0*fPtr[4]-72.0*fBPtr[4]+51.96152422706631*fPtr[2]+51.96152422706631*fBPtr[2])*hRecovB[4])*dy^2+((17.32050807568877*fPtr[3]+17.32050807568877*fBPtr[3]+10.0*fPtr[1]-10.0*fBPtr[1])*hRecovB[7]*dx+60.0*hRecovB[3]*fPtr[4]+hRecovB[3]*(60.0*fBPtr[4]+34.64101615137754*fPtr[2]-34.64101615137754*fBPtr[2]))*dy+(17.32050807568877*fPtr[3]-17.32050807568877*fBPtr[3])*hRecovB[6]*dx+60.0*hRecovB[2]*fPtr[4]-60.0*hRecovB[2]*fBPtr[4]))/(dx*dy^2) 
hyfRecovB[8] = (0.04166666666666666*(((36.37306695894642*fPtr[3]+36.37306695894642*fBPtr[3]+15.0*fPtr[1]-15.0*fBPtr[1])*hRecovB[8]*dx+(126.0*fPtr[4]+126.0*fBPtr[4]+51.96152422706631*fPtr[2]-51.96152422706631*fBPtr[2])*hRecovB[4])*dy^2+((34.64101615137754*fPtr[3]-34.64101615137754*fBPtr[3])*hRecovB[7]*dx+120.0*hRecovB[3]*fPtr[4]-120.0*hRecovB[3]*fBPtr[4])*dy+(51.96152422706631*fPtr[3]+51.96152422706631*fBPtr[3]-30.0*fPtr[1]+30.0*fBPtr[1])*hRecovB[6]*dx+180.0*hRecovB[2]*fPtr[4]+hRecovB[2]*(180.0*fBPtr[4]-103.9230484541326*fPtr[2]+103.9230484541326*fBPtr[2])))/(dx*dy^3) 

if isTopEdge then
    fOutPtr[1] = 1.0*hyfRecovB[1]*dt+fOutPtr[1] 
    fOutPtr[2] = 0.2886751345948129*hyfRecovB[5]*dt*dx+fOutPtr[2] 
    fOutPtr[3] = fOutPtr[3]-1.732050807568877*hyfRecovB[1]*dt 
    fOutPtr[4] = fOutPtr[4]-0.5*hyfRecovB[5]*dt*dx 
elseif isBotEdge then
    fOutPtr[1] = fOutPtr[1]-hyfRecovT[1]*dt 
    fOutPtr[2] = fOutPtr[2]-0.2886751345948129*hyfRecovT[5]*dt*dx 
    fOutPtr[3] = fOutPtr[3]-1.732050807568877*hyfRecovT[1]*dt 
    fOutPtr[4] = fOutPtr[4]-0.5*hyfRecovT[5]*dt*dx 
else
    fOutPtr[1] = fOutPtr[1]-(hyfRecovT[1]-1.0*hyfRecovB[1])*dt 
    fOutPtr[2] = fOutPtr[2]-dt*(0.2886751345948129*hyfRecovT[5]*dx-0.2886751345948129*hyfRecovB[5]*dx) 
    fOutPtr[3] = fOutPtr[3]-(1.732050807568877*hyfRecovT[1]+1.732050807568877*hyfRecovB[1])*dt 
    fOutPtr[4] = fOutPtr[4]-dt*(0.5*hyfRecovT[5]*dx+0.5*hyfRecovB[5]*dx) 
end

fOutPtr[1] = fOutPtr[1] 
fOutPtr[2] = fOutPtr[2] 
fOutPtr[3] = dt*((1.732050807568877*fPtr[4]*hRecovT[5]*dx)/dy^2+(fPtr[2]*hRecovT[5]*dx)/dy^2+(1.732050807568877*fPtr[4]*hRecovB[5]*dx)/dy^2-(1.0*fPtr[2]*hRecovB[5]*dx)/dy^2-(6.0*hPtr[2]*fPtr[4])/dy^2+(6.0*hRecovT[1]*fPtr[3])/dy^2+(6.0*hRecovB[1]*fPtr[3])/dy^2-(6.0*hPtr[1]*fPtr[3])/dy^2+(3.464101615137754*fPtr[1]*hRecovT[1])/dy^2-(3.464101615137754*fPtr[1]*hRecovB[1])/dy^2)+fOutPtr[3] 
fOutPtr[4] = dt*((1.732050807568877*fPtr[3]*hRecovT[5]*dx)/dy^2+(fPtr[1]*hRecovT[5]*dx)/dy^2+(1.732050807568877*fPtr[3]*hRecovB[5]*dx)/dy^2-(1.0*fPtr[1]*hRecovB[5]*dx)/dy^2+(6.0*hRecovT[1]*fPtr[4])/dy^2+(6.0*hRecovB[1]*fPtr[4])/dy^2-(6.0*hPtr[1]*fPtr[4])/dy^2-(6.0*hPtr[2]*fPtr[3])/dy^2+(3.464101615137754*hRecovT[1]*fPtr[2])/dy^2-(3.464101615137754*hRecovB[1]*fPtr[2])/dy^2)+fOutPtr[4] 

 


end

return dragStencilFunc
