function diffStencilFunc(dt, dx, dy,
			 fPtr, fLPtr, fRPtr, fTPtr, fBPtr,
			 gPtr, gLPtr, gRPtr, gTPtr, gBPtr,
			 isTopEdge, isBotEdge, isLeftEdge, isRightEdge,
			 fOutPtr) 
   local gRecovL = {} 
   gRecovL[1] = 0.08333333333333333*((-3.464101615137754*gPtr[2])+3.464101615137754*gLPtr[2]+3.0*gPtr[1]+3.0*gLPtr[1]) 
   gRecovL[2] = (0.5*((-2.0*gPtr[4])+2.0*gLPtr[4]+1.732050807568877*(gPtr[3]+gLPtr[3])))/dy 
   gRecovL[3] = -(0.125*(8.660254037844386*gPtr[2]+8.660254037844386*gLPtr[2]-9.0*gPtr[1]+9.0*gLPtr[1]))/dx 
   gRecovL[4] = -(0.25*(15.0*gPtr[4]+15.0*gLPtr[4]+1.732050807568877*(9.0*gLPtr[3]-9.0*gPtr[3])))/(dx*dy) 
   gRecovL[5] = -(0.5*(1.732050807568877*gLPtr[2]-1.732050807568877*gPtr[2]))/dx^2 
   gRecovL[6] = -(1.0*(3.0*gLPtr[4]-3.0*gPtr[4]))/(dx^2*dy) 
   gRecovL[7] = (0.25*(8.660254037844386*gPtr[2]+8.660254037844386*gLPtr[2]-5.0*gPtr[1]+5.0*gLPtr[1]))/dx^3 
   gRecovL[8] = (0.5*(15.0*gPtr[4]+15.0*gLPtr[4]+1.732050807568877*(5.0*gLPtr[3]-5.0*gPtr[3])))/(dx^3*dy) 

   local gRecovR = {} 
   gRecovR[1] = 0.08333333333333333*((-3.464101615137754*gRPtr[2])+3.464101615137754*gPtr[2]+3.0*gRPtr[1]+3.0*gPtr[1]) 
   gRecovR[2] = (0.5*((-2.0*gRPtr[4])+2.0*gPtr[4]+1.732050807568877*(gRPtr[3]+gPtr[3])))/dy 
   gRecovR[3] = -(0.125*(8.660254037844386*gRPtr[2]+8.660254037844386*gPtr[2]-9.0*gRPtr[1]+9.0*gPtr[1]))/dx 
   gRecovR[4] = -(0.25*(15.0*gRPtr[4]+15.0*gPtr[4]+1.732050807568877*(9.0*gPtr[3]-9.0*gRPtr[3])))/(dx*dy) 
   gRecovR[5] = -(0.5*(1.732050807568877*gPtr[2]-1.732050807568877*gRPtr[2]))/dx^2 
   gRecovR[6] = -(1.0*(3.0*gPtr[4]-3.0*gRPtr[4]))/(dx^2*dy) 
   gRecovR[7] = (0.25*(8.660254037844386*gRPtr[2]+8.660254037844386*gPtr[2]-5.0*gRPtr[1]+5.0*gPtr[1]))/dx^3 
   gRecovR[8] = (0.5*(15.0*gRPtr[4]+15.0*gPtr[4]+1.732050807568877*(5.0*gPtr[3]-5.0*gRPtr[3])))/(dx^3*dy) 

   local gxxfRecovL = {} 
   gxxfRecovL[1] = -0.04166666666666666*(((3.0*fPtr[4]+3.0*fLPtr[4]-1.732050807568877*fPtr[3]+1.732050807568877*fLPtr[3])*gRecovL[8]*dx+(4.0*fPtr[4]-4.0*fLPtr[4]-3.464101615137754*fPtr[3]-3.464101615137754*fLPtr[3])*gRecovL[6])*dy+(10.39230484541326*fPtr[2]+10.39230484541326*fLPtr[2]-6.0*fPtr[1]+6.0*fLPtr[1])*gRecovL[7]*dx+(13.85640646055102*fPtr[2]-13.85640646055102*fLPtr[2]-12.0*fPtr[1]-12.0*fLPtr[1])*gRecovL[5]) 
   gxxfRecovL[2] = -(0.08333333333333333*(((5.196152422706631*fPtr[2]+5.196152422706631*fLPtr[2]-3.0*fPtr[1]+3.0*fLPtr[1])*gRecovL[8]*dx+(6.928203230275509*fPtr[2]-6.928203230275509*fLPtr[2]-6.0*fPtr[1]-6.0*fLPtr[1])*gRecovL[6])*dy+(18.0*fPtr[4]+18.0*fLPtr[4]-10.39230484541326*fPtr[3]+10.39230484541326*fLPtr[3])*gRecovL[7]*dx+(24.0*fPtr[4]-24.0*fLPtr[4]-20.78460969082652*fPtr[3]-20.78460969082652*fLPtr[3])*gRecovL[5]))/dy 
   gxxfRecovL[3] = (0.125*((((-3.0*fPtr[4])+3.0*fLPtr[4]+3.464101615137754*fPtr[3]+3.464101615137754*fLPtr[3])*gRecovL[8]*dx+((-5.0*fPtr[4])-5.0*fLPtr[4]+5.196152422706631*fPtr[3]-5.196152422706631*fLPtr[3])*gRecovL[6])*dy+((-10.39230484541326*fPtr[2])+10.39230484541326*fLPtr[2]+12.0*fPtr[1]+12.0*fLPtr[1])*gRecovL[7]*dx+((-17.32050807568877*fPtr[2])-17.32050807568877*fLPtr[2]+18.0*fPtr[1]-18.0*fLPtr[1])*gRecovL[5]))/dx 
   gxxfRecovL[4] = (0.25*((((-5.196152422706631*fPtr[2])+5.196152422706631*fLPtr[2]+6.0*fPtr[1]+6.0*fLPtr[1])*gRecovL[8]*dx+((-8.660254037844386*fPtr[2])-8.660254037844386*fLPtr[2]+9.0*fPtr[1]-9.0*fLPtr[1])*gRecovL[6])*dy+((-18.0*fPtr[4])+18.0*fLPtr[4]+20.78460969082652*fPtr[3]+20.78460969082652*fLPtr[3])*gRecovL[7]*dx+((-30.0*fPtr[4])-30.0*fLPtr[4]+31.17691453623978*fPtr[3]-31.17691453623978*fLPtr[3])*gRecovL[5]))/(dx*dy) 
   
   local gxxfRecovR = {} 
   gxxfRecovR[1] = -0.04166666666666666*(((3.0*fRPtr[4]+3.0*fPtr[4]-1.732050807568877*fRPtr[3]+1.732050807568877*fPtr[3])*gRecovR[8]*dx+(4.0*fRPtr[4]-4.0*fPtr[4]-3.464101615137754*fRPtr[3]-3.464101615137754*fPtr[3])*gRecovR[6])*dy+(10.39230484541326*fRPtr[2]+10.39230484541326*fPtr[2]-6.0*fRPtr[1]+6.0*fPtr[1])*gRecovR[7]*dx+(13.85640646055102*fRPtr[2]-13.85640646055102*fPtr[2]-12.0*fRPtr[1]-12.0*fPtr[1])*gRecovR[5]) 
   gxxfRecovR[2] = -(0.08333333333333333*(((5.196152422706631*fRPtr[2]+5.196152422706631*fPtr[2]-3.0*fRPtr[1]+3.0*fPtr[1])*gRecovR[8]*dx+(6.928203230275509*fRPtr[2]-6.928203230275509*fPtr[2]-6.0*fRPtr[1]-6.0*fPtr[1])*gRecovR[6])*dy+(18.0*fRPtr[4]+18.0*fPtr[4]-10.39230484541326*fRPtr[3]+10.39230484541326*fPtr[3])*gRecovR[7]*dx+(24.0*fRPtr[4]-24.0*fPtr[4]-20.78460969082652*fRPtr[3]-20.78460969082652*fPtr[3])*gRecovR[5]))/dy 
   gxxfRecovR[3] = (0.125*((((-3.0*fRPtr[4])+3.0*fPtr[4]+3.464101615137754*fRPtr[3]+3.464101615137754*fPtr[3])*gRecovR[8]*dx+((-5.0*fRPtr[4])-5.0*fPtr[4]+5.196152422706631*fRPtr[3]-5.196152422706631*fPtr[3])*gRecovR[6])*dy+((-10.39230484541326*fRPtr[2])+10.39230484541326*fPtr[2]+12.0*fRPtr[1]+12.0*fPtr[1])*gRecovR[7]*dx+((-17.32050807568877*fRPtr[2])-17.32050807568877*fPtr[2]+18.0*fRPtr[1]-18.0*fPtr[1])*gRecovR[5]))/dx 
   gxxfRecovR[4] = (0.25*((((-5.196152422706631*fRPtr[2])+5.196152422706631*fPtr[2]+6.0*fRPtr[1]+6.0*fPtr[1])*gRecovR[8]*dx+((-8.660254037844386*fRPtr[2])-8.660254037844386*fPtr[2]+9.0*fRPtr[1]-9.0*fPtr[1])*gRecovR[6])*dy+((-18.0*fRPtr[4])+18.0*fPtr[4]+20.78460969082652*fRPtr[3]+20.78460969082652*fPtr[3])*gRecovR[7]*dx+((-30.0*fRPtr[4])-30.0*fPtr[4]+31.17691453623978*fRPtr[3]-31.17691453623978*fPtr[3])*gRecovR[5]))/(dx*dy) 
   
   if isLeftEdge then
      fOutPtr[1] = 0.5*gxxfRecovR[3]*dt+fOutPtr[1] 
      fOutPtr[2] = 0.5*dt*(1.732050807568877*gxxfRecovR[3]-(3.464101615137754*gxxfRecovR[1])/dx)+fOutPtr[2] 
      fOutPtr[3] = 0.1443375672974065*gxxfRecovR[4]*dt*dy+fOutPtr[3] 
      fOutPtr[4] = 0.5*dt*(0.5*gxxfRecovR[4]*dy-(1.0*gxxfRecovR[2]*dy)/dx)+fOutPtr[4] 
   elseif isRightEdge then
      fOutPtr[1] = fOutPtr[1]-0.5*gxxfRecovL[3]*dt 
      fOutPtr[2] = 0.5*dt*((3.464101615137754*gxxfRecovL[1])/dx+1.732050807568877*gxxfRecovL[3])+fOutPtr[2] 
      fOutPtr[3] = fOutPtr[3]-0.1443375672974065*gxxfRecovL[4]*dt*dy 
    fOutPtr[4] = 0.5*dt*((gxxfRecovL[2]*dy)/dx+0.5*gxxfRecovL[4]*dy)+fOutPtr[4] 
   else
      fOutPtr[1] = 0.5*(gxxfRecovR[3]-1.0*gxxfRecovL[3])*dt+fOutPtr[1] 
      fOutPtr[2] = 0.5*dt*((-(3.464101615137754*gxxfRecovR[1])/dx)+(3.464101615137754*gxxfRecovL[1])/dx+1.732050807568877*gxxfRecovR[3]+1.732050807568877*gxxfRecovL[3])+fOutPtr[2] 
      fOutPtr[3] = 0.5*dt*(0.2886751345948129*gxxfRecovR[4]*dy-0.2886751345948129*gxxfRecovL[4]*dy)+fOutPtr[3] 
      fOutPtr[4] = 0.5*dt*((-(1.0*gxxfRecovR[2]*dy)/dx)+(gxxfRecovL[2]*dy)/dx+0.5*gxxfRecovR[4]*dy+0.5*gxxfRecovL[4]*dy)+fOutPtr[4] 
   end
   
   local gRecovT = {} 
   gRecovT[1] = 0.08333333333333333*(1.732050807568877*(2.0*gPtr[3]-2.0*gTPtr[3])+3.0*gTPtr[1]+3.0*gPtr[1]) 
   gRecovT[2] = -(0.125*(1.732050807568877*(5.0*gTPtr[3]+5.0*gPtr[3])-9.0*gTPtr[1]+9.0*gPtr[1]))/dy 
   gRecovT[3] = -(0.8660254037844386*(gPtr[3]-1.0*gTPtr[3]))/dy^2 
   gRecovT[4] = (0.25*(1.732050807568877*(5.0*gTPtr[3]+5.0*gPtr[3])-5.0*gTPtr[1]+5.0*gPtr[1]))/dy^3 
   gRecovT[5] = (0.5*((-2.0*gTPtr[4])+2.0*gPtr[4]+1.732050807568877*gTPtr[2]+1.732050807568877*gPtr[2]))/dx 
   gRecovT[6] = -(0.25*(15.0*gTPtr[4]+15.0*gPtr[4]-15.58845726811989*gTPtr[2]+15.58845726811989*gPtr[2]))/(dx*dy) 
   gRecovT[7] = -(1.0*(3.0*gPtr[4]-3.0*gTPtr[4]))/(dx*dy^2) 
   gRecovT[8] = (0.5*(15.0*gTPtr[4]+15.0*gPtr[4]-8.660254037844386*gTPtr[2]+8.660254037844386*gPtr[2]))/(dx*dy^3) 
   
   local gRecovB = {} 
   gRecovB[1] = 0.08333333333333333*(1.732050807568877*(2.0*gBPtr[3]-2.0*gPtr[3])+3.0*gPtr[1]+3.0*gBPtr[1]) 
   gRecovB[2] = -(0.125*(1.732050807568877*(5.0*gPtr[3]+5.0*gBPtr[3])-9.0*gPtr[1]+9.0*gBPtr[1]))/dy 
   gRecovB[3] = -(0.8660254037844386*(gBPtr[3]-1.0*gPtr[3]))/dy^2 
   gRecovB[4] = (0.25*(1.732050807568877*(5.0*gPtr[3]+5.0*gBPtr[3])-5.0*gPtr[1]+5.0*gBPtr[1]))/dy^3 
   gRecovB[5] = (0.5*((-2.0*gPtr[4])+2.0*gBPtr[4]+1.732050807568877*gPtr[2]+1.732050807568877*gBPtr[2]))/dx 
   gRecovB[6] = -(0.25*(15.0*gPtr[4]+15.0*gBPtr[4]-15.58845726811989*gPtr[2]+15.58845726811989*gBPtr[2]))/(dx*dy) 
   gRecovB[7] = -(1.0*(3.0*gBPtr[4]-3.0*gPtr[4]))/(dx*dy^2) 
   gRecovB[8] = (0.5*(15.0*gPtr[4]+15.0*gBPtr[4]-8.660254037844386*gPtr[2]+8.660254037844386*gBPtr[2]))/(dx*dy^3) 
   
   local gyyfRecovT = {} 
   gyyfRecovT[1] = -0.04166666666666666*(((3.0*fTPtr[4]+3.0*fPtr[4]-1.732050807568877*fTPtr[2]+1.732050807568877*fPtr[2])*gRecovT[8]*dx+(10.39230484541326*fTPtr[3]+10.39230484541326*fPtr[3]-6.0*fTPtr[1]+6.0*fPtr[1])*gRecovT[4])*dy+(4.0*fTPtr[4]-4.0*fPtr[4]-3.464101615137754*fTPtr[2]-3.464101615137754*fPtr[2])*gRecovT[7]*dx+(13.85640646055102*fTPtr[3]-13.85640646055102*fPtr[3]-12.0*fTPtr[1]-12.0*fPtr[1])*gRecovT[3]) 
   gyyfRecovT[2] = -(0.125*(((3.0*fTPtr[4]-3.0*fPtr[4]-3.464101615137754*fTPtr[2]-3.464101615137754*fPtr[2])*gRecovT[8]*dx+(10.39230484541326*fTPtr[3]-10.39230484541326*fPtr[3]-12.0*fTPtr[1]-12.0*fPtr[1])*gRecovT[4])*dy+(5.0*fTPtr[4]+5.0*fPtr[4]-5.196152422706631*fTPtr[2]+5.196152422706631*fPtr[2])*gRecovT[7]*dx+(17.32050807568877*fTPtr[3]+17.32050807568877*fPtr[3]-18.0*fTPtr[1]+18.0*fPtr[1])*gRecovT[3]))/dy 
   gyyfRecovT[3] = (0.25*(((3.0*fTPtr[4]+3.0*fPtr[4]+1.732050807568877*fTPtr[2]-1.732050807568877*fPtr[2])*gRecovT[8]*dx+(10.39230484541326*fTPtr[3]+10.39230484541326*fPtr[3]+6.0*fTPtr[1]-6.0*fPtr[1])*gRecovT[4])*dy+(2.0*fTPtr[4]-2.0*fPtr[4])*gRecovT[7]*dx+(6.928203230275509*fTPtr[3]-6.928203230275509*fPtr[3])*gRecovT[3]))/dy^2 
   gyyfRecovT[4] = (0.08333333333333333*(((15.0*fTPtr[4]-15.0*fPtr[4])*gRecovT[8]*dx+(51.96152422706631*fTPtr[3]-51.96152422706631*fPtr[3])*gRecovT[4])*dy+(15.0*fTPtr[4]+15.0*fPtr[4]-8.660254037844386*fTPtr[2]+8.660254037844386*fPtr[2])*gRecovT[7]*dx+(51.96152422706631*fTPtr[3]+51.96152422706631*fPtr[3]-30.0*fTPtr[1]+30.0*fPtr[1])*gRecovT[3]))/dy^3 
   gyyfRecovT[5] = -(0.08333333333333333*(((5.196152422706631*fTPtr[3]+5.196152422706631*fPtr[3]-3.0*fTPtr[1]+3.0*fPtr[1])*gRecovT[8]*dx+(18.0*fTPtr[4]+18.0*fPtr[4]-10.39230484541326*fTPtr[2]+10.39230484541326*fPtr[2])*gRecovT[4])*dy+(6.928203230275509*fTPtr[3]-6.928203230275509*fPtr[3]-6.0*fTPtr[1]-6.0*fPtr[1])*gRecovT[7]*dx+24.0*gRecovT[3]*fTPtr[4]+gRecovT[3]*((-24.0*fPtr[4])-20.78460969082652*fTPtr[2]-20.78460969082652*fPtr[2])))/dx 
   gyyfRecovT[6] = -(0.25*(((5.196152422706631*fTPtr[3]-5.196152422706631*fPtr[3]-6.0*fTPtr[1]-6.0*fPtr[1])*gRecovT[8]*dx+(18.0*fTPtr[4]-18.0*fPtr[4]-20.78460969082652*fTPtr[2]-20.78460969082652*fPtr[2])*gRecovT[4])*dy+(8.660254037844386*fTPtr[3]+8.660254037844386*fPtr[3]-9.0*fTPtr[1]+9.0*fPtr[1])*gRecovT[7]*dx+30.0*gRecovT[3]*fTPtr[4]+gRecovT[3]*(30.0*fPtr[4]-31.17691453623978*fTPtr[2]+31.17691453623978*fPtr[2])))/(dx*dy) 
   gyyfRecovT[7] = (0.5*(((5.196152422706631*fTPtr[3]+5.196152422706631*fPtr[3]+3.0*fTPtr[1]-3.0*fPtr[1])*gRecovT[8]*dx+(18.0*fTPtr[4]+18.0*fPtr[4]+10.39230484541326*fTPtr[2]-10.39230484541326*fPtr[2])*gRecovT[4])*dy+(3.464101615137754*fTPtr[3]-3.464101615137754*fPtr[3])*gRecovT[7]*dx+12.0*gRecovT[3]*fTPtr[4]-12.0*gRecovT[3]*fPtr[4]))/(dx*dy^2) 
   gyyfRecovT[8] = (0.5*(((8.660254037844386*fTPtr[3]-8.660254037844386*fPtr[3])*gRecovT[8]*dx+(30.0*fTPtr[4]-30.0*fPtr[4])*gRecovT[4])*dy+(8.660254037844386*fTPtr[3]+8.660254037844386*fPtr[3]-5.0*fTPtr[1]+5.0*fPtr[1])*gRecovT[7]*dx+30.0*gRecovT[3]*fTPtr[4]+gRecovT[3]*(30.0*fPtr[4]-17.32050807568877*fTPtr[2]+17.32050807568877*fPtr[2])))/(dx*dy^3) 
   
   local gyyfRecovB = {} 
   gyyfRecovB[1] = -0.04166666666666666*(((3.0*fPtr[4]+3.0*fBPtr[4]-1.732050807568877*fPtr[2]+1.732050807568877*fBPtr[2])*gRecovB[8]*dx+(10.39230484541326*fPtr[3]+10.39230484541326*fBPtr[3]-6.0*fPtr[1]+6.0*fBPtr[1])*gRecovB[4])*dy+(4.0*fPtr[4]-4.0*fBPtr[4]-3.464101615137754*fPtr[2]-3.464101615137754*fBPtr[2])*gRecovB[7]*dx+(13.85640646055102*fPtr[3]-13.85640646055102*fBPtr[3]-12.0*fPtr[1]-12.0*fBPtr[1])*gRecovB[3]) 
   gyyfRecovB[2] = -(0.125*(((3.0*fPtr[4]-3.0*fBPtr[4]-3.464101615137754*fPtr[2]-3.464101615137754*fBPtr[2])*gRecovB[8]*dx+(10.39230484541326*fPtr[3]-10.39230484541326*fBPtr[3]-12.0*fPtr[1]-12.0*fBPtr[1])*gRecovB[4])*dy+(5.0*fPtr[4]+5.0*fBPtr[4]-5.196152422706631*fPtr[2]+5.196152422706631*fBPtr[2])*gRecovB[7]*dx+(17.32050807568877*fPtr[3]+17.32050807568877*fBPtr[3]-18.0*fPtr[1]+18.0*fBPtr[1])*gRecovB[3]))/dy 
   gyyfRecovB[3] = (0.25*(((3.0*fPtr[4]+3.0*fBPtr[4]+1.732050807568877*fPtr[2]-1.732050807568877*fBPtr[2])*gRecovB[8]*dx+(10.39230484541326*fPtr[3]+10.39230484541326*fBPtr[3]+6.0*fPtr[1]-6.0*fBPtr[1])*gRecovB[4])*dy+(2.0*fPtr[4]-2.0*fBPtr[4])*gRecovB[7]*dx+(6.928203230275509*fPtr[3]-6.928203230275509*fBPtr[3])*gRecovB[3]))/dy^2 
   gyyfRecovB[4] = (0.08333333333333333*(((15.0*fPtr[4]-15.0*fBPtr[4])*gRecovB[8]*dx+(51.96152422706631*fPtr[3]-51.96152422706631*fBPtr[3])*gRecovB[4])*dy+(15.0*fPtr[4]+15.0*fBPtr[4]-8.660254037844386*fPtr[2]+8.660254037844386*fBPtr[2])*gRecovB[7]*dx+(51.96152422706631*fPtr[3]+51.96152422706631*fBPtr[3]-30.0*fPtr[1]+30.0*fBPtr[1])*gRecovB[3]))/dy^3 
   gyyfRecovB[5] = -(0.08333333333333333*(((5.196152422706631*fPtr[3]+5.196152422706631*fBPtr[3]-3.0*fPtr[1]+3.0*fBPtr[1])*gRecovB[8]*dx+(18.0*fPtr[4]+18.0*fBPtr[4]-10.39230484541326*fPtr[2]+10.39230484541326*fBPtr[2])*gRecovB[4])*dy+(6.928203230275509*fPtr[3]-6.928203230275509*fBPtr[3]-6.0*fPtr[1]-6.0*fBPtr[1])*gRecovB[7]*dx+24.0*gRecovB[3]*fPtr[4]+gRecovB[3]*((-24.0*fBPtr[4])-20.78460969082652*fPtr[2]-20.78460969082652*fBPtr[2])))/dx 
   gyyfRecovB[6] = -(0.25*(((5.196152422706631*fPtr[3]-5.196152422706631*fBPtr[3]-6.0*fPtr[1]-6.0*fBPtr[1])*gRecovB[8]*dx+(18.0*fPtr[4]-18.0*fBPtr[4]-20.78460969082652*fPtr[2]-20.78460969082652*fBPtr[2])*gRecovB[4])*dy+(8.660254037844386*fPtr[3]+8.660254037844386*fBPtr[3]-9.0*fPtr[1]+9.0*fBPtr[1])*gRecovB[7]*dx+30.0*gRecovB[3]*fPtr[4]+gRecovB[3]*(30.0*fBPtr[4]-31.17691453623978*fPtr[2]+31.17691453623978*fBPtr[2])))/(dx*dy) 
   gyyfRecovB[7] = (0.5*(((5.196152422706631*fPtr[3]+5.196152422706631*fBPtr[3]+3.0*fPtr[1]-3.0*fBPtr[1])*gRecovB[8]*dx+(18.0*fPtr[4]+18.0*fBPtr[4]+10.39230484541326*fPtr[2]-10.39230484541326*fBPtr[2])*gRecovB[4])*dy+(3.464101615137754*fPtr[3]-3.464101615137754*fBPtr[3])*gRecovB[7]*dx+12.0*gRecovB[3]*fPtr[4]-12.0*gRecovB[3]*fBPtr[4]))/(dx*dy^2) 
   gyyfRecovB[8] = (0.5*(((8.660254037844386*fPtr[3]-8.660254037844386*fBPtr[3])*gRecovB[8]*dx+(30.0*fPtr[4]-30.0*fBPtr[4])*gRecovB[4])*dy+(8.660254037844386*fPtr[3]+8.660254037844386*fBPtr[3]-5.0*fPtr[1]+5.0*fBPtr[1])*gRecovB[7]*dx+30.0*gRecovB[3]*fPtr[4]+gRecovB[3]*(30.0*fBPtr[4]-17.32050807568877*fPtr[2]+17.32050807568877*fBPtr[2])))/(dx*dy^3) 
   
   if isTopEdge then
      fOutPtr[1] = fOutPtr[1]-0.5*gyyfRecovB[2]*dt 
      fOutPtr[2] = fOutPtr[2]-0.1443375672974065*gyyfRecovB[6]*dt*dx 
      fOutPtr[3] = 0.5*dt*((3.464101615137754*gyyfRecovB[1])/dy+1.732050807568877*gyyfRecovB[2])+fOutPtr[3] 
      fOutPtr[4] = 0.5*dt*((gyyfRecovB[5]*dx)/dy+0.5*gyyfRecovB[6]*dx)+fOutPtr[4] 
   elseif isBotEdge then
      fOutPtr[1] = 0.5*gyyfRecovT[2]*dt+fOutPtr[1] 
      fOutPtr[2] = 0.1443375672974065*gyyfRecovT[6]*dt*dx+fOutPtr[2] 
      fOutPtr[3] = 0.5*dt*(1.732050807568877*gyyfRecovT[2]-(3.464101615137754*gyyfRecovT[1])/dy)+fOutPtr[3] 
      fOutPtr[4] = 0.5*dt*(0.5*gyyfRecovT[6]*dx-(1.0*gyyfRecovT[5]*dx)/dy)+fOutPtr[4] 
   else
      fOutPtr[1] = 0.5*(gyyfRecovT[2]-1.0*gyyfRecovB[2])*dt+fOutPtr[1] 
      fOutPtr[2] = 0.5*dt*(0.2886751345948129*gyyfRecovT[6]*dx-0.2886751345948129*gyyfRecovB[6]*dx)+fOutPtr[2] 
      fOutPtr[3] = 0.5*dt*((-(3.464101615137754*gyyfRecovT[1])/dy)+(3.464101615137754*gyyfRecovB[1])/dy+1.732050807568877*gyyfRecovT[2]+1.732050807568877*gyyfRecovB[2])+fOutPtr[3] 
      fOutPtr[4] = 0.5*dt*((-(1.0*gyyfRecovT[5]*dx)/dy)+(gyyfRecovB[5]*dx)/dy+0.5*gyyfRecovT[6]*dx+0.5*gyyfRecovB[6]*dx)+fOutPtr[4] 
   end   
end

return diffStencilFunc
   
