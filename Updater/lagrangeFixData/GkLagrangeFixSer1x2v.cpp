#include <math.h> 
#include <GkLagrangeFixDecl.h> 

void GkLagrangeFixSer1x2v1p(double *dm0, double *dm1, double *dm2, double *B, double mass, double *L, double *Nv, double *vc, double *f) {
  f[0] = f[0] +  -1.0*((((60.0*dm0[0]*gkyl_ipow(L[0],8)-720.0*dm2[0]*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],8)+(720.0*dm2[0]*gkyl_ipow(L[0],6)-60.0*dm0[0]*gkyl_ipow(L[0],8))*gkyl_ipow(Nv[0],4))*gkyl_ipow(vc[0],2)+vc[0]*((-48.0*dm1[0]*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],8))+96.0*dm1[0]*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],4)-48.0*dm1[0]*gkyl_ipow(L[0],8))+(60.0*dm2[0]*gkyl_ipow(L[0],8)-9.0*dm0[0]*gkyl_ipow(L[0],10))*gkyl_ipow(Nv[0],8)+(5.0*dm0[0]*gkyl_ipow(L[0],10)-60.0*dm2[0]*gkyl_ipow(L[0],8))*gkyl_ipow(Nv[0],6)+(13.0*dm0[0]*gkyl_ipow(L[0],10)-60.0*dm2[0]*gkyl_ipow(L[0],8))*gkyl_ipow(Nv[0],4)+(60.0*dm2[0]*gkyl_ipow(L[0],8)-5.0*dm0[0]*gkyl_ipow(L[0],10))*gkyl_ipow(Nv[0],2)-4.0*dm0[0]*gkyl_ipow(L[0],10))*gkyl_ipow(mass,5)+((L[1]*(509.1168824543146*B[1]*dm0[1]+509.1168824543146*B[0]*dm0[0])*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],8)+L[1]*((-509.1168824543146*B[1]*dm0[1])-509.1168824543146*B[0]*dm0[0])*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],4))*gkyl_ipow(vc[0],2)+vc[1]*(((84.85281374238575*B[1]*dm0[1]+84.85281374238575*B[0]*dm0[0])*gkyl_ipow(L[0],8)+((-1018.23376490863*B[1]*dm2[1])-1018.23376490863*B[0]*dm2[0])*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],8)+(((-84.85281374238575*B[1]*dm0[1])-84.85281374238575*B[0]*dm0[0])*gkyl_ipow(L[0],8)+(1018.23376490863*B[1]*dm2[1]+1018.23376490863*B[0]*dm2[0])*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],4))+L[1]*(((-84.85281374238575*B[1]*dm0[1])-84.85281374238575*B[0]*dm0[0])*gkyl_ipow(L[0],8)+(509.1168824543146*B[1]*dm2[1]+509.1168824543146*B[0]*dm2[0])*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],8)+L[1]*(42.42640687119286*B[1]*dm0[1]+42.42640687119286*B[0]*dm0[0])*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],6)+L[1]*((84.85281374238575*B[1]*dm0[1]+84.85281374238575*B[0]*dm0[0])*gkyl_ipow(L[0],8)+((-509.1168824543146*B[1]*dm2[1])-509.1168824543146*B[0]*dm2[0])*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],4)+L[1]*((-42.42640687119286*B[1]*dm0[1])-42.42640687119286*B[0]*dm0[0])*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],2))*gkyl_ipow(mass,4)+(vc[0]*(gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],4)*(2880.0*dm1[0]*gkyl_ipow(B[1],2)+2880.0*dm1[0]*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2)+gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],8)*((-2880.0*dm1[0]*gkyl_ipow(B[1],2))-2880.0*dm1[0]*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2))+gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],4)*(1800.0*dm2[0]*gkyl_ipow(B[1],2)+1800.0*dm2[0]*gkyl_ipow(B[0],2)-3600.0*B[0]*B[1]*dm2[1])+gkyl_ipow(L[0],6)*((-750.0*dm0[0]*gkyl_ipow(B[1],2))-750.0*dm0[0]*gkyl_ipow(B[0],2)-420.0*B[0]*B[1]*dm0[1]))*gkyl_ipow(L[1],2)+gkyl_ipow(Nv[0],6)*(gkyl_ipow(L[0],4)*((-1800.0*dm2[0]*gkyl_ipow(B[1],2))-1800.0*dm2[0]*gkyl_ipow(B[0],2)+3600.0*B[0]*B[1]*dm2[1])+gkyl_ipow(L[0],6)*(150.0*dm0[0]*gkyl_ipow(B[1],2)+150.0*dm0[0]*gkyl_ipow(B[0],2)-300.0*B[0]*B[1]*dm0[1]))*gkyl_ipow(L[1],2)+gkyl_ipow(Nv[0],8)*gkyl_ipow(vc[0],2)*(gkyl_ipow(L[0],2)*((-21600.0*dm2[0]*gkyl_ipow(B[1],2))-21600.0*dm2[0]*gkyl_ipow(B[0],2)+43200.0*B[0]*B[1]*dm2[1])+gkyl_ipow(L[0],4)*(1800.0*dm0[0]*gkyl_ipow(B[1],2)+1800.0*dm0[0]*gkyl_ipow(B[0],2)-3600.0*B[0]*B[1]*dm0[1]))*gkyl_ipow(L[1],2)+gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],4)*(600.0*dm0[0]*gkyl_ipow(B[1],2)+600.0*dm0[0]*gkyl_ipow(B[0],2)+720.0*B[0]*B[1]*dm0[1])*gkyl_ipow(L[1],2)+vc[1]*(L[1]*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],8)*(720.0*dm0[0]*gkyl_ipow(B[1],2)+720.0*dm0[0]*gkyl_ipow(B[0],2)+1440.0*B[0]*B[1]*dm0[1])+L[1]*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],4)*((-720.0*dm0[0]*gkyl_ipow(B[1],2))-720.0*dm0[0]*gkyl_ipow(B[0],2)-1440.0*B[0]*B[1]*dm0[1])))*gkyl_ipow(mass,3)+(gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],2)*(dm2[1]*(15273.50647362944*gkyl_ipow(B[1],3)-15273.50647362944*B[1]*gkyl_ipow(B[0],2))-15273.50647362944*B[0]*dm2[0]*gkyl_ipow(B[1],2)+15273.50647362944*dm2[0]*gkyl_ipow(B[0],3))+gkyl_ipow(L[0],4)*(dm0[1]*(2545.584412271573*B[1]*gkyl_ipow(B[0],2)-2545.584412271573*gkyl_ipow(B[1],3))+2545.584412271573*B[0]*dm0[0]*gkyl_ipow(B[1],2)-2545.584412271573*dm0[0]*gkyl_ipow(B[0],3)))*gkyl_ipow(L[1],3)+gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],8)*gkyl_ipow(vc[0],2)*(dm0[1]*(15273.50647362944*gkyl_ipow(B[1],3)-15273.50647362944*B[1]*gkyl_ipow(B[0],2))-15273.50647362944*B[0]*dm0[0]*gkyl_ipow(B[1],2)+15273.50647362944*dm0[0]*gkyl_ipow(B[0],3))*gkyl_ipow(L[1],3)+gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],6)*(dm0[1]*(1272.792206135786*gkyl_ipow(B[1],3)-1272.792206135786*B[1]*gkyl_ipow(B[0],2))-1272.792206135786*B[0]*dm0[0]*gkyl_ipow(B[1],2)+1272.792206135786*dm0[0]*gkyl_ipow(B[0],3))*gkyl_ipow(L[1],3)+vc[1]*gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],4)*(dm0[1]*(2545.584412271573*gkyl_ipow(B[1],3)-2545.584412271573*B[1]*gkyl_ipow(B[0],2))-2545.584412271573*B[0]*dm0[0]*gkyl_ipow(B[1],2)+2545.584412271573*dm0[0]*gkyl_ipow(B[0],3))+gkyl_ipow(L[0],2)*(dm2[1]*(30547.01294725889*B[1]*gkyl_ipow(B[0],2)-30547.01294725889*gkyl_ipow(B[1],3))+30547.01294725889*B[0]*dm2[0]*gkyl_ipow(B[1],2)-30547.01294725889*dm2[0]*gkyl_ipow(B[0],3)))*gkyl_ipow(L[1],2))*gkyl_ipow(mass,2)+(vc[0]*gkyl_ipow(Nv[0],8)*((-43200.0*dm1[0]*gkyl_ipow(B[1],4))+86400.0*dm1[0]*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)-43200.0*dm1[0]*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],4)+gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],8)*((-14400.0*dm0[0]*gkyl_ipow(B[1],4))+28800.0*dm0[0]*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)-14400.0*dm0[0]*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],4)+vc[1]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],8)*(21600.0*dm0[0]*gkyl_ipow(B[1],4)-43200.0*dm0[0]*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+21600.0*dm0[0]*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],3))*mass)*gkyl_ipow((12.56637061435917*L[1]*gkyl_ipow(L[0],11)*gkyl_ipow(Nv[0],8)-25.13274122871834*L[1]*gkyl_ipow(L[0],11)*gkyl_ipow(Nv[0],4)+12.56637061435917*L[1]*gkyl_ipow(L[0],11))*gkyl_ipow(mass,4)+(gkyl_ipow(L[0],7)*gkyl_ipow(Nv[0],8)*(753.9822368615503*gkyl_ipow(B[1],2)+753.9822368615503*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3)+gkyl_ipow(L[0],7)*gkyl_ipow(Nv[0],4)*((-753.9822368615503*gkyl_ipow(B[1],2))-753.9822368615503*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3))*gkyl_ipow(mass,2)+gkyl_ipow(L[0],3)*gkyl_ipow(Nv[0],8)*(11309.73355292326*gkyl_ipow(B[1],4)-22619.46710584651*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+11309.73355292326*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],5),-1);

  f[1] = f[1] +  -1.0*((((60.0*dm0[1]*gkyl_ipow(L[0],8)-720.0*dm2[1]*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],8)+(720.0*dm2[1]*gkyl_ipow(L[0],6)-60.0*dm0[1]*gkyl_ipow(L[0],8))*gkyl_ipow(Nv[0],4))*gkyl_ipow(vc[0],2)+vc[0]*((-48.0*dm1[1]*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],8))+96.0*dm1[1]*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],4)-48.0*dm1[1]*gkyl_ipow(L[0],8))+(60.0*dm2[1]*gkyl_ipow(L[0],8)-9.0*dm0[1]*gkyl_ipow(L[0],10))*gkyl_ipow(Nv[0],8)+(5.0*dm0[1]*gkyl_ipow(L[0],10)-60.0*dm2[1]*gkyl_ipow(L[0],8))*gkyl_ipow(Nv[0],6)+(13.0*dm0[1]*gkyl_ipow(L[0],10)-60.0*dm2[1]*gkyl_ipow(L[0],8))*gkyl_ipow(Nv[0],4)+(60.0*dm2[1]*gkyl_ipow(L[0],8)-5.0*dm0[1]*gkyl_ipow(L[0],10))*gkyl_ipow(Nv[0],2)-4.0*dm0[1]*gkyl_ipow(L[0],10))*gkyl_ipow(mass,5)+((L[1]*(509.1168824543146*B[0]*dm0[1]+509.1168824543146*dm0[0]*B[1])*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],8)+L[1]*((-509.1168824543146*B[0]*dm0[1])-509.1168824543146*dm0[0]*B[1])*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],4))*gkyl_ipow(vc[0],2)+vc[1]*(((84.85281374238575*B[0]*dm0[1]+84.85281374238575*dm0[0]*B[1])*gkyl_ipow(L[0],8)+((-1018.23376490863*B[0]*dm2[1])-1018.23376490863*dm2[0]*B[1])*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],8)+(((-84.85281374238575*B[0]*dm0[1])-84.85281374238575*dm0[0]*B[1])*gkyl_ipow(L[0],8)+(1018.23376490863*B[0]*dm2[1]+1018.23376490863*dm2[0]*B[1])*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],4))+L[1]*(((-84.85281374238575*B[0]*dm0[1])-84.85281374238575*dm0[0]*B[1])*gkyl_ipow(L[0],8)+(509.1168824543146*B[0]*dm2[1]+509.1168824543146*dm2[0]*B[1])*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],8)+L[1]*(42.42640687119286*B[0]*dm0[1]+42.42640687119286*dm0[0]*B[1])*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],6)+L[1]*((84.85281374238575*B[0]*dm0[1]+84.85281374238575*dm0[0]*B[1])*gkyl_ipow(L[0],8)+((-509.1168824543146*B[0]*dm2[1])-509.1168824543146*dm2[0]*B[1])*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],4)+L[1]*((-42.42640687119286*B[0]*dm0[1])-42.42640687119286*dm0[0]*B[1])*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],2))*gkyl_ipow(mass,4)+(vc[0]*(dm1[1]*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],4)*(2880.0*gkyl_ipow(B[1],2)+2880.0*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2)+dm1[1]*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],8)*((-2880.0*gkyl_ipow(B[1],2))-2880.0*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2))+gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],4)*(dm2[1]*(1800.0*gkyl_ipow(B[1],2)+1800.0*gkyl_ipow(B[0],2))-3600.0*B[0]*dm2[0]*B[1])+gkyl_ipow(L[0],6)*(dm0[1]*((-750.0*gkyl_ipow(B[1],2))-750.0*gkyl_ipow(B[0],2))-420.0*B[0]*dm0[0]*B[1]))*gkyl_ipow(L[1],2)+gkyl_ipow(Nv[0],8)*gkyl_ipow(vc[0],2)*(gkyl_ipow(L[0],4)*(dm0[1]*(1800.0*gkyl_ipow(B[1],2)+1800.0*gkyl_ipow(B[0],2))-3600.0*B[0]*dm0[0]*B[1])+gkyl_ipow(L[0],2)*(dm2[1]*((-21600.0*gkyl_ipow(B[1],2))-21600.0*gkyl_ipow(B[0],2))+43200.0*B[0]*dm2[0]*B[1]))*gkyl_ipow(L[1],2)+gkyl_ipow(Nv[0],6)*(gkyl_ipow(L[0],6)*(dm0[1]*(150.0*gkyl_ipow(B[1],2)+150.0*gkyl_ipow(B[0],2))-300.0*B[0]*dm0[0]*B[1])+gkyl_ipow(L[0],4)*(dm2[1]*((-1800.0*gkyl_ipow(B[1],2))-1800.0*gkyl_ipow(B[0],2))+3600.0*B[0]*dm2[0]*B[1]))*gkyl_ipow(L[1],2)+gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],4)*(dm0[1]*(600.0*gkyl_ipow(B[1],2)+600.0*gkyl_ipow(B[0],2))+720.0*B[0]*dm0[0]*B[1])*gkyl_ipow(L[1],2)+vc[1]*(L[1]*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],8)*(dm0[1]*(720.0*gkyl_ipow(B[1],2)+720.0*gkyl_ipow(B[0],2))+1440.0*B[0]*dm0[0]*B[1])+L[1]*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],4)*(dm0[1]*((-720.0*gkyl_ipow(B[1],2))-720.0*gkyl_ipow(B[0],2))-1440.0*B[0]*dm0[0]*B[1])))*gkyl_ipow(mass,3)+(gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],2)*(15273.50647362944*dm2[0]*gkyl_ipow(B[1],3)+dm2[1]*(15273.50647362944*gkyl_ipow(B[0],3)-15273.50647362944*B[0]*gkyl_ipow(B[1],2))-15273.50647362944*dm2[0]*B[1]*gkyl_ipow(B[0],2))+gkyl_ipow(L[0],4)*((-2545.584412271573*dm0[0]*gkyl_ipow(B[1],3))+dm0[1]*(2545.584412271573*B[0]*gkyl_ipow(B[1],2)-2545.584412271573*gkyl_ipow(B[0],3))+2545.584412271573*dm0[0]*B[1]*gkyl_ipow(B[0],2)))*gkyl_ipow(L[1],3)+gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],8)*gkyl_ipow(vc[0],2)*(15273.50647362944*dm0[0]*gkyl_ipow(B[1],3)+dm0[1]*(15273.50647362944*gkyl_ipow(B[0],3)-15273.50647362944*B[0]*gkyl_ipow(B[1],2))-15273.50647362944*dm0[0]*B[1]*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3)+gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],6)*(1272.792206135786*dm0[0]*gkyl_ipow(B[1],3)+dm0[1]*(1272.792206135786*gkyl_ipow(B[0],3)-1272.792206135786*B[0]*gkyl_ipow(B[1],2))-1272.792206135786*dm0[0]*B[1]*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3)+vc[1]*gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],2)*((-30547.01294725889*dm2[0]*gkyl_ipow(B[1],3))+dm2[1]*(30547.01294725889*B[0]*gkyl_ipow(B[1],2)-30547.01294725889*gkyl_ipow(B[0],3))+30547.01294725889*dm2[0]*B[1]*gkyl_ipow(B[0],2))+gkyl_ipow(L[0],4)*(2545.584412271573*dm0[0]*gkyl_ipow(B[1],3)+dm0[1]*(2545.584412271573*gkyl_ipow(B[0],3)-2545.584412271573*B[0]*gkyl_ipow(B[1],2))-2545.584412271573*dm0[0]*B[1]*gkyl_ipow(B[0],2)))*gkyl_ipow(L[1],2))*gkyl_ipow(mass,2)+(dm0[1]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],8)*((-14400.0*gkyl_ipow(B[1],4))+28800.0*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)-14400.0*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],4)+vc[0]*dm1[1]*gkyl_ipow(Nv[0],8)*((-43200.0*gkyl_ipow(B[1],4))+86400.0*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)-43200.0*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],4)+dm0[1]*vc[1]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],8)*(21600.0*gkyl_ipow(B[1],4)-43200.0*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+21600.0*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],3))*mass)*gkyl_ipow((12.56637061435917*L[1]*gkyl_ipow(L[0],11)*gkyl_ipow(Nv[0],8)-25.13274122871834*L[1]*gkyl_ipow(L[0],11)*gkyl_ipow(Nv[0],4)+12.56637061435917*L[1]*gkyl_ipow(L[0],11))*gkyl_ipow(mass,4)+(gkyl_ipow(L[0],7)*gkyl_ipow(Nv[0],8)*(753.9822368615503*gkyl_ipow(B[1],2)+753.9822368615503*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3)+gkyl_ipow(L[0],7)*gkyl_ipow(Nv[0],4)*((-753.9822368615503*gkyl_ipow(B[1],2))-753.9822368615503*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3))*gkyl_ipow(mass,2)+gkyl_ipow(L[0],3)*gkyl_ipow(Nv[0],8)*(11309.73355292326*gkyl_ipow(B[1],4)-22619.46710584651*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+11309.73355292326*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],5),-1);

  f[2] = f[2] +  -1.0*((vc[0]*((15.0*dm0[0]*gkyl_ipow(L[0],8)-180.0*dm2[0]*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],8)+(180.0*dm2[0]*gkyl_ipow(L[0],6)-15.0*dm0[0]*gkyl_ipow(L[0],8))*gkyl_ipow(Nv[0],4))-6.0*dm1[0]*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],8)+12.0*dm1[0]*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],4)-6.0*dm1[0]*gkyl_ipow(L[0],8))*gkyl_ipow(mass,5)+vc[0]*(L[1]*(127.2792206135786*B[1]*dm0[1]+127.2792206135786*B[0]*dm0[0])*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],8)+L[1]*((-127.2792206135786*B[1]*dm0[1])-127.2792206135786*B[0]*dm0[0])*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],4))*gkyl_ipow(mass,4)+(vc[0]*gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],2)*((-5400.0*dm2[0]*gkyl_ipow(B[1],2))-5400.0*dm2[0]*gkyl_ipow(B[0],2)+10800.0*B[0]*B[1]*dm2[1])+gkyl_ipow(L[0],4)*(450.0*dm0[0]*gkyl_ipow(B[1],2)+450.0*dm0[0]*gkyl_ipow(B[0],2)-900.0*B[0]*B[1]*dm0[1]))*gkyl_ipow(L[1],2)+gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],4)*(360.0*dm1[0]*gkyl_ipow(B[1],2)+360.0*dm1[0]*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2)+gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],8)*((-360.0*dm1[0]*gkyl_ipow(B[1],2))-360.0*dm1[0]*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2))*gkyl_ipow(mass,3)+vc[0]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],8)*(dm0[1]*(3818.376618407359*gkyl_ipow(B[1],3)-3818.376618407359*B[1]*gkyl_ipow(B[0],2))-3818.376618407359*B[0]*dm0[0]*gkyl_ipow(B[1],2)+3818.376618407359*dm0[0]*gkyl_ipow(B[0],3))*gkyl_ipow(L[1],3)*gkyl_ipow(mass,2)+gkyl_ipow(Nv[0],8)*((-5400.0*dm1[0]*gkyl_ipow(B[1],4))+10800.0*dm1[0]*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)-5400.0*dm1[0]*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],4)*mass)*gkyl_ipow((5.441398092702653*L[1]*gkyl_ipow(L[0],10)*gkyl_ipow(Nv[0],9)-10.88279618540531*L[1]*gkyl_ipow(L[0],10)*gkyl_ipow(Nv[0],5)+5.441398092702653*Nv[0]*L[1]*gkyl_ipow(L[0],10))*gkyl_ipow(mass,4)+(gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],9)*(326.4838855621592*gkyl_ipow(B[1],2)+326.4838855621592*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3)+gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],5)*((-326.4838855621592*gkyl_ipow(B[1],2))-326.4838855621592*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3))*gkyl_ipow(mass,2)+gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],9)*(4897.258283432387*gkyl_ipow(B[1],4)-9794.516566864773*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+4897.258283432387*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],5),-1);

  f[3] = f[3] +  -1.0*((((12.24744871391589*B[1]*dm0[1]+12.24744871391589*B[0]*dm0[0])*gkyl_ipow(L[0],6)+((-146.9693845669908*B[1]*dm2[1])-146.9693845669908*B[0]*dm2[0])*gkyl_ipow(L[0],4))*gkyl_ipow(Nv[0],8)+(((-12.24744871391589*B[1]*dm0[1])-12.24744871391589*B[0]*dm0[0])*gkyl_ipow(L[0],6)+(146.9693845669908*B[1]*dm2[1]+146.9693845669908*B[0]*dm2[0])*gkyl_ipow(L[0],4))*gkyl_ipow(Nv[0],4))*gkyl_ipow(mass,4)+(L[1]*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],8)*(103.9230484541326*dm0[0]*gkyl_ipow(B[1],2)+103.9230484541326*dm0[0]*gkyl_ipow(B[0],2)+207.8460969082653*B[0]*B[1]*dm0[1])+L[1]*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],4)*((-103.9230484541326*dm0[0]*gkyl_ipow(B[1],2))-103.9230484541326*dm0[0]*gkyl_ipow(B[0],2)-207.8460969082653*B[0]*B[1]*dm0[1]))*gkyl_ipow(mass,3)+gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],2)*(dm0[1]*(367.4234614174767*gkyl_ipow(B[1],3)-367.4234614174767*B[1]*gkyl_ipow(B[0],2))-367.4234614174767*B[0]*dm0[0]*gkyl_ipow(B[1],2)+367.4234614174767*dm0[0]*gkyl_ipow(B[0],3))+dm2[1]*(4409.081537009723*B[1]*gkyl_ipow(B[0],2)-4409.081537009723*gkyl_ipow(B[1],3))+4409.081537009723*B[0]*dm2[0]*gkyl_ipow(B[1],2)-4409.081537009723*dm2[0]*gkyl_ipow(B[0],3))*gkyl_ipow(L[1],2)*gkyl_ipow(mass,2)+gkyl_ipow(Nv[0],8)*(3117.691453623978*dm0[0]*gkyl_ipow(B[1],4)-6235.382907247957*dm0[0]*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+3117.691453623978*dm0[0]*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],3)*mass)*gkyl_ipow(Nv[1]*(6.283185307179586*gkyl_ipow(L[0],9)*gkyl_ipow(Nv[0],8)-12.56637061435917*gkyl_ipow(L[0],9)*gkyl_ipow(Nv[0],4)+6.283185307179586*gkyl_ipow(L[0],9))*gkyl_ipow(mass,4)+Nv[1]*(gkyl_ipow(L[0],5)*gkyl_ipow(Nv[0],8)*(376.9911184307751*gkyl_ipow(B[1],2)+376.9911184307751*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2)+gkyl_ipow(L[0],5)*gkyl_ipow(Nv[0],4)*((-376.9911184307751*gkyl_ipow(B[1],2))-376.9911184307751*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2))*gkyl_ipow(mass,2)+L[0]*Nv[1]*gkyl_ipow(Nv[0],8)*(5654.866776461628*gkyl_ipow(B[1],4)-11309.73355292326*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+5654.866776461628*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],4),-1);

  f[4] = f[4] +  -1.0*((vc[0]*((15.0*dm0[1]*gkyl_ipow(L[0],8)-180.0*dm2[1]*gkyl_ipow(L[0],6))*gkyl_ipow(Nv[0],8)+(180.0*dm2[1]*gkyl_ipow(L[0],6)-15.0*dm0[1]*gkyl_ipow(L[0],8))*gkyl_ipow(Nv[0],4))-6.0*dm1[1]*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],8)+12.0*dm1[1]*gkyl_ipow(L[0],8)*gkyl_ipow(Nv[0],4)-6.0*dm1[1]*gkyl_ipow(L[0],8))*gkyl_ipow(mass,5)+vc[0]*(L[1]*(127.2792206135786*B[0]*dm0[1]+127.2792206135786*dm0[0]*B[1])*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],8)+L[1]*((-127.2792206135786*B[0]*dm0[1])-127.2792206135786*dm0[0]*B[1])*gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],4))*gkyl_ipow(mass,4)+(vc[0]*gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],4)*(dm0[1]*(450.0*gkyl_ipow(B[1],2)+450.0*gkyl_ipow(B[0],2))-900.0*B[0]*dm0[0]*B[1])+gkyl_ipow(L[0],2)*(dm2[1]*((-5400.0*gkyl_ipow(B[1],2))-5400.0*gkyl_ipow(B[0],2))+10800.0*B[0]*dm2[0]*B[1]))*gkyl_ipow(L[1],2)+dm1[1]*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],4)*(360.0*gkyl_ipow(B[1],2)+360.0*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2)+dm1[1]*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],8)*((-360.0*gkyl_ipow(B[1],2))-360.0*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2))*gkyl_ipow(mass,3)+vc[0]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],8)*(3818.376618407359*dm0[0]*gkyl_ipow(B[1],3)+dm0[1]*(3818.376618407359*gkyl_ipow(B[0],3)-3818.376618407359*B[0]*gkyl_ipow(B[1],2))-3818.376618407359*dm0[0]*B[1]*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3)*gkyl_ipow(mass,2)+dm1[1]*gkyl_ipow(Nv[0],8)*((-5400.0*gkyl_ipow(B[1],4))+10800.0*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)-5400.0*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],4)*mass)*gkyl_ipow((5.441398092702653*L[1]*gkyl_ipow(L[0],10)*gkyl_ipow(Nv[0],9)-10.88279618540531*L[1]*gkyl_ipow(L[0],10)*gkyl_ipow(Nv[0],5)+5.441398092702653*Nv[0]*L[1]*gkyl_ipow(L[0],10))*gkyl_ipow(mass,4)+(gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],9)*(326.4838855621592*gkyl_ipow(B[1],2)+326.4838855621592*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3)+gkyl_ipow(L[0],6)*gkyl_ipow(Nv[0],5)*((-326.4838855621592*gkyl_ipow(B[1],2))-326.4838855621592*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],3))*gkyl_ipow(mass,2)+gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],9)*(4897.258283432387*gkyl_ipow(B[1],4)-9794.516566864773*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+4897.258283432387*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],5),-1);

  f[5] = f[5] +  -1.0*((((12.24744871391589*B[0]*dm0[1]+12.24744871391589*dm0[0]*B[1])*gkyl_ipow(L[0],6)+((-146.9693845669908*B[0]*dm2[1])-146.9693845669908*dm2[0]*B[1])*gkyl_ipow(L[0],4))*gkyl_ipow(Nv[0],8)+(((-12.24744871391589*B[0]*dm0[1])-12.24744871391589*dm0[0]*B[1])*gkyl_ipow(L[0],6)+(146.9693845669908*B[0]*dm2[1]+146.9693845669908*dm2[0]*B[1])*gkyl_ipow(L[0],4))*gkyl_ipow(Nv[0],4))*gkyl_ipow(mass,4)+(L[1]*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],8)*(dm0[1]*(103.9230484541326*gkyl_ipow(B[1],2)+103.9230484541326*gkyl_ipow(B[0],2))+207.8460969082653*B[0]*dm0[0]*B[1])+L[1]*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],4)*(dm0[1]*((-103.9230484541326*gkyl_ipow(B[1],2))-103.9230484541326*gkyl_ipow(B[0],2))-207.8460969082653*B[0]*dm0[0]*B[1]))*gkyl_ipow(mass,3)+gkyl_ipow(Nv[0],8)*(gkyl_ipow(L[0],2)*(367.4234614174767*dm0[0]*gkyl_ipow(B[1],3)+dm0[1]*(367.4234614174767*gkyl_ipow(B[0],3)-367.4234614174767*B[0]*gkyl_ipow(B[1],2))-367.4234614174767*dm0[0]*B[1]*gkyl_ipow(B[0],2))-4409.081537009723*dm2[0]*gkyl_ipow(B[1],3)+dm2[1]*(4409.081537009723*B[0]*gkyl_ipow(B[1],2)-4409.081537009723*gkyl_ipow(B[0],3))+4409.081537009723*dm2[0]*B[1]*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2)*gkyl_ipow(mass,2)+dm0[1]*gkyl_ipow(Nv[0],8)*(3117.691453623978*gkyl_ipow(B[1],4)-6235.382907247957*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+3117.691453623978*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],3)*mass)*gkyl_ipow(Nv[1]*(6.283185307179586*gkyl_ipow(L[0],9)*gkyl_ipow(Nv[0],8)-12.56637061435917*gkyl_ipow(L[0],9)*gkyl_ipow(Nv[0],4)+6.283185307179586*gkyl_ipow(L[0],9))*gkyl_ipow(mass,4)+Nv[1]*(gkyl_ipow(L[0],5)*gkyl_ipow(Nv[0],8)*(376.9911184307751*gkyl_ipow(B[1],2)+376.9911184307751*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2)+gkyl_ipow(L[0],5)*gkyl_ipow(Nv[0],4)*((-376.9911184307751*gkyl_ipow(B[1],2))-376.9911184307751*gkyl_ipow(B[0],2))*gkyl_ipow(L[1],2))*gkyl_ipow(mass,2)+L[0]*Nv[1]*gkyl_ipow(Nv[0],8)*(5654.866776461628*gkyl_ipow(B[1],4)-11309.73355292326*gkyl_ipow(B[0],2)*gkyl_ipow(B[1],2)+5654.866776461628*gkyl_ipow(B[0],4))*gkyl_ipow(L[1],4),-1);

  f[6] = f[6] +  0.0;

  f[7] = f[7] +  0.0;

}

