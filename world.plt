#2016〜
#print ARGC
if (ARGC<=0) {
set size ratio -1
if (exist("WORLDMAP")==0) {
#http://www.gnuplotting.org/plotting-the-world-revisited/
#https://github.com/Gnuplotting/blog-entries/tree/master/040-plotting_the_world_revisited
  WORLDMAP='world_110m.txt'
}

#主要参考文献（他はその場で記述します）
#pp1453 https://pubs.er.usgs.gov/publication/pp1453
#pp1395 https://pubs.er.usgs.gov/publication/pp1395
#政春(2011) https://id.ndl.go.jp/bib/000011273603
#Special thanks to
#http://www.quadibloc.com/maps/mapint.htm


#略記号 Abbrev.
#N:north polar 方位図法で極中心のイメージ（正軸）
#E:equatorial 方位図法で赤道中心のイメージ（横軸）、円筒図法で正軸（省略可）
#T:transverse, second transverse 円筒図法で赤道中心で図法軸が赤道方向を向いているイメージ（横軸）かつXY平面でも図法軸がX方向となっている
#P:polar, first transverse 円筒図法で北極中心のイメージ（横軸）
#  direct   Az  C   Cy
#       \
#NP-view    N       P  / PT
#              (D)
#Eq-view    E      (E) / T(90W-90E axis)
#
#SP-view    S?
#Az:azimuthal, Cy:cylindrical
#Ed:equidistant, Ea:equal-area, C:conformal
#C:conic
#HS:hemisphere in a square, WS: world in a square, WT: world in a triangle, WE: world in an ellipse
#G:generalized (ex. Lagrange(): by Lambert, LagrangeG(): generalized version by Lagrange)

#変数 Variables
#C: 緯度経度, xy: 変換後平面座標, C1/C2/z/...: 一時変数
#lam: lambda, longitude, phi: latitude, rho: radius in xy-plane

#基本ツール Base
#(double)pi<pi
I={0,1}
DPR=180/pi
RPD=pi/180
goldenratio=(1+5**0.5)/2
silverratio=1+2**0.5
#Distance(C1,C2)=acos(sin(imag(C1))*sin(imag(C2))+cos(imag(C1))*cos(imag(C2))*cos(real(C2)-real(C1)))
#pp1395 (5-3)...(5-3a)
Distance1(C1,C2)=2*asin( (sin((imag(C1)-imag(C2))/2)**2 +cos(imag(C1))*cos(imag(C2))*sin((real(C2)-real(C1))/2)**2)**0.5 )
Distance2(C1,C2)=pi-Distance1(C1,real(C2)+pi+I*-imag(C2))
Distance(C1,C2)=(sin(imag(C1))*sin(imag(C2))+cos(imag(C1))*cos(imag(C2))*cos(real(C2)-real(C1))>=0?Distance1(C1,C2):Distance2(C1,C2))
#政春(2011) 最後に経度の範囲を調整
Eulerzyz(C,lam0,phi0,psi)=(lam=real(C),phi=imag(C),alp=asin(sin(phi0)*sin(phi)+cos(phi0)*cos(phi)*cos(lam-lam0)),bet=acos((sin(phi0)*sin(alp)-sin(phi))/(cos(phi0)*cos(alp)))*sgn(sin(lam-lam0)),bet-psi-floor((bet-psi)/2/pi+0.5)*2*pi+I*alp)
EachScale(C,sr,si)=real(C)*sr+I*imag(C)*si
conj(C)=real(C)-I*imag(C)
HarMean(C1,C2)=abs(C1)<1e-308?0:abs(C2)<1e-308?0:2*(C1**-1+C2**-1)**-1
GeoMean(C1,C2)=(abs(arg(C1**0.5)-arg(C2**0.5))<pi/2?1:-1)*C1**0.5*C2**0.5

#円筒図法 Cylindrical projections
CyEd(C)=C
Equirectangular(C,phi1)=EachScale(C,1.0,1.0/cos(phi1))
#Mercator(C)=abs(imag(C)*DPR)>=86+4?NaN:real(C)+I*(log(tan(imag(C)/2+pi/4)))
#tan(pi/4) don't meet 1
#Mercator(C)=abs(imag(C)*DPR)>=86+4?NaN:real(C)+I*(asinh(tan(imag(C))))
#asinh(x) is bad at x<0
Mercator(C)=abs(imag(C)*DPR)>=86+4?NaN:real(C)+I*(sgn(imag(C))*asinh(tan(abs(imag(C)))))
#Mercator(C)=2*atan(EStereo(C)/2)
CentralCy(C)=abs(imag(C)*DPR)>=86+4?NaN:real(C)+I*(tan(imag(C)))
LambertCyEa(C)=real(C)+I*(sin(imag(C)))
GallPeters(C)=EachScale(LambertCyEa(C),2**-0.5,2**0.5)
MillerCy(C)=EachScale(Mercator(EachScale(C,1,0.8)),1,1.25)
BraunStereo(C)=real(C)+I*(2*tan(imag(C)/2))
#BraunStereo(C)=LambertCyEa(InvMercator(Mercator(C)*0.5))*2
InvCyEd(xy)=xy
#InvMercator(xy)=real(xy)+I*(atan(exp(imag(xy)))*2-pi/2)
InvMercator(xy)=real(xy)+I*(atan(sinh(imag(xy))))
InvCentralCy(xy)=real(xy)+I*(atan(imag(xy)))
InvLambertCyEa(xy)=real(xy)+I*(asin(imag(xy)))
TCyEd(C)=(lam=real(C),phi=imag(C),B=cos(phi)*sin(lam),asin(B)+I*(atan2(tan(phi),cos(lam))))
TMercator(C)=(lam=real(C),phi=imag(C),B=cos(phi)*sin(lam),B*B==1?NaN:0.5*log((1+B)/(1-B))+I*(atan2(tan(phi),cos(lam))))
TLambertCyEa(C)=(lam=real(C),phi=imag(C),B=cos(phi)*sin(lam),B+I*(atan2(tan(phi),cos(lam))))
InvTCyEd(xy)=(x=real(xy),y=imag(xy),phi=asin(sin(y)*cos(x)),lam=atan2(tan(x),cos(y)),lam+I*phi)
#ArdenCloseCy -> into #arithmetical

#方位図法 Azimuthal projections
NOrtho(C)=(lam=real(C),phi=imag(C),rho=cos(phi),phi<0?NaN:(rho*sin(lam))+I*(-rho*cos(lam)))
NStereo(C)=(lam=real(C),phi=imag(C),rho=2*tan(pi/4-phi/2),sin(phi)<-0.54-0.46?NaN:(rho*sin(lam))+I*(-rho*cos(lam)))
NGnomo(C)=(lam=real(C),phi=imag(C),rho=1/tan(phi),phi<=0?NaN:(rho*sin(lam))+I*(-rho*cos(lam)))
NAzEd(C)=(lam=real(C),phi=imag(C),rho=pi/2-phi,(rho*sin(lam))+I*(-rho*cos(lam)))
NAzEa(C)=(lam=real(C),phi=imag(C),rho=2*sin(pi/4-phi/2),(rho*sin(lam))+I*(-rho*cos(lam)))
EOrtho(C)=(lam=real(C),phi=imag(C),cosz=cos(phi)*cos(lam),cosz<0?NaN:(cos(phi)*sin(lam))+I*(sin(phi)))  # K=1
EStereo(C)=(lam=real(C),phi=imag(C),cosz=cos(phi)*cos(lam),K=(cosz<=-1?0:2/(1+cosz)),cosz<=-0.88-0.12?NaN:(K*cos(phi)*sin(lam))+I*(K*sin(phi)))
EGnomo(C)=(lam=real(C),phi=imag(C),cosz=cos(phi)*cos(lam),cosz<=0?NaN:(tan(lam))+I*(tan(phi)/cos(lam)))  # K=1/cosz
EAzEd(C)=(lam=real(C),phi=imag(C),z=acos(cos(phi)*cos(lam)),K=(z==0?1:z/sin(z)),z>=pi?NaN:(K*cos(phi)*sin(lam))+I*(K*sin(phi)))
EAzEa(C)=(lam=real(C),phi=imag(C),cosz=cos(phi)*cos(lam),K=(2/(1+cosz))**0.5,cosz<=-1?NaN:(K*cos(phi)*sin(lam))+I*(K*sin(phi)))
SStereo(C)=conj(NStereo(conj(C)))
EStereoTwoHemisphere(C)=(signl=(real(C)<0?-1:1),EStereo(C-signl*pi/2)+signl*2)
EAzEdTwoHemisphere(C)=(signl=(real(C)<0?-1:1),EAzEd(C-signl*pi/2)+signl*pi/2)
InvNStereo(xy)=(x=real(xy),y=imag(xy),rho=(x*x+y*y)**0.5,phi=pi/2-2*atan(rho/2),lam=(rho==0?0:atan2(x,-y)),lam+I*phi)
InvEStereo(xy)=(x=real(xy),y=imag(xy),rho=(x*x+y*y)**0.5,c=2*atan(rho/2),phi=asin(y*(rho==0?1:sin(c)/rho)),lam=atan2(x*(rho==0?1:sin(c)/rho),cos(c)),lam+I*phi)
InvEOrtho(xy)=(x=real(xy),y=imag(xy),phi=asin(y),lam=asin(x/cos(phi)),lam+I*phi)
InvEAzEa(xy)=(x=real(xy),y=imag(xy),rho=(x*x+y*y)**0.5,c=2*asin(rho/2),phi=asin(y*(rho==0?1:sin(c)/rho)),lam=atan2(x*(rho==0?1:sin(c)/rho),cos(c)),lam+I*phi)

#擬円筒図法 Pseudocylindrical projections
Sinusoidal(C)=real(C)*cos(imag(C))+I*(imag(C))
MollTh(th,phi)=2*th+sin(2*th)-pi*sin(phi)
Newton_MollTh(th,phi)=abs(MollTh(th,phi))<1e-12?th:Newton_MollTh(th-MollTh(th,phi)/((MollTh(th+1e-8,phi)-MollTh(th-1e-8,phi))/2/1e-8),phi)
Mollweide(C)=(th=Newton_MollTh(0,imag(C)),2**1.5/pi*real(C)*cos(th)+I*(2**0.5*sin(th)))
#phi=asin((sin(2*th)+2*th)/pi)
#pi**2/8*(cos(phi)/cos(th))**2==1
# -> th=0.5705354790951345
# -> phi=0.7109888814838442
ph40d44m=40.73666218975136*RPD
ph50d28m=acos(2/pi)
BromleyMollweide(C)=EachScale(Mollweide(C),pi/2**1.5,2**1.5/pi)
CrasterParabolic(C)=(3/pi)**0.5*real(C)*(2*cos(imag(2*C/3))-1)+I*(3*pi)**0.5*sin(imag(C)/3)
QuarticAuthalic(C)=real(C)*cos(imag(C))/cos(imag(C)/2)+I*(2*sin(imag(C)/2))
InvSinusoidal(xy)=(x=real(xy),y=imag(xy),x/cos(y)+I*y)
InvMollweide(xy)=(x=real(xy),y=imag(xy),th=asin(y/2**0.5),pi/2**1.5*x/cos(th)+I*asin((2*th+sin(2*th))/pi))
TSinusoidal(C)=(lam=real(C),phi=imag(C),B=cos(phi)*sin(lam),asin(B)+I*(cos(asin(B))*atan2(tan(phi),cos(lam))))
#WagnerI ~ WagnerVI -> into #Aitoff...
#Hufnagel projection family is a further develeopment of the parametric angle
#https://www.mapthematics.com/interactive/hufnagel.html

#円錐図法 Conic projections
#標準緯線1本の式のみ提示 only tangent projections
EdC(C,phi1)=(lam=real(C),phi=imag(C),n=sin(phi1),rho=cos(phi1)/n+phi1-phi,rho*sin(n*lam)+I*(-rho*cos(n*lam)))
LambertCC(C,phi1)=(lam=real(C),phi=imag(C),n=sin(phi1),rho=cos(phi1)*tan(pi/4+phi1/2)**n/n/tan(pi/4+phi/2)**n,rho*sin(n*lam)+I*(-rho*cos(n*lam)))
#LambertCC(C,phi1)=exp((Mercator(C)-Mercator(I*phi1))*I*sin(phi1))*-I/tan(phi1)
AlbersEaC(C,phi1)=(lam=real(C),phi=imag(C),n=sin(phi1),rho=(cos(phi1)**2+2*n*sin(phi1)-2*n*sin(phi))**0.5/n,rho*sin(n*lam)+I*(-rho*cos(n*lam)))
LambertEaC(C,phi1)=(lam=real(C),phi=imag(C),n=(1+sin(phi1))/2,rho=2*sin(pi/4-phi/2)/n**0.5,rho*sin(n*lam)+I*(-rho*cos(n*lam)))

#擬円錐図法 Pseudoconic projections
Bonne(C,phi1)=(lam=real(C),phi=imag(C),rho=1/tan(phi1)+phi1-phi,n=cos(phi)/rho,rho*sin(n*lam)+I*(-rho*cos(n*lam)))
Werner(C)=(lam=real(C),phi=imag(C),rho=pi/2-phi,n=cos(phi)/rho,rho*sin(n*lam)+I*(-rho*cos(n*lam)))

#多円錐図法 Polyconic projections
AmericanPolyconic(C)=(lam=real(C),phi=imag(C),E=lam*sin(phi),rho=1/tan(phi),phi==0?lam:rho*sin(E)+I*(phi+rho*(1-cos(E))))
RectangularPolyconic(C)=(lam=real(C),phi=imag(C),E=2*atan(lam/2*sin(phi)),rho=1/tan(phi),phi==0?lam:rho*sin(E)+I*(phi+rho*(1-cos(E))))

#Ruler-and-compass constructions, or such conventional projections
NWiechel(C)=(lam=real(C),phi=imag(C),rho=1,(-rho*cos(lam)+sin(lam+phi))+I*(-rho*sin(lam)-cos(lam+phi)))
#Exercise: Prove that this projection is equal-area projection, in several ways.
#memo https://commons.wikimedia.org/wiki/File:Tsubaki_crescent_shaped_plastic_chain.png
WiechelPseudoconic(C,phi1)=(lam=real(C),phi=imag(C),n=sin(phi1),rho=(n==0?0:1/n),n==0?lam+cos(phi)-1+I*sin(phi):(rho*sin(n*lam-phi1)+cos(n*lam+phi-phi1))+I*(-rho*cos(n*lam-phi1)+sin(n*lam+phi-phi1)+tan(pi/2-phi1)))
#https://www.mapthematics.com/ProjectionsList.php?Projection=214
NIsoperimetricPseudoazimuthal(C)=(imag(C)<0?NaN:NWiechel(EachScale(C,1,2)-I*pi/2)/2)
Bacon(C)=(lam=real(C),phi=imag(C),y=pi/2*sin(phi),F=(lam==0?0:((pi/2)**2/abs(lam)+abs(lam))/2),lam==0?I*y:sgn(lam)*(abs(lam)-F+(F**2-y**2)**0.5)+I*y)
Ortelius(C)=(lam=real(C),phi=imag(C),y=phi,F=(lam==0?0:((pi/2)**2/abs(lam)+abs(lam))/2),lam==0?I*y:abs(lam)<pi/2?sgn(lam)*(abs(lam)-F+(F**2-y**2)**0.5)+I*y:sgn(lam)*(abs(lam)-pi/2+((pi/2)**2-y**2)**0.5)+I*y)
Nicolosi(C)=(lam=real(C),phi=imag(C),phi==0?lam:abs(lam)==pi/2?lam*cos(phi)+I*pi/2*sin(phi):(b=pi/2/lam-2*lam/pi,d=(1-(2*phi/pi)**2)/(sin(phi)-(2*phi/pi)),M=(b*sin(phi)/d-b/2)/(1+b**2/d**2),N=(d**2*sin(phi)/b**2+d/2)/(1+d**2/b**2),pi/2*(M+sgn(lam)*(M**2+cos(phi)**2/(1+b**2/d**2))**0.5)+I*pi/2*(N-sgn(phi*b*lam)*(N**2-(d**2*sin(phi)**2/b**2+d*sin(phi)-1)/(1+d**2/b**2))**0.5)))

#Realizing the unusual properties
#Eisenlohr: scale factor 0.5 @ center (differ from pp1453), 2.9141 @ peripheral
Eisenlohr(C)=(lam=real(C),phi=imag(C),C1=cos(lam/2),T=sin(phi/2)/(cos(phi/2)+(2*cos(phi))**0.5*C1),C2=(2/(1+T**2))**0.5,V=((cos(phi/2)+(cos(phi)/2)**0.5*(C1+sin(lam/2)))/(cos(phi/2)+(cos(phi)/2)**0.5*(C1-sin(lam/2))))**0.5,(1.5+2**0.5)*((-2*log(V)+C2*(V-V**-1))+I*(-2*atan(T)+C2*T*(V+V**-1))))
Littrow(C)=(lam=real(C),phi=imag(C),cos(phi)*cos(lam)<=0?NaN:sin(lam)/cos(phi)+I*(tan(phi)*cos(lam)))
#pp1453 says that Littrow is one of configuration of the generalized Lagrange
#but no explanation why it becomes retroazimuthal
#Littrow(C)=(z=EStereo(C)/2,abs(z)>=1?NaN:z/(1+z**2)*2)
#Littrow(C)=(z=EStereo(C)/2,abs(z)>=1?NaN:sin(2*atan(z)))
#Littrow(C)=(z=EStereo(C)/2,abs(z)>=1?NaN:tanh(2*atanh(z)))
LittrowNHemisphere(C)=(lam=real(C),phi=imag(C),sin(phi)<=0?NaN:sin(lam)/cos(phi)+I*(tan(phi)*cos(lam)))
LittrowSHemisphere(C)=(lam=real(C),phi=imag(C),sin(phi)>=0?NaN:sin(lam)/cos(phi)+I*(tan(phi)*cos(lam)))
CraigRetroazimuthal(C,phi1)=(lam=real(C),phi=imag(C),lam==0?I*(sin(phi)-cos(phi)*tan(phi1)):(lam-0)+I*(lam-0)*(sin(phi)*cos(lam-0)-cos(phi)*tan(phi1))/sin(lam-0))
Loximuthal(C,phi1)=(lam=real(C),phi=imag(C),phi==phi1?lam*cos(phi1):lam*(phi-phi1)/log(tan(pi/4+phi/2)/tan(pi/4+phi1/2))+I*(phi-phi1))
#TwoPointAz
#TwoPointEd

#Lagrange projection and its modifications
#Lagrange: scale factor 0.5 @ center (same as pp1453)
#http://www.quadibloc.com/maps/mcf0701.htm
Lagrange(C)=(imag(C)==pi/2?2*I:imag(C)==-pi/2?-2*I:EStereo(InvMercator(Mercator(C)*0.5)))
LagrangeG(C,W,phi1)=(lam=real(C),phi=imag(C),A1=((1+sin(phi1))/(1-sin(phi1)))**(0.5/W),V=((1+sin(phi))/(1-sin(phi)))**(0.5/W)/A1,LC=(V+1/V)/2+cos(lam/W),real(2*sin(lam/W)/LC)+I*((V-1/V)/LC))
#LagrangeG(C,W,phi1)=(imag(C)==pi/2?2*I:imag(C)==-pi/2?-2*I:EStereo(InvMercator((Mercator(C)-I*imag(Mercator(I*phi1)))/W)))
#Lagrange(C)=(C==pi?2:C==-pi?-2:(z=EStereo(C)/2/I,z/(1+(1-z)**0.5*(1+z)**0.5)*2*I))
#Lagrange(C)=(lam=real(C),phi=imag(C),c1=(1-tan(phi/2)**2)**0.5,c=1+c1*cos(lam/2),2*sin(lam/2)*c1/c+I*2*tan(phi/2)/c)
Lagrange240(C)=(imag(C)==pi/2?2*I:imag(C)==-pi/2?-2*I:EStereo(InvMercator(Mercator(C)/1.5)))
Lagrange270(C)=(imag(C)==pi/2?2*I:imag(C)==-pi/2?-2*I:EStereo(InvMercator(Mercator(C)*0.75)))
#https://web.archive.org/web/20180929084109/http://progonos.com:80/furuti/MapProj/Normal/ProjConf/projConf.html
EOrthoGilbertDoubleWorld(C)=EOrtho(InvMercator(Mercator(C)*0.5))
InvLagrange(xy)=InvMercator(Mercator(InvEStereo(xy))*2)
TLagrange(C)=(C==pi/2?2:C==-pi/2?-2:EStereo(InvMercator(TMercator(C)/2/I))*I)
TLagrange2UHalfPlane(C)=imag(C)<0?TLagrange(C)**-1:TLagrange(C)
TLagrange2LHalfPlane(C)=imag(C)>0?TLagrange(C)**-1:TLagrange(C)

#Aitoff projection and its modifications
#Meridian duplication ~ Das Umbeziffern
#Strebe(2018) https://research.tableau.com/paper/bevy-area-preserving-transforms-map-projection-designers
#https://map-projections.net/wagner-umbeziffern.php
#or Hatano's study for making the equal-area projections derived from the Mollweide
#https://www.jstage.jst.go.jp/article/grj1925/39/4/39_4_229/_article
#
Aitoff(C)=EachScale(EAzEd(EachScale(C,0.5,1)),2,1)
Hammer(C)=EachScale(EAzEa(EachScale(C,0.5,1)),2,1)
EckertGreifendorff(C)=EachScale(EAzEa(EachScale(C,0.25,1)),4,1)
CircularHammer(C)=EachScale(EAzEa(EachScale(C,0.5,1)),2**0.5,2**0.5)
#https://en.wikipedia.org/wiki/Bottomley_projection
Bottomley(C,phi1)=EachScale(Werner(EachScale(C,sin(phi1),1)),1/sin(phi1),1)
#WagnerIV: case phi=pi/2 -> make th=2*pi/3 -> truncated Mollweide 2:3**0.5/2 -> make 2:1 and retrieve the equivalence
WagnerIV(C)=EachScale(Mollweide(InvLambertCyEa(EachScale(LambertCyEa(C),1,2.0/3+3**0.5/2/pi))),(0.86310)/(2**1.5/pi),1.56548/(2**0.5))
#http://www.quadibloc.com/maps/mmi0902.htm
WagnerIX(C)=EachScale(EAzEa(EachScale(C,5/18.0,7/9.0)),18.0/5,9.0/7)
#神谷(2012) https://www.jstage.jst.go.jp/article/jjca/50/3/50_3_1/_article/-char/ja/
#エイトフ変換という言い方は以前から存在している
#https://www.jstage.jst.go.jp/article/jjca1963/3/1/3_1_33/_article/-char/ja
#英語だとAitoff's transformation methodとなるが、用例が少ない
#http://www.cartography.oregonstate.edu/pdf/2014_Savric_Jenny_ANewPseudocylindricalEqual-areaProjection.pdf
AitoffTrans_Stereo(C)=EachScale(EStereo(EachScale(C,0.5,1)),2,1)
#https://www.mapthematics.com/forums/viewtopic.php?f=8&t=559
OrthographicAitoff(C)=EachScale(EOrtho(EachScale(C,0.5,1)),2,1)
#https://arxiv.org/abs/astro-ph/0608500
GottEaElliptical(C)=EachScale(BromleyMollweide(InvEStereo(EStereo(EachScale(C,0.5,1))/I))*I,2,1)
InvHammer(xy)=EachScale(InvEAzEa(EachScale(xy,0.5,1)),2,1)

#Oblique projections
#http://www.quadibloc.com/maps/mps0404.htm
Atlantis(C)=Mollweide(Eulerzyz(C,-120*RPD,0,135*RPD))
#http://www.quadibloc.com/maps/meq0801.htm
EBriesemeister(C)=EachScale(Hammer(C),0.875**0.5,0.875**-0.5)
Briesemeister(C)=EBriesemeister(Eulerzyz(C,10*RPD,135*RPD,0))

#Arithmetical
#Using complex polynomial, or mixing a few projections
COblation(C,K,Q)=K*(C+Q/12.0*C**3)
MillerOblatedStereo(C)=COblation(NStereo(Eulerzyz(C,20*RPD,18*RPD,0)),0.9245,0.2522)
#pp1453 eqs. 146, 147, 148, 149.
MillerOblateG(C,da,db)=(Ta=(1-cos(da))/(1+cos(da)),Tb=(1-cos(db))/(1+cos(db)),Q=(Tb-Ta)/(Ta+Ta**2+Tb+Tb**2),COblation(C,1/((1+Ta)*(1+Q*Ta))**0.5,Q))
MillerOblateG_SC(C,da,db)=(Ta=(1-cos(da))/(1+cos(da)),Tb=(1-cos(db))/(1+cos(db)),Q=(Tb-Ta)/(Ta+Ta**2+Tb+Tb**2),COblation(C,1,Q))
#Ta=tan(da/2)**2, Tb=tan(db/2)**2
#EStereo(da)=2*tan(da/2), EStereo(db)=2*tan(db/2)
#d COblation(C)/d C = K*(1+Q/4.0*C**2)
#d EStereo(da)/d(da) = 1/cos(da/2)**2 = 1 + tan(da/2)**2
#d COblation(EStereo(C),K=1,Q)/d C = (1 + tan(C/2)**2) * (1+Q/4.0*EStereo(C)**2) = (1 + tan(C/2)**2) * (1+Q*tan(C/2)**2)
#(1 + tan(da/2)**2) * (1+Q*tan(da/2)**2) = (1 + tan(db/2)**2) * (1-Q*tan(db/2)**2)
#(1 + Ta) * (1+Q*Ta) = (1 + Tb) * (1-Q*Tb)
#(1 + Ta) + Q*(Ta + Ta**2) = (1 + Tb) - Q*(Tb + Tb**2)
#Q*(Ta + Ta**2 + Tb + Tb**2) = Tb - Ta
#Q = (Tb - Ta)/(Ta + Ta**2 + Tb + Tb**2)
#https://www.researchgate.net/publication/315860317_Combining_World_Map_Projections
ArdenCloseCy(C)=(abs(imag(C))>85*RPD?1/0:(Mercator(C)+LambertCyEa(C))/2)
#https://web.archive.org/web/20181004024945/http://www.progonos.com/furuti/MapProj/Normal/ProjOth/projOth.html
ArdenCloseNovelty(C)=(cos(real(C))<0?NaN:(LambertCyEa(C)+TLambertCyEa(C))/2)
#AugustEpicycloidal: scale factor: 0.5 @ center (differ from pp1453), max 4
#http://www.quadibloc.com/maps/mcf0702.htm
AugustEpicycloidal(C)=COblation(Lagrange(C),1,1)
TAugustEpicycloidal(C)=COblation(TLagrange(C),1,-1)
#FalseEisenlohr: scale factor: 0.5 @ center, max 3.1464
FalseEisenlohr(C)=sin(Lagrange(C)*pi/4/I)/pi*4*I
InvFalseEisenlohr(xy)=InvLagrange(asin(xy*pi/4/I)/pi*4*I)
ACN_August(C)=(AugustEpicycloidal(C)+TAugustEpicycloidal(C))/2
ACN_Lagrange(C)=(Lagrange(C)+TLagrange(C))/2
#http://www.quadibloc.com/maps/mmi0901.htm
WinkelTripel(C)=(Aitoff(C)+Equirectangular(C,acos(2/pi)))/2
#http://www.quadibloc.com/maps/mmi0903.htm
LaskowskiTriOptimal(C)=(lam=real(C),phi=imag(C),0.975534*lam-0.119161*lam*phi**2-0.0547009*lam*phi**4-0.0143059*lam**3*phi**2+I*(1.00384*phi+0.0802894*phi*lam**2+0.000199025*phi*lam**4+0.0998909*phi**3-0.02855*phi**3*lam**2-0.0491032*phi**5))
#https://www.mapthematics.com/ProjectionsList.php?Projection=320#lateral%20equidistant
LatTLat(C)=(cos(real(C))<0?NaN:real(TCyEd(C))+I*imag(CyEd(C)))
LateralEd(C)=EachScale(LatTLat(EachScale(C,0.5,1)),2,1)
#https://www.mapthematics.com/forums/viewtopic.php?f=8&t=633&sid=f8ccb8db8fed3fdb4b19a2480e4281bb#p1618
RoundedMercator(C)=abs(imag(C))>pi/4?(signp=(imag(C)<0?-1:1),EachScale(EAzEdTwoHemisphere(EachScale(C-signp*I*pi/4,1,2)),1,0.5**0.5)+I*imag(Mercator(signp*I*pi/4))):Mercator(C)
#SemicircularlyRoundedMercator(C)=abs(imag(C))>pi/4?(signp=(imag(C)<0?-1:1),EAzEdTwoHemisphere(InvMercator(Mercator(C)-signp*I*imag(Mercator(I*pi/4))))+I*imag(Mercator(signp*I*pi/4))):Mercator(C)
#https://www.mapthematics.com/forums/viewtopic.php?f=8&t=856
MiloDihedral(C)=(C1=LatTLat(Eulerzyz(C,90*RPD,45*RPD,-90*RPD))*(1-I)/2,asin(sin(real(C1))*sqrt(2))+I*asin(sin(imag(C1))*sqrt(2)))
#FalseMiloDihedral(C)=(C1=ApproxGuyou(C),InvApproxGuyou(real(C1))+I*InvApproxGuyou(imag(C1)))

#自作 My own work
#FalseCyEd(C)=real(C)+I*2/(1/imag(Mercator(C))+1/imag(LambertCyEa(C)))
FalseMercator(C)=real(C)+I/(2/imag(C)-1/imag(LambertCyEa(C)))
FalseArdenCloseCy(C)=real(C)+I*imag(LambertCyEa(C))/(2-imag(C)/imag(LambertCyEa(C)))
TwoLeavesWerner(C)=(real(C)>=0?Werner(C-pi/2)*I:Werner(C+pi/2)/I)
AugustTwice(C)=AugustEpicycloidal(InvEStereo(AugustEpicycloidal(C)*1.5))
NoName1(C)=EStereo(InvEStereo(real(C)*4/pi)+I*(imag(C)+0.5*(cos(real(C)/1.5)-0.5)*sin(imag(C))*cos(imag(C))))*pi/4
StereoEa(C)=(lam=real(C),phi=imag(C),(cos(2*Newton_MollTh(0,phi))+1+I*sin(2*Newton_MollTh(0,phi)))*sqrt(lam*2/pi+2))
ESquashedStereo(C)=(C1=EStereo(C),InvEStereo(real(C1))+I*InvEStereo(imag(C1)))
#
#Combined
EAzEaOval(C)=(lam=real(C),phi=imag(C),abs(lam)<=pi/2?EAzEa(C):((lam-pi/2*sgn(lam))*0.5+cos(phi)*sgn(lam))*2**0.5+I*sin(phi)*2**0.5)
EAzEaWiechel(C)=(lam=real(C),phi=imag(C),abs(lam)<=pi/2?EAzEa(C):lam>=0?(phi>=0?NWiechel(lam*0.5+pi/4+I*phi)+I:conj(NWiechel(lam*0.5+pi/4-I*phi))-I)*2**0.5:(phi>=0?-conj(NWiechel(-lam*0.5+pi/4+I*phi))+I:-NWiechel(-lam*0.5+pi/4-I*phi)-I)*2**0.5)
EAzEdCapsule(C)=(abs(real(C))<=pi/2?CyEd(C):EAzEd(C-pi/2*sgn(real(C)))+pi/2*sgn(real(C)))
#
#Directional path offset
#Invphi_Hammer(phi,lam,y)=abs(y-imag(Hammer(lam+I*phi)))<1e-12?phi:Invphi_Hammer(phi+y-imag(Hammer(lam+I*phi)),lam,y)
#Invphi_Hammer(phi,lam,y)=y==0?0:sgn(y)*2*atan((((cos(0.5*lam)**2*y**4-8*y**2+16)**0.5+y**2-4)/(cos(0.5*lam)-1)/y**2)**0.5)
#ComplementaryEckertGreifendorff(C)=(signl=(real(C)<0?-1:1),z=Hammer(EachScale(C,0.5,1)+signl*pi/2),EachScale(z-real(Hammer(signl*pi/2+I*Invphi_Hammer(0,signl*pi/2,imag(z)))),2,1))
lamy2x_EAzEa(lam,y)=2**-0.5*((cos(lam)**2-2)*y**2-cos(lam)*(cos(lam)**2*y**4-8*y**2+16)**0.5+4)**0.5*sgn(lam)*sgn((cos(lam)**2*y**4-8*y**2+16)**0.5-2*cos(lam))
ComplementaryEckertGreifendorff(C)=(signl=(real(C)<0?-1:1),z=EAzEa(EachScale(C,0.25,1)+signl*pi/4),EachScale(z-lamy2x_EAzEa(signl*pi/4,imag(z)),4,1))
ComplementaryEckertGreifendorffG(C,lams)=(W=(pi/2-lams)/pi,signl=(real(C)<0?-1:1),z=EAzEa(EachScale(C,W,1)+signl*lams),EachScale(z-lamy2x_EAzEa(signl*lams,imag(z)),W**-1,1))
#
#Strebe's homotopy
HamCyEa(C,n)=(n==0?Hammer(C):LambertCyEa(InvHammer(Hammer(C)*n))/n)
HamTCyEa(C,n)=(n==0?Hammer(C):EachScale(LambertCyEa(InvEAzEa(EAzEa(EachScale(C,0.5,1))*n/I))/n*I,2,1))
HammerPeters(C,n)=(ymax=(n==0?2**0.5:sin(asin(sin(pi/4)*n)*2)/n),EachScale(HamCyEa(C,n),(2**0.5/ymax)**-1,2**0.5/ymax))
#
#ApproxQuincuncial(C)=(C1=ACN_Lagrange(C),C2=GeoMean(Lagrange(C),TLagrange(C)),z=GeoMean(C1,C2),z*(1-r1-r2)+C1*r1+C2*r2+z**5*r3/100.0)
#f(x)=real(ApproxQuincuncial((x+1)*pi/2)/1.3110287771460599)
#fit [x=0:1] f(x) '0to1-0-1.dat' using 1:3 via r1,r2,r3
#r1=-0.580510, r2=0.119212, r3=0.155830
#g(x)=real(ApproxQuincuncial((x+1)*pi/2)-1.3110287771460599)/real(ApproxQuincuncial((x+1)*pi/2)-ApproxQuincuncial((x+1)*pi/2+I*pi/2/10000))
#within +-0.18 km
ApproxQuincuncial(C)=(C1=ACN_Lagrange(C),C2=GeoMean(Lagrange(C),TLagrange(C)),z=GeoMean(C1,C2),z*1.461298+C1*-0.580510+C2*0.119212+z**5*0.00155830)
N2E(C)=InvEStereo(NStereo(C))
NApproxQuincuncial(C)=ApproxQuincuncial(N2E(C))
ApproxAdamsHS(C)=(abs(real(C))>pi/2?NaN:ApproxQuincuncial(C))
ApproxGuyou(C)=(abs(real(C))>pi/2?NaN:ApproxQuincuncial(Eulerzyz(C,90*RPD,45*RPD,-90*RPD))*(1-I))
#
#ApproxEisenlohr(C)=(Q=-4.0/(2+2**3*r1/10+2**5*r2/100+2**7*r3/1000)**2,z=Lagrange(C)/I,z=COblation((z+z**3*r1/10+z**5*r2/100+z**7*r3/1000),1,Q),(z+z**3*r4/10)*I)
#f(x)=abs(ApproxEisenlohr(pi+I*(x)*pi/2-I*0.1e-6)-ApproxEisenlohr(pi+I*(x)*pi/2-I*1.1e-6))*1e6
#fit [x=0:1] f(1-(1-x)**2)*2/(3+8**0.5) '0to1-0-1.dat'using 1:3 via r1,r2,r3,r4
#r1=-0.167040, r2=0.0408032, r3=-0.00796779, r4=-0.00751721
#scale factor: 0.5 @ center, max 2.9147
#ApproxEisenlohr(C)=(Q=-4.0/1.87840514688**2,z=Lagrange(C)/I,z=COblation((z+z**3*-0.0167040+z**5*0.000408032+z**7*-0.00000796779),1,Q),(z+z**3*-0.000751721)*I)
#ApproxEisenlohr(C)=(Q=(2+2**3*r1/10+2**5*r3/100)/(1+2**2*r2/10+2**4*r4/100),z=Lagrange(C)/I,z=(z+z**3*r1/10+z**5*r3/100)/(1+z**2*r2/10+z**4*r4/100),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)*I)
#r1=-0.0505242, r2=0.00344757, r3=-0.0641668, r4=-0.0903435
#scale factor: 0.5 @ center, max 2.9143
ApproxEisenlohr(C)=(Q=(1.939047264)/(0.986924068),z=Lagrange(C)/I,z=(z+z**3*-0.00505242+z**5*-0.000641668)/(1+z**2*0.000344757+z**4*-0.000903435),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)*I)
Inflation(C,n)=C==pi?C:C==-pi?C:InvEStereo(EStereo(C)*n)
#ApproxEisenlohr90(C)=(n=EStereo(pi/2)/EStereo(pi/4),Q=(2+2**3*r1/10+2**5*r3/100)/(1+2**2*r2/10+2**4*r4/100),z=Lagrange(Inflation(C,n))/I,z=(z+z**3*r1/10+z**5*r3/100)/(1+z**2*r2/10+z**4*r4/100),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)/n*I)
#f(x)=abs(ApproxEisenlohr90(pi+I*(x)*pi*3/4)-ApproxEisenlohr90(pi+I*(x)*pi*3/4-I*1e-6))*1e6
#fit [x=0:1] f(1-(1-x)**2)/r5 '0to1-0-1.dat'using 1:3 via r1,r2,r3,r4,r5
#r1=0.961831, r2=1.27208, r3=0.126022, r4=0.31818, r5=1.11993
#scale factor: 0.5 @ center, min 0.46389 @ 0°N 66°E/W, max 1.1203
ApproxEisenlohr90(C)=(n=EStereo(pi/2)/EStereo(pi/4),Q=(2.80979184)/(1.5597408),z=Lagrange(Inflation(C,n))/I,z=(z+z**3*0.0961831+z**5*0.00126022)/(1+z**2*0.127208+z**4*0.0031818),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)/n*I)
#ApproxEisenlohr270(C)=(n=EStereo(pi/2)/EStereo(pi*3/4),Q=(2+2**3*r1/10+2**5*r3/100)/(1+2**2*r2/10+2**4*r4/100),z=Lagrange(Inflation(C,n))/I,z=(z+z**3*r1/10+z**5*r3/100)/(1+z**2*r2/10+z**4*r4/100),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)/n*I)
#f(x)=abs(ApproxEisenlohr270(pi+I*(x)*pi/4)-ApproxEisenlohr270(pi+I*(x)*pi/4-I*1e-6))*1e6
#something is bad numerically
#r1=-0.975065, r2=-1.04292, r3=0.170893, r4=0.220625, r5=12.6371
#scale factor: 0.5 @ center, max 12.638
ApproxEisenlohr270(C)=(n=EStereo(pi/2)/EStereo(pi*3/4),Q=(1.27463376)/(0.618132),z=Lagrange(Inflation(C,n))/I,z=(z+z**3*-0.0975065+z**5*0.00170893)/(1+z**2*-0.104292+z**4*0.00220625),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)/n*I)
#
#conformal version of RoundedMercator
NHemisphere2Sphere(C)=InvEStereo(LittrowNHemisphere(C)*2)
SHemisphere2Sphere(C)=InvEStereo(LittrowSHemisphere(C)*2)
Sphere2NHemisphere(C)=InvEStereo(TLagrange2UHalfPlane(C))
Sphere2SHemisphere(C)=InvEStereo(TLagrange2LHalfPlane(C))
#ApproxEisenlohr_CRM(C)=(Q=(2+2**3*r1/10+2**5*r3/100)/(1+2**2*r2/10+2**4*r4/100),z=Lagrange(C)/I,z=(z+z**3*r1/10+z**5*r3/100)/(1+z**2*r2/10+z**4*r4/100),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)*I)
#ApproxConformalRoundedMercator(C)=Mercator(imag(C)<0?Sphere2SHemisphere(InvEStereo(ApproxEisenlohr_CRM(SHemisphere2Sphere(C))/ApproxEisenlohr_CRM(pi/2)*2)):Sphere2NHemisphere(InvEStereo(ApproxEisenlohr_CRM(NHemisphere2Sphere(C))/ApproxEisenlohr_CRM(pi/2)*2)))
#CRM1: constant scale along interruption
#f(x)=abs(ConformalRoundedMercator(I*(1+x)*pi/4+I*1.1e-4)-ConformalRoundedMercator(I*(1+x)*pi/4+I*0.1e-4))*1e4
#fit [x=0:1] f(x**2)/r5 '0to1-0-1.dat'using 1:3 via r1,r2,r3,r4,r5
#r1=-0.164059, r2=0.0997807, r3=-0.0275333, r4=-0.0994837, r5=2.62206
#scale factor on CRM1: 1.086 @ center, 0.927 @ 90°W/E, 2.62198~2.62214 @ peripheral
ApproxEisenlohr_CRM1(C)=(Q=1.81635881760378,z=Lagrange(C)/I,z=(z+z**3*-0.0164059+z**5*-0.000275333)/(1+z**2*0.00997807+z**4*-0.000994837),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)*I)
ApproxConformalRoundedMercator1(C)=Mercator(imag(C)<0?Sphere2SHemisphere(InvEStereo(ApproxEisenlohr_CRM1(SHemisphere2Sphere(C))/ApproxEisenlohr_CRM1(pi/2)*2)):Sphere2NHemisphere(InvEStereo(ApproxEisenlohr_CRM1(NHemisphere2Sphere(C))/ApproxEisenlohr_CRM1(pi/2)*2)))
#CRM2: makes interruption semicircles
#f(x)=abs(ConformalRoundedMercator(I*(1+x)*pi/4)-(pi/2+I*imag(ConformalRoundedMercator(I*pi/4))))
#fit [x=0:1] f(x**2)/(pi/2) '0to1-0-1.dat'using 1:3 via r1,r2,r3,r4
#r1=0.0444918, r2=0.568859, r3=-0.0725106, r4=-0.151152
#radius: 1.57076~1.57082
#scale factor on CRM2: 1.046 @ center, 0.959 @ 90°W/E, max 4.15 @ pole
ApproxEisenlohr_CRM2(C)=(Q=1.67231024137696,z=Lagrange(C)/I,z=(z+z**3*0.00444918+z**5*-0.000725106)/(1+z**2*0.0568859+z**4*-0.00151152),z=sin(z*pi/2/Q)/pi*2*Q,(z)/(1)*I)
ApproxConformalRoundedMercator2(C)=Mercator(imag(C)<0?Sphere2SHemisphere(InvEStereo(ApproxEisenlohr_CRM2(SHemisphere2Sphere(C))/ApproxEisenlohr_CRM2(pi/2)*2)):Sphere2NHemisphere(InvEStereo(ApproxEisenlohr_CRM2(NHemisphere2Sphere(C))/ApproxEisenlohr_CRM2(pi/2)*2)))
#
#Generalized Circular Hammer
EaOblation_TRemap2R(C,n)=(C1=TCyEd(C),imag(C1)>(1-n)*pi?InvTCyEd(EachScale(C1-I*pi,1,0.5/n)+I*pi):imag(C1)>-(1-n)*pi?InvTCyEd(EachScale(C1,1,0.5/(1-n))):InvTCyEd(EachScale(C1+I*pi,1,0.5/n)-I*pi))
EaOblation_RemapR2(C,n)=(C1=CyEd(C),real(C1)>0.5*pi?InvCyEd(EachScale(C1-pi,n/0.5,1)+pi):real(C1)>-0.5*pi?InvCyEd(EachScale(C1,(1-n)/0.5,1)):InvCyEd(EachScale(C1+pi,n/0.5,1)-pi))
EaOblation(C,n)=EaOblation_RemapR2(EaOblation_TRemap2R(C,n),n)
GeneralizedCircularHammer(C,n)=(n==0?EachScale(EAzEa(InvTCyEd(EachScale(TCyEd(C),1,0.5))),2**0.5,2**0.5):CircularHammer(EaOblation(C,n)))
#Critical Circular Equal-area, under development
CriticalCircularEa_Test0(C)=GeneralizedCircularHammer(C,pi/(2+pi))
EaOblation_TRemapR2(C,n)=(C1=TCyEd(C),imag(C1)>0.5*pi?InvTCyEd(EachScale(C1-I*pi,1,n/0.5)+I*pi):imag(C1)>-0.5*pi?InvTCyEd(EachScale(C1,1,(1-n)/0.5)):InvTCyEd(EachScale(C1+I*pi,1,n/0.5)-I*pi))
EaOblation_Remap2R(C,n)=(C1=CyEd(C),real(C1)>(1-n)*pi?InvCyEd(EachScale(C1-pi,0.5/n,1)+pi):real(C1)>-(1-n)*pi?InvCyEd(EachScale(C1,0.5/(1-n),1)):InvCyEd(EachScale(C1+pi,0.5/n,1)-pi))
TEaOblation(C,n)=EaOblation_TRemapR2(EaOblation_Remap2R(C,n),n)
#plot [0:1]imag(Hoge(I*x*pi*2/(2+pi)))-x*2
#fit [x=0:1] imag(Hoge(I*(x*0.5-0.5+1)*pi*2/(2+pi)))-(x*0.5-0.5+1)*2 '0to1-0-1.dat' using 1:2 via a1,a2,a3,a4,a5,a6,a7,a8
#pi-(pi-1.534268/(2-(1/a1)))/a1
CriticalCircularEa_Test1(C)=(iota=pi/(2+pi),a1=1.1138,a2=1.08237,a3=1.07127,a4=1.06021,a5=1.05384,a6=1.04878,a7=1.0483,a8=1.07467,kappa=(2-a1)*(2-(1/a1))*(2-a2)*(2-(1/a2))*(2-a3)*(2-(1/a3))*(2-a4)*(2-(1/a4))*(2-a5)*(2-(1/a5))*(2-a6)*(2-(1/a6))*(2-a7)*(2-(1/a7))*(2-a8)*(2-(1/a8)),C1=EaOblation(C,iota),abs(real(C1))*0.5*kappa/(1-iota)>=pi?CircularHammer(C1):CircularHammer(EachScale( TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation(TEaOblation( EachScale(C1,0.5*kappa/(1-iota),1) ,0.5/a1),0.5/a2),0.5/a3),0.5/a4),0.5/a5),0.5/a6),0.5/a7),0.5/a8),0.5*a8),0.5*a7),0.5*a6),0.5*a5),0.5*a4),0.5*a3),0.5*a2),0.5*a1) ,(1-iota)/0.5/kappa,1) ))
PCyEd(C)=(lam=real(C),phi=imag(C),B=-cos(phi)*cos(lam),(atan2(sin(lam),tan(phi))+I*asin(B)))
InvPCyEd(xy)=(x=real(xy),y=imag(xy),phi=asin(cos(y)*cos(x)),lam=atan2(sin(x),-tan(y)),lam+I*phi)
EaOblation_QRemapR2(C,n)=(C1=PCyEd(C),signl=(real(C1)>pi/2?+1:real(C1)>-pi/2?0:-1),InvPCyEd(EachScale(PCyEd(EaOblation_RemapR2(InvPCyEd(EachScale(C1-pi/2*signl,2,1)),n)),0.5,1)+pi/2*signl))
DQEaOblation(C,n)=EaOblation_QRemapR2(EaOblation_TRemap2R(C,n),n)
EaOblation_QRemap2R(C,n)=(C1=PCyEd(C),signl=(real(C1)>pi/2?+1:real(C1)>-pi/2?0:-1),InvPCyEd(EachScale(PCyEd(EaOblation_RemapR2(InvPCyEd(EachScale(C1-pi/2*signl,2,1)),n)),0.5,1)+pi/2*signl))
DREaOblation(C,n)=(m=(n>0.5?(n+0.5)/2:(n+0.5)/2),EaOblation_RemapR2(EaOblation_QRemap2R(EaOblation_QRemapR2(EaOblation_TRemap2R(C,n),n),m),m))
#
#TwoCurlyHemisphere
CurlyCurve(t,amp)=InvMercator(t+I*cos(t*2)*amp)
#set parametric
#set trange [-pi:pi]
#set samples 3601
#course A: minimum curvature, nearest to Antarctica, pass through the east side of Tenerife.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.226),-12*RPD,13*RPD,97*RPD)),real(xy)):(imag(xy)) w l
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-97*RPD,167*RPD,12*RPD)),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-97*RPD,167*RPD,12*RPD)),real(xy)):(imag(xy)) w l, '+' using (xy=-NAzEa(Eulerzyz(Eulerzyz(CurlyCurve(t,0.226),-12*RPD,13*RPD,97*RPD),-97*RPD,167*RPD,12*RPD)),real(xy)):(imag(xy)) w l
#course B: pass through the west side of Tenerife and La Gomera.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.232),-11.5*RPD,13.5*RPD,97.5*RPD)),real(xy)):(imag(xy)) w l
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-97.5*RPD,167.5*RPD,11.5*RPD)),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-97.5*RPD,167.5*RPD,11.5*RPD)),real(xy)):(imag(xy)) w l, '+' using (xy=-NAzEa(Eulerzyz(Eulerzyz(CurlyCurve(t,0.232),-11.5*RPD,13.5*RPD,97.5*RPD),-97.5*RPD,167.5*RPD,11.5*RPD)),real(xy)):(imag(xy)) w l
#course C: pass through between Viti Levu and Yasawa Islands. (need to check about Mamanuca)
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.254),-10.5*RPD,14.5*RPD,97.5*RPD)),real(xy)):(imag(xy)) w l
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-97.5*RPD,165.5*RPD,10.5*RPD)),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-97.5*RPD,165.5*RPD,10.5*RPD)),real(xy)):(imag(xy)) w l, '+' using (xy=-NAzEa(Eulerzyz(Eulerzyz(CurlyCurve(t,0.254),-10.5*RPD,14.5*RPD,97.5*RPD),-97.5*RPD,165.5*RPD,10.5*RPD)),real(xy)):(imag(xy)) w l
#course D: pass through Koro Sea and the west side of Yacata Island.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.268),-8.5*RPD,15.0*RPD,96.0*RPD)),real(xy)):(imag(xy)) w l
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-96.0*RPD,165.0*RPD,8.5*RPD)),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-96.0*RPD,165.0*RPD,8.5*RPD)),real(xy)):(imag(xy)) w l, '+' using (xy=-NAzEa(Eulerzyz(Eulerzyz(CurlyCurve(t,0.268),-8.5*RPD,15.0*RPD,96.0*RPD),-96.0*RPD,165.0*RPD,8.5*RPD)),real(xy)):(imag(xy)) w l
#course E: pass through between Ogea Lev and Vatoa.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.300),-7.0*RPD,16.0*RPD,95.0*RPD)),real(xy)):(imag(xy)) w l
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-95.0*RPD,164.0*RPD,7.0*RPD)),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=-NAzEa(Eulerzyz(C,-95.0*RPD,164.0*RPD,7.0*RPD)),real(xy)):(imag(xy)) w l, '+' using (xy=-NAzEa(Eulerzyz(Eulerzyz(CurlyCurve(t,0.300),-7.0*RPD,16.0*RPD,95.0*RPD),-95.0*RPD,164.0*RPD,7.0*RPD)),real(xy)):(imag(xy)) w l
#course F: a bit tilted.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.263),-3.0*RPD,16.0*RPD,93.0*RPD)),real(xy)):(imag(xy)) w l
#course G: cut North Island, New Zealand.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.296),-0.0*RPD,17.5*RPD,91.0*RPD)),real(xy)):(imag(xy)) w l
#course H: cut Iceland.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.376),-0.0*RPD,21.5*RPD,93.0*RPD)),real(xy)):(imag(xy)) w l
#can't through between Greenland and Iceland, without tilt.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.536),-0.0*RPD,29.0*RPD,96.0*RPD)),real(xy)):(imag(xy)) w l
#course I: cut Chukchi Peninsula, Russia.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.250),-0.0*RPD,14.5*RPD,94.0*RPD)),real(xy)):(imag(xy)) w l
#course J: cut Greenland.
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=CyEd(C),real(xy)):(imag(xy)) w l, '+' using (xy=CyEd(Eulerzyz(CurlyCurve(t,0.165),-0.0*RPD,15.0*RPD,97.0*RPD)),real(xy)):(imag(xy)) w l
ECurlyHemisphere_Ham(C,amp)=(gam=(pi/2-imag(CurlyCurve(0,amp)))/(pi/2+imag(CurlyCurve(0,amp))),EachScale(EAzEa(EachScale(C,gam,1)),1/gam,1))
OnEastCurlyHemisphere(C,amp)=imag(Mercator(C))>amp*cos(2*real(C))
TwoCurlyHemisphereA(C)=(gamma=(pi/2-imag(CurlyCurve(0,0.226)))/(pi/2+imag(CurlyCurve(0,0.226))),C1=Eulerzyz(C,-97*RPD,167*RPD,-168*RPD),OnEastCurlyHemisphere(C1,0.226)?ECurlyHemisphere_Ham(N2E(C1),0.226)+EAzEa(pi/2-imag(CurlyCurve(0,0.226)))/gamma:ECurlyHemisphere_Ham(N2E(-C1+pi/2),0.226)*-I-EAzEa(pi/2-imag(CurlyCurve(0,0.226))))
TwoCurlyHemisphereA_NU(C)=(gamma=(pi/2-imag(CurlyCurve(0,0.226)))/(pi/2+imag(CurlyCurve(0,0.226))),C1=Eulerzyz(C,-97*RPD,167*RPD,-168*RPD),OnEastCurlyHemisphere(C1,0.226)?ECurlyHemisphere_Ham(N2E(C1),0.226)*exp(I*12*RPD*gamma)+EAzEa(pi/2-imag(CurlyCurve(0,0.226)))/gamma:ECurlyHemisphere_Ham(N2E(-C1+pi/2),0.226)*-I*exp(I*-12*RPD)-EAzEa(pi/2-imag(CurlyCurve(0,0.226))))
#set trange [-2*pi:2*pi]
#set samples 3601
#plot '+' using (C=Eulerzyz(CurlyCurve(t+pi/2,0.226)+(t<0?-1e-6*I:1e-6*I),-12*RPD,13*RPD,97*RPD),xy=TwoCurlyHemisphereA(C),real(xy)):(imag(xy)) w l, 'world_10m.txt' using (C=($1+I*$2)*RPD,xy=TwoCurlyHemisphereA(C),real(xy)):(imag(xy)) w l, '+' using (C=t,xy=TwoCurlyHemisphereA(C),real(xy)):(imag(xy)) w l
TwoCurlyHemisphereE(C)=(gamma=(pi/2-imag(CurlyCurve(0,0.300)))/(pi/2+imag(CurlyCurve(0,0.300))),C1=Eulerzyz(C,-95*RPD,164*RPD,-173*RPD),OnEastCurlyHemisphere(C1,0.300)?ECurlyHemisphere_Ham(N2E(C1),0.300)+EAzEa(pi/2-imag(CurlyCurve(0,0.300)))/gamma:ECurlyHemisphere_Ham(N2E(-C1+pi/2),0.300)*-I-EAzEa(pi/2-imag(CurlyCurve(0,0.300))))
TwoCurlyHemisphereF(C)=(gamma=(pi/2-imag(CurlyCurve(0,0.263)))/(pi/2+imag(CurlyCurve(0,0.263))),C1=Eulerzyz(C,-93*RPD,164*RPD,-177*RPD),OnEastCurlyHemisphere(C1,0.263)?ECurlyHemisphere_Ham(N2E(C1),0.263)+EAzEa(pi/2-imag(CurlyCurve(0,0.263)))/gamma:ECurlyHemisphere_Ham(N2E(-C1+pi/2),0.263)*-I-EAzEa(pi/2-imag(CurlyCurve(0,0.263))))
TwoCurlyHemisphereG(C)=(gamma=(pi/2-imag(CurlyCurve(0,0.296)))/(pi/2+imag(CurlyCurve(0,0.296))),C1=Eulerzyz(C,-91*RPD,162.5*RPD,-180*RPD),OnEastCurlyHemisphere(C1,0.296)?ECurlyHemisphere_Ham(N2E(C1),0.296)+EAzEa(pi/2-imag(CurlyCurve(0,0.296)))/gamma:ECurlyHemisphere_Ham(N2E(-C1+pi/2),0.296)*-I-EAzEa(pi/2-imag(CurlyCurve(0,0.296))))
TwoCurlyHemisphereH(C)=(gamma=(pi/2-imag(CurlyCurve(0,0.376)))/(pi/2+imag(CurlyCurve(0,0.376))),C1=Eulerzyz(C,-93*RPD,158.5*RPD,-180*RPD),OnEastCurlyHemisphere(C1,0.376)?ECurlyHemisphere_Ham(N2E(C1),0.376)+EAzEa(pi/2-imag(CurlyCurve(0,0.376)))/gamma:ECurlyHemisphere_Ham(N2E(-C1+pi/2),0.376)*-I-EAzEa(pi/2-imag(CurlyCurve(0,0.376))))
TwoCurlyHemisphereI(C)=(gamma=(pi/2-imag(CurlyCurve(0,0.250)))/(pi/2+imag(CurlyCurve(0,0.250))),C1=Eulerzyz(C,-94*RPD,165.5*RPD,-180*RPD),OnEastCurlyHemisphere(C1,0.250)?ECurlyHemisphere_Ham(N2E(C1),0.250)+EAzEa(pi/2-imag(CurlyCurve(0,0.250)))/gamma:ECurlyHemisphere_Ham(N2E(-C1+pi/2),0.250)*-I-EAzEa(pi/2-imag(CurlyCurve(0,0.250))))
#search for another projection, analytically controllable
ECurlyHemisphere_Ait(C,amp)=(gam=(pi/2-imag(CurlyCurve(0,amp)))/(pi/2+imag(CurlyCurve(0,amp))),EachScale(EAzEd(EachScale(C,gam,1)),1/gam,1))
ECurlyHemisphere_PsS(C,amp)=(gam=(pi/2-imag(CurlyCurve(0,amp)))/(pi/2+imag(CurlyCurve(0,amp))),EachScale(EStereo(EachScale(C,gam,1)),1/gam,1))
ECurlyHemisphere_MOS(C,amp)=MillerOblateG_SC(EStereo(C),pi/2+imag(CurlyCurve(0,amp)),pi/2-imag(CurlyCurve(0,amp)))
#
#pointed pole, caring about resolution-efficiency
#https://www.mapthematics.com/forums/viewtopic.php?f=8&t=841
parallelarc(a,s)=a*(sin(s/a)+I-I*cos(s/a))
#AmericanPolyconic(C)=(lam=real(C),phi=imag(C),phi==0?lam:parallelarc(1/tan(phi),lam*cos(phi))+I*phi)
StereoPolyconic(C)=(lam=real(C),phi=imag(C),phi==0?lam:parallelarc(2/tan(phi),lam*cos(phi))+I*2*tan(phi/2))
catenary_s(a,s)=(a*asinh(s/a))+I*(sqrt(a*a+s*s)*sgn(a))
catenarym1_s(a,s)=catenary_s(a,s)-I*a
Catesoidal(C)=(lam=real(C),phi=imag(C),phi==0?lam:catenarym1_s(1/tan(phi),lam*cos(phi))+I*phi)
Catereo(C)=(lam=real(C),phi=imag(C),phi==0?lam:catenarym1_s(2/tan(phi),lam*cos(phi))+I*2*tan(phi/2))
RoundedPolyconic(C)=(lam=real(C),phi=imag(C),phi==0?lam:parallelarc(2/tan(phi),lam*(1-(phi/pi*2)**2))+I*2*tan(phi/2))
RoundyPolyconic(C)=(lam=real(C),phi=imag(C),phi==0?lam:parallelarc(2/tan(phi),lam*(2*cos(phi*2/3)-1))+I*2*tan(phi/2))
#EAzEdCapsule: area=17.621, min-minor=1.000, resolution-efficiency=0.844
#AmericanPolyconic: area=25.114, min-minor=0.860, resolution-efficiency=0.609
#StereoPolyconic: area=22.220, min-minor=0.956, resolution-efficiency=0.719
#Catereo: area=21.476, min-minor=0.815, resolution-efficiency=0.624
#RoundedPolyconic: area=23.757, min-minor=1.000, resolution-efficiency=0.727
#RoundyPolyconic: area=23.098, min-minor=0.995, resolution-efficiency=0.734
#
#polyhedral
RadialModLatTLat(C)=(C1=InvNStereo(EStereo(C)),LatTLat(N2E(real(C1)+I*(pi/2-asin(sin((pi/2-imag(C1))/2)*sqrt(2))))))
Trilinear(C)=imag(C)<0||real(C)<0||real(C)>pi/2?NaN:(d1=Distance(C,I*pi/2),d2=Distance(C,0),d3=Distance(C,pi/2),d=(d1+d2+d3-pi)/3,(d2+pi/2-d3)/2+I*(pi/2-d1+d)*sin(pi/3))
Hemisphere2Octant(C)=InvMercator(Mercator(InvEStereo(EStereo(InvMercator(Mercator(InvEStereo(EStereo(InvNStereo(NStereo(C)*EStereo(pi/2)/EStereo(pi/1.5)))*I))/2-pi/4))/I))/2+pi/4)
TrilinearHT(C)=Trilinear(Hemisphere2Octant(C))-pi/4-I*pi/2*sin(pi/3)/3
TrilinearHex(C)=abs(real(C))>pi/2?NaN:(d1=Distance(C,I*pi/2)-pi/2,d2=Distance(C,-pi/2-I*pi/6)-pi/2,d3=Distance(C,pi/2-I*pi/6)-pi/2,d=(d1+d2+d3)/3,(d2-d3)/2+I*(-d1+d)*sin(pi/3))
#consider resolution-efficiency: map (45°E,0°)-(45°E,35.26438968°N) to (pi/4,0)-(pi/4,pi/(4√3)) linearly
TrilinearC0_ModDist(x)=(b=atan(2**-0.5),d0=pi/4,db=pi/2-b,x<b?x-Distance(0,pi/4+I*x)+d0+x*(db-d0)/b:x<pi/2-b?x:x+Distance(0,pi/4+I*(pi/2-x))-d0+(x-pi/2)*(db-d0)/b)
#plot [0:pi/2]TrilinearC0_ModDist(x)
TrilinearC0(C)=imag(C)<0||real(C)<0||real(C)>pi/2?NaN:(d1=TrilinearC0_ModDist(Distance(C,I*pi/2)),d2=TrilinearC0_ModDist(Distance(C,0)),d3=TrilinearC0_ModDist(Distance(C,pi/2)),d=(d1+d2+d3-pi)/3,(d2+pi/2-d3)/2+I*(pi/2-d1+d)*sin(pi/3))
TrilinearC1_ModDist(x)=(b=atan(2**-0.5),d0=pi/4,db=pi/2-b,g1=1.0412134,g2=(pi/4-g1*b)/(pi/4-b),x<b?x*g1+(-Distance(0,pi/4+I*x)+d0+x*(db-d0)/b)*g2:x<pi/2-b?(x-pi/4)*g2+pi/4:pi/2-(pi/2-x)*g1+(+Distance(0,pi/4+I*(pi/2-x))-d0+(x-pi/2)*(db-d0)/b)*g2)
#plot [db-0.00001:db+0.00001](TrilinearC1_ModDist(x+1e-8)-TrilinearC1_ModDist(x-1e-8))/2e-8
TrilinearC1(C)=imag(C)<0||real(C)<0||real(C)>pi/2?NaN:(d1=TrilinearC1_ModDist(Distance(C,I*pi/2)),d2=TrilinearC1_ModDist(Distance(C,0)),d3=TrilinearC1_ModDist(Distance(C,pi/2)),d=(d1+d2+d3-pi)/3,(d2+pi/2-d3)/2+I*(pi/2-d1+d)*sin(pi/3))
TrilinearC2_Mod3(x)=(g3=0.7225,x+(x-b)*(x-pi/4)*(x-(pi/2-b))*g3)
TrilinearC2_ModDist(x)=(b=atan(2**-0.5),d0=pi/4,db=pi/2-b,g1=1.0506194992,g2=(pi/4-g1*b)/(pi/4-b),x<b?x*g1+(-TrilinearC2_Mod3(Distance(0,pi/4+I*x))+d0+x*(db-d0)/b)*g2:x<pi/2-b?(TrilinearC2_Mod3(x)-pi/4)*g2+pi/4:pi/2-(pi/2-x)*g1+(+TrilinearC2_Mod3(Distance(0,pi/4+I*(pi/2-x)))-d0+(x-pi/2)*(db-d0)/b)*g2)
#plot [db-0.001:db+0.001](TrilinearC2_ModDist(x+1e-5)-2*TrilinearC2_ModDist(x)+TrilinearC2_ModDist(x-1e-5))/4e-5
TrilinearC2(C)=imag(C)<0||real(C)<0||real(C)>pi/2?NaN:(d1=TrilinearC2_ModDist(Distance(C,I*pi/2)),d2=TrilinearC2_ModDist(Distance(C,0)),d3=TrilinearC2_ModDist(Distance(C,pi/2)),d=(d1+d2+d3-pi)/3,(d2+pi/2-d3)/2+I*(pi/2-d1+d)*sin(pi/3))
#call 'world.plt' defstrainfunc TrilinearC2
#splot 'upright-1.dat' using (C=($1+I*$2)*RPD,xy=TrilinearC2(C),real(xy)):(imag(xy)):((aplusb(C)-aminusb(C))/2<0.73671?0:1) with pm3d
TrilinearC3_Mod5(x)=(g3=0.862867,g5=-3.5,x+(x-b)*(x-pi/4)*(x-(pi/2-b))*g3+(x-b)**2*(x-pi/4)*(x-(pi/2-b))**2*g5)
TrilinearC3_ModDist(x)=(b=atan(2**-0.5),d0=pi/4,db=pi/2-b,g1=1.05236018082,g2=(pi/4-g1*b)/(pi/4-b),x<b?x*g1+(-TrilinearC3_Mod5(Distance(0,pi/4+I*x))+d0+x*(db-d0)/b)*g2:x<pi/2-b?(TrilinearC3_Mod5(x)-pi/4)*g2+pi/4:pi/2-(pi/2-x)*g1+(+TrilinearC3_Mod5(Distance(0,pi/4+I*(pi/2-x)))-d0+(x-pi/2)*(db-d0)/b)*g2)
TrilinearC3(C)=imag(C)<0||real(C)<0||real(C)>pi/2?NaN:(d1=TrilinearC3_ModDist(Distance(C,I*pi/2)),d2=TrilinearC3_ModDist(Distance(C,0)),d3=TrilinearC3_ModDist(Distance(C,pi/2)),d=(d1+d2+d3-pi)/3,(d2+pi/2-d3)/2+I*(pi/2-d1+d)*sin(pi/3))
#
#


#Baseball projection, alpha version
PMercator(C)=(lam=real(C),phi=imag(C),B=-cos(phi)*cos(lam),(atan2(sin(lam),tan(phi)))+I*0.5*log((1+B)/(1-B)))
PLagrange(C)=(cos(imag(C))*cos(real(C))==-1?2:cos(imag(C))*cos(real(C))==1?-2:EStereo(InvMercator(PMercator(C)*0.5)))
PStercator(C,n)=(n==0?NStereo(C):PMercator(InvNStereo(NStereo(C)*n))/n)
NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4)+real(z**2*1.4+z**6*0.02)))
NBaseballStrip(C)=(z=PStercator(C,1.00896),z*(1-0.0377262)+z/(1+(1-z*0.287150)**0.5*(1+z*0.287150)**0.5)*2*0.0377262+z**5*0.000453101)
BaseballCurve(t)=NBaseballStrip(InvNStereo(NS_BaseballCurve(t)))
#f(x)=abs(BaseballCurve(x*pi/2+1e-6)-BaseballCurve(x*pi/2-1e-6))/Distance(InvNStereo(NS_BaseballCurve(x*pi/2+1e-6)),InvNStereo(NS_BaseballCurve(x*pi/2-1e-6)))
#g(x)=(f(x*0.5)-f(1-x*0.5))*10
#set fit lambda_factor 2
#fit [x=0:1] g(x) '0to1-0-1.dat' using 1:2 via r1,r2
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=NStereo(C),real(xy)):(imag(xy)) w l, 'full-15.dat' using (C=($1+I*$2)*RPD,xy=NBaseballStrip(C),real(xy)):(imag(xy)) w l, real(BaseballCurve(t)),imag(BaseballCurve(t))
#plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=NOrtho(C),real(xy)):(imag(xy)) w l, real(NOrtho(InvNStereo(NS_BaseballCurve(t)))),imag(NOrtho(InvNStereo(NS_BaseballCurve(t)))),cos(t),sin(t)
#plot t,f(t), t,g(t)
#NBaseballStrip(C)=(z=PStercator(C,r1),z*(1-r3)+z/(1+(1-z*r2)**0.5*(1+z*r2)**0.5)*2*r3+z**5*r4)
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4*0.5)+real(z**2*1.5*0.5)))
#r1=0.889157, r2=0.296751, r3=0.264297, r4=0.000134605
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4)+real(z**2*1.5)))
#r1=1.00442, r2=0.317375, r3=0.0212041, r4=0.000311200
#凸に縫い合わせられるかの検討
#BBCheckStitch1(x,l1)=real((BaseballCurve(x)-BaseballCurve(0))/I)-sin(real(BaseballCurve(pi/2-x)-BaseballCurve(pi/2))*l1)/l1
#BBCheckStitch2(x,l1,b0)=BBCheckStitch2(x,l1,b0,b1)=imag((BaseballCurve(x)-BaseballCurve(0))/I)*exp(I*b0)-imag(BaseballCurve(pi/2-x)-BaseballCurve(pi/2))-I*(1-cos(real(BaseballCurve(pi/2-x)-BaseballCurve(pi/2))*l1))/l1
#plot t,real(BBCheckStitch2(t,1.145,0.935)), t,imag(BBCheckStitch2(t,1.145,0.935)), t,BBCheckStitch1(x,1.145)
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4)+real(z**2*1.35)))
#r1=1.01105, r2=0.319759, r3=0.0623970, r4=0.000364165
#plot t,real(BBCheckStitch2(t,1.010,1.218)), t,imag(BBCheckStitch2(t,1.010,1.218)), t,BBCheckStitch1(x,1.010)
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4)+real(z**2*1.6)))
#r1=1.00210, r2=0.325297, r3=0.00501518, r4=0.000282542
#plot t,real(BBCheckStitch2(t,1.227,0.785)), t,imag(BBCheckStitch2(t,1.227,0.785)), t,BBCheckStitch1(x,1.227)
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4)+real(z**2*1.61)))
#r1=1.00192, r2=0.322182, r3=0.00519325, r4=0.000273823
#plot t,real(BBCheckStitch2(t,1.235,0.772)), t,imag(BBCheckStitch2(t,1.235,0.772)), t,BBCheckStitch1(x,1.235)
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4*0.9)+real(z**2*1.5)))
#r1=0.984959, r2=0.322001, r3=0.00467789, r4=0.000224160
#plot t,real(BBCheckStitch2(t,1.238,0.802)), t,imag(BBCheckStitch2(t,1.238,0.802)), t,BBCheckStitch1(x,1.238)
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4)+real(z**2*1.5+z**6*0.01)))
#r1=1.00418, r2=0.302350, r3=0.0186706, r4=0.000337529
#plot t,real(BBCheckStitch2(t,1.222,0.862)), t,imag(BBCheckStitch2(t,1.222,0.862)), t,BBCheckStitch1(x,1.222)
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4/4)+real(z**2*1.4+z**6*0.02)))
#r1=1.01953, r2=0.259496, r3=0.154880, r4=0.000326087
#plot t,real(BBCheckStitch2(t,1.220,0.944)), t,imag(BBCheckStitch2(t,1.220,0.944)), t,BBCheckStitch1(x,1.220)
#写像から外形を求めてみる
#NS_BaseballCurve(t)=(z=exp(I*t),2*z*exp(I*imag(z**4*r2)+real(z**2*r1+z**6*r3+z**10*r4+z**14*r5)))


#以下未整理

AitoffLagrange(C)=EachScale(Lagrange270(EachScale(C,1/1.5,1)),1.5,1)
AitoffTrans_TSinusoidal(C)=EachScale(TSinusoidal(EachScale(C,0.5,1)),2,1)
AitoffTrans_TCyEd(C)=EachScale(TCyEd(EachScale(C,0.5,1)),2,1)
Aitoff1p6(C)=EachScale(EAzEd(EachScale(C,1/1.6,1)),1.6,1)
AitoffExtrap(C)=Sinusoidal(C)*-0.5+Aitoff(C)*1.5
AitoffEmphasized(C)=EachScale(Aitoff(C),1,2)-I*(imag(C)*1)
AitoffTrans2p5_TSinusoidal(C)=EachScale(TSinusoidal(EachScale(C,1/2.5,1)),2.5,1)
AitoffEmp(C,B1,B2)=EachScale(EAzEd(EachScale(C,1.0/B1,1)),B1,B2)-EAzEd(I*imag(C))*(B2-1)
AitoffEmp2(C,B1,B2,n)=(lam=real(C)/B1/n**0.5,phi=imag(C),z=acos(cos(phi)*cos(lam)),K=(z==0?1:z/sin(z)),z>=pi?NaN:(K*cos(phi)*sin(lam))*B1*n**0.5+I*(K*sin(phi))*B2*n-EAzEd(I*imag(C))*(B2*n-1))
ATSEmp2(C,B1,B2,n)=(lam=real(C)/B1/n**0.5,phi=imag(C),cosz=cos(phi)*cos(lam),K=2/(1+cosz),cosz<=-0.88-0.12?NaN:(K*cos(phi)*sin(lam))*B1*n**0.5+I*(K*sin(phi))*B2*n-EAzEd(I*imag(C))*(B2*n-1))
HammerEmp2(C,B1,B2,n)=(lam=real(C)/B1/n**0.5,phi=imag(C),cosz=cos(phi)*cos(lam),K=(2/(1+cosz))**0.5,cosz<=-0.88-0.12?NaN:(K*cos(phi)*sin(lam))*B1*n**0.5+I*(K*sin(phi))*B2*n-EAzEa(I*imag(C))*(B2*n-1))
AitoffTEmp(C,B1,B2,n)=(lam=real(C)/B1,phi=imag(C)/I,z=acos(cos(phi)*cos(lam)),K=(z==0?1:z/sin(z)),z>=pi?NaN:(K*cos(phi)*sin(lam))*B1*B2+(real(C))*(1-B2)+I*(K*sin(phi))*I)

#マップ端で緯線と経線が直交するパラメータを探す
#右上象限でarg()等を使いやすいように経線南向きと緯線東向きとする
#aem(lat,B1,B2)=(AitoffEmp(pi+I*(lat-1e-6),B1,B2)-AitoffEmp(pi+I*(lat+1e-6),B1,B2))/2e-6
#aep(lat,B1,B2)=(AitoffEmp(pi+1e-6+I*lat,B1,B2)-AitoffEmp(pi-1e-6+I*lat,B1,B2))/2e-6
#aea(lat,B1,B2)=arg(aep(lat,B1,B2))-arg(aem(lat,B1,B2))-pi/2
#aepm(B1,B2)=sin(pi/B1)*pi/2*B1+I*(-cos(pi/B1)*B2+B2-1)
#aepp(B1,B2)=cos(pi/B1)*pi/2*B1+I*(sin(pi/B1)*B2)
#aepa(B1,B2)=arg(aepp(B1,B2))-arg(aepm(B1,B2))-pi/2
#set contour
#unset surface
#splot [1.7:1.725][1.68:1.7] aea(89.9*RPD,x,y),aea(70*RPD,x,y),aea(50*RPD,x,y),aea(30*RPD,x,y),aea(10*RPD,x,y)
#→W1=1.71~1.715,W2=1.68~1.69（一点には決まらなかった）

plot 'world.dat' w l
pause -1
plot WORLDMAP w l
pause -1
plot 'full-15.dat' using (C=($1+I*$2)*RPD,xy=Mercator(C),real(xy)):(imag(xy)) w l, WORLDMAP using (C=($1+I*$2)*RPD,xy=Mercator(C),real(xy)):(imag(xy)) w l
pause -1

#open-mouse sinusoidal I
#plot [-3:3][-2:2]'full-15.dat' using (C=($1+I*$2)*RPD,xy=Sinusoidal(C)+sgn($1)*(((pi/2)**2-imag(C)**2)**0.5-pi/2*cos(imag(C))),real(xy)):(imag(xy)) w l, 'full-15.dat' using (C=($1+I*$2)*RPD,xy=TSinusoidal(C),xy=xy+I*sgn($2)*(((pi/2)**2-real(xy)**2)**0.5-pi/2*cos(real(xy))),real(xy)):(imag(xy)) w l,((pi/2)**2-x**2)**0.5
#open-mouse sinusoidal II
#plot [-3:3][-2:2]'full-15.dat' using (C=($1+I*$2)*RPD,xy=Sinusoidal(C),xy=xy+I*sgn($2)*(((pi/2)**2-real(xy)**2)**0.5-acos(abs(real(xy))/pi*2)),(cos(real(C))<0?xy=NaN:0),real(xy)):(imag(xy)) w l, 'full-15.dat' using (C=($1+I*$2)*RPD,xy=TSinusoidal(C),xy=xy+sgn($1)*(((pi/2)**2-imag(xy)**2)**0.5-acos(abs(imag(xy))/pi*2)),(cos(real(C))<0?xy=NaN:0),real(xy)):(imag(xy)) w l,((pi/2)**2-x**2)**0.5
} else {
if (ARG1 eq 'defstrainfunc') {
str='vmer(C)=('.ARG2.'(C+{0,-1e-6})-'.ARG2.'(C+{0,1e-6}))/2e-6'
eval str
str='vpar(C)=('.ARG2.'(C+{1e-6,0})-'.ARG2.'(C+{-1e-6,0}))/2e-6'
eval str
merscale(C)=abs(vmer(C))
parscale(C)=abs(vpar(C))/cos(imag(C))
areascale(C)=abs(vpar(C))/cos(imag(C))*abs(vmer(C))*sin(arg(vpar(C)/vmer(C)))
a2b2(C)=merscale(C)**2+parscale(C)**2
angledeform(C)=2*asin(((a2b2-2*areascale)/(a2b2+2*areascale))**0.5)
#cf. https://en.wikipedia.org/wiki/Tissot's_indicatrix
#merscale=h,parscale=k,areascale=s=a*b=h*k*sin(thprime),a2b2=a**2+b**2=h**2+k**2
aplusb(C)=(a2b2(C)+2*areascale(C))**0.5
aminusb(C)=(a2b2(C)-2*areascale(C))**0.5
#a/b+b/a=a2b2/areascale
#a/b+b/a-2=(a2b2-2*areascale)/areascale=((a-b)/areascale**0.5)**2 せん断ひずみエネルギーとして採用
#areascale(C)+1/areascale(C)-2 第1項はピストン外の大気圧、第2項はピストン内の等圧変化 面積ひずみエネルギーとして採用
#神谷(2012)「ただし、各図法は、Esnを最小とするように、拡大または縮小を行っている。」
#全体の面積を正積図法と合わせるようにすると面積ひずみエネルギーが最小に近くなる
#NStereo(C)/2**0.5*1.04,NAzEd(C)/pi*2**1.5*1.004
#→理想のひずみエネルギー定義はこれで最小になったりしてほしい気分だが深入りしない（対数ひずみかな？）
#神谷(2012)「一般には、ハンメル図法はモルワイデ図法より優れていると評価されている」
#今回のせん断ひずみエネルギー定義だとハンメル図法は3.162、モルワイデ図法は2.974
#しかしハンメル図法はx方向に引き伸ばされすぎのきらいがあり
#Briesemeisterの比率(1.75/2)を採用すると2.872、0.77の比率で2.777
#モルワイデ図法の方を調整すると0.77の比率で2.607（1°メッシュで特異点付近の近似が十分か？）
#モルワイデ図法のせん断ひずみが発散するところが嫌われているのはわかるが、積分すると収束する程度の特異点である
#エイトフ図法は比率0.79、エイトフ変換平射は比率0.82

#set parametric
#set contour surface
#set isosamples 180,90
#set urange [0.5*RPD:179.5*RPD]
#set vrange [0.5*RPD:89.5*RPD]
#set dgrid3d
#splot real(HogeProjection(u+I*v)),imag(HogeProjection(u+I*v)),areascale(u+I*v)
#pause -1
unset parametric
unset contour
unset dgrid3d
set terminal wxt 11
ene=0
splot 'upright-1.dat' using (C=($1+I*$2)*RPD,xy=@ARG2(C),real(xy)):(imag(xy)):(ene=ene+(areascale(C)+1/areascale(C)-2)*cos(imag(C)),areascale(C)) with pm3d
print 'areaenergy(1/4) =', ene/180/90*pi/2
#splot 'upright-1.dat' using (C=($1+I*$2)*RPD,xy=@ARG2(C),real(xy)):(imag(xy)):(ene=ene+(areascale(C)+1/areascale(C)-2),areascale(C)) with pm3d
#print 'areaenergy(pole-heavy)(1/4) =', ene/180/90
set terminal wxt 12
ene=0
splot 'upright-1.dat' using (C=($1+I*$2)*RPD,xy=@ARG2(C),real(xy)):(imag(xy)):(ene=ene+(a2b2(C)/areascale(C)-2)*cos(imag(C)),a2b2(C)/areascale(C)-2) with pm3d
print 'shearenergy(1/4)=', ene/180/90*pi/2
#splot 'upright-1.dat' using (C=($1+I*$2)*RPD,xy=@ARG2(C),real(xy)):(imag(xy)):(ene=ene+(a2b2(C)/areascale(C)-2),a2b2(C)/areascale(C)-2) with pm3d
#print 'shearenergy(pole-heavy)(1/4)=', ene/180/90
set terminal wxt 13
#ene=0
#splot 'upright-1.dat' using (C=($1+I*$2)*RPD,xy=@ARG2(C),real(xy)):(imag(xy)):(ene=ene+(areascale(C))*cos(imag(C)),areascale(C)) with pm3d
#print 'area(1/4) =', ene/180/90*pi*pi/2
min(a,b)=(b<a?b:a)  #bがNaNなら更新しない
ene=1e300
splot 'upright-1.dat' using (C=($1+I*$2)*RPD,xy=@ARG2(C),real(xy)):(imag(xy)):(ene=min(ene,(aplusb(C)-aminusb(C))/2),(aplusb(C)-aminusb(C))/2) with pm3d
print 'min minor-axis =', ene
}
if (ARG1 eq 'survey') {
print 'AitoffEmphasized'
do for [i=10:20] {
  do for [j=10:20] {
    print 'B1=',i/10.0,', B2=',j/10.0
    ae(C)=EachScale(AitoffEmp(C,i/10.0,j/10.0),0.8,1)
    call 'world.plt' defstrainfunc ae
  }
}
#せん断ひずみエネルギー最小は(B1を決めてB2を振った最小)B1=1.6,B2=2.4→B1=1.7,B2=2.7→B1=1.9,B2=3.4→B1=2.0,B2=3.8
#比率0.8にしてせん断ひずみエネルギー最小はB1=1.6,B2=1.8→B1=1.8,B2=2.3→B1=2.0,B2=2.8

#B1=2.0,B2=2.0から外形をほぼ維持して極限を取る B1=2.0*n,B2=2.0*n**2,n→∞
#ap(C)=(lam=real(C),phi=imag(C),lam*cos(phi)*(phi==0?1:phi/sin(phi))+I*(acos(cos(phi)*cos(lam/32))/sin(acos(cos(phi)*cos(lam/32)))*sin(phi))*512-EAzEd(I*phi)*511)
#ap(C)=(lam=real(C),phi=imag(C),lam*cos(phi)*(phi==0?1:phi/sin(phi))+I*(acos(cos(phi)*(1-(lam/32)**2/2))/sin(acos(cos(phi)*(1-(lam/32)**2/2)))*sin(phi))*512-I*phi*511)
#ap(C)=(lam=real(C),phi=imag(C),lam*cos(phi)*(phi==0?1:phi/sin(phi))+I*((acos(cos(phi))/sin(acos(cos(phi)))+(0+(lam/32)**2/6)*cos(phi)/(sin(phi)/phi)/cos(phi/2)**2)*sin(phi))*512-I*phi*511)
#ap(C)=(lam=real(C),phi=imag(C),lam*cos(phi)*(phi==0?1:phi/sin(phi))+I*((acos(cos(phi))/sin(acos(cos(phi)))+(0+(lam/32)**2/2)*cos(phi)/(sin(phi))*(sin(phi)-phi*cos(phi))/sin(phi)**2)*sin(phi))*512-I*phi*511)
#ap(C)=(lam=real(C),phi=imag(C),lam*cos(phi)*(phi==0?1:phi/sin(phi))+I*(phi==0?0:phi+lam**2/tan(phi)*(1-phi/tan(phi))/4))
#2次の微小量を取り出すとやはり式が煩雑に
#いっそ放物型を通り越して双曲型にいったらどうだろう
#経度を虚数にして→AitoffEmp2
#AitoffEmp2(C,1,0,-1)→正積ではない、緯線が直線、せん断ひずみはMollweideより小さい
#AitoffEmp2(C,1,0.5,-1)→半球が円に近い
#AitoffEmp2(C,2,2,-1)→
#ap(C)→X座標は経度に比例
#AitoffEmp2(C,2,2,1)→X座標はAitoff
#北極のエッジの接線
#f(x,B1,B2)=pi/2+((1-cos(pi/B1))*B2-1)/(pi/2*sin(pi/B1)*B1)*x
#plot [-1:1][-2:2]1/(1-cos(pi/(x**-0.5))),0.5/x
#halfpoleB2(B1)=1/(1-cos(pi/B1))
#plot [0:pi]abs(vpar(x+I*20*RPD)),abs(vpar(x+I*40*RPD)),abs(vpar(x+I*60*RPD)),abs(vpar(x+I*80*RPD))
#エイトフ変換平射に応用→ATSEmp2
#北極のエッジの接線
#f(x,B1,B2)=2+((1-cos(pi/B1))*B2-1)/(sin(pi/B1)*B1)*x
#halfpoleB2(B1)=同上
#zerojacobianB2(B1)=1/(1/-cos(pi/B1)+1)
}
if (ARG1 eq 'makegraticule') {
pause -1

set print 'full-15.dat'
do for [k=-1:0] {
do for [l=-2:1] {

do for [i=0:6] {
  lat=i*15+k*90
  if (i==0) {lat=lat+1e-8}
  if (i==6) {lat=lat-1e-8}
  do for [j=0:90] {
    lon=j+l*90
    if (j==0) {
      print lon+1e-8,lat
      lon=lon+1e-3
    }
    if (j==90) {
      print lon-1e-3,lat
      lon=lon-1e-8
    }
    print lon,lat
  }
  print ''
}
do for [i=0:6] {
  lon=i*15+l*90
  if (i==0) {lon=lon+1e-8}
  if (i==6) {lon=lon-1e-8}
  do for [j=0:90] {
    lat=j+k*90
    if (j==0) {
      print lon,lat+1e-8
      lat=lat+1e-3
    }
    if (j==90) {
      print lon,lat-1e-3
      lat=lat-1e-8
    }
    print lon,lat
  }
  print ''
}

}  # do for l
}  # do for k
unset print

set print 'upright-1.dat'
do for [i=0:179] {
  do for [j=0:89] {
    print i+0.5,j+0.5
  }
  print ''
}
unset print

set print '0to1-0-1.dat'
do for [i=0:20] {
  print i/20.0,0.0,1.0
}
unset print

}
}
#    EOF
