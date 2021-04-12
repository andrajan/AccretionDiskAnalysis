!==============================================================================!
MODULE asym_gauss
use asym_cstes
!==============================================================================!
!
! Module � int�grer pour le calcul d'int�grales par la m�thode de Gauss
! Contient les fonctions : theta, aalind, lindpi0, lindpi2, lindpi4
!
!==============================================================================!
  IMPLICIT NONE
!==============================================================================!

CONTAINS


!==============================================================================!
      SUBROUTINE GSET(AX,BX,Ndim,Z,W)
!==============================================================================!
!     N-POINT GAUSS-Legendre ZEROS AND WEIGHTS FOR THE INTERVAL (AX,BX) ARE
!           STORED IN  ARRAYS Z AND W RESPECTIVELY.
!
!==============================================================================!
    IMPLICIT NONE

    REAL(pr), INTENT(in) :: ax,bx

    INTEGER :: n,ndim
    REAL(pr), INTENT(out) :: z(ndim),w(ndim)
    REAL(pr) , DIMENSION(ndim) :: a,x
    !INTEGER, DIMENSION(ndimmax) :: ktab
    REAL(pr) :: alpha,beta,delta,wtemp
    INTEGER :: j,jmid,jp,jtab,k,m,zn
 !
!
!-----TABLE OF INITIAL SUBSCRIPTS FOR N=2(1)16(4)96
!!$      DATA KTAB(2)/1/
!!$      DATA KTAB(3)/2/
!!$      DATA KTAB(4)/4/
!!$      DATA KTAB(5)/6/
!!$      DATA KTAB(6)/9/
!!$      DATA KTAB(7)/12/
!!$      DATA KTAB(8)/16/
!!$      DATA KTAB(9)/20/
!!$      DATA KTAB(10)/25/
!!$      DATA KTAB(11)/30/
!!$      DATA KTAB(12)/36/
!!$      DATA KTAB(13)/42/
!!$      DATA KTAB(14)/49/
!!$      DATA KTAB(15)/56/
!!$      DATA KTAB(16)/64/
!!$      DATA KTAB(20)/72/
!!$      DATA KTAB(24)/82/
!!$      DATA KTAB(28)/82/
!!$      DATA KTAB(32)/94/
!!$      DATA KTAB(36)/94/
!!$      DATA KTAB(40)/110/
!!$      DATA KTAB(44)/110/
!!$      DATA KTAB(48)/130/
!!$      DATA KTAB(52)/130/
!!$      DATA KTAB(56)/130/
!!$      DATA KTAB(60)/130/
!!$      DATA KTAB(64)/154/
!!$      DATA KTAB(68)/154/
!!$      DATA KTAB(72)/154/
!!$      DATA KTAB(76)/154/
!!$      DATA KTAB(80)/186/
!!$      DATA KTAB(84)/186/
!!$      DATA KTAB(88)/186/
!!$      DATA KTAB(92)/186/
!!$      DATA KTAB(96)/226/
!!$!
!!$!-----TABLE OF ABSCISSAE (X) AND WEIGHTS (A) FOR INTERVAL (-1,+1).
!!$!
!!$!**** N=2
!!$      DATA X(1)  /0.577350269189626_pr/, A(1)  /1.000000000000000_pr/
!!$!**** N=3
!!$      DATA X(2)  /0.774596669241483_pr/, A(2)  /0.555555555555556_pr/
!!$      DATA X(3)  /0.000000000000000_pr/, A(3)  /0.888888888888889_pr/
!!$!**** N=4
!!$      DATA X(4)  /0.861136311594053_pr/, A(4)  /0.347854845137454_pr/
!!$      DATA X(5)  /0.339981043584856_pr/, A(5)  /0.652145154862546_pr/
!!$!**** N=5
!!$      DATA X(6)  /0.906179845938664_pr/, A(6)  /0.236926885056189_pr/
!!$      DATA X(7)  /0.538469310105683_pr/, A(7)  /0.478628670499366_pr/
!!$      DATA X(8)  /0.000000000000000_pr/, A(8)  /0.568888888888889_pr/
!!$!**** N=6
!!$      DATA X(9)  /0.932469514203152_pr/, A(9)  /0.171324492379170_pr/
!!$      DATA X(10) /0.661209386466265_pr/, A(10) /0.360761573048139_pr/
!!$      DATA X(11) /0.238619186083197_pr/, A(11) /0.467913934572691_pr/
!!$!**** N=7
!!$      DATA X(12) /0.949107912342759_pr/, A(12) /0.129484966168870_pr/
!!$      DATA X(13) /0.741531185599394_pr/, A(13) /0.279705391489277_pr/
!!$      DATA X(14) /0.405845151377397_pr/, A(14) /0.381830050505119_pr/
!!$      DATA X(15) /0.000000000000000_pr/, A(15) /0.417959183673469_pr/
!!$!**** N=8
!!$      DATA X(16) /0.960289856497536_pr/, A(16) /0.101228536290376_pr/
!!$      DATA X(17) /0.796666477413627_pr/, A(17) /0.222381034453374_pr/
!!$      DATA X(18) /0.525532409916329_pr/, A(18) /0.313706645877887_pr/
!!$      DATA X(19) /0.183434642495650_pr/, A(19) /0.362683783378362_pr/
!!$!**** N=9
!!$      DATA X(20) /0.968160239507626_pr/, A(20) /0.081274388361574_pr/
!!$      DATA X(21) /0.836031107326636_pr/, A(21) /0.180648160694857_pr/
!!$      DATA X(22) /0.613371432700590_pr/, A(22) /0.260610696402935_pr/
!!$      DATA X(23) /0.324253423403809_pr/, A(23) /0.312347077040003_pr/
!!$      DATA X(24) /0.000000000000000_pr/, A(24) /0.330239355001260_pr/
!!$!**** N=10
!!$      DATA X(25) /0.973906528517172_pr/, A(25) /0.066671344308688_pr/
!!$      DATA X(26) /0.865063366688985_pr/, A(26) /0.149451349150581_pr/
!!$      DATA X(27) /0.679409568299024_pr/, A(27) /0.219086362515982_pr/
!!$      DATA X(28) /0.433395394129247_pr/, A(28) /0.269266719309996_pr/
!!$      DATA X(29) /0.148874338981631_pr/, A(29) /0.295524224714753_pr/
!!$!**** N=11
!!$      DATA X(30) /0.978228658146057_pr/, A(30) /0.055668567116174_pr/
!!$      DATA X(31) /0.887062599768095_pr/, A(31) /0.125580369464905_pr/
!!$      DATA X(32) /0.730152005574049_pr/, A(32) /0.186290210927734_pr/
!!$      DATA X(33) /0.519096129206812_pr/, A(33) /0.233193764591990_pr/
!!$      DATA X(34) /0.269543155952345_pr/, A(34) /0.262804544510247_pr/
!!$      DATA X(35) /0.000000000000000_pr/, A(35) /0.272925086777901_pr/
!!$!**** N=12
!!$      DATA X(36) /0.981560634246719_pr/, A(36) /0.047175336386512_pr/
!!$      DATA X(37) /0.904117256370475_pr/, A(37) /0.106939325995318_pr/
!!$      DATA X(38) /0.769902674194305_pr/, A(38) /0.160078328543346_pr/
!!$      DATA X(39) /0.587317954286617_pr/, A(39) /0.203167426723066_pr/
!!$      DATA X(40) /0.367831498998180_pr/, A(40) /0.233492536538355_pr/
!!$      DATA X(41) /0.125233408511469_pr/, A(41) /0.249147045813403_pr/
!!$!**** N=13
!!$      DATA X(42) /0.984183054718588_pr/, A(42) /0.040484004765316_pr/
!!$      DATA X(43) /0.917598399222978_pr/, A(43) /0.092121499837728_pr/
!!$      DATA X(44) /0.801578090733310_pr/, A(44) /0.138873510219787_pr/
!!$      DATA X(45) /0.642349339440340_pr/, A(45) /0.178145980761946_pr/
!!$      DATA X(46) /0.448492751036447_pr/, A(46) /0.207816047536889_pr/
!!$      DATA X(47) /0.230458315955135_pr/, A(47) /0.226283180262897_pr/
!!$      DATA X(48) /0.000000000000000_pr/, A(48) /0.232551553230874_pr/
!!$!**** N=14
!!$      DATA X(49) /0.986283808696812_pr/, A(49) /0.035119460331752_pr/
!!$      DATA X(50) /0.928434883663574_pr/, A(50) /0.080158087159760_pr/
!!$      DATA X(51) /0.827201315069765_pr/, A(51) /0.121518570687903_pr/
!!$      DATA X(52) /0.687292904811685_pr/, A(52) /0.157203167158194_pr/
!!$      DATA X(53) /0.515248636358154_pr/, A(53) /0.185538397477938_pr/
!!$      DATA X(54) /0.319112368927890_pr/, A(54) /0.205198463721296_pr/
!!$      DATA X(55) /0.108054948707344_pr/, A(55) /0.215263853463158_pr/
!!$!**** N=15
!!$      DATA X(56) /0.987992518020485_pr/, A(56) /0.030753241996117_pr/
!!$      DATA X(57) /0.937273392400706_pr/, A(57) /0.070366047488108_pr/
!!$      DATA X(58) /0.848206583410427_pr/, A(58) /0.107159220467172_pr/
!!$      DATA X(59) /0.724417731360170_pr/, A(59) /0.139570677926154_pr/
!!$      DATA X(60) /0.570972172608539_pr/, A(60) /0.166269205816994_pr/
!!$      DATA X(61) /0.394151347077563_pr/, A(61) /0.186161000015562_pr/
!!$      DATA X(62) /0.201194093997435_pr/, A(62) /0.198431485327111_pr/
!!$      DATA X(63) /0.000000000000000_pr/, A(63) /0.202578241925561_pr/
!!$!**** N=16
!!$      DATA X(64) /0.989400934991650_pr/, A(64) /0.027152459411754_pr/
!!$      DATA X(65) /0.944575023073233_pr/, A(65) /0.062253523938648_pr/
!!$      DATA X(66) /0.865631202387832_pr/, A(66) /0.095158511682493_pr/
!!$      DATA X(67) /0.755404408355003_pr/, A(67) /0.124628971255534_pr/
!!$      DATA X(68) /0.617876244402644_pr/, A(68) /0.149595988816577_pr/
!!$      DATA X(69) /0.458016777657227_pr/, A(69) /0.169156519395003_pr/
!!$      DATA X(70) /0.281603550779259_pr/, A(70) /0.182603415044924_pr/
!!$      DATA X(71) /0.095012509837637_pr/, A(71) /0.189450610455069_pr/
!!$!**** N=20
!!$      DATA X(72) /0.993128599185094_pr/, A(72) /0.017614007139152_pr/
!!$      DATA X(73) /0.963971927277913_pr/, A(73) /0.040601429800386_pr/
!!$      DATA X(74) /0.912234428251325_pr/, A(74) /0.062672048334109_pr/
!!$      DATA X(75) /0.839116971822218_pr/, A(75) /0.083276741576704_pr/
!!$      DATA X(76) /0.746331906460150_pr/, A(76) /0.101930119817240_pr/
!!$      DATA X(77) /0.636053680726515_pr/, A(77) /0.118194531961518_pr/
!!$      DATA X(78) /0.510867001950827_pr/, A(78) /0.131688638449176_pr/
!!$      DATA X(79) /0.373706088715419_pr/, A(79) /0.142096109318382_pr/
!!$      DATA X(80) /0.227785851141645_pr/, A(80) /0.149172986472603_pr/
!!$      DATA X(81) /0.076526521133497_pr/, A(81) /0.152753387130725_pr/
!!$!**** N=24
!!$      DATA X(82) /0.995187219997021_pr/, A(82) /0.012341229799987_pr/
!!$      DATA X(83) /0.974728555971309_pr/, A(83) /0.028531388628933_pr/
!!$      DATA X(84) /0.938274552002732_pr/, A(84) /0.044277438817419_pr/
!!$      DATA X(85) /0.886415527004401_pr/, A(85) /0.059298584915436_pr/
!!$      DATA X(86) /0.820001985973902_pr/, A(86) /0.073346481411080_pr/
!!$      DATA X(87) /0.740124191578554_pr/, A(87) /0.086190161531953_pr/
!!$      DATA X(88) /0.648093651936975_pr/, A(88) /0.097618652104113_pr/
!!$      DATA X(89) /0.545421471388839_pr/, A(89) /0.107444270115965_pr/
!!$      DATA X(90) /0.433793507626045_pr/, A(90) /0.115505668053725_pr/
!!$      DATA X(91) /0.315042679696163_pr/, A(91) /0.121670472927803_pr/
!!$      DATA X(92) /0.191118867473616_pr/, A(92) /0.125837456346828_pr/
!!$      DATA X(93) /0.064056892862605_pr/, A(93) /0.127938195346752_pr/
!!$!**** N=32
!!$      DATA X(94) /0.997263861849481_pr/, A(94) /0.007018610009470_pr/
!!$      DATA X(95) /0.985611511545268_pr/, A(95) /0.016274394730905_pr/
!!$      DATA X(96) /0.964762255587506_pr/, A(96) /0.025392065309262_pr/
!!$      DATA X(97) /0.934906075937739_pr/, A(97) /0.034273862913021_pr/
!!$      DATA X(98) /0.896321155766052_pr/, A(98) /0.042835898022226_pr/
!!$      DATA X(99) /0.849367613732569_pr/, A(99) /0.050998059262376_pr/
!!$      DATA X(100)/0.794483795967942_pr/, A(100)/0.058684093478535_pr/
!!$      DATA X(101)/0.732182118740289_pr/, A(101)/0.065822222776361_pr/
!!$      DATA X(102)/0.663044266930215_pr/, A(102)/0.072345794108848_pr/
!!$      DATA X(103)/0.587715757240762_pr/, A(103)/0.078193895787070_pr/
!!$      DATA X(104)/0.506899908932229_pr/, A(104)/0.083311924226946_pr/
!!$      DATA X(105)/0.421351276130635_pr/, A(105)/0.087652093004403_pr/
!!$      DATA X(106)/0.331868602282127_pr/, A(106)/0.091173878695763_pr/
!!$      DATA X(107)/0.239287362252137_pr/, A(107)/0.093844399080804_pr/
!!$      DATA X(108)/0.144471961582796_pr/, A(108)/0.095638720079274_pr/
!!$      DATA X(109)/0.048307665687738_pr/, A(109)/0.096540088514727_pr/
!!$!**** N=40
!!$      DATA X(110)/0.998237709710559_pr/, A(110)/0.004521277098533_pr/
!!$      DATA X(111)/0.990726238699457_pr/, A(111)/0.010498284531152_pr/
!!$      DATA X(112)/0.977259949983774_pr/, A(112)/0.016421058381907_pr/
!!$      DATA X(113)/0.957916819213791_pr/, A(113)/0.022245849194166_pr/
!!$      DATA X(114)/0.932812808278676_pr/, A(114)/0.027937006980023_pr/
!!$      DATA X(115)/0.902098806968874_pr/, A(115)/0.033460195282547_pr/
!!$      DATA X(116)/0.865959503212259_pr/, A(116)/0.038782167974472_pr/
!!$      DATA X(117)/0.824612230833311_pr/, A(117)/0.043870908185673_pr/
!!$      DATA X(118)/0.778305651426519_pr/, A(118)/0.048695807635072_pr/
!!$      DATA X(119)/0.727318255189927_pr/, A(119)/0.053227846983936_pr/
!!$      DATA X(120)/0.671956684614179_pr/, A(120)/0.057439769099391_pr/
!!$      DATA X(121)/0.612553889667980_pr/, A(121)/0.061306242492928_pr/
!!$      DATA X(122)/0.549467125095128_pr/, A(122)/0.064804013456601_pr/
!!$      DATA X(123)/0.483075801686178_pr/, A(123)/0.067912045815233_pr/
!!$      DATA X(124)/0.413779204371605_pr/, A(124)/0.070611647391286_pr/
!!$      DATA X(125)/0.341994090825758_pr/, A(125)/0.072886582395804_pr/
!!$      DATA X(126)/0.268152185007253_pr/, A(126)/0.074723169057968_pr/
!!$      DATA X(127)/0.192697580701371_pr/, A(127)/0.076110361900626_pr/
!!$      DATA X(128)/0.116084070675255_pr/, A(128)/0.077039818164247_pr/
!!$      DATA X(129)/0.038772417506050_pr/, A(129)/0.077505947978424_pr/
!!$!**** N=48
!!$      DATA X(130)/0.998771007252426_pr/, A(130)/0.003153346052305_pr/
!!$      DATA X(131)/0.993530172266350_pr/, A(131)/0.007327553901276_pr/
!!$      DATA X(132)/0.984124583722826_pr/, A(132)/0.011477234579234_pr/
!!$      DATA X(133)/0.970591592546247_pr/, A(133)/0.015579315722943_pr/
!!$      DATA X(134)/0.952987703160430_pr/, A(134)/0.019616160457355_pr/
!!$      DATA X(135)/0.931386690706554_pr/, A(135)/0.023570760839324_pr/
!!$      DATA X(136)/0.905879136715569_pr/, A(136)/0.027426509708356_pr/
!!$      DATA X(137)/0.876572020274247_pr/, A(137)/0.031167227832798_pr/
!!$      DATA X(138)/0.843588261624393_pr/, A(138)/0.034777222564770_pr/
!!$      DATA X(139)/0.807066204029442_pr/, A(139)/0.038241351065830_pr/
!!$      DATA X(140)/0.767159032515740_pr/, A(140)/0.041545082943464_pr/
!!$      DATA X(141)/0.724034130923814_pr/, A(141)/0.044674560856694_pr/
!!$      DATA X(142)/0.677872379632663_pr/, A(142)/0.047616658492490_pr/
!!$      DATA X(143)/0.628867396776513_pr/, A(143)/0.050359035553854_pr/
!!$      DATA X(144)/0.577224726083972_pr/, A(144)/0.052890189485193_pr/
!!$      DATA X(145)/0.523160974722233_pr/, A(145)/0.055199503699984_pr/
!!$      DATA X(146)/0.466902904750958_pr/, A(146)/0.057277292100403_pr/
!!$      DATA X(147)/0.408686481990716_pr/, A(147)/0.059114839698395_pr/
!!$      DATA X(148)/0.348755886292160_pr/, A(148)/0.060704439165893_pr/
!!$      DATA X(149)/0.287362487355455_pr/, A(149)/0.062039423159892_pr/
!!$      DATA X(150)/0.224763790394689_pr/, A(150)/0.063114192286254_pr/
!!$      DATA X(151)/0.161222356068891_pr/, A(151)/0.063924238584648_pr/
!!$      DATA X(152)/0.097004699209462_pr/, A(152)/0.064466164435950_pr/
!!$      DATA X(153)/0.032380170962869_pr/, A(153)/0.064737696812683_pr/
!!$!**** N=64
!!$      DATA X(154)/0.999305041735772_pr/, A(154)/0.001783280721696_pr/
!!$      DATA X(155)/0.996340116771955_pr/, A(155)/0.004147033260562_pr/
!!$      DATA X(156)/0.991013371476744_pr/, A(156)/0.006504457968978_pr/
!!$      DATA X(157)/0.983336253884625_pr/, A(157)/0.008846759826363_pr/
!!$      DATA X(158)/0.973326827789910_pr/, A(158)/0.011168139460131_pr/
!!$      DATA X(159)/0.961008799652053_pr/, A(159)/0.013463047896718_pr/
!!$      DATA X(160)/0.946411374858402_pr/, A(160)/0.015726030476024_pr/
!!$      DATA X(161)/0.929569172131939_pr/, A(161)/0.017951715775697_pr/
!!$      DATA X(162)/0.910522137078502_pr/, A(162)/0.020134823153530_pr/
!!$      DATA X(163)/0.889315445995114_pr/, A(163)/0.022270173808383_pr/
!!$      DATA X(164)/0.865999398154092_pr/, A(164)/0.024352702568710_pr/
!!$      DATA X(165)/0.840629296252580_pr/, A(165)/0.026377469715054_pr/
!!$      DATA X(166)/0.813265315122797_pr/, A(166)/0.028339672614259_pr/
!!$      DATA X(167)/0.783972358943341_pr/, A(167)/0.030234657072402_pr/
!!$      DATA X(168)/0.752819907260531_pr/, A(168)/0.032057928354851_pr/
!!$      DATA X(169)/0.719881850171610_pr/, A(169)/0.033805161837141_pr/
!!$      DATA X(170)/0.685236313054233_pr/, A(170)/0.035472213256882_pr/
!!$      DATA X(171)/0.648965471254657_pr/, A(171)/0.037055128540240_pr/
!!$      DATA X(172)/0.611155355172393_pr/, A(172)/0.038550153178615_pr/
!!$      DATA X(173)/0.571895646202634_pr/, A(173)/0.039953741132720_pr/
!!$      DATA X(174)/0.531279464019894_pr/, A(174)/0.041262563242623_pr/
!!$      DATA X(175)/0.489403145707052_pr/, A(175)/0.042473515123653_pr/
!!$      DATA X(176)/0.446366017253464_pr/, A(176)/0.043583724529323_pr/
!!$      DATA X(177)/0.402270157963991_pr/, A(177)/0.044590558163756_pr/
!!$      DATA X(178)/0.357220158337668_pr/, A(178)/0.045491627927418_pr/
!!$      DATA X(179)/0.311322871990210_pr/, A(179)/0.046284796581314_pr/
!!$      DATA X(180)/0.264687162208767_pr/, A(180)/0.046968182816210_pr/
!!$      DATA X(181)/0.217423643740007_pr/, A(181)/0.047540165714830_pr/
!!$      DATA X(182)/0.169644420423992_pr/, A(182)/0.047999388596458_pr/
!!$      DATA X(183)/0.121462819296120_pr/, A(183)/0.048344762234802_pr/
!!$      DATA X(184)/0.072993121787799_pr/, A(184)/0.048575467441503_pr/
!!$      DATA X(185)/0.024350292663424_pr/, A(185)/0.048690957009139_pr/
!!$!**** N=80
!!$      DATA X(186)/0.999553822651630_pr/, A(186)/0.001144950003186_pr/
!!$      DATA X(187)/0.997649864398237_pr/, A(187)/0.002663533589512_pr/
!!$      DATA X(188)/0.994227540965688_pr/, A(188)/0.004180313124694_pr/
!!$      DATA X(189)/0.989291302499755_pr/, A(189)/0.005690922451403_pr/
!!$      DATA X(190)/0.982848572738629_pr/, A(190)/0.007192904768117_pr/
!!$      DATA X(191)/0.974909140585727_pr/, A(191)/0.008683945269260_pr/
!!$      DATA X(192)/0.965485089043799_pr/, A(192)/0.010161766041103_pr/
!!$      DATA X(193)/0.954590766343634_pr/, A(193)/0.011624114120797_pr/
!!$      DATA X(194)/0.942242761309872_pr/, A(194)/0.013068761592401_pr/
!!$      DATA X(195)/0.928459877172445_pr/, A(195)/0.014493508040509_pr/
!!$      DATA X(196)/0.913263102571757_pr/, A(196)/0.015896183583725_pr/
!!$      DATA X(197)/0.896675579438770_pr/, A(197)/0.017274652056269_pr/
!!$      DATA X(198)/0.878722567678213_pr/, A(198)/0.018626814208299_pr/
!!$      DATA X(199)/0.859431406663111_pr/, A(199)/0.019950610878141_pr/
!!$      DATA X(200)/0.838831473580255_pr/, A(200)/0.021244026115782_pr/
!!$      DATA X(201)/0.816954138681463_pr/, A(201)/0.022505090246332_pr/
!!$      DATA X(202)/0.793832717504605_pr/, A(202)/0.023731882865930_pr/
!!$      DATA X(203)/0.769502420135041_pr/, A(203)/0.024922535764115_pr/
!!$      DATA X(204)/0.744000297583597_pr/, A(204)/0.026075235767565_pr/
!!$      DATA X(205)/0.717365185362099_pr/, A(205)/0.027188227500486_pr/
!!$      DATA X(206)/0.689637644342027_pr/, A(206)/0.028259816057276_pr/
!!$      DATA X(207)/0.660859898986119_pr/, A(207)/0.029288369583267_pr/
!!$      DATA X(208)/0.631075773046871_pr/, A(208)/0.030272321759557_pr/
!!$      DATA X(209)/0.600330622829751_pr/, A(209)/0.031210174188114_pr/
!!$      DATA X(210)/0.568671268122709_pr/, A(210)/0.032100498673487_pr/
!!$      DATA X(211)/0.536145920897131_pr/, A(211)/0.032941939397645_pr/
!!$      DATA X(212)/0.502804111888784_pr/, A(212)/0.033733214984611_pr/
!!$      DATA X(213)/0.468696615170544_pr/, A(213)/0.034473120451753_pr/
!!$      DATA X(214)/0.433875370831756_pr/, A(214)/0.035160529044747_pr/
!!$      DATA X(215)/0.398393405881969_pr/, A(215)/0.035794393953416_pr/
!!$      DATA X(216)/0.362304753499487_pr/, A(216)/0.036373749905835_pr/
!!$      DATA X(217)/0.325664370747701_pr/, A(217)/0.036897714638276_pr/
!!$      DATA X(218)/0.288528054884511_pr/, A(218)/0.037365490238730_pr/
!!$      DATA X(219)/0.250952358392272_pr/, A(219)/0.037776364362001_pr/
!!$      DATA X(220)/0.212994502857666_pr/, A(220)/0.038129711314477_pr/
!!$      DATA X(221)/0.174712291832646_pr/, A(221)/0.038424993006959_pr/
!!$      DATA X(222)/0.136164022809143_pr/, A(222)/0.038661759774076_pr/
!!$      DATA X(223)/0.097408398441584_pr/, A(223)/0.038839651059051_pr/
!!$      DATA X(224)/0.058504437152420_pr/, A(224)/0.038958395962769_pr/
!!$      DATA X(225)/0.019511383256793_pr/, A(225)/0.039017813656306_pr/
!!$!**** N=96
!!$      DATA X(226)/0.999689503883230_pr/, A(226)/0.000796792065552_pr/
!!$      DATA X(227)/0.998364375863181_pr/, A(227)/0.001853960788946_pr/
!!$      DATA X(228)/0.995981842987209_pr/, A(228)/0.002910731817934_pr/
!!$      DATA X(229)/0.992543900323762_pr/, A(229)/0.003964554338444_pr/
!!$      DATA X(230)/0.988054126329623_pr/, A(230)/0.005014202742927_pr/
!!$      DATA X(231)/0.982517263563014_pr/, A(231)/0.006058545504235_pr/
!!$      DATA X(232)/0.975939174585136_pr/, A(232)/0.007096470791153_pr/
!!$      DATA X(233)/0.968326828463264_pr/, A(233)/0.008126876925698_pr/
!!$      DATA X(234)/0.959688291448742_pr/, A(234)/0.009148671230783_pr/
!!$      DATA X(235)/0.950032717784437_pr/, A(235)/0.010160770535008_pr/
!!$      DATA X(236)/0.939370339752755_pr/, A(236)/0.011162102099838_pr/
!!$      DATA X(237)/0.927712456722308_pr/, A(237)/0.012151604671088_pr/
!!$      DATA X(238)/0.915071423120898_pr/, A(238)/0.013128229566961_pr/
!!$      DATA X(239)/0.901460635315852_pr/, A(239)/0.014090941772314_pr/
!!$      DATA X(240)/0.886894517402420_pr/, A(240)/0.015038721026994_pr/
!!$      DATA X(241)/0.871388505909296_pr/, A(241)/0.015970562902562_pr/
!!$      DATA X(242)/0.854959033434601_pr/, A(242)/0.016885479864245_pr/
!!$      DATA X(243)/0.837623511228187_pr/, A(243)/0.017782502316045_pr/
!!$      DATA X(244)/0.819400310737931_pr/, A(244)/0.018660679627411_pr/
!!$      DATA X(245)/0.800308744139140_pr/, A(245)/0.019519081140145_pr/
!!$      DATA X(246)/0.780369043867433_pr/, A(246)/0.020356797154333_pr/
!!$      DATA X(247)/0.759602341176647_pr/, A(247)/0.021172939892191_pr/
!!$      DATA X(248)/0.738030643744400_pr/, A(248)/0.021966644438744_pr/
!!$      DATA X(249)/0.715676812348967_pr/, A(249)/0.022737069658329_pr/
!!$      DATA X(250)/0.692564536642171_pr/, A(250)/0.023483399085926_pr/
!!$      DATA X(251)/0.668718310043916_pr/, A(251)/0.024204841792364_pr/
!!$      DATA X(252)/0.644163403784967_pr/, A(252)/0.024900633222483_pr/
!!$      DATA X(253)/0.618925840125468_pr/, A(253)/0.025570036005349_pr/
!!$      DATA X(254)/0.593032364777572_pr/, A(254)/0.026212340735672_pr/
!!$      DATA X(255)/0.566510418561397_pr/, A(255)/0.026826866725591_pr/
!!$      DATA X(256)/0.539388108324357_pr/, A(256)/0.027412962726029_pr/
!!$      DATA X(257)/0.511694177154667_pr/, A(257)/0.027970007616848_pr/
!!$      DATA X(258)/0.483457973920596_pr/, A(258)/0.028497411065085_pr/
!!$      DATA X(259)/0.454709422167743_pr/, A(259)/0.028994614150555_pr/
!!$      DATA X(260)/0.425478988407300_pr/, A(260)/0.029461089958167_pr/
!!$      DATA X(261)/0.395797649828908_pr/, A(261)/0.029896344136328_pr/
!!$      DATA X(262)/0.365696861472313_pr/, A(262)/0.030299915420827_pr/
!!$      DATA X(263)/0.335208522892625_pr/, A(263)/0.030671376123669_pr/
!!$      DATA X(264)/0.304364944354496_pr/, A(264)/0.031010332586313_pr/
!!$      DATA X(265)/0.273198812591049_pr/, A(265)/0.031316425596861_pr/
!!$      DATA X(266)/0.241743156163840_pr/, A(266)/0.031589330770727_pr/
!!$      DATA X(267)/0.210031310460567_pr/, A(267)/0.031828758894411_pr/
!!$      DATA X(268)/0.178096882367618_pr/, A(268)/0.032034456231992_pr/
!!$      DATA X(269)/0.145973714654896_pr/, A(269)/0.032206204794030_pr/
!!$      DATA X(270)/0.113695850110665_pr/, A(270)/0.032343822568575_pr/
!!$      DATA X(271)/0.081297495464425_pr/, A(271)/0.032447163714064_pr/
!!$      DATA X(272)/0.048812985136049_pr/, A(272)/0.032516118713868_pr/
!!$      DATA X(273)/0.016276744849602_pr/, A(273)/0.032550614492363_pr/
!
      Z(:)=0.0_pr
      W(:)=0.0_pr  
      
      x(:)=glp(:)
      a(:)=glw(:)

      n=ndim

!-----
      if(ax.eq.bx)then
       Z(:)=0.0_pr
       W(:)=0.0_pr
       return
       write(*,*)'uscito?'
      endif

!----- TEST N
      ALPHA=0.5_pr*(AX+BX)
      BETA=0.5_pr*(BX-AX)
      IF( N.LT.1 ) GO TO 100
      IF(N.NE.1) GO TO 1
      Z(1)=ALPHA
      W(1)=BX-AX
      RETURN
!
    1 IF (N.LE.16) GO TO 3
      IF (N.GT.24) GO TO 4
      N=4*(N/4)
      GO TO 3
    4 IF (N.GT.48) GO TO 5
      N=8*(N/8)
      GO TO 3
    5 N=16*(N/16)
!
!----- SET K EQUAL TO INITIAL SUBSCRIPT AND STORE RESULTS
    3 CONTINUE
      K=1!KTAB(N)
      M=N/2
      DO J=1,M
      JTAB=J!K-1+J 
      WTEMP=BETA*A(JTAB)
      DELTA=BETA*X(JTAB)
      Z(J)=ALPHA-DELTA
      W(J)=WTEMP
      JP=N+1-J
      Z(JP)=ALPHA+DELTA
      W(JP)=WTEMP
      ENDDO
      IF((N-M-M).EQ.0) RETURN
      Z(M+1)=ALPHA
      JMID=K+M
      W(M+1)=BETA*A(JMID)
      RETURN
!
  100 ZN=N
      WRITE(6,*) ZN,AX,BX
!==============================================================================!
  END SUBROUTINE gset


!
!!
!
!==============================================================================!
!
      SUBROUTINE LEGPV(PJ,PJM1,X,J,N)
!
!        SUBROUTINE LEGP  COMPUTES THE LEGENDRE POLYNOMS
!
    IMPLICIT NONE
!
    INTEGER, INTENT(in) :: j,n
    REAL(pr),DIMENSION(n), INTENT(in) :: x
    REAL(pr),DIMENSION(n), INTENT(out) :: pj,pjm1
    INTEGER :: i,k
    REAL(pr) :: a,b
!
!
!        COMPUTE LEGENDRE PLOYNOM FOR J EQUALS ZERO
!
      DO  10  I = 1 , N
          PJM1(I) = 1.0_pr
   10 CONTINUE
!
      IF (.NOT. J .GT. 0) THEN
         DO  20  I = 1 , N
             PJ(I) = 1.0_pr
   20    CONTINUE
      ELSE
!
!        COMPUTE LEGENDRE POLYNOM FOR J EQUALS ONE
!
         DO  30  I = 1 , N
             PJ(I) = X(I)
   30    CONTINUE
!
!        COMPUTE LEGENDRE POLYNOM FOR J GREATER OR EQUAL TWO
!
         IF (J .GT. 1) THEN
            DO  50  K = 1 , N
                DO  40  I = 2 , J
                    A       = X(K)*PJ(K)
                    B       = A-PJM1(K)
                    PJM1(K) = PJ(K)
                    PJ(K)   = -B/ FLOAT(I)+B+A
   40           CONTINUE
   50       CONTINUE
         END IF
      END IF
!
!==============================================================================!
  END SUBROUTINE legpv
!==============================================================================!

      SUBROUTINE LEGP (PJ,PJM1,X,J)
!
!        SUBROUTINE LEGP   COMPUTES THE LEGENDRE POLYNOMS
!
    IMPLICIT NONE

    INTEGER, INTENT(in) :: j
    REAL(pr),INTENT(in) :: x

    REAL(pr),INTENT(out) :: pj,pjm1
    INTEGER :: i
    REAL(pr) :: a,b
!
!        COMPUTE LEGENDRE POLYNOM FOR J EQUALS ZERO
!
      PJM1=1.0_pr
      IF (J.GT.0) GO TO 1
      PJ=1.0_pr
      RETURN
!
!        COMPUTE LEGENDRE POLYNOMS FOR J EQUALS ONE
!
    1 PJ=X
      IF (J.EQ.1) RETURN
!
!        COMPUTE LEGENDRE POLYNOM FOR J GREATER OR EQUAL TWO
!
      DO  I=2,J
      A=X*PJ
      B=A-PJM1
      PJM1=PJ
      PJ=-B/ FLOAT(I)+B+A
      ENDDO
!
!==============================================================================!
  END SUBROUTINE legp
!==============================================================================!


!==============================================================================!
END module asym_gauss
!==============================================================================!













