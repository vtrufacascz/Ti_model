      SUBROUTINE KERRTI(MIRREQ,DOUT)
C////////////////////////////////coefficients - error Ti model part//////////////////////
      REAL DOUT(4,3,81)
      INTEGER MIRREQ(81),I,J,K
      REAL DERRTI(4,3,81)
C     350km equinox
      DATA (DERRTI(1,1,J),J=1,81)/ 1.18564E+02,
     &        3.18751E-14, 3.59164E+01, 2.79981E-14, 5.95898E+01,
     &       -1.11370E-13, 2.79364E+01, 6.57277E-14,-8.23019E+00,
     &       -1.80917E+01, 3.42510E-16, 9.48300E-01, 5.16058E-17,
     &       -7.09536E-01,-1.94405E-15,-6.69293E-01, 6.12856E-16,
     &       -1.58039E+01, 2.41275E-15, 3.72558E+00, 4.94600E-16,
     &        6.66242E+00,-1.15642E-15, 2.93621E+00,-1.89934E-15,
     &       -1.44546E+00,-1.53196E-15, 1.13341E+00,-2.04757E-16,
     &       -1.33517E+00, 2.64274E-16,-6.99407E-01, 3.15330E+00,
     &       -3.60806E-17,-1.82733E+00,-5.09619E-17, 1.57457E-01,
     &        1.04994E-16, 2.53836E-01, 3.20500E+00, 1.64755E-16,
     &       -1.10286E-01, 2.50994E-16, 5.44585E-01,-3.77891E-16,
     &        3.20647E-01, 5.28060E-17,-7.51849E-03, 2.31451E-16,
     &        5.51323E-02, 6.16019E-17, 3.09390E+00,-2.32050E-16,
     &       -1.24613E+00, 2.18079E-17, 1.03237E-01, 1.16614E+00,
     &        3.15034E-16, 1.72207E+00, 5.61521E-16, 2.05215E-01,
     &        4.10359E-01, 1.02359E-16, 1.28407E-01, 1.95844E-16,
     &       -1.56232E+00, 1.39189E-16,-4.61026E-01,-2.06581E-16,
     &       -1.32577E+00,-2.65229E-16, 3.34807E-01, 1.94837E+00,
     &       -3.12076E-17, 1.13317E-01, 3.31384E-01,-1.91889E-16,
     &        1.62496E+00, 1.71324E-16, 6.52200E-01,-1.07607E+00/
C     350km June solstice
      DATA (DERRTI(1,2,J),J=1,81)/ 1.07574E+02,
     &       -7.72258E-01, 2.32614E+01, 1.24898E+01, 5.31720E+01,
     &       -1.46889E+01,-1.70653E+00, 4.96241E+00, 1.43538E+01,
     &       -1.72787E+01, 6.74264E+00, 2.83005E+00,-3.54154E+00,
     &       -3.11546E+00,-7.98558E-01,-6.71766E-01,-4.56043E-01,
     &       -1.11616E+01, 5.69570E+00, 3.30846E+00,-9.64915E-01,
     &        4.99295E+00,-9.19128E-01, 2.70966E+00, 3.43650E-01,
     &       -7.19993E+00,-5.17564E+00,-1.39794E+00,-4.42421E-01,
     &       -7.03531E-01,-1.58199E-01,-1.73569E-01, 9.24434E+00,
     &        6.13910E-01, 1.88814E-01, 7.68338E-01, 2.55643E-01,
     &       -2.39615E-01, 2.09864E-01, 3.48734E+00,-1.67166E-01,
     &       -5.56985E-01,-6.40148E-02, 7.50126E-02, 1.86245E-01,
     &        1.36568E+00,-9.78133E-01, 6.97798E-01,-3.80973E-01,
     &       -1.67794E-01,-1.57859E-02, 4.74883E+00,-5.68822E-01,
     &       -1.10808E+00,-2.18929E-01, 1.47298E-02,-5.26333E+00,
     &       -6.13654E-01,-1.49004E-01,-2.81596E-01,-1.42796E-01,
     &        1.95396E-01, 1.19278E-01,-9.44040E-02,-9.41098E-02,
     &       -1.13718E+00, 3.56602E-01, 3.99941E-02, 2.26460E-01,
     &       -8.74859E-01,-2.50866E-02, 1.46838E-01, 2.19310E+00,
     &       -4.75867E-02, 4.77529E-01, 4.26530E-01, 6.51821E-01,
     &       -1.23607E+00,-1.36776E-01, 8.14853E-01,-1.08404E+00/
C     430km equinox
      DATA (DERRTI(2,1,J),J=1,81)/ 1.21647E+02,
     &       -2.06692E-14, 4.70328E+01, 3.32685E-14, 5.99539E+01,
     &       -7.06774E-14, 1.04145E+01, 9.40271E-14, 1.62029E+01,
     &       -9.31278E+00,-2.64716E-16, 1.82682E-01,-2.67932E-16,
     &       -6.03486E-01,-1.99704E-15, 9.05546E-01, 1.07211E-15,
     &        5.30366E-01, 1.23598E-15, 3.46088E+00,-1.40831E-15,
     &        7.34456E+00,-1.26526E-15, 3.99194E+00, 2.28787E-16,
     &       -4.63520E+00,-2.25803E-15, 3.23866E+00, 7.62787E-16,
     &       -1.24390E+00,-1.35829E-16,-5.65174E-01,-2.34324E+00,
     &        3.84218E-16,-3.81166E+00,-7.86551E-16, 8.45964E-01,
     &        8.26804E-16, 1.10637E-01, 1.18986E+00, 8.29562E-16,
     &       -2.20685E+00, 2.94483E-16, 3.02831E-01,-3.90945E-16,
     &       -2.67467E+00,-5.33985E-16, 1.48632E+00, 3.46196E-16,
     &       -3.94140E-01, 2.08699E-16, 2.38084E+00, 5.35550E-16,
     &       -5.04958E-01,-1.13483E-16, 1.12453E-01, 4.30663E+00,
     &        4.78772E-16, 4.66434E-01, 1.72272E-16, 8.04668E-02,
     &        5.87160E-01,-2.41875E-16, 1.89272E-01, 1.74341E-16,
     &       -1.80391E+00,-4.01537E-16,-2.43615E-01,-1.85997E-16,
     &        2.43836E+00,-5.05613E-16, 4.77041E-01,-1.16041E+00,
     &       -2.09567E-16, 4.62420E-01,-2.02920E+00,-6.73698E-17,
     &        2.38431E+00,-1.26168E-17,-1.12194E+00,-4.86171E-01/
C     430km June solstice
      DATA (DERRTI(2,2,J),J=1,81)/ 1.11251E+02,
     &        1.25995E+00, 2.30368E+01, 6.94590E-01, 5.34230E+01,
     &       -1.29723E+01,-2.70927E+01,-8.22855E+00, 2.35017E+01,
     &       -1.16840E+01, 4.50792E+00, 3.14085E+00,-4.29104E+00,
     &       -3.87568E+00,-3.12458E+00, 4.93011E-01,-7.49196E-01,
     &        2.52959E-01, 5.97611E+00,-1.92354E+00,-1.03596E+00,
     &        3.13092E+00, 1.11003E+00,-6.40338E-01,-5.37781E-02,
     &       -8.93062E+00,-4.23534E-01, 1.35469E+00, 7.45449E-01,
     &       -9.89864E-01, 4.86793E-01, 3.54719E-01, 4.85628E+00,
     &        1.07172E+00, 1.19581E+00, 9.64873E-01,-4.06925E-01,
     &        4.38580E-03,-1.43846E-01, 4.76234E+00, 1.77481E-01,
     &       -6.75945E-01, 2.24885E-01,-6.17006E-03,-2.18917E-01,
     &       -5.90814E+00,-1.01195E+00, 2.37143E-01,-1.70297E-01,
     &        6.87353E-03, 4.13783E-02, 6.50062E+00, 8.41616E-01,
     &       -8.91077E-01, 1.08171E-01, 1.01300E-01, 1.34543E+00,
     &       -6.39061E-01,-2.76283E-02,-1.51646E-01,-2.37588E-02,
     &        1.21243E+00, 6.50141E-01, 5.80030E-01,-4.25991E-02,
     &        1.75884E+00,-1.17163E+00,-6.50553E-03,-3.26797E-01,
     &       -2.98915E+00, 1.04993E+00,-3.90346E-02,-2.54316E+00,
     &        1.81616E-01,-9.52828E-02,-6.55020E-01,-3.24229E-01,
     &       -1.93405E+00,-1.01143E-01, 7.90640E-01, 7.64876E-01/
C     600km equinox
      DATA (DERRTI(3,1,J),J=1,81)/ 1.23383E+02,
     &        1.12497E-14, 9.61563E+01,-3.46605E-14, 5.81021E+01,
     &       -4.29513E-14,-2.61454E+01, 7.19187E-14,-1.47740E+01,
     &        2.32403E+00,-6.95805E-15, 7.11185E+00, 1.46171E-15,
     &       -3.19486E+00,-4.23749E-16, 3.61459E+00,-1.05103E-15,
     &        9.43126E+00,-3.20193E-16, 1.80498E+01,-3.36473E-15,
     &        1.16200E+01,-7.38766E-15, 1.54538E+01, 5.96812E-15,
     &       -1.67483E+01,-1.73913E-15,-3.32677E+00,-6.27806E-16,
     &       -1.76122E+00,-1.58462E-17,-2.39731E+00,-2.48407E+00,
     &       -6.13267E-16, 7.64519E-01, 4.70514E-16, 5.17077E-01,
     &        3.62839E-16, 6.26572E-02,-3.54952E+00, 1.48570E-16,
     &        2.58294E-01, 1.89742E-16, 1.75728E-01,-1.00770E-16,
     &       -1.29695E+01, 1.33766E-15,-2.43337E+00, 7.59815E-17,
     &        3.02448E-01, 4.15933E-20, 1.19717E+01, 5.13024E-16,
     &        4.85522E-01,-2.90220E-17,-4.29641E-02,-4.60473E-01,
     &       -3.14148E-16, 5.07662E-01, 1.48825E-16, 7.32770E-02,
     &        5.80716E+00, 5.03390E-16, 7.22441E-02, 3.81509E-17,
     &        9.80827E+00, 5.03265E-16, 6.28970E-01,-5.94293E-17,
     &       -6.45706E+00,-1.25291E-16, 3.96573E-02, 2.48721E+00,
     &        1.75540E-16,-1.39359E-01,-5.87244E+00,-4.58380E-16,
     &       -5.01084E+00,-3.08818E-16, 7.01727E-01,-3.73470E+00/
C     600km June solstice
      DATA (DERRTI(3,2,J),J=1,81)/ 1.29884E+02,
     &       -2.16283E+01, 9.11480E+01,-1.43655E+01, 8.05628E+00,
     &       -8.87552E+00,-6.34942E+01, 1.02339E+01,-3.26207E+01,
     &       -4.29686E+00, 5.85997E-02, 6.26688E+00,-4.85498E+00,
     &       -1.81359E+00, 5.08726E-01, 3.67758E+00, 1.88498E+00,
     &        7.37253E+00, 1.20254E+00, 8.32220E+00,-6.07879E-01,
     &        5.90250E+00,-9.50781E-01, 4.21618E+00,-7.32862E-02,
     &       -1.57466E+01, 1.90796E+00,-1.07067E+00,-1.68803E-01,
     &       -9.58362E-01,-6.87828E-02,-5.42266E-01,-3.99709E+00,
     &        4.73353E+00, 2.32839E-01, 2.42368E-01, 7.29817E-01,
     &        1.17218E-01, 7.68827E-01,-5.62172E+00,-2.81949E+00,
     &       -9.47799E-01,-6.15302E-01,-5.31118E-01,-3.52999E-01,
     &       -7.92378E+00,-1.79083E+00,-1.23697E+00, 4.57899E-01,
     &       -2.64944E-02, 2.23595E-01, 7.19864E+00,-1.03424E+00,
     &       -3.86267E-01,-4.72060E-01,-1.07357E-01,-7.59322E-01,
     &       -2.49914E+00,-3.42013E-01,-7.82734E-02,-1.58324E-01,
     &        9.80662E+00, 1.02683E+00, 2.28538E-01, 7.92079E-02,
     &        6.61293E+00, 3.70353E-01,-2.38804E-01,-4.93482E-02,
     &       -1.79546E+00, 1.10522E+00, 4.58724E-01, 3.42003E+00,
     &        1.22014E+00, 1.17760E-01,-6.10683E+00, 1.98447E-02,
     &       -2.02123E+00, 5.34722E-01,-3.21204E+00,-3.53669E+00/
C     850km equinox
      DATA (DERRTI(4,1,J),J=1,81)/ 1.79767E+02,
     &        4.99863E-14, 8.53012E+01,-4.27398E-14, 4.41491E+01,
     &       -9.78417E-14,-3.00528E+01, 1.27719E-13, 1.41935E+01,
     &       -9.37537E+00, 4.98015E-17, 1.92402E+01, 2.06054E-15,
     &       -6.42973E-01,-2.61885E-15,-2.65640E-01,-1.14981E-15,
     &        1.80054E+01, 9.96690E-15, 1.34061E+01, 2.61056E-15,
     &        1.02332E+00, 6.37023E-16,-3.95954E+00,-8.76311E-16,
     &       -3.52988E+01,-2.37946E-15, 8.54724E+00, 7.94609E-16,
     &        7.48440E-02, 1.24207E-15,-2.88888E+00,-2.00441E+01,
     &        2.41794E-16,-3.09010E+00,-1.02331E-15,-1.12850E+00,
     &       -4.44576E-16,-2.49107E-01,-5.58694E+00,-1.51204E-16,
     &       -6.44809E+00,-3.31683E-16,-1.11938E+00,-1.19668E-15,
     &       -2.38515E+01,-1.88703E-15, 5.88608E+00, 1.79372E-15,
     &        1.43047E+00, 4.94514E-16, 3.68967E+01,-2.10857E-16,
     &        2.72530E+00, 2.12777E-16, 1.77923E-01,-1.70715E+00,
     &        7.00875E-16,-5.71091E-01,-5.90993E-16,-1.01214E-01,
     &        7.25814E+00, 5.25188E-16,-8.80783E-01, 7.07762E-16,
     &        2.31478E+01,-1.10374E-15, 2.72138E+00,-2.78230E-16,
     &       -6.22856E+00,-1.78175E-15, 2.36877E+00,-2.84380E+00,
     &        8.36714E-16,-1.90972E+00,-2.39276E+00,-2.11964E-16,
     &        3.80104E+00,-3.43599E-16,-6.72237E+00,-5.99081E+00/
C     850km June solstice
      DATA (DERRTI(4,2,J),J=1,81)/ 1.82446E+02,
     &       -4.58919E+01, 5.20246E+01,-3.47508E+01, 2.02343E+01,
     &       -2.02712E+01,-2.87610E+01, 2.58023E+01, 3.92093E+01,
     &        5.02578E+00, 1.52096E+01, 2.52272E+01,-2.18689E+00,
     &       -9.23062E+00,-1.10061E+00,-3.17171E+00,-2.93270E+00,
     &        1.31091E+01,-8.80226E+00, 9.11917E+00, 4.86888E+00,
     &       -1.60386E+00,-1.92173E+00,-5.03796E+00, 5.53181E-02,
     &       -3.40267E+01, 1.77352E+00, 5.83784E+00,-2.78806E-01,
     &        5.37019E-01, 5.97907E-01,-2.06253E+00,-1.88425E+01,
     &        5.52838E+00,-2.26824E+00, 2.02123E+00,-1.38536E+00,
     &        5.58939E-02,-2.71219E-01,-7.54871E+00, 1.07884E-01,
     &       -4.43096E+00,-1.88325E-01,-1.03097E-01,-1.04328E-01,
     &       -1.69987E+01,-5.08520E+00, 4.76611E+00,-1.27428E+00,
     &        8.05233E-01,-3.67458E-02, 2.88803E+01,-1.34976E+00,
     &       -8.58966E-01,-3.37431E-01, 4.51486E-01,-6.18559E+00,
     &       -2.36339E-01,-8.59574E-01, 5.88599E-01, 2.56310E-02,
     &        1.17924E+01, 1.09272E+00,-8.47421E-01,-1.87615E-01,
     &        1.15823E+01,-1.53139E+00, 2.79530E-01,-2.08045E-02,
     &       -5.62334E+00,-6.12214E-01, 1.45218E+00, 4.89297E+00,
     &        1.54671E+00,-7.30855E-01,-6.65479E+00, 8.58245E-01,
     &        4.80930E+00, 1.13300E+00,-1.03658E+01, 4.79901E-01/
      DO 10 I=1,81
       DERRTI(1,3,I)=DERRTI(1,2,I)*MIRREQ(I)
       DERRTI(2,3,I)=DERRTI(2,2,I)*MIRREQ(I)
       DERRTI(3,3,I)=DERRTI(3,2,I)*MIRREQ(I)
10     DERRTI(4,3,I)=DERRTI(4,2,I)*MIRREQ(I)
      DO 40 K=1,81
       DO 30 J=1,3
        DO 20 I=1,4
         DOUT(I,J,K)=DERRTI(I,J,K)
20      CONTINUE
30     CONTINUE
40    CONTINUE
C////////////////////////////////////////////////////////////////////////////////////
      RETURN
      END
