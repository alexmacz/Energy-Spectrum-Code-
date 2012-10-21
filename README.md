subroutine spectra( N, XMIN, XMAX, XMIX, Zz)
   Implicit none	
         Real*8 L, Pi, r, No, nw, nh
         REAL*8 RANF
         EXTERNAL RANF
         Real*8 f_sum, maxf, minf, v, pa
         Real*8 f_avg, f_rms, frms_good
         Real*8 fmean_good, frms_ratio
	 REAL*8 XMIN, XMAX, argg, XMIX
	 Integer m, b, N, rank
         Real*8, Dimension (:), ALLOCATABLE :: k, xn
	 Real*8 Zz(*)
	 Complex, Dimension (:), ALLOCATABLE ::E ,P, S 
         Real*8 Ke, Kd
	 Complex i_complex
	 OPEN(41, FILE='mfield.dat')
  	   L = XMAX - XMIN
	   nw = N - 1
	   nh = (nw/2) + 1
	   rank = N
	   WRITE (*, *) L, N, XMIN, XMAX
         Allocate ( K(rank), xn(rank))
         Allocate ( E(rank), P(rank), S(rank))
	 Pi = 4*ATAN(1.0)
         Ke = 2 * PI / XMIX  
      	 Kd = 600 
	   i_complex = (0, 1)
!	   write (*,*) Ke, Kd, i_complex
	 Do m = 1, N
            k(m) = 2*Pi*( m-1 )/L
            v = (((k(m)/Ke)**4)/(1.0+(k(m)/Ke)**2)**(17/6))
	    pa = exp(-2.25*((k(m)/Kd)**(4/3)))
	    E(m) = v * pa
           write (41,*) k(m), real(E(m))	
         End Do
         Do m = 1, nh
      	    r = RANF () 
	    argg = 2 * Pi * ( r - 0.5 ) 
    	    P(m) = sqrt(E(m))* EXP(i_complex * argg)
!	    write (41,*) real(P(m))	
         End Do
C        symmetry ....
	   Do m = nh+1, N
            b = nw + 2 - m
	      P(m) = conjg(P(b))
!	      write (41,*) real(P(m))
	   End Do
!        compute f(x) 			   
	   Do b = 1, N
  	      S(b) = 0
	      xn(b) = (b-1)*L/(N-1)
	   Do m = 1, N
	      k(m) = (m-1)*2*Pi/L
		S(b) = S(b) + real(P(m)*exp(i_complex*k(m)*xn(b)))
	   End Do
	      S(b) = real(S(b)) / N
!		write (41,*) xn(b), real(S(b))
	   End Do
         f_sum = 0  
         Do b = 1, N
            f_sum = f_sum+real(S(b))
         End do			  
!	   write(*,*) f_sum
         f_avg = f_sum/N
         frms_good = 0.5
         fmean_good = 0.16 
         Do b = 1, N
            f_rms = f_rms + (S(b)*S(b))
         End Do
	   f_rms = sqrt(f_rms/N )
	   frms_ratio = frms_good / f_rms
         Do b = 1, N
		Zz(b) = fmean_good + frms_ratio * S(b)
!		fscale(n) = MIN (0.99999,fscale(n))
!               fscale(n) = MAX (0.00001,fscale(n))
               IF (Zz(b).GT.1.) THEN
	          Zz(b) = 1
	       Else IF (Zz(b).LT.0.) THEN
	               Zz(b) = 0
	       End IF
!	       WRITE (*, *) Zz(b)
!		write (41,*) Zz(b)
	 End Do
	 DEALLOCATE ( k, xn, E, P, S)
	 RETURN