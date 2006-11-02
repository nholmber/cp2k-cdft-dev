      SUBROUTINE LSFBTR(F,G,L,RHOMIN,KAPMIN,H,N,SAVED,XA,TA,            
     1  TRBE,WW,NDIM,DISC)                                              
C                                                                       
C  THIS PROGRAM CALCULATES THE SPHERICAL BESSEL TRANSFORM OF ORDER      
C  L FOR A FUNCTION WHOSE VALUES ARE GIVEN IN THE ARRAY F FOR 2**N      
C  VALUES OF R.  THE VALUES OF R ARE GIVEN BY R(I)=EXP(RHOMIN+(I-1)*H)  
C  WHERE I=1,2,...,NT AND NT=2**N.  THE VALUES OF THE TRANSFORM         
C  ARE RETURNED IN THE ARRAY G FOR NT K VALUES.  THE VALUES OF K        
C  ARE EXP(KAPMIN+(I-1)*H), I=1,2,...,NT.                               
C                                                                       
C  TWO DIFFERENT SETS OF RESULTS ARE CALCULATED.  ONE SET IS ACCU-      
C  RATE FOR LARGE VALUES OF K AND THE OTHER IS ACCURATE FOR K CLOSE     
C  TO ZERO.  THE TWO SETS ARE JOINED AT THE K VALUE FOR WHICH THE       
C  DIFFERENCE BETWEEN THEM IS LEAST.  FOR L=0 OR L=1 THE SMALL K        
C  RESULTS ARE CALCULATED USING SIEGMANS METHOD.                        
C                                                                       
C  IF THE LOGICAL VARIABLE SAVED IS .FALSE. AUXILIARY DATA IN THE       
C  ARRAYS TA, TRBE, AND WW ARE CALCULATED ANEW.  THIS MUST BE DONE      
C  IF THE MESH PARAMETERS N, H, RHOMIN OR KAPMIN HAVE CHANGED FROM      
C  THE PREVIOUS CALL TO THE SUBROUTINE AND ON THE FIRST CALL TO         
C  THE SUBROUTINE.  OTHERWISE A SUBSTANTIAL IMPROVEMENT IN SPEED        
C  CAN BE OBTAINED BY TAKING SAVED=.TRUE..                              
C                                                                       
C  THE VARIABLE DISC IS THE DIFFERENCE BETWEEN THE SMALL K AND          
C  LARGE K VALUES FOR THE RESULT AT THE MATCHING POINT OF MINIMUM       
C  DISCREPANCY.  IT PROVIDES A ROUGH GUIDE TO THE ABSOLUTE ACC-         
C  URACY OF THE RESULTS AT THE MAXIMUM VALUE OF THE TRANSFORM.          
C  IT IS, HOWEVER, A GREAT OVERESTIMATE OF THE ABSOLUTE ERROR AT        
C  LARGE K VALUES, AND SMALL K VALUES IN THE CASE OF NON-ZERO L.        
C                                                                       
C  XA, TA, TRBE AND WW ARE COMPLEX WORKING AREA ARRAYS DIMENSIONED      
C  AS XA(NDIM), TA(NDIM), TRBE(NDIM,2) AND WW(NDIM) WHERE NDIM          
C  MUST BE AT LEAST AS LARGE AS 2**N.  AS NOTED, IT MAY BE USEFUL       
C  TO PRESERVE THE CONTENTS OF TA, TRBE AND WW BUT THE CONTENTS OF      
C  XA ARE NOT OF VALUE.                                                 
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XA(NDIM),TA(NDIM),TRBE(NDIM,2),WW(NDIM),                  
     1  Y1,Y2,W1,W2                                                     
      DIMENSION F(NDIM),G(NDIM)                                         
      REAL*8 KAPMIN                                                       
      LOGICAL SAVED                                                     
      NT=2**N                                                           
      IF (NT.GT.NDIM) GO TO 131                                         
      N2=2*NT                                                           
      NH=NT/2                                                           
      NHP=NH+1                                                          
      PI=3.14159265358979323846264D0
      AN=NT                                                             
      DT=2.0*PI/(H*AN)                                                  
      CM=SQRT(2.0*PI)/AN                                                
      CL=0.25*H/AN                                                      
      IF (.NOT.SAVED) GO TO 100                                         
C                                                                       
C  CALCULATE RESULTS ACCURATE AT SMALL K VALUES                         
C                                                                       
    1 IF (L.LT.2) GO TO 11                                              
      AA=EXP((L+1.5)*RHOMIN)                                            
      BB=EXP((L+1.5)*H)                                                 
      DO 2 I=1,NT                                                       
      XA(I)=F(I)*AA                                                     
    2 AA=AA*BB                                                          
      CALL NLOGN(N,XA,WW,NDIM)                                          
      DO 3 I=1,NH                                                       
    3 XA(I)=TA(I)*XA(I)                                                 
      DO 4 JJ=1,L                                                       
      AA=2*JJ-0.5                                                       
      T=0.0                                                             
      DO 5 I=1,NH                                                       
      XA(I)=XA(I)/CMPLX(AA,T)                                           
    5 T=T+DT                                                            
    4 CONTINUE                                                          
      DO 6 I=NHP,NT                                                     
    6 XA(I)=0.0                                                         
      CALL NLOGN(N,XA,WW,NDIM)                                          
      AA=EXP((L-1.5)*KAPMIN)*CM                                         
      BB=EXP((L-1.5)*H)                                                 
      DO 7 I=1,NT                                                       
      G(I)=AA*XA(I)                                                     
    7 AA=AA*BB                                                          
      GO TO 21                                                          
C                                                                       
C  FOR L=0 OR 1 CALCULATE RESULTS ACCURATE AT SMALL K VALUES USING      
C  SIEGMANS METHOD                                                      
C                                                                       
   11 LP=L+1                                                            
      DO 12 I=NH,NT                                                     
   12 XA(I)=0.0                                                         
      AA=EXP(3.0*RHOMIN)                                                
      BB=EXP(3.0*H)                                                     
      DO 13 I=1,NH                                                      
      XX=AA*F(2*I-1)                                                    
      AA=AA*BB                                                          
      YY=AA*F(2*I)                                                      
      XA(I)=CMPLX(XX,YY)                                                
   13 AA=AA*BB                                                          
      CALL NLOGN(N,XA,WW,NDIM)                                          
      Y1=TRBE(1,LP)*XA(1)                                               
      Y2=TRBE(1,LP)*CONJG(XA(1))                                        
      XA(1)=2.0*(Y1+Y2+CONJG(Y2-Y1))                                    
      XA(NHP)=4.0*CONJG(XA(NHP)*TRBE(NHP,LP))                           
      DO 14 I=2,NH                                                      
      IC=NT-I+2                                                         
      Y1=XA(I)                                                          
      Y2=CONJG(XA(IC))                                                  
      W1=Y1+Y2                                                          
      W2=WW(I+NH)*(Y1-Y2)                                               
      Y1=(W1-W2)*TRBE(I,LP)                                             
      Y2=(W1+W2)*CONJG(TRBE(IC,LP))                                     
      W1=Y1+Y2                                                          
      W2=WW(I+NH)*(Y1-Y2)                                               
      XA(I)=W1+W2                                                       
   14 XA(IC)=CONJG(W1-W2)                                               
      CALL NLOGN(N,XA,WW,NDIM)                                          
      DO 15 I=1,NH                                                      
      G(2*I-1)=CL*REAL(XA(I))                                           
   15 G(2*I)=CL*AIMAG(XA(I))                                            
C                                                                       
C  CALCULATE RESULTS ACCURATE AT LARGE K VALUES                         
C                                                                       
   21 AA=EXP(1.5*RHOMIN)                                                
      BB=EXP(1.5*H)                                                     
      DO 22 I=1,NT                                                      
      XA(I)=AA*F(I)                                                     
   22 AA=AA*BB                                                          
      CALL NLOGN(N,XA,WW,NDIM)                                          
      IJ=MOD(L,2)                                                       
      IJK=MOD(L,4)                                                      
      DO 23 I=1,NH                                                      
      Y1=XA(I)*TA(I+IJ*NH)                                              
      IF (IJK.GT.1) Y1=-Y1                                              
   23 XA(I)=Y1                                                          
      IF (L.EQ.0) GO TO 24                                              
      DO 25 JJ=1,L                                                      
      AA=2*JJ-L-0.5                                                     
      BB=JJ-0.5                                                         
      T=0.0                                                             
      DO 26 I=1,NH                                                      
      XA(I)=XA(I)*CMPLX(BB,-T)/CMPLX(AA,T)                              
   26 T=T+DT                                                            
   25 CONTINUE                                                          
   24 DO 27 I=NHP,NT                                                    
   27 XA(I)=0.0                                                         
      CALL NLOGN(N,XA,WW,NDIM)                                          
      AA=EXP(-1.5*KAPMIN)*CM                                            
      BB=EXP(-1.5*H)                                                    
      DO 28 I=1,NT                                                      
      XA(I)=AA*XA(I)                                                    
   28 AA=AA*BB                                                          
C                                                                       
C  FIND POINT OF MINIMUM DISCREPANCY BETWEEN SMALL K AND LARGE K VALUES 
C  OF TRANSFORM AND CONSTRUCT BEST APPROXIMATE RESULT.                  
C                                                                       
      D1=ABS(G(1)-REAL(XA(1)))                                          
      D2=ABS(G(2)-REAL(XA(2)))                                          
      TE=D1+D2                                                          
      II=2                                                              
      D1=D2                                                             
      DO 31 I=3,NT                                                      
      D2=ABS(G(I)-REAL(XA(I)))                                          
      AA=D1+D2                                                          
      IF (AA.GE.TE) GO TO 32                                            
      II=I                                                              
      TE=AA                                                             
   32 D1=D2                                                             
   31 CONTINUE                                                          
      DO 33 I=II,NT                                                     
   33 G(I)=XA(I)                                                        
      DISC=TE                                                           
      RETURN                                                            
C                                                                       
C  INITIALIZE AUXILIARY DATA FOR TRANSFORM                              
C                                                                       
C  THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY TA USED IN CALCULATING    
C  THE TRANSFORM AT LARGE K VALUES AND SMALL K VALUES FOR L GREATER     
C  THAN 1.                                                              
C                                                                       
  100 CONTINUE                                                          
      Y1=1.0                                                            
      AA=(RHOMIN+KAPMIN)*DT                                             
      Y2=CMPLX(COS(AA),SIN(AA))                                         
      DO 102 I=1,NH                                                     
      T=(I-1)*DT                                                        
      S=0.0                                                             
      RR=SQRT(110.25+T*T)                                               
      PHI=ATAN(T/10.5)                                                  
      DO 103 NA=1,10                                                    
  103 S=S+ATAN(T/(NA-0.5))                                              
      S=S-T*LOG(RR)+T-10.0*PHI+SIN(PHI)/(12.0*RR)                      
      S=S-SIN(3.0*PHI)/(360.0*RR**3)+SIN(5.0*PHI)/(1260.0*RR**5)        
     1   -SIN(7.0*PHI)/(1680.0*RR**7)                                   
      XX=EXP(PI*T)                                                      
      XX=ATAN((XX-1.0)/(XX+1.0))                                        
      CC=S-XX                                                           
      TA(I)=Y1*CMPLX(COS(CC),SIN(CC))                                   
      CC=S+XX                                                           
      TA(I+NH)=Y1*CMPLX(COS(CC),SIN(CC))                                
  102 Y1=Y1*Y2                                                          
      TA(1)=TA(1)/2.0                                                   
      TA(1+NH)=TA(1+NH)/2.0                                             
      TA(NH)=TA(NH)/2.0                                                 
      TA(NT)=TA(NT)/2.0                                                 
C                                                                       
C  THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY WW USED BY                
C  THE NLOGN SUBROUTINE.  THE ELEMENTS IN THE SECOND HALF OF THE        
C  ARRAY WW ARE USED IN THE IMPLEMENTATION OF SIEGMANS METHOD.          
C                                                                       
      DO 111 I=1,NH                                                     
      XX=(I-1)*PI/AN                                                    
      WW(I+NH)=CMPLX(-SIN(XX),COS(XX))                                  
      XX=2.0*XX                                                         
  111 WW(I)=CMPLX(COS(XX),SIN(XX))                                      
C                                                                       
C  THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY TRBE USED IN THE          
C  IMPLEMENTATION OF SIEGMANS METHOD.                                   
C                                                                       
      DO 112 I=1,NT                                                     
  112 XA(I)=0.0                                                         
      CF=EXP(H)                                                         
      XX=EXP(RHOMIN+KAPMIN)                                             
      DO 113 I=1,NT                                                     
      AA=SIN(XX)/XX                                                     
      XX=CF*XX                                                          
      BB=SIN(XX)/XX                                                     
      XA(I)=CMPLX(AA,BB)                                                
      IF (XX.GT.1.0E8) GO TO 114                                        
  113 XX=XX*CF                                                          
  114 CALL NLOGN(N,XA,WW,NDIM)                                          
      TRBE(1,1)=XA(1)                                                   
      TRBE(NHP,1)=CONJG(XA(NHP))                                        
      DO 115 I=2,NH                                                     
      IC=NT-I+2                                                         
      Y1=XA(I)                                                          
      Y2=CONJG(XA(IC))                                                  
      W1=Y1+Y2                                                          
      W2=WW(NH+I)*(Y1-Y2)                                               
      Y1=W1-W2                                                          
      Y2=CONJG(W1+W2)                                                   
      TRBE(I,1)=0.5*CONJG(Y1)                                           
  115 TRBE(IC,1)=0.5*CONJG(Y2)                                          
      DO 120 I=1,NT                                                     
  120 XA(I)=0.0                                                         
      XX=EXP(RHOMIN+KAPMIN)                                             
      DO 121 I=1,NT                                                     
      IF (XX.LT.0.1) GO TO 122                                          
      AA=(SIN(XX)/XX-COS(XX))/XX                                        
      XX=CF*XX                                                          
      BB=(SIN(XX)/XX-COS(XX))/XX                                        
      XA(I)=CMPLX(AA,BB)                                                
      IF (XX.GT.1.0E8) GO TO 123                                        
      GO TO 121                                                         
  122 CC=XX*XX/2.0                                                      
      CD=1.0-CC/5.0+CC*CC/70.0-CC*CC*CC/1890.0+CC**4/83160.0            
      AA=XX*CD/3.0                                                      
      XX=CF*XX                                                          
      CC=XX*XX/2.0                                                      
      CD=1.0-CC/5.0+CC*CC/70.0-CC*CC*CC/1890.0+CC**4/83160.0            
      BB=XX*CD/3.0                                                      
      XA(I)=CMPLX(AA,BB)                                                
  121 XX=XX*CF                                                          
  123 CALL NLOGN(N,XA,WW,NDIM)                                          
      TRBE(1,2)=XA(1)                                                   
      TRBE(NHP,2)=CONJG(XA(NHP))                                        
      DO 124 I=2,NH                                                     
      IC=NT-I+2                                                         
      Y1=XA(I)                                                          
      Y2=CONJG(XA(IC))                                                  
      W1=Y1+Y2                                                          
      W2=WW(NH+I)*(Y1-Y2)                                               
      Y1=W1-W2                                                          
      Y2=CONJG(W1+W2)                                                   
      TRBE(I,2)=0.5*CONJG(Y1)                                           
  124 TRBE(IC,2)=0.5*CONJG(Y2)                                          
      GO TO 1                                                           
  131 WRITE (6,132)                                                     
  132 FORMAT(56H SUBROUTINE LSFBTR NOT EXECUTED-DIMENSION NDIM TOO SMALL
     1)                                                                 
      RETURN                                                            
      END                                                               
      SUBROUTINE NLOGN(N,X,WW,NDIM)                                     
C                                                                       
C   FAST FOURIER TRANSFORM ROUTINE                                      
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NDIM),MM(15)                                          
      COMPLEX*16 WW(1),X,WK,HOLD,Q                                         
      DO 1 I=1,N                                                        
    1 MM(I)=2**(N-I)                                                    
      LX=2*MM(1)                                                        
      DO 4 L=1,N                                                        
      NBLOCK=2**(L-1)                                                   
      LBLOCK=LX/NBLOCK                                                  
      LBHALF=LBLOCK/2                                                   
      K=0                                                               
      DO 4 IBLOCK=1,NBLOCK                                              
      WK=WW(K+1)                                                        
      ISTART=LBLOCK*(IBLOCK-1)                                          
      DO 2 I=1,LBHALF                                                   
      J=ISTART+I                                                        
      JH=J+LBHALF                                                       
      Q=X(JH)*WK                                                        
      X(JH)=X(J)-Q                                                      
      X(J)=X(J)+Q                                                       
   2  CONTINUE                                                          
      DO 3 I=2,N                                                        
      II=I                                                              
      IF (K.LT.MM(I)) GO TO 4                                           
   3  K=K-MM(I)                                                         
   4  K=K+MM(II)                                                        
      K=0                                                               
      DO 7 J=1,LX                                                       
      IF (K.LT.J) GO TO 5                                               
      HOLD= X(J)                                                        
      X(J)=X(K+1)                                                       
      X(K+1)= HOLD                                                      
   5  DO 6 I=1,N                                                        
      II=I                                                              
      IF (K.LT.MM(I)) GO TO 7                                           
   6  K=K-MM(I)                                                         
   7  K=K+MM(II)                                                        
      RETURN                                                            
      END                                                               
                                                                                
