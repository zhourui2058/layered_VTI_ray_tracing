function thomsen=C2thomsen(Cij)  %Tsvankin book pdf 13, equation 1.44
C11=Cij(1,1);
C13=Cij(1,3);
C33=Cij(3,3);
C55=Cij(5,5);
C66=Cij(6,6);

thomsen(1)=sqrt(C33);
thomsen(2)=sqrt(C55);
thomsen(3)=(C11-C33)/(2*C33);                           %epsilon
thomsen(4)=((C13+C55)^2-(C33-C55)^2)/(2*C33*(C33-C55)); %delta     
thomsen(5)=(C66-C55)/(2*C55);                           %gamma
end
