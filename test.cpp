#include"KineCal.h"
#include"TaggedLambda_DIS.h"
#include"piIMParton.h"

#include"KineCal.cpp"
#include"TaggedLambda_DIS.cpp"
#include"piIMParton.cpp"


void test(){

	TaggedLambda_DIS dis;

	dis.SetQ2max(60);
	dis.SetQ2min(1);
	dis.SetTmax(2);
	dis.SetTmin(0.001);
	dis.SetxLmax(0.995);
	dis.SetxLmin(0.2);

	char filename[50] = "TaggedLambda-DIS-EicC.root";
	dis.SetOutputFileName(filename); 

	dis.SetElecBeamEnergy(3.5);
	dis.SetProtBeamEnergy(20);
	//dis.SetBeamCrossAngle(0.05);  //// 50 mrad
	dis.SetBeamCrossAngle(0);   //// 0 mrad

	dis.Generate(20000);

}


