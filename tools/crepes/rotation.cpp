//app for rotation about a line
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <cstring>

using namespace std;

#define pi 3.14159265

//define variables
	double rot1[3], rot2[3];
	double ux,uy,uz;
	double uxn,uxd;
	string line;
	string first;
	string second;
	string strxrot1, stryrot1, strzrot1;
	string strxrot2, stryrot2, strzrot2;
	int deg;
	int atoms=0;
	double transmat[3][3] = {{ 1,1,1},{1,1,1},{1,1,1}};
	int rotid=0;
	
// args are #of atoms to rotate
// degree of rotation
// input file
double **createcoordarrays (int);
string *createnamearray (int);
void assigntransmat(double* , double*, int );
void readfile();
void printdata(int, string*, double* , double*,double**);




int main (int argc, char **argv){
//	cout << "program initialized \n";		
	atoms = atoi(argv[1]);
//	cout << "input atoms : " << atoms << endl; 

	//maintain precision
	cout.precision(8);
	cout << fixed;

//	cout << "SET PRECISION\n";

	// setup some variables	
	double **coord = createcoordarrays(atoms);
//	cout << "COORD ARRAY CREATED\n";
	
	double **transf = createcoordarrays(atoms);
//	cout << "TRANSF ARRAY CREATED\n";

	int allatoms = atoms + 2;
	
//	cout << "allatoms: " << allatoms << endl;
	
	string *name = createnamearray(allatoms);
//	cout << "NAME " << endl;


	
	// open file
//	cout << "OPENING FILE" << endl;
	ifstream coordsfile;
	coordsfile.open(argv[3]);
	//error if file does not exist
	if (!coordsfile){
		cout << " File does not exist\n";
		exit (1);
	}

//	cout << "FILE IS OPEN " << endl;
	if (coordsfile.is_open()){


	int counter = 0;
	while (getline(coordsfile,line)){
		counter++;
//		cout << "\n\t LINE READING NUMBER " << counter << endl;
		size_t i = line.find("\t");
		 
//		cout << "\t\t searched for tab found at i: " << i << endl; 
		 
		if (i != string::npos){
			size_t y=0;
			if (!line.empty()){
				string first="";
				string second="";
				while (y!=i){
					first += line[y++];
				}
				//store name
				name[rotid] = first;
//				cout << "\n\t\t " << first << endl;
				y++;
				
				while (y!=line.length()){
					second += line[y++];
				}
//				cout << "\n\t\t " << second << endl;
				//Parse numeric half
				size_t r=second.find(",");
				size_t s=second.rfind(",");
//				cout << "\t\t r : " << r << endl;
//				cout << "\t\t s : " << s << endl;
				
				int z=0;
				
//				cout << "\n\t\tbeforeif " << endl;
				if (!second.empty()){
					strxrot1="";
					stryrot1="";
					strzrot1="";
					
//					cout << "\t\t\t z < r "<< endl;
					
					while (z<r){
						if (rotid==0) {
							strxrot1 += second[z++];
							rot1[0] = atof(strxrot1.c_str());
//							cout << "x1";
						}
						if (rotid==1) {
							strxrot2 += second[z++];
							rot2[0] =atof(strxrot2.c_str());
//							cout << "x2";
						}
						if (rotid>=2){
//							cout << "x";
							strxrot1 += second[z++];
							int coordid=rotid-2;
							coord[0][coordid]=atof(strxrot1.c_str());
//							cout << "x3";
						}
					}
//					cout << endl;
					z++;
					
//					cout << "\t\t\t z < r && z > r"<< endl;
					while (z<s and z>r){
					
						if (rotid==0) {
							stryrot1 += second[z++];
							rot1[1] =atof (stryrot1.c_str());
//							cout << "y1";
						}
						if (rotid==1) {
							stryrot2 += second[z++];
							rot2[1] =atof(stryrot2.c_str());
//							cout << "y2";
						}
						if (rotid>=2){
							stryrot1 += second[z++];
							int coordid=rotid-2;
							coord[1][coordid]=atof(stryrot1.c_str());
//							cout << "y3";
						}
					}
					z++;
//					cout << endl;
//					cout << "\t\t\t z != second.length "<< endl;
					while (z!=second.length()){

//						cout << "\n\t\t\t if 1 " << endl;					
						if (rotid==0) {
							strzrot1 += second[z++];
							rot1[2] =atof (strzrot1.c_str());
//							cout << "z1";
						}
//						cout << "\n\t\t\t if 2 " << endl;					
						if (rotid==1) {
							strzrot2 += second[z++];
							rot2[2] =atof (strzrot2.c_str());
							
						}
//						cout << "\n\t\t\t if 3 " << endl;					
						if (rotid>=2){
							strzrot1 += second[z++];
							int coordid=rotid-2;
							coord[2][coordid]=atof(strzrot1.c_str());
//							cout << "z3" << endl;
						}
//						cout << "\n\t\t\t end of three ifs " << endl;
					}
				}
//				cout << endl << "before : " <<	 rotid;
				rotid++;
//				cout << " after : " << rotid << endl;
			}
		}
		else {
			string first=line;
			string second="";			
		}		
//		cout << "\n\t LINE READ NUMBER " << counter << endl << endl;
	}
}


//	cout << "CLOSING FILE\n";

	coordsfile.close();
		
 	
 	deg = atoi(argv[2]);
// 	cout << "SET the ANGLE OF ROTATION " << deg << "\n";
	assigntransmat(rot1, rot2, deg);	


//translate coords
	for(int i=0; i < 3; i++){
		for (int j=0; j < atoms; j++){
			coord[i][j] = coord[i][j] - rot1[i];
		}
	}


		
	// do matrix transformation (rotation) about line
	
//	cout << "PERFORMING MATRIX TRANSFORMATION\n";
	for (int i=0; i < 3; i++){
		for (int j=0; j < atoms; j++){
			for (int k=0; k < 3 ; k++){
				transf[i][j] += transmat[i][k]*coord[k][j];
			}
		}			
	}
	
//translate coords back
	for(int i=0; i <3; i++){
		for (int j=0; j < atoms; j++){
			transf[i][j] = transf[i][j] + rot1[i];
		}
	}

//	printdata(allatoms ,name, rot1, rot2, transf);				
 	for (int i=0; i<allatoms; i++){
		if (i==0){
			cout << name[i];
			for (int r1=0; r1<3; r1++){
				printf("  %*.*f",15,8,rot1[r1] );
			}
			cout << endl;
		}
		if (i==1){
			cout << name[i];
			for (int r2=0; r2<3; r2++){
				printf("  %*.*f",15,8,rot2[r2] );
			}
			cout << endl;
		}
		if (i>1){
			cout << name[i];
			for (int r3=0; r3<3; r3++){
				printf("  %*.*f",15,8,transf[r3][i-2] );
			}
			cout << endl;
		}
	}
	/*
	//cleanup

	delete[] name;
	delete[] coord;
	delete[] transf;
	
	//cout << "\nDone\n";
	*/
	return 0;
}

double **createcoordarrays (int atoms){

	//cout << endl << "-------- start createcoordarrays (atoms: " <<  atoms << ") " << endl;

	double **j;
	j = new double*[3];
	for (int i=0;i<3;i++){
		j[i]=new double[atoms];
	}
	
	//cout << endl << "-------- end createcoordarrays " << endl;
	return j;
}

string *createnamearray (int allatoms){
	// create name array
	string *name;
	name = new string[allatoms];
	for (int i=0;i<allatoms;i++){
		name[i]="Some Literal string  ";
	}
	return name;
}


void assigntransmat(double rot1[], double rot2[],int deg){
	double u=(rot1[0]-rot2[0]);
	double v=(rot1[1]-rot2[1]);
	double w=(rot1[2]-rot2[2]);
	
	double l=sqrt(u*u+v*v+w*w);
	
	transmat[0][0]=(u*u+(v*v+w*w)*(cos(deg*pi/180)))/(l*l);
	transmat[0][1]=(u*v*(1-cos(deg*pi/180))-w*l*sin(deg*pi/180))/(l*l);
	transmat[0][2]=(u*w*(1-cos(deg*pi/180))+v*l*sin(deg*pi/180))/(l*l);
	transmat[1][0]=(u*v*(1-cos(deg*pi/180))+w*l*sin(deg*pi/180))/(l*l);
	transmat[1][1]=(v*v+(u*u+w*w)*(cos(deg*pi/180)))/(l*l);
	transmat[1][2]=(v*w*(1-cos(deg*pi/180))-u*l*sin(deg*pi/180))/(l*l);
	transmat[2][0]=(u*w*(1-cos(deg*pi/180))-v*l*sin(deg*pi/180))/(l*l);
	transmat[2][1]=(v*w*(1-cos(deg*pi/180))+u*l*sin(deg*pi/180))/(l*l);
	transmat[2][2]=(w*w+(u*u+v*v)*(cos(deg*pi/180)))/(l*l);
}

/*
void readfile(){
// read in the file
if (coordsfile.is_open()){
	while (getline(coordsfile,line)){
		size_t i = line.find("\t");
		if (i != string::npos){
			size_t y=0;
			if (!line.empty()){
				string first="";
				string second="";
				while (y!=i){
					first += line[y++];
				}
				//store name
				name[rotid] = first;
				
				y++;
				
				while (y!=line.length()){
					second += line[y++];
				}
				//Parse numeric half
				size_t r=second.find(",");
				size_t s=second.rfind(",");
				
				if  (i!=string::npos){
					size_t y=0;
					if (!second.empty()){
						strxrot1="";
						stryrot1="";
						strzrot1="";
						while (y<r){
							if (rotid==0) {
								strxrot1 += second[y++];
								rot1[0] = atof(strxrot1.c_str());
							}
							if (rotid==1) {
								strxrot2 += second[y++];
								rot2[0] =atof(strxrot2.c_str());
							}
							if (rotid>=2){
								strxrot1 += second[y++];
								int coordid=rotid-2;
								coord[0][coordid]=atof(strxrot1.c_str());
							}
						}
						
						y++;
						while (y<s){
							if (rotid==0) {
								stryrot1 += second[y++];
								rot1[1] =atof (stryrot1.c_str());
							}
							if (rotid==1) {
								stryrot2 += second[y++];
								rot2[1] =atof(stryrot2.c_str());
							}
							if (rotid>=2){
								stryrot1 += second[y++];
								int coordid=rotid-2;
								coord[1][coordid]=atof(stryrot1.c_str());
							}
						}
						y++;
						while (y!=second.length()){
							if (rotid==0) {
								strzrot1 += second[y++];
								rot1[2] =atof (strzrot1.c_str());
							}
							if (rotid==1) {
								strzrot2 += second[y++];
								rot2[2] =atof (strzrot2.c_str());
							}
							if (rotid>=2){
								strzrot1 += second[y++];
								int coordid=rotid-2;
								coord[2][coordid]=atof(strzrot1.c_str());									
							}
							
						}
					}
				}
				rotid++;
				
			}
		}
		else {
		string first=line;
		string second="";
		}		
	}
}
coordsfile.close();

}
*/

/*void printdata(int allatoms ,string* name[], double rot1[], double rot2[], double** transf[][]){
	for (int i=0; i<allatoms; i++){
		if (i==0){
			cout << name[i];
			for (int r1=0; r1<3; r1++){
				printf("  %*.*f",15,8,rot1[r1] );
			}
			cout << endl;
		}
		if (i==1){
			cout << name[i];
			for (int r2=0; r2<3; r2++){
				printf("  %*.*f",15,8,rot2[r2] );
			}
			cout << endl;
		}
		if (i>1){
			cout << name[i];
			for (int r3=0; r3<3; r3++){
				printf("  %*.*f",15,8,transf[r3][i-2] );
			}
			cout << endl;
		}
	}
}
*/
