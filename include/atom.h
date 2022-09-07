#ifndef atom_h
#define atom_h
#include <iostream>
#include <list>
#include "polarconfig.h"
#include <string>
typedef struct Atom{
	double position[3];
	double charge[3];
	char type;
	int tick;
}atom;
void sort(double* input,int dim);
void unitcell(std::string file,int cell);
double* distance(atom *a,atom *b,double *p);
double norm(double* p,int dim);
/*no periodical boudary condition*/
double* distance(double*a ,double* b);
double far(atom* a,atom* b,double* p);
double far(double* a,double* b,double* p);
double* displace_average_A(atom *A,atom *oxygen,double* p,int cell);
double* displace_average_B(atom *B,atom *oxygen,double* p,int cell);
double* polar_average(atom *A,atom *B,atom *oxygen,double* p,int cell);
void displace_A_unit(atom* A,atom* oxygen,double* p,int cell);
void displace_B_unit(atom* B,atom* oxygen,double* p,int cell);
double displace_average_B_scalar(atom *B,atom *oxygen,double* p,int cell);
double* displace_average_Asite(atom* A,atom* oxygen,double* p,int cell,char type_id);
double* displace_average_Bsite(atom* B,atom* oxygen,double* p,int cell,char type_id);
double displace_average_Asite_scalar(atom* A,atom* oxygen,double *p,int cell,char type_id);
double displace_average_Bsite_scalar(atom* B,atom* oxygen,double* p,int cell,char type_id);
void sum_together(double*,double*,int);
int* changeindex(int index,int cell);
int  changeback(int x,int y,int z,int cell);
double average(std::list<double> &input);
double variance(std::list<double> &input);
int* neighbor_o_forA(int index,int cell);
int* neighbor_o_forB(int index,int cell);
int* neighbor_A_forB(int index,int cell);
double tiltangle(atom* a, atom* b,atom* c,double* p);
void analyzepolar(atom* A,atom* B,atom* oxygen, double* period, int cell);
void analyzeposition_variance(atom* A,atom* B,atom* oxygen,double* period,int cell,size_t signal);
void outpolar();
void calculate_local_die(int cell,double local_volume,double temperature);
double dielectric(double polarvar,double volume,double temp);
void calculate_local_variance(int cell,double temperature);
#endif
