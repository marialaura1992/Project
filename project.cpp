#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <iomanip> //to use setw(5) i.e. set width to 5
#include <cstdlib>
#include <ctime>
using namespace std;
using Matrix = vector<vector<int>>;
int up = 1;
int down = -1;
Matrix* p1;
Matrix* p2;
Matrix* p3;
Matrix* p5;
Matrix* p6;
Matrix* p7;
Matrix* p8;
Matrix* p9;
Matrix* p10;
vector<Matrix*>* p4;
double* a;
double* b;
double damage= .8;
vector<Matrix>* possible_inputs;
Matrix* output;
vector<Matrix>* second_possible_inputs;
Matrix* second_output;

//defining a nicer print function (all the -1 are neglected from the grid)

void print_true (const Matrix & m ) {
  
  for (int i=0; i< 50; i++) {
   cout << '-' ;
  }
  cout<< '\n'<< '\n';
  
  for (unsigned row =0; row < m.size(); row++) {
    for (unsigned col=0; col< m[row].size(); col++) {
      if (m[row][col]== -1){
	cout << setw(5) << " ";
      }
      else {
	cout << setw(5) << m[row][col];
      }
    }
    cout << '\n';

  }
  
  cout<< '\n';
  for (int i =0; i< 50; i++) {
    cout << '-' ;
  }
  cout<< '\n';
  
}

//Defining a normal print function

void print (const Matrix & m ) {
  for (unsigned row =0; row < m.size(); row++) {
    for (unsigned col=0; col< m[row].size(); col++)
      cout << setw(5) << m[row][col];
    cout << '\n';
  }
}


//Defining a function spin(m,n,matrix) that gives the value of the spin matrix[m,n]

int spin_value (int m, int n, const Matrix & grid) {
  //cout<< grid[m][n];
  return grid[m][n];
}

    
//Defining the spin identifier: given the row and column for a spin calculate a unique integer number between 1 and 100/
   
int I(int m, int n) {
  int res =  10 * (m-1) + n;
  cout<< res<< '\n';
  return res;
}

    
//Given an integer between 1 and 100 find the unique value of m and n 

double res2;
vector<int> res_vec (2);
vector<int>  find_m_n(int i ) {
 double res1= i /10.0;
 double res3= modf(res1,&res2); //take the modulus of res1 and put in res2 which is still a float
 int intpart= (int) res2;
  if (res3==0) {
    res_vec[0] = intpart-1;
    res_vec[1]= 9;
  }
  else {
    res_vec[0]= intpart;
    res_vec[1]= (i-1) - (10 * intpart);
  }
  //cout << "m is equal to " << res_vec[0]<< " and n is equal to "<< res_vec[1]<< '\n';
  return res_vec;
}


//Defining the Hamming distance: gives how well defined and far from each other are the energy minima associated to each pattern

double hamm_dist(const Matrix &m, const Matrix &n) {
  double result=0, tot;
  for (int i = 1; i<= 100; i++ ) {
    int s_m= spin_value(find_m_n(i)[0], find_m_n(i)[1], m);
    int s_n= spin_value(find_m_n(i)[0], find_m_n(i)[1], n);
    result += (s_m-s_n)*(s_m-s_n);
  }
  tot= 1/100. * result;
  cout<< tot << '\n';
  return tot;
}

//Checking the orthogonality between two patterns
double ortho(const Matrix&m, const Matrix&n) {
  double result=0;
  for (int i =1; i<=100; i++){
    int s_m= spin_value(find_m_n(i)[0], find_m_n(i)[1], m);
    int s_n= spin_value(find_m_n(i)[0], find_m_n(i)[1], n);
    result += s_m * s_n;
  }
  cout<< result << '\n';
  return result;
}

//General definition of the interaction energies (for a one  pattern)

double J_ij(int i1, int i2, const Matrix &m) {
  double value;
  int s1= spin_value(find_m_n(i1)[0], find_m_n(i1)[1], m);
  int s2= spin_value(find_m_n(i2)[0], find_m_n(i2)[1], m);
  value= s1*s2;
  return value;
}

//Defining the interaction energies for a set of patterns

double J(int i1, int i2, const vector<Matrix*> &pattern) {
  double value=0, total;
  for  (unsigned j = 0; j < pattern.size(); j++) {
    value+=J_ij(i1,i2, *pattern[j]);
  }
  total= value/ pattern.size();
  return total;
}


//Define the energies associated to A, B and C but with some of these modified with probability prob: suppose to damage some of the interaction energies stored, deleting randmomnly some J_ij

double J_rand_damage  (int i1, int i2, const vector< Matrix*>  &pattern, double prob ) {
  double value=0, fin_value;
  for (unsigned j = 0; j < pattern.size(); j++) {
    int s1=spin_value(find_m_n(i1)[0],find_m_n(i1)[1], *pattern[j]);
    int s2=spin_value(find_m_n(i2)[0],find_m_n(i2)[1], *pattern[j]);
    value += s1*s2;
  }
  fin_value= value/pattern.size();
  // cout<< "before damaging"<< fin_value<< '\n';
  float r =  ((float) rand())/ (float) RAND_MAX;
  if (r < prob) {
    fin_value = 0;
  }
  return fin_value;
}


//Learing a new pattern: starting from a J_ij the memory can with time learn a new pattern and forget the old one, these proccesses are regulated by the functions a and b 

double new_J_ij(  int i1, int i2, const Matrix &m, double* a, double* b ) {
  int s_i= spin_value(find_m_n(i1)[0],find_m_n(i1)[1],m);
  int s_j= spin_value(find_m_n(i2)[0],find_m_n(i2)[1],m);
  double new_j= ((*b) * J(i1, i2, *p4) + (*a) * s_i *s_j);
  return new_j;
}

//summing over all the possible spins

double new_J(const Matrix &m) {
  double newj=0;
  double oldj=0;
  for (int i = 1; i <= 100 ; i++ ) {
    for (int j=1; j<= 100; j++ ) {
      newj += new_J_ij( i, j, m, a, b);
      oldj += J(i,j,*p4);
    }
  }
  cout<< "the old one is "<< oldj<< '\n';
  cout<<"the new one is "<<  newj<<'\n';
  return newj;
}

//Asymmetry in the energies

double asym_J_ij(int i_in, int i_out, Matrix output, vector<Matrix> possible_inputs, Matrix second_output, vector<Matrix> second_possible_inputs) {
  double J_ij_asym= 0;
  double J_ij_asym_second=0;
  int s_out= spin_value(find_m_n(i_out)[0],find_m_n(i_out)[1], output);
  int s_out_2=spin_value(find_m_n(i_out)[0],find_m_n(i_out)[1], second_output);
  for (unsigned i =0; i< possible_inputs.size(); i++){ 
    int  s_in= spin_value(find_m_n(i_in)[0], find_m_n(i_in)[1], possible_inputs[i]);
    J_ij_asym += s_in * s_out;
  }
  for (unsigned i =0; i< second_possible_inputs.size(); i++){ 
    int  s_in= spin_value(find_m_n(i_in)[0], find_m_n(i_in)[1], second_possible_inputs[i]);
    J_ij_asym_second += s_in * s_out_2;
    //cout<<J_ij_asym<<'\n';
  }
  double tot_asym_J= 1/3. * J_ij_asym+ 1/3. * J_ij_asym_second;
  //cout<< tot_asym_J<<'\n';
  return tot_asym_J;
}
  
    

//Defining the energy of a spin in a random matrix m 

double E_i ( int i, const Matrix &m){
  double sumJ=0, tot_sum;
  int s_i= spin_value(find_m_n(i)[0],find_m_n(i)[1], m);
  for (int j=1; j <= 100; j++) { //loop over all the spins for a given pattern
    int s_j=  spin_value(find_m_n(j)[0],find_m_n(j)[1], m);
    sumJ+=  J(i, j, (*p4) ) * s_j;
  }
  tot_sum= - sumJ * s_i;
  //cout<<  tot_sum<< '\n';
  return tot_sum;
}

//Defining the total energy of a configuration

double Energy(const Matrix &m){
  double sumE=0;
  for ( int i=1; i<=100; i++){
    sumE+= E_i(i, m);
  }
  cout << sumE<< '\n';
  return sumE;
}

//Defining the energy difference between two patterns

double delta_energy (const Matrix &m , const Matrix &n){
  double res= Energy(m)-Energy(n);
  cout << res<< '\n';
  return res;
}


//Defining the asymmetric energies
double asym_E_i ( int i, const Matrix m){
  double sumJ=0, tot_sum;
  int s_i= spin_value(find_m_n(i)[0],find_m_n(i)[1], m);
  for (int j=1; j <= 100; j++) { //loop over all the spins for a given pattern
    int s_j=  spin_value(find_m_n(j)[0],find_m_n(j)[1], m);
    sumJ+=  asym_J_ij(i, j, *output, *possible_inputs, *second_output, *second_possible_inputs ) * s_j;
  }
  tot_sum= - sumJ * s_i;
  //cout<<  tot_sum<< '\n';
  return tot_sum;
}



//Defining the new energies during the learning process

double new_E_i ( int i, const Matrix &m){
  double sumJ=0, tot_sum;
  int s_i= spin_value(find_m_n(i)[0],find_m_n(i)[1], m);
  for (int j=1; j <= 100; j++) { //loop over all the spins for a given pattern
    int s_j=  spin_value(find_m_n(j)[0],find_m_n(j)[1], m);
    sumJ+= new_J_ij(i, j, *p3, a, b ) * s_j;
  }
  tot_sum= - sumJ * s_i;
  // cout<< tot_sum<< '\n';
  return tot_sum;
}



//Energies with damage

double damaged_E_i ( int i, const Matrix &m, float damage){
  double sumJ=0, tot_sum;
  int s_i= spin_value(find_m_n(i)[0],find_m_n(i)[1], m);
  for (int j=1; j <= 100; j++) { //loop over all the spins for a given pattern
    int s_j=  spin_value(find_m_n(j)[0],find_m_n(j)[1], m);
    sumJ+=  J_rand_damage(i, j, (*p4), damage) * s_j; //we assign by hand the percentage of damaging
  }
  tot_sum= - sumJ * s_i;
  //cout<< tot_sum<< '\n';
  return tot_sum;
}


//MonteCarlo method at T=0

Matrix MC (Matrix mat) {
  Matrix new_mat=mat;
  for ( int i = 1; i <=100; i++) {
    if ( ( -1 * new_E_i(i, mat)  - new_E_i(i, mat) ) < 0) {
      int r= find_m_n(i)[0];
      int c= find_m_n(i)[1];
      new_mat[r][c]= (-1) * new_mat[r][c];
    }
  }
  cout<< '\n';
  // cout<<"After one MC sweep: "<< '\n';
  print_true(new_mat);
  return new_mat;
}


//MonteCarlo method at T!=0

Matrix MC_T (Matrix &m, double t) {
  for ( int i = 1; i <=100; i++) {
    double E_flip=( -1 * E_i(i, m ) - E_i(i,m) );
    int r= find_m_n(i)[0];
    int c= find_m_n(i)[1];
    if ( E_flip < 0) {
      m[r][c]= (-1) * m[r][c];
    }
    else {
      double boltz_fac= exp(- E_flip/ t);
      float r = ((float) rand())/ (float) RAND_MAX;
      if (r > boltz_fac) {
	m[r][c]=(-1)* m[r][c];
      }
    }
  }
  cout<< '\n';
  print_true(m);
  return m;
}



//Generate a matrix starting from one of the stored with a percentage 'prob' of spins flipped

Matrix random_matrix (Matrix &m, float prob) {
  float p=prob;
  Matrix mat=m;
  for (int i=1; i<=100; i++){  
    float r =  ((float) rand())/ (float) RAND_MAX;
    // cout<< r<<'\n';
    if (r < prob) {
      int row = find_m_n(i)[0];
      int col = find_m_n(i)[1];
      mat[row][col]= (-1)* mat[row][col];
    }
  }
  float percentage= p* 100;
  cout<< "With "<< percentage << " % of spins flipped" << '\n' << '\n';
  print_true(mat);
  return mat;
}




int main () {

  
//Create the pattern with the letter A

  Matrix A (10, vector<int> (10)  );
  for (unsigned row =0; row < A.size(); row++) {
    for (unsigned col=0; col< A[row].size(); col++) 
      A[row][col]= -1;
      }
  
 
  A[0][4] = up;  A[0][5] = up;  A[1][3] = up;  A[1][4] = up;  A[1][5] = up;  A[1][6] = up;  A[2][2] = up;
  A[2][3] = up;  A[2][6] = up;  A[2][7] = up;  A[3][1] = up;  A[3][2] = up;  A[3][7] = up;  A[3][8] = up;
  A[4][1] = up;  A[4][2] = up;  A[4][3] = up;  A[4][5] = up;  A[4][6] = up;  A[4][7] = up;  A[4][8] = up;
  A[5][1] = up;  A[5][2] = up;  A[5][3] = up;  A[5][4] = up;  A[5][5] = up;  A[5][6] = up;  A[5][7] = up;
  A[5][8] = up;  A[5][9] = up;  A[6][0] = up;  A[6][1] = up;  A[6][2] = up;  A[6][7] = up;  A[6][8] = up;
  A[6][9] = up;  A[7][0] = up;  A[7][1] = up;  A[7][8] = up;  A[7][9] = up;  A[8][0] = up;  A[8][1] = up;
  A[8][8] = up;  A[8][9] = up;  A[9][0] = up;  A[9][1] = up;  A[9][8] = up;  A[9][9] = up;  A[5][0] = up;
  A[4][4]=up;
  p1=&A;
  //cout<< "The pattern A is : " << '\n';  
  //print_true (A);
  


  //First Changing to A font
  
   Matrix A_1 (10, vector<int> (10)  );
  for (unsigned row =0; row < A_1.size(); row++) {
    for (unsigned col=0; col< A_1[row].size(); col++) 
      A_1[row][col]= -1;
      }

  for(int i =3; i<=6; i++){
    A_1[0][i]=up;
  }
   for(int i =2; i<=9; i++){
    A_1[i][1]=up;
    A_1[i][8]=up;
  }
   for(int i =1; i<=8; i++){
    A_1[5][i]=up;
  }
   A_1[1][2]=up; A_1[1][7]=up;
 
  p7=&A_1;
  //cout<< "The pattern A_1 is : " << '\n';  
  // print_true (A_1);

   

   // Second Changing of A font
  
   Matrix A_2 (10, vector<int> (10)  );
  for (unsigned row =0; row < A_2.size(); row++) {
    for (unsigned col=0; col< A_2[row].size(); col++) 
      A_2[row][col]= -1;
      }

  for(int i =3; i<=5; i++){
    A_2[0][i]=up;
  }
   for(int i =4; i<=5; i++){
    A_2[1][i]=up;
  }
   for(int i =3; i<=6; i++){
     A_2[6][i]=up;
  }
   for(int i =3; i<=9; i++){
    A_2[i][2]=up;
    A_2[i][7]=up;
  }
   for(int i =1; i<=3; i++){
     A_2[9][i]=up;
  }
   for(int i =6; i<=8; i++){
     A_2[9][i]=up;
   }

   A_2[2][3]=up; A_2[2][6]=up;
 
  p8=&A_2;
  //cout<< "The pattern A_2 is : " << '\n';  
  //print_true (A_2);
 


//Create the pattern with the letter B

  Matrix B ((10), vector<int> (10)  );
  for (unsigned row =0; row < B.size(); row++) {
    for (unsigned col=0; col< B[row].size(); col++) 
      B[row][col]= -1;
  }
 


  
  B[0][2] = up;  B[0][3] = up;  B[0][4] = up;  B[0][5] = up;  B[0][6] = up;  B[1][2] = up;  B[1][3] = up;
  B[1][6] = up;  B[1][7] = up;  B[2][2] = up;  B[2][3] = up;  B[2][6] = up;  B[2][7] = up;  B[3][2] = up;
  B[3][3] = up;  B[3][6] = up;  B[3][7] = up;  B[4][2] = up;  B[4][3] = up;  B[4][4] = up;  B[4][5] = up;
  B[4][6] = up;  B[5][2] = up;  B[5][3] = up;  B[5][4] = up;  B[5][5] = up;  B[5][6] = up;  B[6][2] = up;
  B[6][3] = up;  B[6][6] = up;  B[6][7] = up;  B[7][2] = up;  B[7][3] = up;  B[7][6] = up;  B[7][7] = up;
  B[8][2] = up;  B[8][3] = up;  B[8][6] = up;  B[8][7] = up;  B[9][2] = up;  B[9][3] = up;  B[9][4] = up;
  B[9][5] = up;  B[9][6] = up;
  p2= &B;
  
  // cout<< "The pattern B is : " << '\n';
  //print_true (B);
  

  //First new font for B

  
  Matrix B_1 ((10), vector<int> (10)  );
  for (unsigned row =0; row < B_1.size(); row++) {
    for (unsigned col=0; col< B_1[row].size(); col++) 
      B_1[row][col]= -1;
  }
 
  for (int i =2; i<=6; i++){
    B_1[0][i]=up;
  }
  for (int i =1; i<=9; i++){
    B_1[i][3]=up;
  }
  for (int i =2; i<=6; i++){
    B_1[9][i]=up;
  }
  for (int i =4; i<=6; i++){
    B_1[5][i]=up;
  }
  for (int i =2; i<=4; i++){
    B_1[i][7]=up;
  }
  for (int i =6; i<=8; i++){
    B_1[i][7]=up;
  }
    
   p9=&B_1;
  //cout<< "The pattern B_1 is : " << '\n';  
  //print_true (B_1);
 
  

    //Second new font for B

  
  Matrix B_2 ((10), vector<int> (10)  );
  for (unsigned row =0; row < B_2.size(); row++) {
    for (unsigned col=0; col< B_2[row].size(); col++) 
      B_2[row][col]= -1;
  }
 
  for (int i =3; i<=6; i++){
    B_2[0][i]=up;
  }
  for (int i =1; i<=9; i++){
    B_2[i][3]=up;
  }
  for (int i =3; i<=6; i++){
    B_2[9][i]=up;
  }
  for (int i =4; i<=6; i++){
    B_2[3][i]=up;
  }
  for (int i =4; i<=8; i++){
    B_2[i][7]=up;
  }
  B_2[2][7]=up;
    
   p10=&B_2;
  //cout<< "The pattern B_2 is : " << '\n';  
  //print_true (B_2);
 
  
  

//Create the pattern with the letter C

  Matrix C (10, vector<int> (10)  );
  for (unsigned row =0; row < C.size(); row++) {
    for (unsigned col=0; col< C[row].size(); col++) 
      C[row][col]= -1;
  }
 

  C[0][3] = up;  C[0][4] = up;  C[0][5] = up;  C[0][6] = up;  C[0][7] = up;  C[0][8] = up;  C[1][2] = up;
  C[1][3] = up;  C[1][4] = up;  C[1][5] = up;  C[1][6] = up;  C[1][7] = up;  C[1][8] = up;  C[2][2] = up;
  C[2][3] = up;  C[2][4] = up;  C[3][2] = up;  C[3][3] = up;  C[4][2] = up;  C[4][3] = up;  C[5][2] = up;
  C[5][2] = up;  C[5][3] = up;  C[6][2] = up;  C[6][3] = up;  C[7][2] = up;  C[7][3] = up;  C[8][2] = up;
  C[8][3] = up;  C[8][4] = up;  C[8][5] = up;  C[8][6] = up;  C[8][7] = up;  C[8][8] = up;  C[9][3] = up;
  C[9][4] = up;  C[9][5] = up;  C[9][6] = up;  C[9][7] = up;  C[9][8] = up;  C[7][4] = up;
  p3= &C;
  //cout<< "The pattern C is : " << '\n';
  // print_true (C);


  
//Create the pattern with the letter E

  Matrix E (10, vector<int> (10)  );
  for (unsigned row =0; row < E.size(); row++) {
    for (unsigned col=0; col< E[row].size(); col++) 
      E[row][col]= -1;
  }
 
  for(int i=2; i<= 7; i++) {
    E[1][i]= up;
  }
  for(int i=2; i<= 7; i++) {
    E[0][i]= up;
  }
  for(int i=2; i<= 5; i++) {
    E[4][i]= up;
  }
  for(int i=2; i<= 5; i++) {
    E[5][i]= up;
  }
  for(int i=2; i<= 7; i++) {
    E[9][i]= up;
  }
  for(int i=2; i<= 7; i++) {
    E[8][i]= up;
  }
   
  E[3][2] = up;  E[3][3] = up;  E[6][2] = up;  E[6][3] = up;  E[2][2] = up;  E[2][3] = up;  E[7][2] = up;  E[7][3] = up;
  p5= &E;
  
  //cout<< "The pattern E is : " << '\n';
  // print_true (E);

  
//Create the pattern with the letter D

  Matrix D (10, vector<int> (10)  );
  for (unsigned row =0; row < D.size(); row++) {
    for (unsigned col=0; col< D[row].size(); col++) 
      D[row][col]= -1;
  }
 
  for(int i=0; i<= 9; i++) {
    D[i][2]= up;
    D[i][3]= up;
  }
  for(int i= 4; i<= 7; i++) {
    D[0][i]= up;
  }
  for(int i=4; i<= 8; i++) {
    D[1][i]= up;
  }
  for(int i=7; i<= 8; i++) {
    D[2][i]= up;
    D[3][i]= up;
    D[4][i]= up;
    D[5][i]= up;
    D[6][i]= up;
    D[7][i]= up;
  }
  for(int i=4; i<= 8; i++) {
    D[8][i]= up;
  }
  for(int i=4; i<= 7; i++) {
    D[9][i]= up;
  }
  
  p6= &D;
  
  //cout<< "The pattern D is : " << '\n';
  //  print_true (D);



  
  //Defining the stored patterns
  
  vector<Matrix*> pattern {p1,p2};
  p4= &pattern;

  vector<Matrix*> configurations { p1,p2,p3,p5};

  //Matrix rand_matrix=random_matrix(B,.3);
  // MC(rand_matrix);

  //----------------------------------

  //ASYMMETRIC ENERGIES

  
  /* 
  vector<Matrix> inputs;
  vector<Matrix> second_inputs;
  output= &(*pattern[0]);
  second_output=&(*pattern[1]);
  //print_true(*output);
  
  inputs.push_back(*p1);
  inputs.push_back(*p7);
  inputs.push_back(*p8);

  second_inputs.push_back(*p2);
  second_inputs.push_back(*p9);
  second_inputs.push_back(*p10);
 
  possible_inputs= &inputs;
  second_possible_inputs=&second_inputs;
  
  //input_matrix=&inputs[1];
  // cout<< "the matrix given to MonteCarlo is"<< '\n';
  // print_true(inputs[2]);

  Matrix rand_matrix=random_matrix(*p7, .3);
   MC(rand_matrix);
   Matrix rand_matrix1=random_matrix(*p7, .4);
   MC(rand_matrix1);
  */
  //----------------------------------------
  //DAMAGED ENERGIES
  /*

  float p; 
  
  cout << "Enter a number between 0 and 1: ";
  cin >> p;
   */  

  // Matrix rand_matrix= random_matrix( B , .01);
  // MC(rand_matrix); 
  // delta_energy(rand_matrix, A);

  //hamm_dist(rand_matrix,A);

  // for (int i =0; i < 10; i++) {
  // rand_matrix = MC(rand_matrix);
  //}

  // --------------------------------------------
 
  //---------------------------------------------------
  //LEARNING PROCESS
 /*
  double a_= .1;
  double b_= 1;
  a= &a_;
  b= &b_;
  
  // Matrix rand_matrix= random_matrix(C, .3);
  MC(C);
  
 */ 
  
}

