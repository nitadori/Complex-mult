#include <complex>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using cmplx = std::complex<double>;

struct ComplexMulLhs{
	double sum;
	double diff;
	double imag;

	ComplexMulLhs(){}

	ComplexMulLhs(const cmplx &z){
		sum  = z.real() + z.imag();
		diff = z.real() - z.imag();
		imag = z.imag();
	}
};

struct ComplexMulRhs{
	double real;
	double imag;
	double diff;

	ComplexMulRhs(){}
	
	ComplexMulRhs(const cmplx &z){
		real = z.real();
		imag = z.imag();
		diff = real - imag;
	}
};

struct ComplexMulAccum{
	double real;
	double imag;
	double common;

	operator cmplx() const {
		return {real+common, imag+common};
	}

	void mult(const ComplexMulLhs &lhs, const ComplexMulRhs &rhs){
		real = lhs.diff * rhs.real;
		imag = lhs.sum  * rhs.imag;
		common = lhs.imag * rhs.diff;
	}
	void madd(const ComplexMulLhs &lhs, const ComplexMulRhs &rhs){
		real   += lhs.diff * rhs.real;
		imag   += lhs.sum  * rhs.imag;
		common += lhs.imag * rhs.diff;
	}
	// conj(lhs) * rhs
	void mult_lconj(const ComplexMulLhs &lhs, const ComplexMulRhs &rhs){
		real = lhs.sum  * rhs.real;
		imag = lhs.diff * rhs.imag;
		common = -lhs.imag * rhs.diff;
	}
	void madd_lconj(const ComplexMulLhs &lhs, const ComplexMulRhs &rhs){
		real   += lhs.sum  * rhs.real;
		imag   += lhs.diff * rhs.imag;
		common -= lhs.imag * rhs.diff;
	}
};

extern "C"
__attribute__((noinline))
void matmul1(cmplx b[2][3], const cmplx A[3][3], const cmplx x[2][3]){
	for(int s=0; s<2; s++){
		for(int i=0; i<3; i++){
			b[s][i] = A[i][0] * x[s][0] 
			        + A[i][1] * x[s][1] 
					+ A[i][2] * x[s][2] ;
		}
	}
}
extern "C"
__attribute__((noinline))
void matmul_dag1(cmplx b[2][3], const cmplx A[3][3], const cmplx x[2][3]){
	for(int s=0; s<2; s++){
		for(int i=0; i<3; i++){
			b[s][i] = conj(A[0][i]) * x[s][0] 
			        + conj(A[1][i]) * x[s][1] 
					+ conj(A[2][i]) * x[s][2] ;
		}
	}
}
extern "C"
__attribute__((noinline))
void matmul2(cmplx b[2][3], const cmplx A[3][3], const cmplx x[2][3]){
	ComplexMulLhs Al[3][3];
	for(int a=0; a<3; a++) for(int b=0; b<3; b++)  Al[a][b] = ComplexMulLhs(A[a][b]);
	ComplexMulRhs xr[2][3];
	for(int s=0; s<2; s++) for(int c=0; c<3; c++)  xr[s][c] = ComplexMulRhs(x[s][c]);
	ComplexMulAccum ba[2][3];
	for(int s=0; s<2; s++){
		for(int i=0; i<3; i++){
			ba[s][i].mult(Al[i][0], xr[s][0]);
			ba[s][i].madd(Al[i][1], xr[s][1]);
			ba[s][i].madd(Al[i][2], xr[s][2]);

			b[s][i] = cmplx(ba[s][i]);
		}
	}
}
extern "C"
__attribute__((noinline))
void matmul_dag2(cmplx b[2][3], const cmplx A[3][3], const cmplx x[2][3]){
	ComplexMulLhs Al[3][3];
	for(int a=0; a<3; a++) for(int b=0; b<3; b++)  Al[a][b] = ComplexMulLhs(A[a][b]);
	ComplexMulRhs xr[2][3];
	for(int s=0; s<2; s++) for(int c=0; c<3; c++)  xr[s][c] = ComplexMulRhs(x[s][c]);
	ComplexMulAccum ba[2][3];
	for(int s=0; s<2; s++){
		for(int i=0; i<3; i++){
			ba[s][i].mult_lconj(Al[0][i], xr[s][0]);
			ba[s][i].madd_lconj(Al[1][i], xr[s][1]);
			ba[s][i].madd_lconj(Al[2][i], xr[s][2]);

			b[s][i] = cmplx(ba[s][i]);
		}
	}
}

void mm_test(){
	cmplx A[3][3], x[2][3], b[2][3];
	for(int a=0; a<3; a++) for(int b=0; b<3; b++)  A[a][b] = {drand48(), drand48()};
	for(int s=0; s<2; s++) for(int c=0; c<3; c++)  x[s][c] = {drand48(), drand48()};

	matmul1(b, A, x);
	std::cout << b[0][0] << b[0][1] << b[0][2] << "\t";
	std::cout << b[1][0] << b[1][1] << b[1][2] << std::endl;

	matmul2(b, A, x);
	std::cout << b[0][0] << b[0][1] << b[0][2] << "\t";
	std::cout << b[1][0] << b[1][1] << b[1][2] << std::endl;

	matmul_dag1(b, A, x);
	std::cout << b[0][0] << b[0][1] << b[0][2] << "\t";
	std::cout << b[1][0] << b[1][1] << b[1][2] << std::endl;

	matmul_dag2(b, A, x);
	std::cout << b[0][0] << b[0][1] << b[0][2] << "\t";
	std::cout << b[1][0] << b[1][1] << b[1][2] << std::endl;
}


int main(){
	srand48(57);

	cmplx lhs{drand48(), drand48()};
	cmplx rhs{drand48(), drand48()};
	cmplx prod1 = lhs * rhs;
	cmplx cprod1 = conj(lhs) * rhs;


	ComplexMulAccum acc;
	acc.mult(ComplexMulLhs(lhs), ComplexMulRhs(rhs));
	cmplx prod2 = acc;
	acc.mult_lconj(ComplexMulLhs(lhs), ComplexMulRhs(rhs));
	cmplx cprod2 = acc;

	std::cout << prod1 << std::endl;
	std::cout << prod2 << std::endl;

	std::cout << cprod1 << std::endl;
	std::cout << cprod2 << std::endl;

	puts("");
	mm_test();

	return 0;

}
