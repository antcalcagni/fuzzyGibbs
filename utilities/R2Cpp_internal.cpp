#include <Rcpp.h>
#include <cmath>        
#include <cstdlib>     
#include <vector>
#include <functional>
#include <cstdint> //`uint32_t`
#include <iostream>
#include <stdexcept>
#include <random>
#include <numeric> // For std::iota
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
using namespace std;
namespace ublas = boost::numeric::ublas;


// ----- Insert your own Cpp functions here ----- //
// ----- This section can be modified at your convenience! ----- //

// ----- head ----- //

class Random {
private:
  static std::mt19937 rng; // Static RNG instance
  
public:
  // Manually set the seed
  static void seed(uint32_t seed) {
    rng.seed(static_cast<std::mt19937::result_type>(seed));
  }
  
  // Auto set.seed() using time
  static void seed() {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = now.time_since_epoch();
    uint64_t nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
  }
  
  // Access the generator
  static std::mt19937& get_rng() {
    return rng;
  }
  
  // Uniform distribution (similar to R's runif)
  static double runif(double min = 0.0, double max = 1.0) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
  }
  
  // Normal distribution (similar to R's rnorm)
  static double rnorm(double mean = 0.0, double sd = 1.0) {
    std::normal_distribution<double> dist(mean, sd);
    return dist(rng);
  }
  
  // Beta distribution (similar to R's rbeta)
  static double rbeta(double alpha, double beta) {
    if (alpha <= 0 || beta <= 0) {
      throw std::invalid_argument("Alpha and Beta parameters must be positive.");
    }
    
    std::gamma_distribution<double> gamma_alpha(alpha, 1.0);
    std::gamma_distribution<double> gamma_beta(beta, 1.0);
    
    double x = gamma_alpha(rng); // Sample from Gamma(alpha, 1)
    double y = gamma_beta(rng);  // Sample from Gamma(beta, 1)
    
    return x / (x + y); // Transform to Beta(alpha, beta)
  }
};

const double pi = 3.14159265358979323846;
const double NaN = std::numeric_limits<double>::quiet_NaN();

template <typename T>
  void print_vector(const std::vector<T>& vec) {
    for (const T& val : vec) {
      std::cout << val << " ";
    }
    std::cout << std::endl;
  }

template <typename T>
  void print(const T& value) {
    std::cout << value << std::endl;
  }
  
template <typename T>
std::vector<T> c(std::initializer_list<T> values) {
    return std::vector<T>(values);
}  

double rescale(double x, const double lb=0.0, const double ub=1.0){
  double xout = std::nan("");
  if(lb<ub){xout = (x-lb)/(ub-lb);}
  return xout;
}

double rescale_inv(double x, const double lb=0.0, const double ub=1.0){
  double xout = std::nan("");
  if(lb<ub){xout = lb+(ub-lb)*x;}
  return xout;
}

double matrixmult1d(const vector<double>& x, const vector<double>& b){
  if (x.size() != b.size()) {
    throw std::invalid_argument("matrixMult1d: Vectors must have the same size");
  }
  
  double out = 0.0;
  for(int j = 0; j < x.size(); j++){
    out += x[j]*b[j];
  }
  return out;
}

// ..the same as above but it works for int types
double matrixmult1d(const vector<double>& x, const vector<int>& b){
  if (x.size() != b.size()) {
    throw std::invalid_argument("matrixMult1d: Vectors must have the same size");
  }
  
  double out = 0.0;
  for(int j = 0; j < x.size(); j++){
    if (b[j] != 0){out += x[j];}
  }
  return out;
}

// Overload for vector / scalar
std::vector<double> operator/(const std::vector<double>& vec, double scalar) {
    if (scalar == 0.0) {
        throw std::invalid_argument("Division by zero is not allowed.");
    }

    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] / scalar;
    }
    return result;
}

// Overload for scalar / vector
std::vector<double> operator/(double scalar, const std::vector<double>& vec) {
    if (scalar == 0.0) {
        throw std::invalid_argument("Division by zero is not allowed.");
    }

    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar / vec[i];
    }
    return result;
}

// Overload for vector * scalar
std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

// Overload for scalar * vector (commutative)
std::vector<double> operator*(double scalar, const std::vector<double>& vec) {
    return vec * scalar;
}

std::vector<double> x_seq(const std::vector<double>& x,int a, int b){
  vector<double> xout;
  if((b-a)==0){a-=1;}
    for(int i = a; i < b; i++){
      xout.push_back(x[i]);
      }
    return xout;
}

// Overload for pow function to work with vectors with a single element only x[0]
double pow(const std::vector<double>& x, int b) {
  return pow(x[0],b);
}

// Function to compute Cartesian product of two integer vectors
std::vector<std::vector<int>> expand_grid(const std::vector<int>& vec1, const std::vector<int>& vec2,bool unique=false) {
  std::vector<std::vector<int>> result;
  
  // Generate Cartesian product
  for (int i : vec1) {
    for (int j : vec2) {
      if(unique){
        if(i>=j){result.push_back({i, j});}
      }else{
        result.push_back({i, j});
      }
    }
  }
  return result;
}

// Convert a vector to a symmetric matrix
std::vector<std::vector<double>> vector_to_symmetric_matrix(const std::vector<double>& x, int JH,bool symm=true) {
  
  int K = (JH * (JH + 1)) / 2;
  if (x.size() != K) {
    throw std::invalid_argument("Vector size does not match the required lower triangular matrix size.");
  }
  
  std::vector<std::vector<double>> X(JH, std::vector<double>(JH, 0.0));
  
  int k = 0; // Index for the input vector
  for (int j = 0; j < JH; ++j) {
    for (int i = 0; i < JH; ++i) { 
      if(i>=j){
        X[i][j] = x[k];
        if(symm){X[j][i] = x[k];} // Make it symmetric 
        ++k;
        }
    }
  }
  return X;
}

double digamma(const double x){
  return boost::math::digamma(x);
}

double trigamma(const double x){
  return boost::math::trigamma(x);
}

// Create a diagonal matrix from a vector
std::vector<std::vector<double>> diag(const std::vector<double>& vec) {
    size_t n = vec.size();
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        matrix[i][i] = vec[i]; // Set diagonal elements
    }
    return matrix;
}

// Overload for vector - scalar
std::vector<double> operator-(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

// Overload for vector + scalar
std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

// Matrix multiplication function
std::vector<std::vector<double>> matrix_multiply(const std::vector<std::vector<double>>& A,
                                                 const std::vector<std::vector<double>>& B) {
  size_t rows_A = A.size();
  size_t cols_A = A[0].size();
  size_t rows_B = B.size();
  size_t cols_B = B[0].size();
  
  // Ensure dimensions are compatible for multiplication
  if (cols_A != rows_B) {
    throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
  }
  
  // Initialize the result matrix with zeros
  std::vector<std::vector<double>> result(rows_A, std::vector<double>(cols_B, 0.0));
  
  // Perform multiplication
  for (size_t i = 0; i < rows_A; ++i) {
    for (size_t j = 0; j < cols_B; ++j) {
      for (size_t k = 0; k < cols_A; ++k) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  
  return result;
}

// Matrix-Vector Multiplication
std::vector<double> matrix_multiply(const std::vector<std::vector<double>>& A,
                                    const std::vector<double>& b) {
  size_t rows_A = A.size();
  size_t cols_A = A[0].size();
  
  // Ensure dimensions are compatible
  if (cols_A != b.size()) {
    throw std::invalid_argument("Matrix and vector dimensions do not match.");
  }
  
  // Result vector
  std::vector<double> result(rows_A, 0.0);
  
  // Matrix-vector multiplication
  for (size_t i = 0; i < rows_A; ++i) {
    for (size_t j = 0; j < cols_A; ++j) {
      result[i] += A[i][j] * b[j];
    }
  }
  
  return result;
}

// Vector-Matrix Multiplication
std::vector<double> matrix_multiply(const std::vector<double>& a,
                                    const std::vector<std::vector<double>>& B) {
  size_t cols_B = B[0].size();
  size_t rows_B = B.size();
  
  // Ensure dimensions are compatible
  if (a.size() != rows_B) {
    throw std::invalid_argument("Vector and matrix dimensions do not match.");
  }
  
  // Result vector
  std::vector<double> result(cols_B, 0.0);
  
  // Vector-matrix multiplication
  for (size_t j = 0; j < cols_B; ++j) {
    for (size_t i = 0; i < rows_B; ++i) {
      result[j] += a[i] * B[i][j];
    }
  }
  
  return result;
}

// Matrix subtraction function
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& A,
                                           const std::vector<std::vector<double>>& B) {
    size_t rows = A.size();
    size_t cols = A[0].size();

    // Ensure matrices have the same dimensions
    if (rows != B.size() || cols != B[0].size()) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction.");
    }

    // Result matrix
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols, 0.0));

    // Perform element-wise subtraction
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }

    return result;
}

// Matrix addition function
std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& A,
                                           const std::vector<std::vector<double>>& B) {
  size_t rows = A.size();
  size_t cols = A[0].size();
  
  // Ensure matrices have the same dimensions
  if (rows != B.size() || cols != B[0].size()) {
    throw std::invalid_argument("Matrix dimensions must match for subtraction.");
  }
  
  // Result matrix
  std::vector<std::vector<double>> result(rows, std::vector<double>(cols, 0.0));
  
  // Perform element-wise subtraction
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      result[i][j] = A[i][j] + B[i][j];
    }
  }
  
  return result;
}

template<typename T>
ublas::matrix<T> invert_matrix(const ublas::matrix<T>& input) {
  using namespace boost::numeric::ublas;
  
  // Create a working copy of the input
  matrix<T> A(input);
  
  // Create a permutation matrix for LU decomposition
  permutation_matrix<std::size_t> pm(A.size1());
  
  // Perform LU factorization
  if (lu_factorize(A, pm) != 0) {
    throw std::runtime_error("Matrix is singular and cannot be inverted.");
  }
  
  // Create an identity matrix of the same size as input
  matrix<T> inverse = identity_matrix<T>(A.size1());
  
  // Perform back substitution to get the inverse
  lu_substitute(A, pm, inverse);
  
  return inverse;
}

// Function to convert std::vector<std::vector<double>> to boost::ublas::matrix
ublas::matrix<double> to_ublas_matrix(const std::vector<std::vector<double>>& std_matrix) {
  size_t rows = std_matrix.size();
  size_t cols = std_matrix.empty() ? 0 : std_matrix[0].size();
  
  // Initialize the ublas matrix with the same dimensions
  ublas::matrix<double> ublas_matrix(rows, cols);
  
  // Copy data from std::vector<std::vector<double>> to ublas::matrix
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      ublas_matrix(i, j) = std_matrix[i][j];
    }
  }
  
  return ublas_matrix;
}

// Function to convert boost::ublas::matrix to std::vector<std::vector<double>>
std::vector<std::vector<double>> to_std_matrix(const ublas::matrix<double>& ublas_matrix) {
  size_t rows = ublas_matrix.size1();
  size_t cols = ublas_matrix.size2();
  
  // Initialize a std::vector<std::vector<double>> with the same dimensions
  std::vector<std::vector<double>> std_matrix(rows, std::vector<double>(cols));
  
  // Copy data from ublas::matrix to std::vector<std::vector<double>>
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      std_matrix[i][j] = ublas_matrix(i, j);
    }
  }
  
  return std_matrix;
}

ublas::matrix<double> vector_to_ublas(const std::vector<double>& flat_vector, size_t rows, size_t cols) {
  if (flat_vector.size() != rows * cols) {
    throw std::invalid_argument("Flat vector size does not match the specified matrix dimensions.");
  }
  
  // Initialize the Boost matrix
  ublas::matrix<double> ublas_matrix(rows, cols);
  
  // Copy the data from the vector to the matrix
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      ublas_matrix(i, j) = flat_vector[i * cols + j];
    }
  }
  
  return ublas_matrix;
}

std::vector<double> ublas_to_vector(const ublas::matrix<double>& ublas_matrix) {
  size_t rows = ublas_matrix.size1();
  size_t cols = ublas_matrix.size2();
  
  // Create a flat vector with the required size
  std::vector<double> flat_vector(rows * cols);
  
  // Copy the data from the matrix to the vector
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      flat_vector[i * cols + j] = ublas_matrix(i, j);
    }
  }
  
  return flat_vector;
}

double dnorm(double x, double mean = 0.0, double sd = 1.0, bool log = false) {
  const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
  double z = (x - mean) / sd;
  double density = inv_sqrt_2pi / sd * std::exp(-0.5 * z * z);

  return log ? std::log(density) : density;
}

std::vector<double> abs(vector<double>& x){
  int n = x.size();
  vector<double> xabs(n);
  for(int i=0;i<n;++i){
    xabs[i] = abs(x[i]);
  }
  return xabs;
}

std::mt19937 Random::rng; //= std::mt19937{std::random_device{}()}; // Default seed

// ----- R2Cpp conversion starts here ----- //
// ----- Please edit manually only if really needed! ----- //

 
// [[Rcpp::export]]
double fy(double y,vector<double> x,vector<double> z,double lb,double ub,vector<double> pars)
{
auto     b =  x_seq(pars,0, x.size()) ;
auto     g =  x_seq(pars,(0 +x.size()), (x.size() + z.size())) ;
auto     eta = matrixmult1d(x, b) ;
auto     theta = exp(eta)/(1 + exp(eta)) ;
auto     phi = exp(matrixmult1d(z, g)) ;
auto     alpha = theta * phi ;
auto     beta = phi - theta * phi ;
return  ( pow( (y-lb),(alpha-1) ) )*( pow( (ub-y),(beta-1) ) )/( pow( (ub-lb),(alpha+beta-1) ) *(tgamma(alpha)*tgamma(beta)/tgamma(alpha+beta))) ;
}

 
// [[Rcpp::export]]
double d1lnfy(double y,vector<double> x,vector<double> z,double lb,double ub,vector<double> pars)
{
auto     e1 = x.size() ;
auto     e3 = exp(matrixmult1d(x,  x_seq(pars,0, e1))) ;
auto     e5 = exp(matrixmult1d(z,  x_seq(pars,(0 +e1), (e1 + z.size())))) ;
auto     e6 = 1 + e3 ;
auto     e7 = (1 - e3/e6) * e5 ;
auto     e9 = e3 * e5/e6 ;
auto     e10 = e7 - 1 ;
auto     e11 = e9 - 1 ;
auto     e12 = ub - y ;
auto     e13 = y - lb ;
return  (e11* pow( e12,e10 ) * pow( e13,(e9-2) ) -e10* pow( e12,(e7-2) ) * pow( e13,e11 ) )*tgamma(1*e5)/(fy(y,x,z,lb,ub,pars)*tgamma(e7)*tgamma(e9)* pow( (ub-lb),(e5-1) ) ) ;
}

 
// [[Rcpp::export]]
double d2lnfy(double y,vector<double> x,vector<double> z,double lb,double ub,vector<double> pars)
{
auto     e1 = x.size() ;
auto     e3 = exp(matrixmult1d(x,  x_seq(pars,0, e1))) ;
auto     e5 = exp(matrixmult1d(z,  x_seq(pars,(0 +e1), (e1 + z.size())))) ;
auto     e6 = 1 + e3 ;
auto     e7 = (1 - e3/e6) * e5 ;
auto     e9 = e3 * e5/e6 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e12 = e7 - 1 ;
auto     e13 = e9 - 1 ;
auto     e14 = e7 - 2 ;
auto     e15 = e9 - 2 ;
auto     e16   =  pow( e10,e14 ) ;
auto     e17   =  pow( e11,e15 ) ;
auto     e18 = e12 * e16 ;
auto     e20   =  pow( e10,e12 ) ;
auto     e21   =  pow( e11,e13 ) ;
auto     e22 = fy(y, x, z, lb, ub, pars) ;
auto     e25   =  tgamma(e7)*tgamma(e9)* pow( (ub-lb),(e5-1) ) ;
auto     e26 = tgamma(e5) ;
return  (((e15*e20* pow( e11,(e9-3) ) -e18*e17)*e13-e12*(e13*e16*e17-e14* pow( e10,(e7-3) ) *e21))*e26/e25- pow( ((e13*e20*e17-e18*e21)*e26/e25),2 ) /e22)/e22 ;
}

// ----- tail ----- //

// [[Rcpp::export]]
std::vector<double> approx_pi_y_theta(double m,double s,const vector<double>& x,const vector<double>& z,const vector<double>& pars, double lambda0=NA_REAL,double sigma0=NA_REAL,const double lb=0.0,const double ub=1.0,const double eps=1E-9,const bool print_out=1){

  double lambda = lambda0; 
  double sigma = sigma0; 
  double y = 1.99; 
  m=rescale(m,lb,ub);
  
  int k = 0; double K1; double K2;
  while((abs(y/lambda-1))>eps){ 
    y = lambda;
    K1 = d1lnfy(rescale_inv(y,lb,ub),x,z,lb,ub,pars) + s*log(m/(1 - m)) + s*(-boost::math::polygamma(0,s*y) + boost::math::polygamma(0,s - s*y));
    K2 = d2lnfy(rescale_inv(y,lb,ub),x,z,lb,ub,pars) - pow(s,2)*(boost::math::polygamma(1,s*y) + boost::math::polygamma(1,s - s*y));
    lambda = (1 + y*(-2 + K1 + sigma - K1*y))/sigma;
    sigma = (1 - (-1 + y)*y*(-2 + K2*(-1 + y)*y))/(lambda - 2*lambda*y + pow(y,2));
    
    if(print_out){std::cout << k << ", " << lambda << ", " << rescale_inv(lambda,lb,ub) << ", " << sigma << std::endl;}
    if(std::isnan(lambda)||std::isnan(sigma)){break;}
    k++;
  }
  
  return {lambda, sigma};
}

// [[Rcpp::export]]
std::vector<double> sampling_y_theta(double m,double s,const vector<double>& x,const vector<double>& z,const vector<double>& pars, double lambda0=NA_REAL,double sigma0=NA_REAL,const double lb=0.0,const double ub=1.0,const int maxiter=150,const double eps=1E-9,const bool print_out=0, const int seedx=-1, const double sigma0delta=20.0){
  Random::seed();
  
  vector<double> xout(5);
  
  if(std::isnan(lambda0)){lambda0 = rescale(m,lb,ub);} //starting point for lambda
  
  // Smart initialization for sigma0 -- Note: starting point for sigma0 is what makes stochastic this algorithm
  vector<double> sigma_k(maxiter); // -1 is a sentinel value for sigma0
  if(sigma0 < 0){ 
    for(int i=0; i<maxiter; i++){sigma_k[i] = s + Random::runif(abs(s-sigma0delta),s+sigma0delta);}
  }else if(std::isnan(sigma0)){
    for(int i=0; i<maxiter; i++){sigma_k[i] = Random::runif(2.99,300.99);} 
  }else{
    for(int i=0; i<maxiter; i++){sigma_k[i] = sigma0;} 
  }
    
  int k = 0;
  int cond = 0;
  
  // Run the 4P-Beta approximation of pi(y|theta) --solutions {lambda,sigma} are searched untile they are feasible
  vector<double> out(2);
  while(cond<=0 && k<=maxiter){
    out = approx_pi_y_theta(m,s,x,z,pars,lambda0,sigma_k[k],lb,ub,eps,print_out);
    cond = 2-( static_cast<int>(std::isnan(out[0])) + static_cast<int>(std::isnan(out[1])) );
    k++;
  }
  
  // Draw a random realization from the (0,1)-Beta distribution and linearly transform into [lb,ub]
  double ystar;
  if(seedx>0){Random::seed(seedx);}else{Random::seed();} // set the seed of the RNG chain only if required
  if(cond>0 && k<=maxiter){ //random samples are obtained only if the approximation converged
    ystar  = rescale_inv(Random::rbeta(out[0]*out[1],out[1]-out[1]*out[0]),lb,ub);
  }else{
    std::cerr << "approx_pi_y_theta failed: No solution found. Using observed values." << std::endl;
    ystar = rescale_inv(Random::rbeta(m*s,s-m*s),lb,ub);
  }
  
  xout = {ystar,out[0],out[1],static_cast<double>(cond),static_cast<double>(k)};
  return xout;
}

// [[Rcpp::export]]
vector<vector<double>> sampling_y_theta_loop(vector<double>& m, vector<double>& s, const vector<vector<double>>& X,vector<vector<double>>& Z,vector<double> pars,double lb = 0.0,double ub = 1.0,double lambda0 = NA_REAL, double sigma0 = NA_REAL, double eps = 1E-9,bool print_out = 0,int maxiter = 150,int seedx = -1,double sigma0delta=20.0){
  int n = X.size(); // Number of rows in X
  vector<vector<double>> out(n, vector<double>(5));
  bool cond; 
  int k;
    
  for (int i = 0; i < n; i++) {
    if(sigma0==-2){ //brute-force way to retrieve feasible approximations (it depends heavily by sigma0)
      cond = true;
      k = 0;
      while(cond || k < 5e2){
        k += 1;
        out[i] = sampling_y_theta(m[i], s[i], X[i], Z[i], pars, lambda0, std::nan(""), lb, ub, maxiter, eps, print_out, seedx, sigma0delta);
        if(out[i][3]>1){cond=false;}
      }
    }else{
      out[i] = sampling_y_theta(m[i], s[i], X[i], Z[i], pars, lambda0, sigma0, lb, ub, maxiter, eps, print_out, seedx);
    }
  }

  return out; // Return the result
}



 
// [[Rcpp::export]]
double fy_internal(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     eta = matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7})) ;
auto     theta = exp(eta)/(1 + exp(eta)) ;
auto     phi = exp(matrixmult1d(z, c({b8}))) ;
auto     alpha = theta * phi ;
auto     beta = phi - theta * phi ;
return  ( pow( (y-lb),(alpha-1) ) )*( pow( (ub-y),(beta-1) ) )/( pow( (ub-lb),(alpha+beta-1) ) *(tgamma(alpha)*tgamma(beta)/tgamma(alpha+beta))) ;
}

 
// [[Rcpp::export]]
auto grad_lnfy_internal1(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e7 = e6 * e4 ;
auto     e9 = e2 * e4/e5 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( (ub-lb),(e4-1) ) ;
auto     e14   =  pow( e10,(e7-1) ) ;
auto     e15   =  pow( e11,(e9-1) ) ;
auto     e16 = tgamma(e7) ;
auto     e17 = tgamma(e9) ;
auto     e19 = e16 * e17 * e13 ;
return  ((log(e11)*e14*e15-log(e10)*e14*e15)/e19-(digamma(e9)-digamma(e7))*e16*e17*e13*e14*e15/ pow( e19,2 ) )*e6*e2*e4*tgamma(1*e4)*matrixmult1d(x,c({1,0,0,0,0,0,0}))/(e5*fy_internal(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)) ;
}

 
// [[Rcpp::export]]
auto grad_lnfy_internal2(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e7 = e6 * e4 ;
auto     e9 = e2 * e4/e5 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( (ub-lb),(e4-1) ) ;
auto     e14   =  pow( e10,(e7-1) ) ;
auto     e15   =  pow( e11,(e9-1) ) ;
auto     e16 = tgamma(e7) ;
auto     e17 = tgamma(e9) ;
auto     e19 = e16 * e17 * e13 ;
return  ((log(e11)*e14*e15-log(e10)*e14*e15)/e19-(digamma(e9)-digamma(e7))*e16*e17*e13*e14*e15/ pow( e19,2 ) )*e6*e2*e4*tgamma(1*e4)*matrixmult1d(x,c({0,1,0,0,0,0,0}))/(e5*fy_internal(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)) ;
}

 
// [[Rcpp::export]]
auto grad_lnfy_internal3(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e7 = e6 * e4 ;
auto     e9 = e2 * e4/e5 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( (ub-lb),(e4-1) ) ;
auto     e14   =  pow( e10,(e7-1) ) ;
auto     e15   =  pow( e11,(e9-1) ) ;
auto     e16 = tgamma(e7) ;
auto     e17 = tgamma(e9) ;
auto     e19 = e16 * e17 * e13 ;
return  ((log(e11)*e14*e15-log(e10)*e14*e15)/e19-(digamma(e9)-digamma(e7))*e16*e17*e13*e14*e15/ pow( e19,2 ) )*e6*e2*e4*tgamma(1*e4)*matrixmult1d(x,c({0,0,1,0,0,0,0}))/(e5*fy_internal(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)) ;
}

 
// [[Rcpp::export]]
auto grad_lnfy_internal4(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e7 = e6 * e4 ;
auto     e9 = e2 * e4/e5 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( (ub-lb),(e4-1) ) ;
auto     e14   =  pow( e10,(e7-1) ) ;
auto     e15   =  pow( e11,(e9-1) ) ;
auto     e16 = tgamma(e7) ;
auto     e17 = tgamma(e9) ;
auto     e19 = e16 * e17 * e13 ;
return  ((log(e11)*e14*e15-log(e10)*e14*e15)/e19-(digamma(e9)-digamma(e7))*e16*e17*e13*e14*e15/ pow( e19,2 ) )*e6*e2*e4*tgamma(1*e4)*matrixmult1d(x,c({0,0,0,1,0,0,0}))/(e5*fy_internal(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)) ;
}

 
// [[Rcpp::export]]
auto grad_lnfy_internal5(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e7 = e6 * e4 ;
auto     e9 = e2 * e4/e5 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( (ub-lb),(e4-1) ) ;
auto     e14   =  pow( e10,(e7-1) ) ;
auto     e15   =  pow( e11,(e9-1) ) ;
auto     e16 = tgamma(e7) ;
auto     e17 = tgamma(e9) ;
auto     e19 = e16 * e17 * e13 ;
return  ((log(e11)*e14*e15-log(e10)*e14*e15)/e19-(digamma(e9)-digamma(e7))*e16*e17*e13*e14*e15/ pow( e19,2 ) )*e6*e2*e4*tgamma(1*e4)*matrixmult1d(x,c({0,0,0,0,1,0,0}))/(e5*fy_internal(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)) ;
}

 
// [[Rcpp::export]]
auto grad_lnfy_internal6(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e7 = e6 * e4 ;
auto     e9 = e2 * e4/e5 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( (ub-lb),(e4-1) ) ;
auto     e14   =  pow( e10,(e7-1) ) ;
auto     e15   =  pow( e11,(e9-1) ) ;
auto     e16 = tgamma(e7) ;
auto     e17 = tgamma(e9) ;
auto     e19 = e16 * e17 * e13 ;
return  ((log(e11)*e14*e15-log(e10)*e14*e15)/e19-(digamma(e9)-digamma(e7))*e16*e17*e13*e14*e15/ pow( e19,2 ) )*e6*e2*e4*tgamma(1*e4)*matrixmult1d(x,c({0,0,0,0,0,1,0}))/(e5*fy_internal(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)) ;
}

 
// [[Rcpp::export]]
auto grad_lnfy_internal7(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e7 = e6 * e4 ;
auto     e9 = e2 * e4/e5 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( (ub-lb),(e4-1) ) ;
auto     e14   =  pow( e10,(e7-1) ) ;
auto     e15   =  pow( e11,(e9-1) ) ;
auto     e16 = tgamma(e7) ;
auto     e17 = tgamma(e9) ;
auto     e19 = e16 * e17 * e13 ;
return  ((log(e11)*e14*e15-log(e10)*e14*e15)/e19-(digamma(e9)-digamma(e7))*e16*e17*e13*e14*e15/ pow( e19,2 ) )*e6*e2*e4*tgamma(1*e4)*matrixmult1d(x,c({0,0,0,0,0,0,1}))/(e5*fy_internal(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)) ;
}

 
// [[Rcpp::export]]
auto grad_lnfy_internal8(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e7 = e6 * e4 ;
auto     e9 = e2 * e4/e5 ;
auto     e10 = ub - lb ;
auto     e11 = ub - y ;
auto     e13   =  pow( e10,(e4-1) ) ;
auto     e14   =  pow( e11,(e7-1) ) ;
auto     e15 = y - lb ;
auto     e16   =  pow( e15,(e9-1) ) ;
auto     e17 = tgamma(e7) ;
auto     e18 = tgamma(e9) ;
auto     e20 = e17 * e18 * e13 ;
return  z*(((e6*log(e11)*e14+digamma(e4)*e14)*e16+e2*log(e15)*e14*e16/e5)/e20-((e6*digamma(e7)+digamma(e9)*e2/e5)*e13+log(e10)*e13)*e17*e18*e14*e16/ pow( e20,2 ) )*e4*tgamma(1*e4)/fy_internal(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8) ;
}


// [[Rcpp::export]]
std::vector<double> grad_lnfy_exec(double y,vector<double> x,vector<double> z,double lb,double ub, double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8 ){

//Note: we cannot use arrays of function pointers here because of the strange behavior of Deriv::Deriv(), which sometimes return vectors instead of scalars.
std::list<std::vector<double>> grad_eval; //temporary container to store gradient elements
std::vector<double> out_grad_eval; //current container for the output

//because the function based on Deriv::Deriv() returns sometimes a vector instead of a scalar! So safely store the output
vector<double> xout1 = {grad_lnfy_internal1(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)}; grad_eval.push_back(xout1);
vector<double> xout2 = {grad_lnfy_internal2(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)}; grad_eval.push_back(xout2);
vector<double> xout3 = {grad_lnfy_internal3(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)}; grad_eval.push_back(xout3);
vector<double> xout4 = {grad_lnfy_internal4(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)}; grad_eval.push_back(xout4);
vector<double> xout5 = {grad_lnfy_internal5(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)}; grad_eval.push_back(xout5);
vector<double> xout6 = {grad_lnfy_internal6(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)}; grad_eval.push_back(xout6);
vector<double> xout7 = {grad_lnfy_internal7(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)}; grad_eval.push_back(xout7);
vector<double> xout8 = {grad_lnfy_internal8(y,x,z,lb,ub,b1,b2,b3,b4,b5,b6,b7,b8)}; grad_eval.push_back(xout8);

for(const auto& vec:grad_eval) {out_grad_eval.insert(out_grad_eval.end(), vec.begin(), vec.end());}

return out_grad_eval;

}


// [[Rcpp::export]]
std::vector<double> grad_lnfy_exec_loop(vector<double>& y,vector<vector<double>>& X,vector<vector<double>>& Z,double lb,double ub, double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8 ){

int n = X.size();
int JH = X[0].size() + Z[0].size();
std::vector<std::vector<double>> Out_grad_eval(n,vector<double>(JH));
std::vector<double> out_grad_eval(JH);
for (int i = 0; i < n; i++) {
Out_grad_eval[i] = grad_lnfy_exec(y[i], X[i], Z[i],lb,ub,
b1,b2,b3,b4,b5,b6,b7,b8
);
for (int j = 0; j < JH; j++) {
out_grad_eval[j] += Out_grad_eval[i][j];
}
}

return out_grad_eval;

}

 
// [[Rcpp::export]]
auto hess_lnfy_internal1(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32* pow( (matrixmult1d(x,c({1,0,0,0,0,0,0}))),2 ) /e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal2(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,1,0,0,0,0,0}))*matrixmult1d(x,c({1,0,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal3(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,1,0,0,0,0}))*matrixmult1d(x,c({1,0,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal4(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,1,0,0,0}))*matrixmult1d(x,c({1,0,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal5(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,1,0,0}))*matrixmult1d(x,c({1,0,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal6(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,1,0}))*matrixmult1d(x,c({1,0,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal7(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,0,1}))*matrixmult1d(x,c({1,0,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal8(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e3 = 1 + e2 ;
auto     e5 = exp(matrixmult1d(z, c({b8}))) ;
auto     e6 = 1 - e2/e3 ;
auto     e7 = e6 * e5 ;
auto     e9 = e2 * e5/e3 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( e10,(e7-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e7) ;
auto     e16 = digamma(e9) ;
auto     e17 = log(e10) ;
auto     e18 = ub - lb ;
auto     e19 = log(e11) ;
auto     e20   =  pow( e18,(e5-1) ) ;
auto     e21 = e16 - e15 ;
auto     e22 = digamma(e5) ;
auto     e23 = tgamma(e7) ;
auto     e24 = tgamma(e9) ;
auto     e25 = log(e18) ;
auto     e28 = e6 * e15 + e16 * e2/e3 + e25 ;
auto     e29 = e6 * e17 ;
auto     e30 = e28 * e5 ;
auto     e31 = e29 + e22 ;
auto     e32 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e34 = e23 * e24 * e20 ;
auto     e35 = tgamma(e5) ;
return  z*((((e31*e19*e14-((e6*(e17-e15)+e22+e2*(e19-e16)/e3-e25)*((e19*e13*e14-e17*e13*e14)/e34-e21*e23*e24*e20*e13*e14/ pow( e34,2 ) )*e35/e32+(e31+e2*e19/e3)*e21)*e14)*e13-e28*(e19*e14-2*(e21*e14))*e13)*e5-((((e30+1)*e21+(e2*trigamma(e9)/e3-e6*trigamma(e7))*e5)*e13-e30*e17*e13)*e14+((e29*e13+e22*e13)*e5+e13)*e17*e14))*e6+((e7*e19*e14-e14)*e2*e13/e3+(e13-e6*e2*e5*e17*e13/e3)*e14)*e19)*e2*e5*e35*matrixmult1d(x,c({1,0,0,0,0,0,0}))/(e3*e32*e23*e24*e20) ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal9(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32* pow( (matrixmult1d(x,c({0,1,0,0,0,0,0}))),2 ) /e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal10(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,1,0,0,0,0}))*matrixmult1d(x,c({0,1,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal11(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,1,0,0,0}))*matrixmult1d(x,c({0,1,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal12(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,1,0,0}))*matrixmult1d(x,c({0,1,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal13(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,1,0}))*matrixmult1d(x,c({0,1,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal14(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,0,1}))*matrixmult1d(x,c({0,1,0,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal15(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e3 = 1 + e2 ;
auto     e5 = exp(matrixmult1d(z, c({b8}))) ;
auto     e6 = 1 - e2/e3 ;
auto     e7 = e6 * e5 ;
auto     e9 = e2 * e5/e3 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( e10,(e7-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e7) ;
auto     e16 = digamma(e9) ;
auto     e17 = log(e10) ;
auto     e18 = ub - lb ;
auto     e19 = log(e11) ;
auto     e20   =  pow( e18,(e5-1) ) ;
auto     e21 = e16 - e15 ;
auto     e22 = digamma(e5) ;
auto     e23 = tgamma(e7) ;
auto     e24 = tgamma(e9) ;
auto     e25 = log(e18) ;
auto     e28 = e6 * e15 + e16 * e2/e3 + e25 ;
auto     e29 = e6 * e17 ;
auto     e30 = e28 * e5 ;
auto     e31 = e29 + e22 ;
auto     e32 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e34 = e23 * e24 * e20 ;
auto     e35 = tgamma(e5) ;
return  z*((((e31*e19*e14-((e6*(e17-e15)+e22+e2*(e19-e16)/e3-e25)*((e19*e13*e14-e17*e13*e14)/e34-e21*e23*e24*e20*e13*e14/ pow( e34,2 ) )*e35/e32+(e31+e2*e19/e3)*e21)*e14)*e13-e28*(e19*e14-2*(e21*e14))*e13)*e5-((((e30+1)*e21+(e2*trigamma(e9)/e3-e6*trigamma(e7))*e5)*e13-e30*e17*e13)*e14+((e29*e13+e22*e13)*e5+e13)*e17*e14))*e6+((e7*e19*e14-e14)*e2*e13/e3+(e13-e6*e2*e5*e17*e13/e3)*e14)*e19)*e2*e5*e35*matrixmult1d(x,c({0,1,0,0,0,0,0}))/(e3*e32*e23*e24*e20) ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal16(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32* pow( (matrixmult1d(x,c({0,0,1,0,0,0,0}))),2 ) /e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal17(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,1,0,0,0}))*matrixmult1d(x,c({0,0,1,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal18(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,1,0,0}))*matrixmult1d(x,c({0,0,1,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal19(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,1,0}))*matrixmult1d(x,c({0,0,1,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal20(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,0,1}))*matrixmult1d(x,c({0,0,1,0,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal21(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e3 = 1 + e2 ;
auto     e5 = exp(matrixmult1d(z, c({b8}))) ;
auto     e6 = 1 - e2/e3 ;
auto     e7 = e6 * e5 ;
auto     e9 = e2 * e5/e3 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( e10,(e7-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e7) ;
auto     e16 = digamma(e9) ;
auto     e17 = log(e10) ;
auto     e18 = ub - lb ;
auto     e19 = log(e11) ;
auto     e20   =  pow( e18,(e5-1) ) ;
auto     e21 = e16 - e15 ;
auto     e22 = digamma(e5) ;
auto     e23 = tgamma(e7) ;
auto     e24 = tgamma(e9) ;
auto     e25 = log(e18) ;
auto     e28 = e6 * e15 + e16 * e2/e3 + e25 ;
auto     e29 = e6 * e17 ;
auto     e30 = e28 * e5 ;
auto     e31 = e29 + e22 ;
auto     e32 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e34 = e23 * e24 * e20 ;
auto     e35 = tgamma(e5) ;
return  z*((((e31*e19*e14-((e6*(e17-e15)+e22+e2*(e19-e16)/e3-e25)*((e19*e13*e14-e17*e13*e14)/e34-e21*e23*e24*e20*e13*e14/ pow( e34,2 ) )*e35/e32+(e31+e2*e19/e3)*e21)*e14)*e13-e28*(e19*e14-2*(e21*e14))*e13)*e5-((((e30+1)*e21+(e2*trigamma(e9)/e3-e6*trigamma(e7))*e5)*e13-e30*e17*e13)*e14+((e29*e13+e22*e13)*e5+e13)*e17*e14))*e6+((e7*e19*e14-e14)*e2*e13/e3+(e13-e6*e2*e5*e17*e13/e3)*e14)*e19)*e2*e5*e35*matrixmult1d(x,c({0,0,1,0,0,0,0}))/(e3*e32*e23*e24*e20) ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal22(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32* pow( (matrixmult1d(x,c({0,0,0,1,0,0,0}))),2 ) /e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal23(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,1,0,0}))*matrixmult1d(x,c({0,0,0,1,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal24(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,1,0}))*matrixmult1d(x,c({0,0,0,1,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal25(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,0,1}))*matrixmult1d(x,c({0,0,0,1,0,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal26(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e3 = 1 + e2 ;
auto     e5 = exp(matrixmult1d(z, c({b8}))) ;
auto     e6 = 1 - e2/e3 ;
auto     e7 = e6 * e5 ;
auto     e9 = e2 * e5/e3 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( e10,(e7-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e7) ;
auto     e16 = digamma(e9) ;
auto     e17 = log(e10) ;
auto     e18 = ub - lb ;
auto     e19 = log(e11) ;
auto     e20   =  pow( e18,(e5-1) ) ;
auto     e21 = e16 - e15 ;
auto     e22 = digamma(e5) ;
auto     e23 = tgamma(e7) ;
auto     e24 = tgamma(e9) ;
auto     e25 = log(e18) ;
auto     e28 = e6 * e15 + e16 * e2/e3 + e25 ;
auto     e29 = e6 * e17 ;
auto     e30 = e28 * e5 ;
auto     e31 = e29 + e22 ;
auto     e32 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e34 = e23 * e24 * e20 ;
auto     e35 = tgamma(e5) ;
return  z*((((e31*e19*e14-((e6*(e17-e15)+e22+e2*(e19-e16)/e3-e25)*((e19*e13*e14-e17*e13*e14)/e34-e21*e23*e24*e20*e13*e14/ pow( e34,2 ) )*e35/e32+(e31+e2*e19/e3)*e21)*e14)*e13-e28*(e19*e14-2*(e21*e14))*e13)*e5-((((e30+1)*e21+(e2*trigamma(e9)/e3-e6*trigamma(e7))*e5)*e13-e30*e17*e13)*e14+((e29*e13+e22*e13)*e5+e13)*e17*e14))*e6+((e7*e19*e14-e14)*e2*e13/e3+(e13-e6*e2*e5*e17*e13/e3)*e14)*e19)*e2*e5*e35*matrixmult1d(x,c({0,0,0,1,0,0,0}))/(e3*e32*e23*e24*e20) ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal27(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32* pow( (matrixmult1d(x,c({0,0,0,0,1,0,0}))),2 ) /e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal28(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,1,0}))*matrixmult1d(x,c({0,0,0,0,1,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal29(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,0,1}))*matrixmult1d(x,c({0,0,0,0,1,0,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal30(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e3 = 1 + e2 ;
auto     e5 = exp(matrixmult1d(z, c({b8}))) ;
auto     e6 = 1 - e2/e3 ;
auto     e7 = e6 * e5 ;
auto     e9 = e2 * e5/e3 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( e10,(e7-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e7) ;
auto     e16 = digamma(e9) ;
auto     e17 = log(e10) ;
auto     e18 = ub - lb ;
auto     e19 = log(e11) ;
auto     e20   =  pow( e18,(e5-1) ) ;
auto     e21 = e16 - e15 ;
auto     e22 = digamma(e5) ;
auto     e23 = tgamma(e7) ;
auto     e24 = tgamma(e9) ;
auto     e25 = log(e18) ;
auto     e28 = e6 * e15 + e16 * e2/e3 + e25 ;
auto     e29 = e6 * e17 ;
auto     e30 = e28 * e5 ;
auto     e31 = e29 + e22 ;
auto     e32 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e34 = e23 * e24 * e20 ;
auto     e35 = tgamma(e5) ;
return  z*((((e31*e19*e14-((e6*(e17-e15)+e22+e2*(e19-e16)/e3-e25)*((e19*e13*e14-e17*e13*e14)/e34-e21*e23*e24*e20*e13*e14/ pow( e34,2 ) )*e35/e32+(e31+e2*e19/e3)*e21)*e14)*e13-e28*(e19*e14-2*(e21*e14))*e13)*e5-((((e30+1)*e21+(e2*trigamma(e9)/e3-e6*trigamma(e7))*e5)*e13-e30*e17*e13)*e14+((e29*e13+e22*e13)*e5+e13)*e17*e14))*e6+((e7*e19*e14-e14)*e2*e13/e3+(e13-e6*e2*e5*e17*e13/e3)*e14)*e19)*e2*e5*e35*matrixmult1d(x,c({0,0,0,0,1,0,0}))/(e3*e32*e23*e24*e20) ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal31(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32* pow( (matrixmult1d(x,c({0,0,0,0,0,1,0}))),2 ) /e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal32(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32*matrixmult1d(x,c({0,0,0,0,0,0,1}))*matrixmult1d(x,c({0,0,0,0,0,1,0}))/e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal33(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e3 = 1 + e2 ;
auto     e5 = exp(matrixmult1d(z, c({b8}))) ;
auto     e6 = 1 - e2/e3 ;
auto     e7 = e6 * e5 ;
auto     e9 = e2 * e5/e3 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( e10,(e7-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e7) ;
auto     e16 = digamma(e9) ;
auto     e17 = log(e10) ;
auto     e18 = ub - lb ;
auto     e19 = log(e11) ;
auto     e20   =  pow( e18,(e5-1) ) ;
auto     e21 = e16 - e15 ;
auto     e22 = digamma(e5) ;
auto     e23 = tgamma(e7) ;
auto     e24 = tgamma(e9) ;
auto     e25 = log(e18) ;
auto     e28 = e6 * e15 + e16 * e2/e3 + e25 ;
auto     e29 = e6 * e17 ;
auto     e30 = e28 * e5 ;
auto     e31 = e29 + e22 ;
auto     e32 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e34 = e23 * e24 * e20 ;
auto     e35 = tgamma(e5) ;
return  z*((((e31*e19*e14-((e6*(e17-e15)+e22+e2*(e19-e16)/e3-e25)*((e19*e13*e14-e17*e13*e14)/e34-e21*e23*e24*e20*e13*e14/ pow( e34,2 ) )*e35/e32+(e31+e2*e19/e3)*e21)*e14)*e13-e28*(e19*e14-2*(e21*e14))*e13)*e5-((((e30+1)*e21+(e2*trigamma(e9)/e3-e6*trigamma(e7))*e5)*e13-e30*e17*e13)*e14+((e29*e13+e22*e13)*e5+e13)*e17*e14))*e6+((e7*e19*e14-e14)*e2*e13/e3+(e13-e6*e2*e5*e17*e13/e3)*e14)*e19)*e2*e5*e35*matrixmult1d(x,c({0,0,0,0,0,1,0}))/(e3*e32*e23*e24*e20) ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal34(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e4 = exp(matrixmult1d(z, c({b8}))) ;
auto     e5 = 1 + e2 ;
auto     e6 = 1 - e2/e5 ;
auto     e8 = e2 * e4/e5 ;
auto     e9 = e6 * e4 ;
auto     e10 = y - lb ;
auto     e11 = ub - y ;
auto     e12   =  pow( e10,(e8-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e9) ;
auto     e16 = digamma(e8) ;
auto     e17 = log(e11) ;
auto     e18 = log(e10) ;
auto     e19 = e16 - e15 ;
auto     e20   =  pow( (ub-lb),(e4-1) ) ;
auto     e21 = tgamma(e9) ;
auto     e22 = tgamma(e8) ;
auto     e24 = e15 + e18 - (e16 + e17) ;
auto     e25 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e27 = e21 * e22 * e20 ;
auto     e28 = e5 * e25 ;
auto     e30 = e24 * e14 * e12 ;
auto     e31 = e19 * e12 ;
auto     e32 = tgamma(e4) ;
auto     e34 = e17 * e14 * e12 ;
auto     e35 = e18 * e12 ;
return  ((((((e35-e31)*e14-e34)*(e18-e17)-((( pow( e19,2 ) +trigamma(e9)+trigamma(e8))*e14-e19*e17*e14)*e12+e19*(e35-2*e31)*e14))*e6*e4-e30)*e2/e5+e30)/e28-(((e18*e14*e12-e34)/e27-e19*e21*e22*e20*e14*e12/ pow( e27,2 ) )*e6*e4*e32+e25)*e24*e2*e14*e12/ pow( e28,2 ) )*e6*e2*e4*e32* pow( (matrixmult1d(x,c({0,0,0,0,0,0,1}))),2 ) /e27 ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal35(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e3 = 1 + e2 ;
auto     e5 = exp(matrixmult1d(z, c({b8}))) ;
auto     e6 = 1 - e2/e3 ;
auto     e7 = e6 * e5 ;
auto     e9 = e2 * e5/e3 ;
auto     e10 = ub - y ;
auto     e11 = y - lb ;
auto     e13   =  pow( e10,(e7-1) ) ;
auto     e14   =  pow( e11,(e9-1) ) ;
auto     e15 = digamma(e7) ;
auto     e16 = digamma(e9) ;
auto     e17 = log(e10) ;
auto     e18 = ub - lb ;
auto     e19 = log(e11) ;
auto     e20   =  pow( e18,(e5-1) ) ;
auto     e21 = e16 - e15 ;
auto     e22 = digamma(e5) ;
auto     e23 = tgamma(e7) ;
auto     e24 = tgamma(e9) ;
auto     e25 = log(e18) ;
auto     e28 = e6 * e15 + e16 * e2/e3 + e25 ;
auto     e29 = e6 * e17 ;
auto     e30 = e28 * e5 ;
auto     e31 = e29 + e22 ;
auto     e32 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e34 = e23 * e24 * e20 ;
auto     e35 = tgamma(e5) ;
return  z*((((e31*e19*e14-((e6*(e17-e15)+e22+e2*(e19-e16)/e3-e25)*((e19*e13*e14-e17*e13*e14)/e34-e21*e23*e24*e20*e13*e14/ pow( e34,2 ) )*e35/e32+(e31+e2*e19/e3)*e21)*e14)*e13-e28*(e19*e14-2*(e21*e14))*e13)*e5-((((e30+1)*e21+(e2*trigamma(e9)/e3-e6*trigamma(e7))*e5)*e13-e30*e17*e13)*e14+((e29*e13+e22*e13)*e5+e13)*e17*e14))*e6+((e7*e19*e14-e14)*e2*e13/e3+(e13-e6*e2*e5*e17*e13/e3)*e14)*e19)*e2*e5*e35*matrixmult1d(x,c({0,0,0,0,0,0,1}))/(e3*e32*e23*e24*e20) ;
}

 
// [[Rcpp::export]]
auto hess_lnfy_internal36(double y,vector<double> x,vector<double> z,double lb,double ub,double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8)
{
auto     e2 = exp(matrixmult1d(x, c({b1, b2, b3, b4, b5, b6, b7}))) ;
auto     e3 = 1 + e2 ;
auto     e5 = exp(matrixmult1d(z, c({b8}))) ;
auto     e6 = 1 - e2/e3 ;
auto     e7 = e6 * e5 ;
auto     e8 = ub - lb ;
auto     e10 = e2 * e5/e3 ;
auto     e11 = ub - y ;
auto     e12 = y - lb ;
auto     e13   =  pow( e8,(e5-1) ) ;
auto     e15   =  pow( e11,(e7-1) ) ;
auto     e16   =  pow( e12,(e10-1) ) ;
auto     e17 = digamma(e7) ;
auto     e18 = digamma(e10) ;
auto     e19 = log(e8) ;
auto     e20 = log(e11) ;
auto     e21 = e6 * e17 ;
auto     e22 = e18 * e2 ;
auto     e23 = log(e12) ;
auto     e24 = e21 + e22/e3 ;
auto     e25 = digamma(e5) ;
auto     e26 = e6 * e20 ;
auto     e27 = e2 * e23 ;
auto     e28 = tgamma(e7) ;
auto     e29 = tgamma(e10) ;
auto     e30 = e19 * e13 ;
auto     e32 = e24 * e13 + e30 ;
auto     e33 = e24 + e19 ;
auto     e34 = e26 + e25 ;
auto     e35 = e26 * e15 ;
auto     e38 = e27 * e15 * e16/e3 ;
auto     e39 = fy_internal(y, x, z, lb, ub, b1, b2, b3, b4, b5, b6, b7, b8) ;
auto     e41 = e28 * e29 * e13 ;
auto     e42 = tgamma(e5) ;
return  pow( z,2 ) *(((e34*e6*e20*e15+trigamma(e5)*e15)*e16+(e34*e2*e23*e16/e3-e32*(e34+e27/e3)*e16/e13)*e15+(e35*e16+e38)*e2*e23/e3-((((e21*e13+e22*e13/e3+e30)*e33+( pow( e6,2 ) *trigamma(e7)+ pow( e2,2 ) *trigamma(e10)/ pow( e3,2 ) )*e13)*e15+e33*e6*e20*e13*e15)*e16+e33*(e27*e16/e3-2*(e32*e16/e13))*e13*e15)/e13)*e5+(e6*(e20-e17)+e25+e2*(e23-e18)/e3-e19)*((e25-(((e35+e25*e15)*e16+e38)/e41-e32*e28*e29*e15*e16/ pow( e41,2 ) )*e42/e39)*e5+1)*e15*e16)*e5*e42/(e39*e28*e29*e13) ;
}


// [[Rcpp::export]]
std::vector<double> hess_lnfy_exec(double y,vector<double> x,vector<double> z,double lb,double ub, double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8 ){

//Note: we cannot use arrays of function pointers here because of the strange behavior of Deriv::Deriv(), which sometimes return vectors instead of scalars.
std::list<std::vector<double>> hess_eval; //temporary container to store gradient elements
std::vector<double> out_hess_eval; //current container for the output

//because the function based on Deriv::Deriv() returns sometimes a vector instead of a scalar! So safely store the output
vector<double> xout1 = {hess_lnfy_internal1(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout1);
vector<double> xout2 = {hess_lnfy_internal2(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout2);
vector<double> xout3 = {hess_lnfy_internal3(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout3);
vector<double> xout4 = {hess_lnfy_internal4(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout4);
vector<double> xout5 = {hess_lnfy_internal5(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout5);
vector<double> xout6 = {hess_lnfy_internal6(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout6);
vector<double> xout7 = {hess_lnfy_internal7(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout7);
vector<double> xout8 = {hess_lnfy_internal8(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout8);
vector<double> xout9 = {hess_lnfy_internal9(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout9);
vector<double> xout10 = {hess_lnfy_internal10(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout10);
vector<double> xout11 = {hess_lnfy_internal11(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout11);
vector<double> xout12 = {hess_lnfy_internal12(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout12);
vector<double> xout13 = {hess_lnfy_internal13(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout13);
vector<double> xout14 = {hess_lnfy_internal14(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout14);
vector<double> xout15 = {hess_lnfy_internal15(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout15);
vector<double> xout16 = {hess_lnfy_internal16(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout16);
vector<double> xout17 = {hess_lnfy_internal17(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout17);
vector<double> xout18 = {hess_lnfy_internal18(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout18);
vector<double> xout19 = {hess_lnfy_internal19(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout19);
vector<double> xout20 = {hess_lnfy_internal20(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout20);
vector<double> xout21 = {hess_lnfy_internal21(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout21);
vector<double> xout22 = {hess_lnfy_internal22(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout22);
vector<double> xout23 = {hess_lnfy_internal23(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout23);
vector<double> xout24 = {hess_lnfy_internal24(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout24);
vector<double> xout25 = {hess_lnfy_internal25(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout25);
vector<double> xout26 = {hess_lnfy_internal26(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout26);
vector<double> xout27 = {hess_lnfy_internal27(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout27);
vector<double> xout28 = {hess_lnfy_internal28(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout28);
vector<double> xout29 = {hess_lnfy_internal29(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout29);
vector<double> xout30 = {hess_lnfy_internal30(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout30);
vector<double> xout31 = {hess_lnfy_internal31(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout31);
vector<double> xout32 = {hess_lnfy_internal32(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout32);
vector<double> xout33 = {hess_lnfy_internal33(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout33);
vector<double> xout34 = {hess_lnfy_internal34(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout34);
vector<double> xout35 = {hess_lnfy_internal35(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout35);
vector<double> xout36 = {hess_lnfy_internal36(y,x,z,lb,ub, b1,b2,b3,b4,b5,b6,b7,b8)}; hess_eval.push_back(xout36);

for(const auto& vec:hess_eval) {out_hess_eval.insert(out_hess_eval.end(), vec.begin(), vec.end());}

return out_hess_eval;

}


// [[Rcpp::export]]
std::vector<double> hess_lnfy_exec_loop(vector<double>& y,vector<vector<double>>& X,vector<vector<double>>& Z,double lb,double ub, double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8 ){

int n = X.size();
int JH = X[0].size() + Z[0].size();

//std::vector<int> l(JH);
//std::iota(l.begin(), l.end(), 1);
//vector<vector<int>> IJ = expand_grid(l,l,true);
//int K = IJ.size();
int K = (JH*(1+JH))/2;

std::vector<std::vector<double>> Out_hess_eval(n,vector<double>(K));
std::vector<double> out_hess_eval(K);
for (int i = 0; i < n; i++) {
Out_hess_eval[i] = hess_lnfy_exec(y[i], X[i], Z[i],lb,ub,b1,b2,b3,b4,b5,b6,b7,b8);
for (int j = 0; j < K; j++) {
out_hess_eval[j] += Out_hess_eval[i][j];
}
}

return out_hess_eval;

}
