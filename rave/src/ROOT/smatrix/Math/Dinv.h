// @(#)root/smatrix:$Id: Dinv.h 26866 2008-12-12 10:50:07Z moneta $
// Authors: T. Glebe, L. Moneta    2005  

#ifndef  ROOT_Math_Dinv
#define  ROOT_Math_Dinv
// ********************************************************************
//
// source:
//
// type:      source code
//
// created:   03. Apr 2001
//
// author:    Thorsten Glebe
//            HERA-B Collaboration
//            Max-Planck-Institut fuer Kernphysik
//            Saupfercheckweg 1
//            69117 Heidelberg
//            Germany
//            E-mail: T.Glebe@mpi-hd.mpg.de
//
// Description: square Matrix inversion
//              Code was taken from CERNLIB::kernlib dfinv function, translated
//              from FORTRAN to C++ and optimized.
//              n:    Order of the square matrix
//              idim: First dimension of array A
//
// changes:
// 03 Apr 2001 (TG) creation
//
// ********************************************************************
#ifdef OLD_IMPL
#include "Math/Dfactir.h"
#include "Math/Dfinv.h"
#include "Math/Dsinv.h"
#endif

 
namespace ROOT { 

  namespace Math { 



/** 
    Matrix Inverter class 
    Class to specialize calls to Dinv. Dinv computes the inverse of a square
    matrix if dimension idim and order n. The content of the matrix will be
    replaced by its inverse. In case the inversion fails, the matrix content is
    destroyed. Invert specializes Dinv by the matrix order. E.g. if the order
    of the matrix is two, the routine Inverter<2> is called which implements
    Cramers rule.

    @author T. Glebe
*/
//==============================================================================
// Inverter class
//==============================================================================
template <unsigned int idim, unsigned int n = idim>
class Inverter {
public:
  ///
  template <class MatrixRep>
  static bool Dinv(MatrixRep& rhs) {

#ifdef XXX
      if (n < 1 || n > idim) {
	return false;
      }
#endif

#ifdef OLD_IMPL

      /* Initialized data */
      static unsigned int work[n];
      for(unsigned int i=0; i<n; ++i) work[i] = 0;

      static typename MatrixRep::value_type det = 0;

      /* Function Body */
      
      /*  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT. */
      if (Dfactir<MatrixRep,n,idim>(rhs,det,work) == false) {
	std::cerr << "Dfactir failed!!" << std::endl;
	return false;
      }
      return Dfinv<MatrixRep,n,idim>(rhs,work);
#else 

      /* Initialized data */
      static unsigned int work[n+1];
      for(unsigned int i=0; i<n+1; ++i) work[i] = 0;

      static typename MatrixRep::value_type det = 0;
      
      if (DfactMatrix(rhs,det,work) != 0) {
	std::cerr << "Dfact_matrix failed!!" << std::endl;
	return false;
      }

      int ifail =  DfinvMatrix(rhs,work); 
      if (ifail == 0) return true; 
      return false; 
#endif
  } // Dinv


  // symmetric function (copy in a general one) 
  template <class T>
  static bool Dinv(MatRepSym<T,idim> & rhs) {
    // not very efficient but need to re-do Dsinv for new storage of 
    // symmetric matrices
#ifdef OLD_IMPL
    MatRepStd<T,idim>  tmp; 
    for (unsigned int i = 0; i< idim*idim; ++i) 
      tmp[i] = rhs[i];
    if (! Inverter<idim>::Dinv(tmp) ) return false;
    // recopy the data
    for (unsigned int i = 0; i< idim*n; ++i) 
      rhs[i] = tmp[i];

    return true; 
#else
    int ifail = 0; 
    InvertBunchKaufman(rhs,ifail); 
    if (ifail == 0) return true; 
    return false; 
#endif
  }


  /**
     Bunch-Kaufman method for inversion of symmetric matrices
   */
  template <class T>
  static int DfactMatrix(MatRepStd<T,idim,n> & rhs, T & det, unsigned int * work); 
  template <class T>
  static int DfinvMatrix(MatRepStd<T,idim,n> & rhs, unsigned int * work); 

  /**
     Bunch-Kaufman method for inversion of symmetric matrices
   */
  template <class T>
  static void InvertBunchKaufman(MatRepSym<T,idim> & rhs, int &ifail); 



}; // class Inverter

// fast inverter class using Cramer inversion 
// by default use other default inversion
/** 
    Fast Matrix Inverter class 
    Class to specialize calls to Dinv. Dinv computes the inverse of a square
    matrix if dimension idim and order n. The content of the matrix will be
    replaced by its inverse. In case the inversion fails, the matrix content is
    destroyed. Invert specializes Dinv by the matrix order. E.g. if the order
    of the matrix is less than 5 , the class implements
    Cramers rule. 
    Be careful that for matrix with high condition the accuracy of the Cramer rule is much poorer

    @author L. Moneta
*/
template <unsigned int idim, unsigned int n = idim>
class FastInverter {
public:
  ///
  template <class MatrixRep>
  static bool Dinv(MatrixRep& rhs) {
     return Inverter<idim,n>::Dinv(rhs); 
  }
  template <class T>
  static bool Dinv(MatRepSym<T,idim> & rhs) {
     return Inverter<idim,n>::Dinv(rhs); 
  }
};


/** Inverter<0>.
    In case of zero order, do nothing.

    @author T. Glebe
*/
//==============================================================================
// Inverter<0>
//==============================================================================
template <>
class Inverter<0> {
public:
  ///
  template <class MatrixRep>
  inline static bool Dinv(MatrixRep& rhs) { return true; }
};
    

/** 
    1x1 matrix inversion \f$a_{11} \to 1/a_{11}\f$

    @author T. Glebe
*/
//==============================================================================
// Inverter<1>
//==============================================================================
template <>
class Inverter<1> {
public:
  ///
  template <class MatrixRep>
  static bool Dinv(MatrixRep& rhs) {
    
    if (rhs[0] == 0.) {
      return false;
    }
    rhs[0] = 1. / rhs[0];
    return true;
  }
};


/** 
    2x2 matrix inversion  using Cramers rule.

    @author T. Glebe
*/
//==============================================================================
// Inverter<2>: Cramers rule
//==============================================================================

template <>
class Inverter<2> {
public:
  ///
  template <class MatrixRep>
  static bool Dinv(MatrixRep& rhs) {

    typename MatrixRep::value_type det = rhs[0] * rhs[3] - rhs[2] * rhs[1];
    
    if (det == 0.) { return false; }

    // use 1.0f to remove warning on Windows
    typename MatrixRep::value_type s = 1.0f / det; 

    typename MatrixRep::value_type c11 = s * rhs[3];


    rhs[2] = -s * rhs[2];
    rhs[1] = -s * rhs[1];
    rhs[3] =  s * rhs[0];
    rhs[0] = c11;


    return true;
  }

  // specialization for the symmetric matrices
  template <class T>
  static bool Dinv(MatRepSym<T,2> & rep) {
    
    T * rhs = rep.Array(); 

    T det = rhs[0] * rhs[2] - rhs[1] * rhs[1];

    
    if (det == 0.) { return false; }

    T s = 1.0f / det;
    T c11 = s * rhs[2];

    rhs[1] = -s * rhs[1];
    rhs[2] =  s * rhs[0];
    rhs[0] = c11;
    return true;
  }

};


/** 
    3x3 direct matrix inversion  suing Cramer Rule
    use only for FastInverter
*/
//==============================================================================
// FastInverter<3>
//==============================================================================

template <>
class FastInverter<3> {
public:
  ///
  // use Cramer Rule
  template <class MatrixRep>
  static bool Dinv(MatrixRep& rhs); 

  template <class T>
  static bool Dinv(MatRepSym<T,3> & rhs);

};

/** 
    4x4 matrix inversion using Cramers rule.
*/
template <>
class FastInverter<4> {
public:
  ///
  template <class MatrixRep>
  static bool Dinv(MatrixRep& rhs); 

  template <class T>
  static bool Dinv(MatRepSym<T,4> & rhs);

};

/** 
    5x5 Matrix inversion using Cramers rule.
*/
template <>
class FastInverter<5> {
public:
  ///
  template <class MatrixRep>
  static bool Dinv(MatrixRep& rhs); 

  template <class T>
  static bool Dinv(MatRepSym<T,5> & rhs);

};



  }  // namespace Math

}  // namespace ROOT
          
#ifndef ROOT_Math_CramerInversion_icc
#include "CramerInversion.icc"
#endif
#ifndef ROOT_Math_CramerInversionSym_icc
#include "CramerInversionSym.icc"
#endif
#ifndef ROOT_Math_MatrixInversion_icc
#include "MatrixInversion.icc"
#endif

#endif  /* ROOT_Math_Dinv */
