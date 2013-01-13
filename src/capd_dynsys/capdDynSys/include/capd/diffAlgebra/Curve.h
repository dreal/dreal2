/////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2000-2005 by the CAPD Group.
#ifndef _CAPD_DIFFALGEBRA_CURVE_H_
#include <stdexcept>

namespace capd{
namespace diffAlgebra{

class Curve
{
public:
  Curve(Real left, Real right, int dimension, int order);
  Curve& operator=(const Curve& c);
  void setOrder(int order);
  const MatrixType* getMatrixCoefficients() const;
   VectorType* getCoefficientsAtCenter();
   MatrixType* getMatrixCoefficients();
   void clearCoefficients();
protected:
  VectorType *m_coefficientsAtCenter;
  MatrixType *m_matrixCoefficients;
  int m_order;

// ----------------- inline definitions ------------------

template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
template<class MatrixT>
}} // namespace capd::diffAlgebra