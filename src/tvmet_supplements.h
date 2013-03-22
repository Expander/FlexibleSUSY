#ifndef tvmet_supplements_h
#define tvmet_supplements_h


namespace tvmet {


// Shorthand notation of built-in trace

template<class T, std::size_t Sz>
inline T
Tr(const Matrix<T, Sz, Sz>& m) {
  return meta::Matrix<Sz, Sz, 0, 0>::trace(m);
}

template<class E, std::size_t Sz>
inline typename E::value_type
Tr(const XprMatrix<E, Sz, Sz>& m) {
  return meta::Matrix<Sz, Sz, 0, 0>::trace(m);
}


// Shorthand notation of built-in transpose

template<class U, std::size_t Rows, std::size_t Cols>
inline
XprMatrix<
  XprMatrixTranspose<
    MatrixConstReference<U, Rows, Cols>
  >,
  Cols, Rows
>
T(const Matrix<U, Rows, Cols>& rhs) {
  typedef XprMatrixTranspose<
    MatrixConstReference<U, Rows, Cols>
  >							expr_type;
  return XprMatrix<expr_type, Cols, Rows>(
    expr_type(rhs.const_ref()));
}

template<class E, std::size_t Rows, std::size_t Cols>
inline
XprMatrix<
  XprMatrixTranspose<
    XprMatrix<E, Rows, Cols>
  >,
  Cols, Rows
>
T(const XprMatrix<E, Rows, Cols>& rhs) {
  typedef XprMatrixTranspose<
    XprMatrix<E, Rows, Cols>
  >							expr_type;
  return XprMatrix<expr_type, Cols, Rows>(expr_type(rhs));
}


#if defined(TVMET_HAVE_COMPLEX)


// Complex conjugate

template<class E, std::size_t Rows, std::size_t Cols>
inline
XprMatrix<
  XprUnOp<
    Fcnl_conj<typename E::value_type>,
    XprMatrix<E, Rows, Cols>
  >,
  Rows, Cols
>
conj(const XprMatrix<E, Rows, Cols>& rhs) {
  typedef XprUnOp<
    Fcnl_conj<typename E::value_type>,
    XprMatrix<E, Rows, Cols>
  > expr_type;
  return XprMatrix<expr_type, Rows, Cols>(expr_type(rhs));
}

// Hermitian conjugate

template<class T, std::size_t Rows, std::size_t Cols>
inline
XprMatrix<
  XprMatrixTranspose<
    XprUnOp<
      Fcnl_conj<std::complex<T> >,
      MatrixConstReference<std::complex<T>, Rows, Cols>
    >
  >,
  Cols, Rows
>
H(const Matrix<std::complex<T>, Rows, Cols>& rhs) {
  typedef XprUnOp<
    Fcnl_conj<std::complex<T> >,
    MatrixConstReference<std::complex<T>, Rows, Cols>
  > expr_type1;
  typedef XprMatrixTranspose<
    expr_type1
  > expr_type2;
  return XprMatrix<expr_type2, Cols, Rows>(
    expr_type2(expr_type1(rhs.const_ref())));
}

template<class E, std::size_t Rows, std::size_t Cols>
inline
XprMatrix<
  XprMatrixTranspose<
    XprUnOp<
      Fcnl_conj<typename E::value_type>,
      XprMatrix<E, Rows, Cols>
    >
  >,
  Cols, Rows
>
H(const XprMatrix<E, Rows, Cols>& rhs) {
  typedef XprUnOp<
    Fcnl_conj<typename E::value_type>,
    XprMatrix<E, Rows, Cols>
  > expr_type1;
  typedef XprMatrixTranspose<
    expr_type1
  > expr_type2;
  return XprMatrix<expr_type2, Cols, Rows>(expr_type2(expr_type1(rhs)));
}


#endif // defined(TVMET_HAVE_COMPLEX)


} // namespace tvmet


#endif // tvmet_supplements_h


// Local Variables:
// mode:C++
// End:
