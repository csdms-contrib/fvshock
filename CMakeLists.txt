cmake_minimum_required (VERSION 2.6)

enable_language (Fortran)

set (MODEL fvshock)
set (BUILD_SHARED_LIBS ON)

#add_subdirectory (lapack)

set (fvshock_srcs
     CFLCON.f
     LEGGI.f
     RECONSTRUCTION.f
     SUBROUT.f
     CONTORNO.f
     MainPRICE2D.f
     REYNOLDS.f
     UPDATE1.f
     INIZIO.f
     NODALIfisseECC.f
     STAMPA.f
     coeff.f
)

set (lapack_srcs
     lapack/daxpy.f
     lapack/dcopy.f
     lapack/ddot.f
     lapack/dgebak.f
     lapack/dgebal.f
     lapack/dgeev.f
     lapack/dgehd2.f
     lapack/dgehrd.f
     lapack/dgemm.f
     lapack/dgemv.f
     lapack/dger.f
     lapack/dgetf2.f
     lapack/dgetrf.f
     lapack/dgetri.f
     lapack/dhseqr.f
     lapack/dlabad.f
     lapack/dlacpy.f
     lapack/dladiv.f
     lapack/dlaexc.f
     lapack/dlahqr.f
     lapack/dlahr2.f
     lapack/dlaln2.f
     lapack/dlamch.f
     lapack/dlange.f
     lapack/dlanv2.f
     lapack/dlapy2.f
     lapack/dlaqr0.f
     lapack/dlaqr1.f
     lapack/dlaqr2.f
     lapack/dlaqr3.f
     lapack/dlaqr4.f
     lapack/dlaqr5.f
     lapack/dlarfb.f
     lapack/dlarf.f
     lapack/dlarfg.f
     lapack/dlarft.f
     lapack/dlarfx.f
     lapack/dlartg.f
     lapack/dlascl.f
     lapack/dlaset.f
     lapack/dlassq.f
     lapack/dlaswp.f
     lapack/dlasy2.f
     lapack/dnrm2.f
     lapack/dorg2r.f
     lapack/dorghr.f
     lapack/dorgqr.f
     lapack/drot.f
     lapack/dscal.f
     lapack/dswap.f
     lapack/dtrevc.f
     lapack/dtrexc.f
     lapack/dtrmm.f
     lapack/dtrmv.f
     lapack/dtrsm.f
     lapack/dtrti2.f
     lapack/dtrtri.f
     lapack/idamax.f
     lapack/ieeeck.f
     lapack/ilaenv.f
     lapack/iparmq.f
     lapack/lsame.f
     lapack/xerbla.f
)

add_executable (${MODEL} ${fvshock_srcs} ${lapack_srcs})

install (TARGETS ${MODEL} RUNTIME DESTINATION bin)
#install (TARGETS ${MODEL} DESTINATION lib)

