/*****************************************************************************
  Copyright (c) 2014, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************
* Contents: Native C interface to LAPACK
* Author: Intel Corporation
* Generated August, 2015
*****************************************************************************/

#ifndef _LAPACKE_H_
#define _LAPACKE_H_

/*
 *  Turn on HAVE_LAPACK_CONFIG_H to redefine C-LAPACK datatypes
 */
#ifdef HAVE_LAPACK_CONFIG_H
#include "lapacke_config.h"
#endif

#include <R_ext/Lapack.h>
#include <stdlib.h>

#ifndef lapack_int
#define lapack_int int
#endif

#ifndef lapack_logical
#define lapack_logical lapack_int
#endif

/* Complex types are structures equivalent to the
 * Fortran complex types COMPLEX(4) and COMPLEX(8).
 *
 * One can also redefine the types with his own types
 * for example by including in the code definitions like
 *
 * #define lapack_complex_float std::complex<float>
 * #define lapack_complex_double std::complex<double>
 *
 * or define these types in the command line:
 *
 * -Dlapack_complex_float="std::complex<float>"
 * -Dlapack_complex_double="std::complex<double>"
 */

#ifndef LAPACK_COMPLEX_CUSTOM

/* Complex type (single precision) */
#ifndef lapack_complex_float
#include <complex.h>
#define lapack_complex_float float _Complex
#endif

#ifndef lapack_complex_float_real
#define lapack_complex_float_real(z) (creal(z))
#endif

#ifndef lapack_complex_float_imag
#define lapack_complex_float_imag(z) (cimag(z))
#endif

lapack_complex_float lapack_make_complex_float(float re, float im);

/* Complex type (double precision) */
#ifndef lapack_complex_double
#include <complex.h>
#define lapack_complex_double double _Complex
#endif

#ifndef lapack_complex_double_real
#define lapack_complex_double_real(z) (creal(z))
#endif

#ifndef lapack_complex_double_imag
#define lapack_complex_double_imag(z) (cimag(z))
#endif

lapack_complex_double lapack_make_complex_double(double re, double im);

#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef LAPACKE_malloc
#define LAPACKE_malloc(size) malloc(size)
#endif
#ifndef LAPACKE_free
#define LAPACKE_free(p) free(p)
#endif

#define LAPACK_C2INT(x) (lapack_int)(*((float*)&x))
#define LAPACK_Z2INT(x) (lapack_int)(*((double*)&x))

#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

#define LAPACK_WORK_MEMORY_ERROR -1010
#define LAPACK_TRANSPOSE_MEMORY_ERROR -1011

/* Callback logical functions of one, two, or three arguments are used
 *  to select eigenvalues to sort to the top left of the Schur form.
 *  The value is selected if function returns TRUE (non-zero). */

typedef lapack_logical (*LAPACK_S_SELECT2)(const float*, const float*);
typedef lapack_logical (*LAPACK_S_SELECT3)(const float*, const float*,
                                           const float*);
typedef lapack_logical (*LAPACK_D_SELECT2)(const double*, const double*);
typedef lapack_logical (*LAPACK_D_SELECT3)(const double*, const double*,
                                           const double*);

typedef lapack_logical (*LAPACK_S_SELECT2_NON_CONST)(float*, float*);
typedef lapack_logical (*LAPACK_S_SELECT3_NON_CONST)(float*, float*, float*);
typedef lapack_logical (*LAPACK_D_SELECT2_NON_CONST)(double*, double*);
typedef lapack_logical (*LAPACK_D_SELECT3_NON_CONST)(double*, double*, double*);

typedef lapack_logical (*LAPACK_C_SELECT1)(const lapack_complex_float*);
typedef lapack_logical (*LAPACK_C_SELECT2)(const lapack_complex_float*,
                                           const lapack_complex_float*);
typedef lapack_logical (*LAPACK_Z_SELECT1)(const lapack_complex_double*);
typedef lapack_logical (*LAPACK_Z_SELECT2)(const lapack_complex_double*,
                                           const lapack_complex_double*);

/*
#include "lapacke_mangling.h"

#define LAPACK_lsame LAPACK_GLOBAL(lsame,LSAME)
lapack_logical LAPACK_lsame( char* ca,  char* cb,
                              lapack_int lca, lapack_int lcb );
*/

#include "lapacke_utils.h"

/* APIs for set/get nancheck flags */
void LAPACKE_set_nancheck(int flag);
int LAPACKE_get_nancheck();

/* Selected Double Precision Lapack Routines
 ========
 */

//* Double precision BiDiagonal and DIagonal matrices  -> DBD & DDI

#include "lapacke_bits/lapacke_dbdsqr_work.h"
#include "lapacke_bits/lapacke_ddisna_work.h"

#include "lapacke_bits/lapacke_dbdsqr.h"
#include "lapacke_bits/lapacke_ddisna.h"

//* Double precision General Banded matrices -> DGB

#include "lapacke_bits/lapacke_dgbbrd_work.h"
#include "lapacke_bits/lapacke_dgbcon_work.h"
#include "lapacke_bits/lapacke_dgbequ_work.h"
#include "lapacke_bits/lapacke_dgbrfs_work.h"
#include "lapacke_bits/lapacke_dgbsv_work.h"
#include "lapacke_bits/lapacke_dgbsvx_work.h"
#include "lapacke_bits/lapacke_dgbtrf_work.h"
#include "lapacke_bits/lapacke_dgbtrs_work.h"

#include "lapacke_bits/lapacke_dgbbrd.h"
#include "lapacke_bits/lapacke_dgbcon.h"
#include "lapacke_bits/lapacke_dgbequ.h"
#include "lapacke_bits/lapacke_dgbrfs.h"
#include "lapacke_bits/lapacke_dgbsv.h"
#include "lapacke_bits/lapacke_dgbsvx.h"
#include "lapacke_bits/lapacke_dgbtrf.h"
#include "lapacke_bits/lapacke_dgbtrs.h"

//* Double precision GEneral matrices -> DGE

#include "lapacke_bits/lapacke_dgebak_work.h"
#include "lapacke_bits/lapacke_dgebal_work.h"
#include "lapacke_bits/lapacke_dgebrd_work.h"
#include "lapacke_bits/lapacke_dgecon_work.h"
#include "lapacke_bits/lapacke_dgeequ_work.h"
#include "lapacke_bits/lapacke_dgees_work.h"
#include "lapacke_bits/lapacke_dgeesx_work.h"
#include "lapacke_bits/lapacke_dgeev_work.h"
#include "lapacke_bits/lapacke_dgeevx_work.h"
#include "lapacke_bits/lapacke_dgehrd_work.h"
#include "lapacke_bits/lapacke_dgelq2_work.h"
#include "lapacke_bits/lapacke_dgelqf_work.h"
#include "lapacke_bits/lapacke_dgels_work.h"
#include "lapacke_bits/lapacke_dgelss_work.h"
#include "lapacke_bits/lapacke_dgelsy_work.h"
#include "lapacke_bits/lapacke_dgeqlf_work.h"
#include "lapacke_bits/lapacke_dgeqp3_work.h"
#include "lapacke_bits/lapacke_dgeqr2_work.h"
#include "lapacke_bits/lapacke_dgeqrf_work.h"
#include "lapacke_bits/lapacke_dgerfs_work.h"
#include "lapacke_bits/lapacke_dgerqf_work.h"
#include "lapacke_bits/lapacke_dgesv_work.h"
#include "lapacke_bits/lapacke_dgesvd_work.h"
#include "lapacke_bits/lapacke_dgesvx_work.h"
#include "lapacke_bits/lapacke_dgetf2_work.h"
#include "lapacke_bits/lapacke_dgetrf_work.h"
#include "lapacke_bits/lapacke_dgetri_work.h"
#include "lapacke_bits/lapacke_dgetrs_work.h"

#include "lapacke_bits/lapacke_dgebak.h"
#include "lapacke_bits/lapacke_dgebal.h"
#include "lapacke_bits/lapacke_dgebrd.h"
#include "lapacke_bits/lapacke_dgecon.h"
#include "lapacke_bits/lapacke_dgeequ.h"
#include "lapacke_bits/lapacke_dgees.h"
#include "lapacke_bits/lapacke_dgeesx.h"
#include "lapacke_bits/lapacke_dgeev.h"
#include "lapacke_bits/lapacke_dgeevx.h"
#include "lapacke_bits/lapacke_dgehrd.h"
#include "lapacke_bits/lapacke_dgelq2.h"
#include "lapacke_bits/lapacke_dgelqf.h"
#include "lapacke_bits/lapacke_dgels.h"
#include "lapacke_bits/lapacke_dgelss.h"
#include "lapacke_bits/lapacke_dgelsy.h"
#include "lapacke_bits/lapacke_dgeqlf.h"
#include "lapacke_bits/lapacke_dgeqp3.h"
#include "lapacke_bits/lapacke_dgeqr2.h"
#include "lapacke_bits/lapacke_dgeqrf.h"
#include "lapacke_bits/lapacke_dgerfs.h"
#include "lapacke_bits/lapacke_dgerqf.h"
#include "lapacke_bits/lapacke_dgesv.h"
#include "lapacke_bits/lapacke_dgesvd.h"
#include "lapacke_bits/lapacke_dgesvx.h"
#include "lapacke_bits/lapacke_dgetf2.h"
#include "lapacke_bits/lapacke_dgetrf.h"
#include "lapacke_bits/lapacke_dgetri.h"
#include "lapacke_bits/lapacke_dgetrs.h"

//* Double precision General matrices Generalized problems -> DGG

#include "lapacke_bits/lapacke_dggbak_work.h"
#include "lapacke_bits/lapacke_dggbal_work.h"
#include "lapacke_bits/lapacke_dgges_work.h"
#include "lapacke_bits/lapacke_dggglm_work.h"
#include "lapacke_bits/lapacke_dgghrd_work.h"
#include "lapacke_bits/lapacke_dgglse_work.h"
#include "lapacke_bits/lapacke_dggqrf_work.h"
#include "lapacke_bits/lapacke_dggrqf_work.h"

#include "lapacke_bits/lapacke_dggbak.h"
#include "lapacke_bits/lapacke_dggbal.h"
#include "lapacke_bits/lapacke_dgges.h"
#include "lapacke_bits/lapacke_dggglm.h"
#include "lapacke_bits/lapacke_dgghrd.h"
#include "lapacke_bits/lapacke_dgglse.h"
#include "lapacke_bits/lapacke_dggqrf.h"
#include "lapacke_bits/lapacke_dggrqf.h"

//* Double precision General Tridiagonal matrices  -> DGT

#include "lapacke_bits/lapacke_dgtcon_work.h"
#include "lapacke_bits/lapacke_dgtrfs_work.h"
#include "lapacke_bits/lapacke_dgtsv_work.h"
#include "lapacke_bits/lapacke_dgtsvx_work.h"
#include "lapacke_bits/lapacke_dgttrf_work.h"
#include "lapacke_bits/lapacke_dgttrs_work.h"

#include "lapacke_bits/lapacke_dgtcon.h"
#include "lapacke_bits/lapacke_dgtrfs.h"
#include "lapacke_bits/lapacke_dgtsv.h"
#include "lapacke_bits/lapacke_dgtsvx.h"
#include "lapacke_bits/lapacke_dgttrf.h"
#include "lapacke_bits/lapacke_dgttrs.h"

#include "lapacke_bits/lapacke_dpotrf_work.h"

#include "lapacke_bits/lapacke_dpotrf.h"

#include "lapacke_bits/lapacke_zgesv_work.h"
#include "lapacke_bits/lapacke_zgetrf_work.h"

#include "lapacke_bits/lapacke_zgesv.h"
#include "lapacke_bits/lapacke_zgetrf.h"

#include "lapacke_bits/lapacke_nancheck.h"

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LAPACKE_H_ */
