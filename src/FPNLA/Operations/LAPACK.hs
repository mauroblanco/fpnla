{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables   #-}
{-# LANGUAGE TypeFamilies          #-}

-- | This module defines all the LAPACK (Linear Algebra PACKage) operations supported by the framework.
-- See <http://www.netlib.org/lapack/> for more information about LAPACK.
module FPNLA.Operations.LAPACK (

    POTRF(..),

) where

import FPNLA.Matrix                (MatrixVector)
import FPNLA.Operations.Parameters (Elt (), ResM (), StratCtx (), TriangType ())



-- Po: Symmetric matrix or Hermitian matrix positive definite
-- Trs: Triangular Factorization

-- | Defines the signature of the LAPACK /potrf/ operation in the framework.
-- This operation takes a symmetric (or hermitian) positive definite (SPD) matrix (flagged with TriangType) and computes the Cholesky factorization of the matrix.
-- The Cholesky decomposition of an SPD matrix /M/ is a lower triangular matrix /L/ where /M = L L*/ being /L*/ the conjugate transpose of /L/.
class (Elt e, MatrixVector m v e) => POTRF s m v e where
    potrf :: StratCtx s -> TriangType (m e) -> ResM s v m e
