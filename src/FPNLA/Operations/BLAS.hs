{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables   #-}
{-# LANGUAGE TypeFamilies          #-}

-- | This module defines all the BLAS (Basic Linear Algebra Subprograms) operations supported by the framework.
-- See <http://www.netlib.org/blas/> for more information about BLAS and <http://www.ugcs.caltech.edu/~srbecker/blasqr_betterFonts.pdf> for a quick description of all BLAS operation signatures and behaviour.
module FPNLA.Operations.BLAS (
    -- *  Level One
    -- | Vector-Vector operations
    DOT(..),
    -- *  Level Two
    -- | Matrix-Vector operations
    GEMV(..),
    -- *  Level Three
    -- | Matrix-Matrix operations
    SYRK(..),
    GEMM(..),
    TRSM(..),
    Elt()

) where

import FPNLA.Matrix                (MatrixVector, Vector)
import FPNLA.Operations.Parameters (Elt (), ResM (), ResV (), ResS(),
                                    StratCtx (), TransType (),
                                    TriangType (), UnitType ())

-- | Defines the signature of the level-1 BLAS /dot/ operation in the framework.
class (Elt e, Vector v e) => DOT s v e where
    dot :: StratCtx s -- ^ The context of the operation
        -> v e -- ^ A vector /x/
        -> v e -- ^ A vector /y/
        -> ResS s e -- ^ The scalar product between /x/ and /y/

-- | Defines the signature of the level-2 BLAS /gemv/ operation in the framework.
class (Elt e, MatrixVector m v e) => GEMV s m v e where
    gemv :: StratCtx s -- ^ The context of the operation
        -> TransType (m e) -- ^ A matrix /A/
        -> v e -- ^ A vector /x/
        -> e -- ^ A scalar /alpha/
        -> e -- ^ A scalar /beta/
        -> v e  -- ^ A vecor /y/
        -> ResV s v e -- ^ @alpha * A * x + beta * y@

-- | Defines the signature of the level-3 BLAS /syrk/ operation in the framework.
class (Elt e, MatrixVector m v e) => SYRK s m v e where
    syrk :: StratCtx s -- ^ The context of the operation
        -> e -- ^ A scalar /alpha/
        -> TransType (m e) -- ^ A matrix /A/
        -> e -- ^ A scalar /beta/
        -> TriangType (m e) -- ^ A triangular matrix /C/
        -> ResM s v m e -- ^ @alpha * A * A' + beta * C@ where @A'@ is the conjugate transposed of /A/.

-- | Defines the signature of the level-3 BLAS /gemm/ operation in the framework.
class (Elt e, MatrixVector m v e) => GEMM s m v e where
    gemm :: StratCtx s -- ^ The context of the operation
        -> TransType (m e) -- ^ A matrix /A/
        -> TransType (m e) -- ^ A matrix /B/
        -> e -- ^ A scalar /alpha/
        -> e -- ^ A scalar /beta/
        -> m e -- ^ A matrix /C/
        -> ResM s v m e -- ^ @alpha * A * B + beta * C@

-- | Defines the signature of the level-3 BLAS /trsm/ operation in the framework.
class (Elt e, MatrixVector m v e) => TRSM s m v e where
    trsm :: StratCtx s -- ^ The context of the operation
        -> e -- ^ A scalar /alpha/
        -> TransType (TriangType (UnitType (m e))) -- ^ A triangular matrix /A/
        -> m e -- ^ A matrix /B/
        -> ResM s v m e -- ^ @alpha * A_inv * B@ where @A_inv@ is the inverse of /A/.
