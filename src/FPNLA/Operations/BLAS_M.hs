{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables   #-}
{-# LANGUAGE TypeFamilies          #-}

module FPNLA.Operations.BLAS_M (

    GEMM_M(..)

) where

import FPNLA.Matrix_M                (Matrix(), Vector())
import FPNLA.Operations.Parameters (Elt (), ResM (),
                                    StratCtx (), TransType ())


class (Elt e, Matrix mon m e, Vector mon v e, Monad mon) => GEMM_M mon s m v e where
    gemm_m :: StratCtx s
        -> TransType (m e)
        -> TransType (m e)
        -> e
        -> e
        -> m e
        -> mon(ResM s v m e)



