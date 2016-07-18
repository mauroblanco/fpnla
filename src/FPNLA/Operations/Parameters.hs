{-# LANGUAGE TypeFamilies #-}

module FPNLA.Operations.Parameters(
    -- * Elements
    Elt(..),
    -- * Strategies and contexts
    StratCtx(), 
    
    -- * Result type
    -- | In BLAS it's common that operations in higher levels use operations in the lower levels, so, an operation in level three that by its signature manipulates matrices only, internally uses level two operations that manipulates vectors. In order to avoid the /show . read/ problem, the type of the vector (or any other internal data type) must appear in the signature of an operation.
    -- To solve the problem we use phantom types to pass the internally used types to the Haskell type system.
    ResM(),
    ResV(),
    ResS(),
    blasResultM,
    blasResultV,
    blasResultS,
    getResultDataM,
    getResultDataV,
    getResultDataS,
    
    -- * Miscellaneous
    TransType(..),
    UnitType(..),
    TriangType(..),
    unTransT, 
    unUnitT, 
    unTriangT, 
    elemTrans_m,
    dimTrans_m,
    elemSymm,
    dimTriang,
    elemUnit_m,
    dimUnit_m,
    elemTransUnit_m,
    dimTransUnit_m,
    transTrans_m
    
) where

import FPNLA.Matrix(Matrix(..))
import Data.Complex (Complex, conjugate)
import Data.Tuple (swap)


-- | This class represents the elements that can be used in the BLAS operations.
-- The elements in BLAS are real or complex numbers, so we provide default instances for the Haskell 'Double', 'Float' and 'Complex' types.
class (Eq e, Floating e) => Elt e where
    -- | Returns the conjugate of a number. For real numbers it's the identity function and for complex numbers it's the common 'Complex.conjugate' function.
    getConjugate :: e -> e
    getConjugate = id

instance Elt Double
instance Elt Float
instance (RealFloat e) => Elt (Complex e) where
    getConjugate = conjugate

-- | This type family is used to represent the /context/ of an operation.
-- A particular implementation is a combination of an algorithm and a parallelism technique, and we call it a /strategy/. A particular strategy may need particular information to execute. For example, an operation that computes the matrix-matrix multiplication by splitting the matrices in blocks must require the size of the blocks.
-- With this context we allows to pass any additional information that the operation needs to execute as parameters, but maintaining a common signature.
-- The /s/ type parameter is the strategy so, there must exist a Haskell data type to represent a particular strategy.
type family StratCtx s :: *


-- | The 'ResM' data type is used as result of level three BLAS operations and returns a matrix /m/ of elements /e/ and contains the strategy /s/ and vector /v/ as phantom types.
data ResM s (v :: * -> *) m e = ResM { unResM :: m e } deriving (Show)
-- | The 'ResV' data type is used as result of level two BLAS operations and returns a vector /v/ of elements /e/ and contains the strategy /s/ as phantom types.
data ResV s v e = ResV { unResV :: v e } deriving (Show)
-- | The 'ResS' data type is used as result of level one BLAS operations and returns an scalar /e/ and contains the strategy /s/ as phantom types.
data ResS s e = ResS { unResS :: e } deriving (Show)


-- | Wrap a matrix into a 'ResM'.
blasResultM :: m e -> ResM s v m e
blasResultM = ResM
-- | Unwrap a matrix from a 'ResM'.
getResultDataM :: ResM s v m e -> m e
getResultDataM = unResM

-- | Wrap a vector into a 'ResV'.
blasResultV :: v e -> ResV s v e
blasResultV = ResV
-- | Unwrap a vector from a 'ResV'.
getResultDataV :: ResV s v e -> v e
getResultDataV = unResV

-- | Wrap a scalar into a 'ResS'.
blasResultS :: e -> ResS s e
blasResultS = ResS
-- | Unwrap a scalar from a 'ResS'.
getResultDataS :: ResS s e -> e
getResultDataS = unResS


-- | Indicates if a matrix must be considered as normal, transposed or transposed conjugated.
-- This is part of the common flags in the BLAS operation signatures and it's useful to work with a transposed matrix without really computing the transposed matrix.
data TransType m = Trans m | NoTrans m | ConjTrans m deriving (Eq, Show)
-- | Indicates if a matrix must be considered as unitary or not. An unitary matrix is a matrix that contains ones in the diagonal.
-- This is part of the common flags in the BLAS operation signatures.
data UnitType m = Unit m | NoUnit m deriving (Eq, Show)
-- | Indicates that a matrix is symmetric and with which triangular part of the matrix the operation is going to work ('Upper' or 'Lower').
-- The operation only will see the indicated part of the matrix and should not try to access the other part.
-- This is part of the common flags in the BLAS operation signatures.
data TriangType m = Lower m | Upper m deriving (Eq, Show)

-- | Given a data type flagged by a TransType, returns a pair containing the TransType constructor and the data type.
unTransT :: TransType a -> (b -> TransType b, a)
unTransT (Trans a) = (Trans, a)
unTransT (NoTrans a) = (NoTrans, a)
unTransT (ConjTrans a) = (ConjTrans, a)

-- | Given a data type flagged by a UnitType, returns a pair containing the UnitType constructor and the data type.
unUnitT :: UnitType a -> (b -> UnitType b, a)
unUnitT (Unit a) = (Unit, a)
unUnitT (NoUnit a) = (NoUnit, a)

-- | Given a data type flagged by a TriangType, returns a pair containing the TriangType constructor and the data type.
unTriangT :: TriangType a -> (b -> TriangType b, a)
unTriangT (Lower a) = (Lower, a)
unTriangT (Upper a) = (Upper, a)


-- | Given an /i,j/ position and a TransType flagged matrix, returns the element in that position without computing the transpose.
elemTrans_m :: (Elt e, Matrix m e) => Int -> Int -> TransType (m e) -> e
elemTrans_m i j (NoTrans m) = elem_m i j m
elemTrans_m i j (Trans m) = elem_m j i m
elemTrans_m i j (ConjTrans m) = getConjugate $ elem_m j i m


-- | Given a TransType flagged matrix, returns the dimension of the matrix without computing the transpose.
dimTrans_m :: (Matrix m e) => TransType (m e) -> (Int, Int)
dimTrans_m (NoTrans m) = dim_m m
dimTrans_m (ConjTrans m) = swap $ dim_m m
dimTrans_m (Trans m) = swap $ dim_m m


-- | Given an /i,j/ position and a TransType flagged matrix, returns the element in that position only accessing the part indicated by the TransType.
elemSymm :: (Matrix m e) => Int -> Int -> TriangType (m e) -> e
elemSymm i j (Upper m)
    | i > j = elem_m j i m
    | otherwise = elem_m i j m
elemSymm i j (Lower m)
    | i > j = elem_m i j m
    | otherwise = elem_m j i m

-- | Given a TransType flagged matrix, returns the dimension of the matrix.
dimTriang :: (Matrix m e) => TriangType (m e) -> (Int, Int)
dimTriang = dim_m . snd . unTriangT

-- | Given an /i,j/ position and a UnitType flagged matrix, returns the element in that position. If the matrix is flagged as Unit and /i == j/ (the element is in the diagonal) returns one.
elemUnit_m :: (Elt e, Matrix m e) => Int -> Int -> UnitType (m e) -> e
elemUnit_m i j (Unit m)
    | i == j = 1
    | otherwise = elem_m i j m
elemUnit_m i j (NoUnit m) = elem_m i j m

-- | Given a UnitType flagged matrix, returns the dimension of the matrix.
dimUnit_m :: (Matrix m e) => UnitType (m e) -> (Int, Int)
dimUnit_m (Unit m) = dim_m m
dimUnit_m (NoUnit m) = dim_m m

-- | Given an /i,j/ position and a TransType-UnitType flagged matrix, returns the element in that position without computing the transpose.
elemTransUnit_m :: (Elt e, Matrix m e) => Int -> Int -> TransType (UnitType (m e)) -> e
elemTransUnit_m i j (NoTrans pmA) = elemUnit_m i j pmA
elemTransUnit_m i j (Trans pmA) = elemUnit_m j i pmA
elemTransUnit_m i j (ConjTrans pmA) = getConjugate $ elemUnit_m j i pmA

-- | Given a TransType-UnitType flagged matrix, returns the dimension of the matrix.
dimTransUnit_m :: Matrix m e => TransType (UnitType (m e)) -> (Int, Int)
dimTransUnit_m = dimUnit_m . snd . unTransT

-- | Given a TransType flagged matrix, computes and returns its transpose.
transTrans_m :: (Elt e, Matrix m e) => TransType (m e) -> m e
transTrans_m (NoTrans m) = m
transTrans_m (ConjTrans m) = map_m getConjugate $ transpose_m m
transTrans_m (Trans m) = transpose_m m
