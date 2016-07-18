{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}

-- | This module defines classes used to handle Matrices and Vectors in a generic way.
module FPNLA.Matrix (
    -- * Vector
    Vector(..),
    -- * Matrix
    Matrix(..),
    -- * MatrixVectror
    MatrixVector(..),
    -- * Miscellaneous
    cantRows_m,
    cantCols_m,
    size_m,
    asColumn_vm,
    diagonalBlock,
    verticalBlock
) where

import FPNLA.Utils (iif, splitSlice)

-- | This class represents a vector structure
-- For any vector data structure used in this framework it must have an instance of this class.
-- A minimal instance of this class must provide an implementation for 'generate_v', 'elem_v' and 'length_v'.
-- The index of the first element stored in the vector must be 0.
-- This is a multi-param type class because usually libraries of Matrix and Vector imposes restrictions over the type of the elements (for example, to have a "Storable" instance).
class Vector v e where
    -- *Constructors. This functions allows to create a new vector.
    
    -- | Creates a new vector from the length of the vector and a function that returns the value for each index.
    generate_v :: Int -- ^ The length of the vector
        -> (Int -> e) -- ^ A function that returns the value for each index
        -> v e
    -- | Creates a new vector from a list, th vector will contain all (and only) the elements of the list, and in the same order.
    fromList_v :: [e] -- ^ A list that contains the elements of the vector. The first element of the list will be the first element of the vector and so on.
        -> v e
    -- | Join multiple vectors in one single vector.
    concat_v :: [v e] -- ^ A list containing the vectors to be joined. If the list is @[<1>, <2,3>, <4,5,6>]@ then the vector will be @<1,2,3,4,5,6>@.
        -> v e
    
    -- * Selectors. This functions returns information from a vector.
    -- | Returns the element from a given index. The index must be between 0 and the length of the vector.
    elem_v :: Int -> v e -> e
    -- | Returns the length of the vector.
    length_v :: v e -> Int
    
    -- * High-order. This functions are the classics /fold/, /map/ and /zipWith/ but are in the class to allow to get an specific implementation from a specific data structure in order to achieve better performance.
    -- | The common /foldr/ function but specific for this instance.
    foldr_v :: (e -> b -> b) -> b -> v e -> b
    -- | The common /map/ function but specific for this instance.
    map_v :: (e -> e) -> v e -> v e
    -- | The common /zipWith/ function but specific for this instance.
    zipWith_v :: (e -> e -> e) -> v e -> v e -> v e

    -- Default implementations:
    fromList_v l = generate_v (length l) (\i ->  l !! i)
    concat_v [] = generate_v 0 undefined
    concat_v [v] = v
    concat_v (v1:v2:vs) = concat_v (join:vs)
        where l1 = length_v v1
              l2 = length_v v2
              gen i | i < l1 = elem_v i v1
                    | otherwise = elem_v (i-l1) v2
              join = generate_v (l1 + l2) gen
    foldr_v f z v = go 0
        where len = length_v v
              go i | i == len   = z
                   | otherwise  = f (elem_v i v) (go (i+1))
    map_v f v = generate_v (length_v v) (\i -> f (elem_v i v))
    zipWith_v f v1 v2 = generate_v (length_v v1) (\i -> f (elem_v i v1) (elem_v i v2))


-- | This class represents any Matrix structure
-- For any matrix data structure used in this framework it must have an instance of this class.
-- A minimal instance of this class must provide an implementation for 'generate_m', 'dim_m' and 'elem_m'.
-- The index of the first element stored in the vector must be (0,0).
-- This is a multi-param type class because usually libraries of Matrix and Vector imposes restrictions over the type of the elements (for example, to have an "Storable" instance).
class Matrix m e where
    -- * Constructors. This functions allows to create a new matrix.
    -- | Creates a new matrix from the dimension (length of rows and columns) of the matrix and a function that returns the value for each index.
    generate_m :: Int -- ^ The length of the rows
        -> Int -- ^ A function that returns the value for each index
        -> (Int -> Int -> e) -- ^ 
        -> m e
    -- | Creates a new matrix from the dimension (length of rows and columns) of and a list.
    fromList_m :: Int -> Int -> [e] -> m e
    -- | Given a matrix, creates a new matrix that will be the transpose of the original.
    transpose_m :: m e -> m e

    -- * Selectors. This functions returns information from a vector.
    -- | Returns the dimension of a given matrix. The first element in the pair is the number of rows and the second the number of columns.
    dim_m :: m e -> (Int, Int)
    -- | Returns the element that is in the position @(i,j)@.
    elem_m :: Int -- ^ The row component of the position (@i@).
        -> Int -- ^ The column component of the position (@j@).
        -> m e -- ^ The matrix.
        -> e

    -- * High-order. This functions are the classics /fold/, /map/ and /zipWith/ but are in the class to allow to get an specific implementation from a specific data structure in order to achieve better performance.
    -- | The common /map/ function but specific for this instance.
    map_m :: (e -> e) -> m e -> m e
    -- | The common /zipWith/ function but specific for this instance.
    zipWith_m :: (e -> e -> e) -> m e -> m e -> m e

    -- * Blocks. This functions are constructors too, but this splits or allow to get a sub-matrix from a given matrix.
    -- | Given an initial position @(i,j)@, a dimension @(r,c)@ and a matrix @m@, creates a new matrix that is a sub matrix of @m@ containing the elements in positions @[(p1, p2) | p1 <- [i..i+r], p2 <- [j..j+c]]@
    -- For example, if @m@ is:
    -- @1 2 3
    -- 4 5 6
    -- 7 8 9@
    -- and @(i,j)@ and @(r,c)@ are @(0,1)@ and @(2,3)@ respectively, then the result will be a new matrix containing:
    -- @2 3
    -- 5 6
    -- 8 9@
    subMatrix_m :: Int -- ^ @i@
        -> Int -- ^ @j@
        -> Int -- ^ @r@
        -> Int -- ^ @c@
        -> m e -- ^ @m@
        -> m e
    -- | Constructs a new matrix from a list of list where the elements are sub-matrices.
    -- This operations is the analogue of 'concat_v' operation from the 'Vector' class.
    -- The \'borders\' of the sub-matrices will be joined.
    -- The dimensions of the sum-matrices must be compatibles.
    fromBlocks_m :: [[m e]] -> m e
    -- | Split a matrix @m@ in blocks of size @(r,c)@ (the dimension of the blocks in the right and the bottom may be smaller).
    toBlocks_m :: Int -- ^ @r@
        -> Int -- ^ @c@
        -> m e -- ^ @m@
        -> [[m e]]

    -- Default implementations:
    fromList_m m n l = generate_m m n (\i j ->  ls !! i !! j)
        where ls = splitSlice n l
    transpose_m m = generate_m (cantCols_m m) (cantRows_m m) (\i j -> elem_m j i m)
    map_m f m = generate_m (cantRows_m m) (cantCols_m m) (\i j -> f (elem_m i j m))
    zipWith_m f m1 m2 = generate_m (cantRows_m m1) (cantCols_m m1) (\i j -> f (elem_m i j m1) (elem_m i j m2))
    subMatrix_m posI posJ cantFilas cantCols m = generate_m cantFilas cantCols (\i j -> elem_m (i + posI) (j + posJ) m)
    fromBlocks_m [] = fromList_m 0 0 []
    fromBlocks_m [[]] = fromList_m 0 0 []
    fromBlocks_m bss = 
        generate_m cr cc $ \i j -> elem_m (iToB i) (jToB j) $ bss !! iToBSS i !! jToBSS j
        where bssm = length bss
              bssn = length $ head bss
              (br, bc) = dim_m . head . head $ bss
              lbr = cantRows_m . head . last $ bss
              lbc = cantCols_m . last . head $ bss
              ilimit = (bssm-1) * br
              jlimit = (bssn-1) * bc
              cr = ilimit + lbr
              cc = jlimit + lbc
              iToBSS i = iif (i < ilimit) (i `div` br) (bssm-1)
              jToBSS j = iif (j < jlimit) (j `div` bc) (bssn-1)
              iToB i = iif (i < ilimit) (i `mod` br) ((i-ilimit) `mod` lbr)
              jToB j = iif (j < jlimit) (j `mod` bc) ((j-jlimit) `mod` lbc)
    toBlocks_m brs bcs m
        | brs == 0 || bcs == 0 = []
        | otherwise = map (rowBlocks 0) [i * brs | i <- [0 .. (rs `div` brs + iif (rs `mod` brs > 0) 1 0) - 1]]
        where
            (rs, cs) = dim_m m
            getBound s b c = iif (s + b > c) (c - s) b
            rowBlocks sJ sI
                    | sJ >= cs = []
                    | otherwise = subMatrix_m sI sJ (getBound sI brs rs) (getBound sJ bcs cs) m : rowBlocks (sJ + bcs) sI


-- | This class allows us to leave separated the two concepts of 'Vector' and 'Matrix'.
-- This class establishes a link between a vector and a matrix structure.
-- We provide a default instance for any 'Vector' and 'Matrix' instances so any vector could be used with any matrix whenever the type of the elements are the same.
class (Vector v e, Matrix m e) => MatrixVector m v e where
    -- | Extracts the @i@-th row of the matrix @m@ and converts it in the vector structure.
    row_vm :: Int -> m e -> v e
    -- | Extracts the @j@-th column of the matrix @m@ and converts it in the vector structure.
    col_vm :: Int -> m e -> v e
    -- | Joins a list of column vectors in a matrix structure.
    -- The length of every vector in the list must be the same.
    fromCols_vm :: [v e] -> m e
    -- | Returns a list containing every column of the matrix converted in the vector structure.
    toCols_vm :: m e -> [v e]

    -- Default implementations:
    row_vm i m = generate_v (cantCols_m m) (\j -> elem_m i j m)
    col_vm j m = generate_v (cantRows_m m) (\i -> elem_m i j m)
    fromCols_vm css | cc == 0 = fromList_m 0 0 []
                  | otherwise = generate_m cr cc (\i j -> elem_v i $ css !! j)
        where cc = length css
              cr = length_v $ head css
    toCols_vm m = map (`col_vm` m) [0 .. cantCols_m m - 1]


instance {-# OVERLAPPABLE #-} (Vector v e, Matrix m e) => MatrixVector m v e

-- | Returns the number of rows of the matrix.
cantRows_m :: (Matrix m e) => m e -> Int
cantRows_m = fst . dim_m

-- | Returns the number of columns of the matrix.
cantCols_m :: (Matrix m e) => m e -> Int
cantCols_m = snd . dim_m

-- | Returns the size of the matrix (number of rows x number of columns).
size_m :: (Matrix m e) => m e -> Int
size_m = uncurry (*) . dim_m

-- | Converts a vector in a matrix of one column.
asColumn_vm :: MatrixVector m v e => v e -> m e
asColumn_vm = fromCols_vm . (:[])

-- | Splits a matrix @m@ into blocks of size @(r,c)@ and returns the block on the position @(i,i)@.
diagonalBlock :: (Matrix m e) =>  (Int, Int) -- ^ @(r,c)@
    -> Int -- ^ @i@
    -> m e -- ^ @m@
    -> m e
diagonalBlock (blockDimM, blockDimN) i m =
    subMatrix_m (blockDimM*i) (blockDimN*i) (min blockDimM (cantRows_m m - (blockDimM * i))) (min blockDimN (cantCols_m m - (blockDimN * i))) m



-- | Splits a matrix @m@ into blocks of size @(s, cantCols_m m)@ and returns the block on the position @i@.
verticalBlock :: (Matrix m e) =>  Int -- ^ @s@
    -> Int -- ^ @i@
    -> m e -- ^ @m@
    -> m e
verticalBlock blockDim i m =
    subMatrix_m (blockDim*i) 0 (min blockDim (cantRows_m m - (blockDim * i))) (cantCols_m m) m

