{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module FPNLA.Matrix_M (
    Vector(..),
    Matrix(..),
    RowMatrixVector(..),
    ColMatrixVector(..),

    cantRows_m,
    cantCols_m,
    size_m
) where

class Monad mon => Vector mon v e where

    generate_v :: Int
        -> (Int -> mon e)
        -> mon (v e)

    fromList_v :: [e]
        -> mon (v e)

    elem_v :: Int -> v e -> mon e

    length_v :: v e -> mon Int

    slice_v :: Int -> Int -> v e -> mon (v e)

    update_v :: (Int -> mon e) -> v e -> mon (v e)

    -- Default implementations:

--    fromList_v l = generate_v (length l) (\i ->  l !! i)
--    generate_v size f = fromList_v [f x | x <- [0 .. size - 1]]

class Monad mon => Matrix mon m e where

    generate_m :: Int
        -> Int
        -> (Int -> Int -> mon e)
        -> mon (m e)

    fromList_m :: Int -> Int -> [e] -> mon (m e)

    dim_m :: m e -> mon (Int, Int)

    elem_m :: Int
        -> Int
        -> m e
        -> mon e

    subMatrix_m :: Int
        -> Int
        -> Int
        -> Int
        -> m e
        -> mon (m e)

    update_m :: (Int -> Int -> mon e) -> m e -> mon(m e)

    -- Default implementations:

--    fromList_m m n l = generate_m m n (\i j ->  ls !! i !! j)
--        where ls = splitSlice n l
--    generate_m cr cc f = fromList_m cr cc [f x y | x <- [0 .. cr], y <- [0 .. cc]]


class (Monad mon, Vector mon v e, Matrix mon m e) => RowMatrixVector mon m v e where

    row_vm :: Int -> m e -> mon (v e)

    toRows_vm :: m e -> mon [v e]

    toRows_vm m =
        do
            cr <- cantRows_m m
            traverse (`row_vm` m) [0 .. cr - 1]


class (Monad mon, Vector mon v e, Matrix mon m e) => ColMatrixVector mon m v e where

    col_vm :: Int -> m e -> mon (v e)

    toCols_vm :: m e -> mon [v e]

    toCols_vm m =
        do
            cc <- cantCols_m m
            traverse (`col_vm` m) [0 .. cc - 1]

-- | Returns the number of rows of the matrix.
cantRows_m :: (Matrix mon m e) => m e -> mon Int
cantRows_m = fmap fst . dim_m

-- | Returns the number of columns of the matrix.
cantCols_m :: (Matrix mon m e) => m e -> mon Int
cantCols_m = fmap snd . dim_m

-- | Returns the size of the matrix (number of rows x number of columns).
size_m :: (Matrix mon m e) => m e -> mon Int
size_m = fmap (uncurry (*)) . dim_m
