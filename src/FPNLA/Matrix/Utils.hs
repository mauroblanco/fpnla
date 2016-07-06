-- | This module defines commonly useful functions that are related specifically with vectors and matrices.
module FPNLA.Matrix.Utils where


import FPNLA.Matrix (Matrix(..))


-- | Prints a matrix to the standard output.
-- This operation requires the elements of the matrix to have an instance of 'Show' but does not requires a 'Show' instance for the matrix data type.
print_m :: (Show e, Matrix m e) => m e -> IO ()
print_m mi = for 0 0
    where (m,n) = dim_m mi
          for i j | i >= m = return ()
                  | j < n = (putStr . show $ elem_m i j mi) >> putStr " " >> for i (j+1)
                  | j >= n = putStrLn "" >> for (i+1) 0
