{-# LANGUAGE FlexibleInstances    #-}
{-# LANGUAGE ScopedTypeVariables  #-}
{-# LANGUAGE UndecidableInstances #-}

-- | This module defines commonly useful functions that are not specifically related with vectors and matrices.
module FPNLA.Utils (
    iif,
    mapPair,
    splitSlice,

    Truncable(..),

    prevEnum,
    nextEnum,
    enumerate

) where

import Data.Complex (Complex ((:+)))

-- | Takes two functions @f@ and @g@, and a pair, and applies @f@ to the first element and @g@ to the second element of the pair.
mapPair :: (a -> b) -- ^ @f@
    -> (c -> d) -- ^ @g@
    -> (a, c) -- ^ Input pair @(a,b)@
    -> (b, d) -- ^ Output pair @(f a, g b)@
mapPair f g (a, c) = (f a,  g c)

-- | It's equivalent to @if cond then e1 else e2@, but defined as a function.
iif :: Bool -- ^ @cond@
    -> a -- ^ @e1@
    -> a -- ^ @e2@
    -> a -- ^ @if cond then e1 else e2@
iif cond expTrue expFalse = if cond then expTrue else expFalse


-- | Splits a list into pieces of size @n@.
-- The last piece could be smaller if the length of the list is not multiple of @n@.
splitSlice :: Int -- ^ @n@
    -> [a] -- ^ The input list
    -> [[a]] -- ^ The input list divided in pieces of size @n@ 
splitSlice 0 _ = error "n must be greater than 0"
splitSlice _ [] = []
splitSlice n ls = sl : splitSlice n sls
    where (sl, sls) = splitAt n ls

-- | This class is used to compare two float point numbers.
-- With float point numbers, two different computations that theoretically returns the same value in fact returns different values. Even if the difference is so small that is ignorable, if we test using equality the test fail.
-- For this, instead to check for equality we define that two float point numbers are the same if the distance between them is less than some epsilon value.
class Truncable a where
    -- | Check if the distance between two float point values is less than some epsilon value.
    inEpsilonRange :: Real e => e -- ^ Epsilon
        -> a -- ^ A float point value
        -> a -- ^ Another float point value
        -> Bool

-- | We provide a default instance of 'Truncable' for any data type that has 'Floating' and 'RealFrac' instances.
instance (Floating n, RealFrac n) => Truncable n where
    inEpsilonRange epsilon n1 n2 = expN1 == expN2 && abs(mantN1 - mantN2) <= eps
        where
            eps = fromRational $ toRational epsilon
            shift n d1 = n * (10 ** fromIntegral d1)
            getExponent r
                    | r == 0 = 0
                    | r < 0 = getExponent (abs r)
                    | r > 10 = 1 + getExponent (r / 10)
                    | r < 1 = (-1) + getExponent (r * 10)
                    | otherwise = 0
            expN1 :: Integer = getExponent n1
            expN2 :: Integer = getExponent n2
            mantN1 = shift n1 (-expN1)
            mantN2 = shift n2 (-expN2)

-- | We define a specific instance of 'Truncable' for the 'Complex' data type.
instance (RealFrac n, Truncable n) => Truncable (Complex n) where
    inEpsilonRange epsilon (r1:+i1) (r2:+i2) = inEpsilonRange epsilon r1 r2 && inEpsilonRange epsilon i1 i2


-- | Returns the predecessor of a enumerated value.
-- The operation is cyclic so, the predecessor of the first value is the last value.
prevEnum :: (Eq a, Bounded a, Enum a) => a -> a
prevEnum a | a == minBound = maxBound
           | otherwise = pred a

-- | Returns the successor of a enumerated value.
-- The operation is cyclic so, the successor of the last value is the first value.
nextEnum :: (Eq a, Bounded a, Enum a) => a -> a
nextEnum a | a == maxBound = minBound
           | otherwise = succ a

-- | Returls a list containing all the values of an enumerated.
enumerate :: forall a. (Enum a, Bounded a) => [a]
enumerate = map toEnum [fromEnum (minBound :: a) .. fromEnum (maxBound :: a)]
