{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}

{-# LANGUAGE TemplateHaskell    #-}
{-# LANGUAGE GADTs              #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE BangPatterns       #-}
{-# LANGUAGE GADTs              #-}
{-# LANGUAGE RankNTypes         #-}
{-# LANGUAGE DeriveDataTypeable #-}

import Data.Random
import Data.Random.Source

import qualified System.Random.MWC as MWC

import Control.Monad.Reader
import Control.Monad.Primitive

-- import Data.Word
-- import Data.Typeable


-- data Prim a where
--     PrimWord8           :: Prim Word8
--     PrimWord16          :: Prim Word16
--     PrimWord32          :: Prim Word32
--     PrimWord64          :: Prim Word64
--     PrimDouble          :: Prim Double
--     PrimNByteInteger    :: !Int -> Prim Integer
--     deriving (Typeable)

-- class Monad m => MonadRandom m where
--     getRandomPrim :: Prim t -> m t
--     getRandomPrim PrimWord8             = getRandomWord8
--     getRandomPrim PrimWord16            = getRandomWord16
--     getRandomPrim PrimWord32            = getRandomWord32
--     getRandomPrim PrimWord64            = getRandomWord64
--     getRandomPrim PrimDouble            = getRandomDouble

$(monadRandom [d|
  instance (PrimMonad m, s ~ PrimState m) => MonadRandom (ReaderT (MWC.Gen s) m) where
    getRandomWord16 = ask >>= lift . MWC.uniform
    getRandomWord32 = ask >>= lift . MWC.uniform
    getRandomWord64 = ask >>= lift . MWC.uniform
    getRandomDouble = ask >>= lift . MWC.uniform
  |])

foo :: MonadRandom m => m Double
foo = do
    let testUniform 0 !x = return x
        testUniform n !x = do
            y <- sample stdUniform
            testUniform (n - 1) (x + y)

    testUniform n 0

n :: Int
n = 10^8

viaIO :: IO Double
viaIO = do
    seed <- MWC.create
    xs <- runReaderT foo seed
    return xs

main :: IO ()
main = do
  x <- viaIO
  print $ x / fromIntegral n
