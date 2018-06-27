{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Main where

import Foreign.Ptr
import Foreign.Storable
import Foreign.C.Types
import Foreign.ForeignPtr
import Foreign.Marshal.Array

import Test.QuickCheck.Monadic
import Test.QuickCheck

import Data.Word8
import Data.ByteString.Unsafe (unsafePackCStringFinalizer)
import ArbBS
import Data.ByteString (ByteString)

import qualified Data.ByteString          as BS
import qualified Data.ByteString.Internal as BS

import qualified Crypto.Hash.SHA256 as H
import qualified Data.Vector.Storable as V

foreign import ccall "WriteG" sha256write :: Ptr Word8 -> CSize -> IO ()
foreign import ccall "FinalizeG" sha256finalize :: Ptr Word8 -> IO ()
foreign import ccall "ResetG" sha256reset :: IO ()
foreign import ccall "GlobalInit" globalInit :: IO ()
foreign import ccall "InitG" sha256init  :: IO ()

write :: ByteString -> IO ()
write bs = do
  let (fptr, off, len) = BS.toForeignPtr bs
  withForeignPtr fptr $ \ptr -> do
    sha256write ptr (fromIntegral len)

finalize :: Ptr Word8 -> IO ByteString
finalize src = do
  sha256finalize src
  BS.create 32 $ \dest -> do
    BS.memcpy dest src 32

prop_pureEquiv buf (fromABS -> bs) = monadicIO $ do
   h <- run (hash buf bs)
   assert (H.hash bs == h)

hash :: Ptr Word8 -> ByteString -> IO ByteString
hash buf bs = do
  sha256reset
  write bs
  finalize buf

main :: IO ()
main = do
  globalInit
  buf <- mallocArray 32
  quickCheckWith stdArgs { maxSuccess = 100000 } (prop_pureEquiv buf)
