module DataUtils (save, load) where                

import Control.Applicative
import Data.Csv
import Data.Char
import qualified Data.ByteString.Lazy as BS       
import qualified Data.Vector as V
import Data.Foldable as Foldable       

save :: (Foldable f, ToRecord a) => FilePath -> f a -> IO ()
save fname xs =
    BS.writeFile fname $ encodeWith opts $ V.fromList $ Foldable.toList xs
  where opts = defaultEncodeOptions { encDelimiter=fromIntegral $ ord '\t' }

load :: (FromRecord a) => FilePath -> IO [a]
load fname = V.toList . either error id . decode False <$> BS.readFile fname