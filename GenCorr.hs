{-# LANGUAGE FlexibleContexts #-}                

module GenCorr where
       
import Control.Applicative                
import Control.Monad                
import System.Random.MWC
import Data.Random
import Data.Random.Distribution.Normal
import Data.Random.Distribution.Categorical
import DataUtils

type Diffusivity = Double
type Time = Double
type Corr = Double

wxy = 1
aspect = 1

components :: [(Double, (Diffusivity, Diffusivity))]     
components = [ (0.5, (10, 1e-4))
             -- , (0.5, (20, 1e-4))
             ]

logNormal :: (Floating a, Distribution Normal a) => a -> a -> RVar a
logNormal mu std = do
    z <- stdNormal
    return $ exp (mu + std*z)

diffusivity :: [(Double, (Diffusivity, Diffusivity))] -> RVar Diffusivity
diffusivity comps = do
    (mu,std) <- categorical comps
    logNormal mu std

acorr :: Diffusivity -> Time -> Corr
acorr d tau = 1 / (1 + tau/tauD) / sqrt (1 + aspect^2 * tau/tauD)
    where tauD = wxy^2 / 4 / d
    
noisify :: (Num a, Distribution Normal a) => (a -> a) -> (a -> RVar a) -> a -> RVar a
noisify std m x = (+) <$> m x <*> normal 0 (std x)

heteroCorr :: Int -> [(Double, (Diffusivity, Diffusivity))] -> RVar (Time -> Corr, [Diffusivity])
heteroCorr n comps = do
    ds <- replicateM n $ diffusivity comps
    return (\tau->sum $ map (flip acorr tau) ds, ds)
           
generateCorr n = withSystemRandom $ asGenIO $ runRVar (heteroCorr n components)

main' = withSystemRandom $ asGenIO $ \mwc->do
    --let g = noisify (const 1e-4) (\tau->flip acorr tau <$> diffusivity components)
    --print =<< runRVar (mapM g [1..10000]) mwc

    --DataUtils.save "hi2" $ map (\tau->(tau, acorr 10 tau)) $ logSpace 1e-6 1e5 200

    (g, ds) <- runRVar (heteroCorr 40 components) mwc
    print ds
    DataUtils.save "hi" $ map (\tau->(tau, g tau)) $ logSpace 1e-6 1e5 200
    
logSpace :: (RealFloat a, Enum a) => a -> a -> Int -> [a]
logSpace a b n = map exp $ [la,la+dx..lb]
  where la = log a
        lb = log b
        dx = (lb - la) / fromIntegral n