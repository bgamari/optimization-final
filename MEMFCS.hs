{-# LANGUAGE DeriveFunctor, DeriveFoldable, DeriveTraversable, RecordWildCards #-}

import           Control.Applicative
import           Data.Csv
import           Data.Foldable as Foldable
import           Data.List hiding (sum)
import           Data.Traversable
import qualified Data.Vector as V
import           DataUtils
import           Debug.Trace
import           GenCorr
import           Linear
import           Numeric.AD
import           Numeric.AD.Types (auto)
import           Optimization.Constrained.ProjectedSubgradient
import           Optimization.LineSearch
import           Optimization.LineSearch.MirrorDescent
import           Optimization.LineSearch.SteepestDescent
import           Prelude hiding (sum, mapM)

type Species a = V.Vector a
     
sumSqResidual :: (Functor f, Foldable f, RealFrac a) => (a -> a) -> f (Obs a) -> a
sumSqResidual f = Foldable.sum . fmap (\(Obs x y s)->(f x - y)^2 / s)

data Obs a = Obs { oX :: !a   -- ^ Abscissa
                 , oY :: !a   -- ^ Ordinate
                 , oE :: !a   -- ^ Variance
                 }
             deriving (Show, Eq, Functor, Foldable, Traversable)

instance ToField a => ToRecord (Obs a) where
    toRecord (Obs {..}) = record $ map toField [oX, oY, oE]

diffToTauD :: GlobalParams Double -> Diffusivity -> Double
diffToTauD (GParams {..}) d = beamWaist^2 / 4 / d

saveFit :: GlobalParams Double -> Species Double -> FilePath -> Species Double -> IO ()
saveFit gp tauDs fname weights =         
    save fname
    $ fmap (\t->(t, diffusionModel gp tauDs weights t))
    $ logSpace 1e-10 1e0 500

main = do
    let tauDs = V.fromList $ logSpace 1e-10 1e0 50
        gp = GParams 1 1
    (g, ds) <- generateCorr 50
    putStrLn $ "tauDs: "++show (map (diffToTauD gp) ds)
    let obs = fmap (\tauD->Obs tauD (g tauD) 0.1) tauDs
    save "generated" obs

    --let iters = take 100 $ testProjSubgrad 1000 gp tauDs 1.01 obs
    --let iters = take 100 $ testBarrier 1e4 gp tauDs obs
    let iters = take 1000 $ testUnregDescent gp tauDs obs
    forM_ iters $ \i->
        putStrLn $ "chi^2= "++show (chiSq obs (diffusionModel gp tauDs i))
            ++"\tS= "++show (entropy i)
            ++"\tnorm1= "++show (sum i)
    putStrLn $ "\n\nFit weights:"
    putStrLn $ intercalate "\n" (map show $ V.toList (V.zip tauDs (last iters)))
    forM_ (zip [0..] iters) $ \(i,weights)->saveFit gp tauDs ("data/fit-"++show i) weights

{-    
testMirrorDescent :: V.Vector (Obs Double) -> [Double]
testMirror obs =
    let psi = sum . map (\a->a * log a) 
        psiStar = log . sum . map exp
        search = wolfeSearch 0.1 0.01 1 1 f
    in mirrorDescent search psi psiStar df 
-}    

-- | Entropy of the mixing distribution      
entropy :: RealFloat a => Species a -> a      
entropy = sum . fmap f
  where f 0 = 0
        f a = a * log a

-- | Derivative of the entropy        
entropy' :: RealFloat a => Species a -> Species a
entropy' = fmap (\a->log a + 1) 
      
-- | Parameters of the model global to all species
data GlobalParams a = GParams { beamWaist, aspectRatio :: a }   
                    deriving (Show, Functor, Foldable, Traversable)
     
-- | Correlation function
diffusionModel :: RealFrac a => GlobalParams a -> Species a -> Species a -> a -> a
diffusionModel (GParams {..}) tauDs weights tau =
    sum $ V.zipWith (\tauD a->a / (1 + tau/tauD) / (1 + aspectRatio^2 * tau/tauD)) tauDs weights
               
chiSq :: RealFloat a => V.Vector (Obs a) -> (a -> a) -> a
chiSq obs f = sum (fmap (\(Obs x y e)->(f x - y)^2 / e^2) obs) / n
  where n = realToFrac $ V.length obs
      
-- | Constraint on chi^2
-- Satisfied when greater than or equal to 0
chiSqConstr :: RealFloat a => a -> GlobalParams a -> Species a -> V.Vector (Obs a) -> Species a -> a
chiSqConstr chiTol gp tauDs obs weights =
    chiTol - chiSq obs (diffusionModel gp tauDs weights)

chiSqConstrSatisfied :: RealFloat a => a -> GlobalParams a -> Species a
                     -> V.Vector (Obs a) -> Species a -> Bool
chiSqConstrSatisfied chiTol gp tauDs obs weights =
    chiSqConstr chiTol gp tauDs obs weights >= 0

-- | Constraint on weight normalization    
-- Satisfied when equal to 0   
normConstr :: RealFloat a => Species a -> a
normConstr weights = 1 - sum weights

normConstrSatisfied :: (Epsilon a, RealFloat a) => Species a -> Bool
normConstrSatisfied weights = nearZero $ normConstr weights

tr x = traceShow x x           
            
normalizeL1 :: (Fractional a, Functor f, Foldable f) => f a -> f a   
normalizeL1 x = x ^/ sum x

testProjSubgrad :: (Show a, Epsilon a, RealFloat a) => a -> GlobalParams a -> Species a 
                -> a -> V.Vector (Obs a) -> [Species a]
testProjSubgrad chiTol gp tauDs alpha obs =
    let stepSched = constStepSched 5e-4
        x0 = fmap (const $ 1/realToFrac (V.length tauDs)) tauDs
        f = chiSqConstr chiTol gp tauDs obs
        df = grad $ chiSqConstr (auto chiTol) (auto <$> gp) (auto <$> tauDs) (fmap (fmap auto) obs)
        proj a = head $ dropWhile (not . chiSqConstrSatisfied chiTol gp tauDs obs)
                 -- $ map (\f->traceShow (chiSqConstr 0 gp tauDs obs f, V.take 2 f) f)
                 $ map (normalizeL1 . fmap (max 0) . (a ^+^))
                 $ take 1000 $ iterate (alpha *^) $ normalizeL1 $ df a
    in projSubgrad stepSched proj df f x0

testBarrier :: (Show a, RealFloat a) => a -> GlobalParams a -> Species a
            -> V.Vector (Obs a) -> [Species a]
testBarrier chiTol gp tauDs obs =
    let obj :: (Show a, RealFloat a) => a -> GlobalParams a -> Species a
            -> V.Vector (Obs a) -> a -> Species a -> a
        obj chiTol gp tauDs obs mu a =
            entropy a - mu * log (chiSqConstr chiTol gp tauDs obs a) - mu * sum (fmap log a) + 1/mu * (sum a - 1)^2
        obj' chiTol mu = grad $ obj (auto chiTol) (auto <$> gp) (auto <$> tauDs) (fmap auto <$> obs) (auto mu)
        x0 = fmap (const $ 1/realToFrac (V.length tauDs)) tauDs
        search chiTol mu = armijoSearch 0.1 100 0.1 (obj chiTol gp tauDs obs mu)
        mu0 = 100
        chiTols = scanl (\a f->f a) 2e5
                  $ take 10000 (cycle $ (/2):replicate 10 id) ++ repeat id
        mus = scanl (\a f->f a) mu0
              $ take 10000 (cycle $ (/2):replicate 1000 id) ++ repeat id
    in runBlock 40 (\(chiTol,mu)->traceShow (chiTol, mu) $ steepestDescent (search chiTol mu) (obj' chiTol mu)) (zip chiTols mus) x0
    
testUnregDescent :: (RealFloat a) => GlobalParams a -> Species a
            -> V.Vector (Obs a) -> [Species a]
testUnregDescent gp tauDs obs =
    let obj :: (RealFloat a) => a -> GlobalParams a -> Species a
            -> V.Vector (Obs a) -> Species a -> a
        obj mu gp tauDs obs a = chiSq obs (diffusionModel gp tauDs a) - mu * sum (fmap log a) + 1/mu * (sum a - 1)^8
        obj' mu = grad $ obj (auto mu) (auto <$> gp) (auto <$> tauDs) (fmap auto <$> obs)
        x0 = fmap (const $ 1/realToFrac (V.length tauDs)) tauDs
        search mu = armijoSearch 0.1 100 0.1 (obj mu gp tauDs obs)
        mu0 = 1
    in runBlock 40 (\mu->steepestDescent (search mu) (obj' mu)) (iterate (0.5*) mu0) x0

runBlock :: Int -> (a -> b -> [b]) -> [a] -> b -> [b]
runBlock n f (x:xs) y0 =
    let ys = f x y0
    in take n ys ++ runBlock n f xs (head $ drop n ys)
