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

main = do
    let tauDs = V.fromList $ logSpace 1e-10 1e0 50
        gp = GParams 1 1
    (g, ds) <- generateCorr 50
    putStrLn $ "tauDs: "++show (map (diffToTauD gp) ds)
    let obs = fmap (\tauD->Obs tauD (g tauD) 0.1) tauDs
    save "generated" obs

    --let iters = take 100 $ testProjSubgrad 1000 gp tauDs 1.01 obs
    let iters = take 100 $ testBarrier 1000 gp tauDs obs
    --let iters = take 100 $ testUnregDescent gp tauDs obs
    forM_ iters $ \i->
        putStrLn $ "chi^2="++show (chiSq obs (diffusionModel gp tauDs i))++"\tS="++show (entropy i)
    putStrLn $ "\n\nFit weights:"
    putStrLn $ intercalate "\n" (map show $ V.toList (V.zip tauDs (last iters)))
    save "fit" $ fmap (\t->(t, diffusionModel gp tauDs (last iters) t)) $ logSpace 1e-10 1e0 500

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
            
testProjSubgrad :: (Show a, Epsilon a, RealFloat a) => a -> GlobalParams a -> Species a 
                -> a -> V.Vector (Obs a) -> [Species a]
testProjSubgrad chiTol gp tauDs alpha obs =
    let stepSched = constStepSched 5e-4
        x0 = fmap (const $ 1/realToFrac (V.length tauDs)) tauDs
        f = chiSqConstr chiTol gp tauDs obs
        df = grad $ chiSqConstr (auto chiTol) (auto <$> gp) (auto <$> tauDs) (fmap (fmap auto) obs)
        proj a = head $ dropWhile (not . chiSqConstrSatisfied chiTol gp tauDs obs)
                 -- $ map (\f->traceShow (chiSqConstr 0 gp tauDs obs f, V.take 2 f) f)
                 $ map (normalize . fmap (max 0) . (a ^+^))
                 $ take 1000 $ iterate (alpha *^) $ normalize $ df a
    in projSubgrad stepSched proj df f x0

testBarrier :: (RealFloat a) => a -> GlobalParams a -> Species a
            -> V.Vector (Obs a) -> [Species a]
testBarrier chiTol gp tauDs obs =
    let obj :: (RealFloat a) => a -> GlobalParams a -> Species a
            -> V.Vector (Obs a) -> a -> Species a -> a
        obj chiTol gp tauDs obs mu a =
            entropy a + mu * (log (chiSqConstr chiTol gp tauDs obs a) + sum (fmap log a) + log (recip $ (sum a)^2))
        obj' mu = grad $ obj (auto chiTol) (auto <$> gp) (auto <$> tauDs) (fmap auto <$> obs) (auto mu)
        x0 = fmap (const $ 1/realToFrac (V.length tauDs)) tauDs
        search mu = armijoSearch 0.5 100 0.1 (obj chiTol gp tauDs obs mu)
        go n mu x0 = let xs = steepestDescent (search mu) (obj' mu) x0
                   in take n xs++go n (2*mu) (head $ drop n xs)
    in go 10 0.1 x0

testUnregDescent :: (RealFloat a) => GlobalParams a -> Species a
            -> V.Vector (Obs a) -> [Species a]
testUnregDescent gp tauDs obs =
    let obj :: (RealFloat a) => a -> GlobalParams a -> Species a
            -> V.Vector (Obs a) -> Species a -> a
        obj mu gp tauDs obs a = chiSq obs (diffusionModel gp tauDs a) + mu * sum (fmap log a)
        obj' mu = grad $ obj (auto mu) (auto <$> gp) (auto <$> tauDs) (fmap auto <$> obs)
        x0 = fmap (const $ 1/realToFrac (V.length tauDs)) tauDs
        search = armijoSearch 0.5 100 0.1 (obj mu gp tauDs obs)
        mu = 100
    in steepestDescent search (obj' mu) x0