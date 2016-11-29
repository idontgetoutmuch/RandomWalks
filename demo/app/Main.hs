{-# LANGUAGE
  TypeFamilies
 #-}

module Main where

import Control.Monad.Bayes.LogDomain
import Control.Monad.Bayes.Primitive
import Control.Monad.Bayes.Class
import Control.Monad.Bayes.Population
import Control.Monad.Bayes.Conditional
import Control.Monad.Bayes.Inference

-- number of particles used in PMMH
n_particles = 10 :: Int

-- model data and constants
h = 0.1
obs = []

-- type of simulation state
data S = S {p :: Double, z :: Double, log_alpha :: Double}

-- prior over model parameters
parameters :: (MonadBayes m, CustomReal m ~ Double) => m (Double, Double)
parameters = do
  mu <- uniform 0 1
  sigma <- uniform 0 5
  return (mu,sigma)

-- initial state of the simulation
initial_state :: (MonadBayes m, CustomReal m ~ Double) => (Double, Double) -> m S
initial_state (mu, sigma) = do
  log_p_init <- normal (log 100) 0.2
  log_z_init <- normal (log  50) 0.1
  log_alpha_init <- normal (log mu) sigma
  return (S (exp log_p_init) (exp log_z_init) log_alpha_init)

-- transition model
transition :: (MonadBayes m, CustomReal m ~ Double) => (Double, Double) -> S -> m S
transition params state = do
  w <- normal 0 (sqrt h)
  -- TODO: ODE solver updates state here
  return state

-- full simulation given parameters, returns the parameters
model :: (MonadBayes m, CustomReal m ~ Double) => (Double, Double) -> m (Double, Double)
model params = foldl step (initial_state params) obs >> return params where
  step old_state p_obs = do
    state <- old_state
    new_state <- transition params state
    observe (Continuous (Normal (log (p state)) 0.1)) (log p_obs)
    return new_state

-- full model with particle filter, returns posterior over model parameters
full_model :: (MonadBayes m, CustomReal m ~ Double) => m (Double, Double)
full_model = parameters >>= (collapse . smc (length obs) n_particles . model)

-- PMMH transition kernel
-- monad-bayes does not currently have truncated normals
pmmh_kernel :: (MonadBayes m, CustomReal m ~ Double) => [Double] -> m [Double]
pmmh_kernel [mu, sigma] = do
  mu' <- normal mu 1
  sigma' <- normal sigma 1
  return [mu', sigma']

-- full PMMH transition step
pmmh_step :: (MonadBayes m, CustomReal m ~ Double) => [Double] -> m [Double]
pmmh_step params = do
  params' <- pmmh_kernel params
  let kernel_density  = unsafeContJointDensity (pmmh_kernel params) params'
  let kernel_density' = unsafeContJointDensity (pmmh_kernel params') params
  pm_density  <- pseudoDensity full_model (map Just params , [])
  pm_density' <- pseudoDensity full_model (map Just params', [])
  let mh_ratio = pm_density' * kernel_density' / (pm_density * kernel_density)
  accept <- bernoulli (min 1 (fromLogDomain mh_ratio))
  return (if accept then params' else params)

main :: IO ()
main = return ()
