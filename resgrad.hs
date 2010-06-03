-- Copyright (C) 2010 Riad S Wahby <rsw@jfet.org>
-- 
-- This file is part of resgrad
--
--  resgrad is free software.  It comes without any warranty, to
--  to the extent permitted by applicable law.  You can redistribute it
--  and/or modify it under the terms of the Do What The Fuck You Want To
--  Public License, Version 2, as published by Sam Hocevar.  See
--  the COPYING file or http://sam.zoy.org/wtfpl/COPYING for more details
--

{-
 - When attempting to make a pair of resistors with a large ratio, it is
 - optimal to devote equivalent area to each resistor.  As such, the most
 - obviously ideal ratios are those which are perfect squares; this allows
 - one string of sqrt(R) in series, and one string of sqrt(R) in parallel,
 - devoting equal area to each resistor and resulting in a ratio R.
 -
 - Unfortunately, laying out such resistors to be gradient insensitive is
 - not as simple as when matching two resistors of equal value: after some
 - thought it should be obvious that common centroiding does not work to
 - cancel the effect of gradients, since one resistor accumulates errors
 - linearly and the other as 1/delta.  As a result, the "correct" layout for
 - a given set of resistors can be very nonintuitive.
 -
 - This software performs 2d optimization of a resistor pack in the presence
 - of linear gradients in X and Y of specified magnitude.  It is assumed
 - that the resistor ratio is equal to nUnits^2, i.e., the resistors will be
 - nUnits in parallel and in series.  The results are written to the present
 - directory in the form of the nOut best configurations.
 -
 - Usage:
 -
 - resgrad <nUnits> <nRows> <deltaY> <deltaX> <nOut>
 -   nUnits : the number of units per resistor
 -   nRows  : the number of rows into which to split
 -     NOTE: nRows must evenly divide 2*nUnits!
 -   deltaY : deltaRho/Rho between top-most and bottom-most row
 -   deltaX : deltaRho/Rho between left-most and right-most row
 -   nOut   : output the nOut best configurations
 -}

module Main where

import Data.List (foldl', sortBy)
import Control.Parallel.Strategies (NFData(..),rnf)
import Control.Parallel (par)
import System (getArgs,getProgName,exitFailure)

-- compute the gain given dR and a specification of resistors by row
resGainL dR nSx = seriesValue / parallelValue
        where lnSx = length nSx
              dRs = take lnSx $ zipWith (\x y -> 1 + x*y) (repeat dR) [0..]
              rMix = \x -> (\n d -> x $ take n $ repeat d)
              seriesValue = seriesRes $ zipWith (rMix seriesRes) (map fst nSx) dRs
              parallelValue = parallelRes $ zipWith (rMix parallelRes) (map snd nSx) dRs
              seriesRes = sum
              parallelRes = (1/) . ( sum . ( map (1/) ) )

-- generate every possible row specification for a certain resistor configuration
genRows _    _    0 = []
genRows tRow tRes 1 = [[(tRes,tRow-tRes)]]
genRows tRow tRes n = concat $ flip map [(max 0 (tRes-(n-1)*tRow))..(min tRow tRes)]
                                        (\x -> mingle (x,tRow-x) $ genRows tRow (tRes-x) (n-1))
    where mingle x = map (x:)

-- choose the optimal (1-d) resistor config for a given resistor configuration and dR
optimalN tRow tRes nRow dR =
         foldl1' minErr (zipGains dR (fromIntegral $ tRes^2) (genRows tRow tRes nRow))

-- create a list of resistor configurations sorted by gradient insensitivity for a given dR
optimalNLs tRow tRes nRow dR =
           sortBy cmpErr (zipGains dR (fromIntegral $ tRes^2) (genRows tRow tRes nRow))

-- create a list of 2-d resistor configurations (up to nx*ny in length)
-- optimizing for combined gradient insensitivity in X and Y directions
-- nx and ny should be on the order of 5ish probably
optimalNLs2dS tRow tRes nRow dRx dRy nP = rnf xs `par` rnf ys `seq` take nP $ combPerpsS xs ys
    where xs = take (2*nP) $ optimalNLs tRow tRes nRow (dRx / (fromIntegral $ nRow - 1))
          ys = take (2*nP) $ optimalNLs nRow tRes tRow (dRy / (fromIntegral $ tRow - 1))

-- zip together config and corresponding gain error, normalizing to dVal (the expected gain)
zipGains dR dVal rows = rnf v1 `par` rnf v2 `seq` zipX rows v1 v2
    where zipX (r:rs) (v1:v1s) v2s = (r,v1) : zipX rs v1s v2s
          zipX (r:rs) []  (v2:v2s) = (r,v2) : zipX rs [] v2s
          zipX _      []  []       = []
          zipX []     _   _        = []
          lR = length rows
          genValues = map $ (1+) . abs . (1-) . (/dVal) . resGainL dR
          (r1,r2) = splitAt (lR`div`2) rows
          v1 = (genValues r1) :: [Double]
          v2 = (genValues r2) :: [Double]

-- combine one spec with a list of perpendicular ones which may or may not be compatible
combPerp _       []           = []
combPerp (rx,dx) ((ry,dy):ys) 
           | perpCompat rx ry = (rx,ry,dx*dy) : combPerp (rx,dx) ys
           | otherwise        = combPerp (rx,dx) ys

-- combine two lists of perpendicular specs
combPerps xs ys = concat $ map (flip combPerp ys) xs

-- as above, but sort the result to get the best joint gradient insensitivity
combPerpsS xs ys = sortBy cmpErrC $ combPerps xs ys

-- select minimum absolute error
minErr d1@(_,a) d2@(_,b) | a-1 < b-1 = d1
                         | otherwise = d2

cmpErr (_,a) (_,b) | a-1 < b-1 = LT
                   | otherwise = GT

cmpErrC (_,_,a) (_,_,b) | a-1 < b-1 = LT
                        | otherwise = GT

-- foldl1, strictly to prevent stack overflow
foldl1' _ []     = error "empty list"
foldl1' f (x:xs) = foldl' f x xs

-- are these two perpendicular specs compatible?
perpCompat []                _   = True
perpCompat pX@(x@(x1,x2):xs) pY@((y1,y2):_)
            | length pX /= y1+y2 = False
            | length pY /= x1+x2 = False
            | otherwise          = perpCompat xs $ reducePerp x pY []
    where reducePerp (0,0)   []           npY = reverse npY
          reducePerp (x1,x2) ((y1,y2):ys) npY 
                     | x2 > x1 && x2 /= 0 && y2 /= 0 = reducePerp (x1,x2-1) ys $ (y1,y2-1):npY
                     | x1 > x2 && x1 /= 0 && y1 /= 0 = reducePerp (x1-1,x2) ys $ (y1-1,y2):npY
                     | x2 /= 0 && y2 /= 0 = reducePerp (x1,x2-1) ys $ (y1,y2-1):npY
                     | x1 /= 0 && y1 /= 0 = reducePerp (x1-1,x2) ys $ (y1-1,y2):npY
                     | otherwise          = [(-1,-1)]

optSizes = [(14,14,2),(10,15,3),(7,14,4),(6,15,5),(5,15,6),(4,14,7),(4,16,8),(4,18,9),(3,15,10)]
optimals dR = map (\(x,y,z) -> optimalN x y z dR) optSizes

{-- SVG STUFF --}

xMargin = 10
yMargin = 10
deltaX = 20 
deltaY = 60
rectX = 18
rectY = 58
fSize = 30
parColor = "purple"
serColor = "silver"

svgHeader fName w h = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
                      ++show w++"\" height=\""++show h++"\"><title>"++fName++"</title>"
svgTail = "</svg>"

svgBox xos yos color = 
   "<rect x=\""++show (xMargin+xos)++"\" y=\""++show (yMargin+yos)++"\" width=\""
   ++show rectX++"\" height=\""++show rectY++"\" style=\"fill: "
   ++color++"; stroke: "++color++";\"/>\n"

svgArray []         x y dR = svgShowdR dR (x+deltaX) (y+deltaY)
svgArray ((0,0):ls) x y dR = svgArray ls (x+deltaX) 0 dR
svgArray ((0,m):ls) x y dR = svgBox x y parColor ++ svgArray ((0,m-1):ls) x (y+deltaY) dR
svgArray ((l,m):ls) x y dR = svgBox x y serColor ++ svgArray ((l-1,m):ls) x (y+deltaY) dR

svgShowdR dR xos yos =
   "<text x=\""++show xos++"\" y=\""++show yos++"\" font-size=\""++show fSize++"\">"
   ++(show $ (fromIntegral (round (dR*1e12)))/1e12)++"</text>"

svgArray2 pX pY dR = svgArrayP pX pY yMargin
    where svgArrayP _ [(-1,-1)] _ = []  -- error signaled from reducePerp, shouldn't happen!
          svgArrayP []     _  yos = svgShowdR dR (2*xMargin) (yos+deltaY)
          svgArrayP (x:xs) pY yos = xLine ++ svgArrayP xs (reverse npY) (yos+deltaY)
              where (xLine,npY) = reducePerp x pY xMargin yos ([],[])
          reducePerp (0,0)   _            _  _   xnp         = xnp
          reducePerp (x1,x2) ((y1,y2):ys) xP yos (xLine,npY)
            | x2 > x1 && x2 /= 0 && y2 /= 0 =
              reducePerp (x1,x2-1) ys (xP+deltaX) yos ((svgBox xP yos parColor)++xLine,(y1,y2-1):npY)
            | x1 > x2 && x1 /= 0 && y1 /= 0 =
              reducePerp (x1-1,x2) ys (xP+deltaX) yos ((svgBox xP yos serColor)++xLine,(y1-1,y2):npY)
            | x2 /= 0 && y2 /= 0 =
              reducePerp (x1,x2-1) ys (xP+deltaX) yos ((svgBox xP yos parColor)++xLine,(y1,y2-1):npY)
            | x1 /= 0 && y1 /= 0 =
              reducePerp (x1-1,x2) ys (xP+deltaX) yos ((svgBox xP yos serColor)++xLine,(y1-1,y2):npY)
            | otherwise          = ("",[(-1,-1)])

svg1d tRow tRes nRow (x,dR) = svgHeader fName (2*xMargin+tRow*deltaX+14*fSize) (2*yMargin+nRow*deltaY) ++
                              svgArray x 0 0 dR ++
                              svgTail
        where fName = show nRow ++ "_" ++ show tRes ++ "_" ++ show dR ++ ".svg"

svg2d tRow tRes nRow (x,y,dR) = svgHeader fName (max (2*xMargin+tRow*deltaX) (14*fSize)) (2*yMargin+nRow*deltaY+2*fSize) ++
                                svgArray2 x y dR ++
                                svgTail
        where fName = show nRow ++ "_" ++ show tRes ++ "_" ++ show dR ++ ".svg"

usage name = "Usage: resgrad <nUnits> <nRows> <deltaY> <deltaX> <nOut> [+RTS -N<x>]\n"
  ++"   nUnits : the number of units per resistor\n"
  ++"   nRows  : the number of rows into which to split\n"
  ++"     NOTE: nRows must evenly divide 2*nUnits!\n"
  ++"   deltaY : deltaRho/Rho between top-most and bottom-most row\n"
  ++"   deltaX : deltaRho/Rho between left-most and right-most row\n"
  ++"   nOut   : output the nOut best configurations\n"
  ++"   x      : number of processor cores on this machine\n"

main = getArgs >>= \argv ->
  if (length argv /= 5)
     then getProgName >>= (\name -> putStr $ usage name) >> exitFailure
     else let nUnits = read $ argv !! 0 :: Int
              nRows  = read $ argv !! 1 :: Int
              deltaX = read $ argv !! 2 :: Double
              deltaY = read $ argv !! 3 :: Double
              nOut   = read $ argv !! 4 :: Int
              fNameB = argv!!1 ++ "_" ++ argv!!0 ++ "_"
              nPRow  = 2*nUnits`div`nRows
          in case ((2*nUnits`mod`nRows /= 0),nRows) of
                  (True ,_) -> putStrLn "Error: nRows must evenly divide 2*nUnits!" >> exitFailure
                  (False,1) -> do mapM_ (\xdR@(_,dR) -> writeFile (fNameB ++ show dR ++ ".svg") $
                                                                  svg1d nPRow nUnits nRows xdR) $
                                        take nOut $ optimalNLs nRows nUnits nPRow deltaY
                  (_    ,_) -> do mapM_ (\xydR@(_,_,dR) -> writeFile (fNameB ++ show dR ++ ".svg") $
                                                                     svg2d nPRow nUnits nRows xydR) $
                                        optimalNLs2dS nPRow nUnits nRows deltaX deltaY nOut
