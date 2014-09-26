import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util

main :: IO ()
main = shakeArgs shakeOptions{shakeVerbosity=Chatty} $ do
  want ["RandomNumberPerformanceExpanded.lhs"]

  "RandomNumberPerformanceExpanded.lhs" *> \out -> do
    need ["RandomNumberPerformance.lhs"]
    cmd "pandoc -s RandomNumberPerformance.lhs --filter=./Include -t markdown+lhs -o" out