(ns rosalind.solutions)

(defn expected-offspring [& couples]
  "Calculating Expected Offspring (http://rosalind.info/problems/iev/)"
  (reduce + (map (partial * 2) couples [1 1 1 0.75 0.5 0])))
