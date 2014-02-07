(ns rosalind.solutions
  (:use [rosalind.core]))

(defn protein-mass [s]
  "Calculating Protein Mass (http://rosalind.info/problems/prtm/)"
  (reduce + (map (comp (partial get monoisotopic-mass-table) str) s)))
