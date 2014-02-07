(ns rosalind.solutions
  (:require [clojure.math.combinatorics :as combo])
  (:require [clojure.string :as string])
  (:require [clojure.set :as sets])
  (:require [clojure.pprint :as p])
  (:use [rosalind.core]))

(defn transition-transversion-ratio [s1 s2]
  "Transitions and Transversions (http://rosalind.info/problems/tran/)"
  (let [checkmutation (fn [c1 c2 exp1 exp2]
                        (or
                         (and (= c1 exp1) (= c2 exp2))
                         (and (= c2 exp1) (= c1 exp2))))
        mutations (map (fn [c1 c2]
                           (cond
                            (= c1 c2) [0 0]
                            (checkmutation c1 c2 \A \G) [1 0]
                            (checkmutation c1 c2 \T \C) [1 0]
                            :else [0 1])) s1 s2)
        transitions (map first mutations)
        transversions (map second mutations)
        total-transitions (reduce + transitions)
        total-transversions (reduce + transversions)]
    (float (/ total-transitions total-transversions))))

(apply transition-transversion-ratio (map second (parse-fasta (slurp "input.txt"))))