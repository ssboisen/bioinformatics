(ns rosalind.solutions
  (:require [clojure.math.combinatorics :as combo])
  (:require [clojure.string :as string])
  (:require [clojure.set :as sets])
  (:require [clojure.pprint :as p])
  (:use [rosalind.core]))

(defn signed-permutations [n]
  "Enumerating Oriented Gene Orderings (http://rosalind.info/problems/sign/)"
  (let [permutations (combo/selections (filter #(not= % 0) (range (- n) (inc n))) n)]
    (filter #(let [xs (map abs %)]
               (= (distinct xs) xs)) permutations)))