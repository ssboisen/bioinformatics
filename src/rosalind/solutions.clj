(ns rosalind.solutions
  (:use [rosalind.core])
  (:require [clojure.string :as string])
  (:require [clojure.math.combinatorics :as comb]))

(defn get-params [file]
  (let [file-content (slurp file)]
    (string/split file-content #"\n")))


(defn partial-right [f & partial-args]
  (fn [& args]
    (apply f (concat args partial-args))))
