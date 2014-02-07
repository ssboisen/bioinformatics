(defn indexes-of [s1 s2]
  (loop [index (.indexOf s1 s2)
         acc []]
    (if (= index -1)
      acc
      (recur (.indexOf s1 s2 (inc index)) (conj acc index)))))

(indexes-of "GATATATGCATATACTT" "ATAT")
