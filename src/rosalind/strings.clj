(ns rosalind.solutions
  (:require [clojure.math.combinatorics :as combo])
  (:require [clojure.string :as string])
  (:require [clojure.set :as sets])
  (:require [clojure.pprint :as p])
  (:use [rosalind.core]))

(defn rna-to-amino-acid [rnastring]
  (->> rnastring
       (partition-string 3)
       (map rna-codon-to-protein)
       (take-while (partial not= "Stop"))
       (apply str)))

(defn reverse-complement
  "Complementing a Strand of DNA (http://rosalind.info/problems/revc/)"
  [s]
  (string/replace (string/reverse s) #"A|C|G|T" {"A" "T" "T" "A" "C" "G" "G" "C"}))

(defn peptide-encoding-problem [dnastring peptide]
  (let [monomers (partition-string 1 peptide)
        codons (map (fn [monomer]
                      (map key (filter #(= monomer (val %)) dna-codon-table))) monomers)
        dnacombs (map (partial apply str) (apply combo/cartesian-product codons))
        dnacombsrc (map reverse-complement dnacombs)
        dnas (concat dnacombs dnacombsrc)]
    (filter #(not= (.indexOf dnastring %) -1) dnas)))

(defn find-candidates [dnastring]
  "Find candidate protein strings for a dnastring"
  (let [partitions (apply vector (partition 3 dnastring))
        start-indexes (keep-indexed #(if (= (seq "ATG") %2) %1) partitions)
        stop-indexes (keep-indexed #(if (= (get dna-codon-table (apply str %2)) "Stop") %1) partitions)
        candidate-indexes (filter (comp not nil? second) (map (fn [i] [i (some #(if (< i %) %) stop-indexes)]) start-indexes))
        candidates (map (fn [[start stop]] (subvec partitions start stop)) candidate-indexes)]
  (map #(map (partial apply str) %) candidates)))


(defn find-all-candidates [dnastring]
  "Find every distinct candidate protein string for a dnastring"
  (for [i (range 3)
        candidate (distinct
                   (concat (find-candidates (subs dnastring i))
                           (find-candidates (subs (reverse-complement dnastring) i))))]
        (apply str (map (partial get dna-codon-table) candidate))))

(defn strings-of-alphabet [alphabet n]
  "Enumerating k-mers Lexicographically (http://rosalind.info/problems/lexf/)"
  (let [strs (comb/selections alphabet n)]
    (map (partial apply str) strs)))

(defn all-strings-of-alphabet [alphabet n]
  "Ordering Strings of Varying Length Lexicographically (http://rosalind.info/problems/lexv/)"
  (let [sortvals (apply hash-map (apply concat (map-indexed #(vector %2 (+ %1 1 (count alphabet))) alphabet)))
        strs (mapcat (partial comb/selections alphabet) (range 1 (inc n)))
        sortval #(apply str (take n (concat (map (partial get sortvals) %) (repeat 0))))]
    (map (partial apply str) (sort #(compare (sortval %1) (sortval %2)) strs))))

(defn finding-motif [s t]
  "Finding a Motif in DNA (http://rosalind.info/problems/subs/)"
  (let [windows (map #(apply str %) (partition (count t) 1 s))
        motifs (keep-indexed #(if (= %2 t) (inc %1)) windows)]
    motifs))

(defn load-data [file]
  (map  #(apply str (drop 1 (string/split % #"\n")))
        (drop 1 (string/split (slurp file) #">"))))

(def input (load-data "/Users/ssboisen/Downloads/rosalind_lcsm.txt"))

(defn consensus-data [dnastrings]
  "Consensus and Profile (http://rosalind.info/problems/cons/)"
  (let [all-freqs (apply map (fn [& args] (merge {\A 0, \C 0, \G 0, \T 0} (frequencies args))) dnastrings)
        sorted-freqs (map (partial sort-by val >) all-freqs)
        consensus (apply str (map (comp first first) sorted-freqs))
        profile (map (fn [[key freqs]] (format "%s: %s" key (string/join " " (map second freqs)))) (group-by first (mapcat identity all-freqs)))]
    (apply str consensus "\n"
           (string/join "\n" profile))))

(defn shared-motif [dnastrings]
  (let [parts (map (fn [dnastring] (mapcat #(partition % 1 dnastring) (range 1 (inc (count dnastring))))) dnastrings)
        sets (map #(apply hash-set %) parts)]
    (apply str (first (sort-by count > (apply sets/intersection sets))))))

(defn- all-substrings [s]
(set (for [i (range (count s))
           j (range (count s))
           :when (<= i j)]
         (subs s i (inc j)))))

(defn- substring? [s x]
    "Is `x` a substring of `s`?"
    (not= -1 (.indexOf s x)))

(defn- filter-substrings [substring-set s]
    "Find the subset of `substring-set` that are substrings of `s`."
    (set (filter (partial substring? s) substring-set)))

(defn common-substrings
    ([ss] (common-substrings (first ss) (rest ss)))
    ([f ss]
        (let [xs (all-substrings f)]
            (reduce filter-substrings xs ss))))

(defn longest-common-substring [ss]
    (apply max-key count (common-substrings ss)))


(defn overlap-graphs [fasta-pairs k]
  "Overlap Graphs (http://rosalind.info/problems/grph/)"
  (let [string-datas (map
                      (fn [[name dnastring]]
                        {:name name, :dna dnastring :prefix (subs dnastring 0 k) :suffix (subs dnastring (- (count dnastring) k))})
                      fasta-pairs)]
    (for [{suffix :suffix name1 :name sudna :dna} string-datas
          :let [adjacencies (filter (fn [{:keys [prefix dna]}] (and (= suffix prefix) (not= sudna dna))) string-datas)]
          adjacency (map (fn [{:keys [name]}] [name1 name]) adjacencies)]
      adjacency)))


(defn indexes-of [s1 s2]
  (loop [index (.indexOf s1 s2)
         acc []]
    (if (= index -1)
      acc
      (recur (.indexOf s1 s2 (inc index)) (conj acc index)))))

(defn- dnastring-overlaps [dnastrings]
  (for [dna1 dnastrings
        dna2 dnastrings
        :when (not= dna1 dna2)
        :let [halfprefix (subs dna2 0 (int (/ (count dna2) 2)))
              index (.lastIndexOf dna1 halfprefix)]
        :when (not= index -1)]
    {:length (- (count dna1) index) :index index :dna1 dna1 :dna2 dna2}))

(defn- combine-overlapping-strings [combination]
  (apply str
         (:dna1 (first combination))
         (map (fn [{:keys [length dna2]}]
                (subs dna2 length)) combination)))

(defn shortest-superstring [dnastrings]
  "Genome Assembly as Shortest Superstring (http://rosalind.info/problems/long/)"
  (let [overlaps (dnastring-overlaps dnastrings)
        combinations (for [overlap overlaps]
                       (loop [remaining (filter #(not= overlap %) overlaps)
                              combined [overlap]]
                         (let [latest (last combined)
                               {[candidate] true others false} (group-by #(= (:dna2 latest) (:dna1 %)) remaining)]
                           (if (nil? candidate)
                             combined
                             (recur others (conj combined candidate))))))
        valid-combinations (filter #(= (count %) (dec (count dnastrings))) combinations)]
   (sort-by count (map combine-overlapping-strings valid-combinations))))


(defn dna-to-protein [fasta-pairs]
  "RNA Splicing (http://rosalind.info/problems/splc/)"
  (let [[dnastring & introns] (map second fasta-pairs)
        dna-transcribed (reduce #(string/replace %1 %2 "") dnastring introns)
        protein (take-while (partial not= "Stop") (map dna-codon-to-protein (partition-dna dna-transcribed)))]
    (string/join "" protein)))

(defn spliced-motif [fasta-pairs]
  "Finding a spliced motif (http://rosalind.info/problems/sseq/)"
  (let [[s t] (map second fasta-pairs)]
    (second (reduce (fn [[ts idxs] [idx c]]
                      (if (= (first ts) c)
                        [(rest ts) (conj idxs (inc idx))]
                        [ts idxs]
                        ))
                    [t []]
                    (map-indexed vector s)))))

(defn most-frequent-words [text k]
   (let [kmers (partition k 1 text)]
     (map first (second (first (sort-by key > (group-by val (frequencies kmers))))))))

(defn t-frequent-words [text k t]
   (let [kmers (partition k 1 text)]
     (map first (filter (comp (partial = t) val) (frequencies kmers)))))

(def input (slurp "E-coli.txt"))

(defn finding-clumps [text k L t]
  (let [windows (partition L 1 text)
        kmers (filter (complement empty?) (mapcat #(t-frequent-words % k t) windows))]
    (string/join " " (map #(string/join "" %) (distinct kmers)))))

(defn cyclic-subpeptides [cyclicpeptide]
  (let [len (count cyclicpeptide)
        cycled (str cyclicpeptide (subs cyclicpeptide 0 (max 0 (- len 2))))]
    (for [[i c] (map-indexed vector cyclicpeptide)
          l (range 1 len)]
      (subs cycled i (+ i l)))))

(defn linear-subpeptides [linearpeptide]
  (let [len (count linearpeptide)]
    (if (= len 1)
      []
      (for [[i c] (map-indexed vector linearpeptide)
            l (range (+ i 1) (+ 1 len))
            :let [sp (subs linearpeptide i l)]
            :when (not= sp linearpeptide)]
        sp))))


(defn peptide-mass [peptide]
  (apply + (map #(get molecule-weight-table (str %)) peptide)))

(defn theoretical-spectrum [subpeptide-selector cyclicpeptide]
  (let [subpeptides (subpeptide-selector cyclicpeptide)]
    (concat
       [0]
       (map peptide-mass subpeptides)
       [(peptide-mass cyclicpeptide)])))

(defn cyclic-subpeptidemasses [cyclicpeptide]
  (let [len (count cyclicpeptide)
        cycled (vec (concat cyclicpeptide (subvec cyclicpeptide 0 (max 0 (- len 2)))))]
    (for [[i c] (map-indexed vector cyclicpeptide)
          l (range 1 len)]
      (subvec cycled i (+ i l)))))

(defn linear-subpeptidemasses [linearpeptidemasses]
  (let [len (count linearpeptidemasses)]
    (if (= len 1)
      []
      (for [[i c] (map-indexed vector linearpeptidemasses)
            l (range (+ i 1) (+ 1 len))
            :let [sp (subvec linearpeptidemasses i l)]
            :when (not= sp linearpeptidemasses)]
        sp))))

(defn theoretical-spectrummasses [subpeptidemasses-selector peptidemasses]
  (let [subpeptidemasses (subpeptidemasses-selector peptidemasses)]
    (concat
       [0]
       (map (partial apply +) subpeptidemasses)
       [(apply + peptidemasses)])))

(def linear-theoretical-spectrum (partial theoretical-spectrummasses linear-subpeptidemasses))
(def cyclic-theoretical-spectrum (partial theoretical-spectrummasses cyclic-subpeptidemasses))

(defn cyclopeptide-sequencing [experimental-spectrum]
  (let [masses (map val molecule-weight-table)
        spectfreq (frequencies experimental-spectrum)]
    (loop [lst [[]]
           peptides []]
      (if (empty? lst)
        (distinct peptides)
        (let [expanded (for [v lst
                             m masses
                             :let [nv (conj v m)
                                   theoretical-spectrum (linear-theoretical-spectrum nv)
                                   freqs (frequencies theoretical-spectrum)]
                             :when (every? (fn [[k v]] (>= (spectfreq k 0) v)) freqs)]
                         nv)
              peptides (if (empty? expanded) lst expanded)]
            (recur expanded peptides))))))

(spit "output.txt" (string/join " " (map #(string/join "-" %) (cyclopeptide-sequencing [0	97	97	99	101	103	196	198	198	200	202 295	297	299	299	301	394	396	398	400	400	497]))))
(def spectrumm (into #{} [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]))

(defn- remove-first [pred coll]
   (let [[n m] (split-with (complement pred) coll)]
     (concat n (rest m))))

(def timerstate (atom []))
(defmacro timer
  "Evaluates expr and conjes the time it took to state. Returns the value of
 expr."
  [state expr]
  `(let [start# (. System (nanoTime))
         state# ~state
         ret# ~expr]
     (swap! state# conj (/ (double (- (. System (nanoTime)) start#)) 1000000.0))
     ret#))

(reset! timerstate [])
(apply + @timerstate)

(defn spectral-convolution-problem [spectrum]
  (let [sorted (sort spectrum)]
    (sort (for [s sorted
                ss (take (dec (count spectrum)) sorted)
                :let [diff (- s ss)]
                :when (> diff 0)]
            diff))))

(defn probability-of-random-string [N A k t]
  "From https://beta.stepic.org/Bioinformatics-Algorithms-2/Detour-Probabilities-of-Patterns-in-a-String-11/#step-6"
  (letfn [(factorial [n] (reduce *' (range 1 (inc n))))
          (binocoef [m k]
                    (/ (factorial m) (* (factorial k) (factorial (- m k)))))
          (exp [x n]
               (reduce * (repeat n x)))]
    (/ (binocoef (- N (* t (- k 1))) t) (exp A (* (- t 1) k)))))

(defn hamming-distance [s1 s2]
  (if (not= (count s1) (count s2))
    (throw (Exception. "Strings must be of equal length"))
    (count (filter (fn [[a b]] (not= a b)) (map vector s1 s2)))))

(defn motif-enumeration [dnas k d]
  (let [all-kmers (strings-of-alphabet ["T" "A" "G" "C"] k)
        kmers (map #(map (partial apply str) (partition k 1 %)) dnas)]
    (sort (distinct
           (for [a (first kmers)
                 a' (filter #(>= d (hamming-distance a %)) all-kmers)
                 :when (every? #(some (comp (partial >= d) (partial hamming-distance a')) %) (rest kmers))]
             a')))))
